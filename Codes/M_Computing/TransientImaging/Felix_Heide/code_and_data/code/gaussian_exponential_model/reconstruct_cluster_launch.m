function []= reconstruct_cluster_launch(job)


% ***Transient image reconstruction for one i and u step***
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%
 
% This function is a launcher function for multiple jobs on MULTIPLE
% MACHINES. Note that the u-step is very computational expensive in this
% version.
%
% PLEASE USE reconstruct_launch.m FOR A FAST APPROXIMATE SOLVE FOR
% JUST THE I-step WHICH RUNS ON A SINGLE MACHINE IN REASONABLE TIME.

% The solver for the u-step depends on Mark Schmidt's MinConf package: 
% http://www.di.ens.fr/~mschmidt/Software/minConf.html
%
%
% Copyright (C) 2013. Felix Heide
% Email: fheide@cs.ubc.ca

%Path
addpath('../');

%Needs minConf !
addpath('./minConf');
addpath('./minConf/minConf');
addpath('./minConf/minFunc');

%Use the files in this folder only if you own the Matlab global
%optimization toolbox.
addpath('./gads');

%% Reconstruction params
num_jobs = 1;  % 128;   %%

%Debug parameters
verbose = 'all'; %Options: 'all', 'brief', 'none'

if nargin < 1
    job = 1;
end
    
%Check job
if job > num_jobs
    fprintf('\n\nJobs exceeded !!! Exiting ( [%d / %d] )\n\n', job, num_jobs)
    return;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Output
output_folder = '.';                                               % 注释 by whm  
%output_folder='../../../../Outputs';               % ../表示相对路径的层次。 % 文件输出 添加 by whm

%Dataset
%img_folder = '../../data/Images/Mirrors/2013-01-14.Mirrors-extracted';
%img_folder = '../../data/Images/Wall/2013-01-12.Room-take1-extracted';
%img_folder = '../../data/Images/DiscoBall/2013-01-12.Discoball-extracted';
img_folder = '../../data/Images/Corner/2013-01-13.Corner-extracted';
%img_folder = '../../data/Images/Bottles/2013-01-14.Bottles-extracted';

img_name = 'capture_take_avg';
img_desc = 'capture_take_desc';

%Calibration
cal_folder = '../../data/Calibration/2013-01-12.Calibration6-results';
cal_name = 'capture_take_massaged_avg-resizeFreq221';

%Debug
fprintf('\n\nLoading img %s/%s.mat with %s.mat\n', img_folder, img_name, img_desc);
fprintf('Loading calib %s/%s.mat\n\n', cal_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load the description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_desc = load( sprintf( '%s/%s.mat', img_folder, img_desc ) );
data_desc = data_desc.data_desc;

%Extract info
imagedims = data_desc.imagedims; %width, height

%Freq
frequencies = data_desc.frequencies; % first step last
frequencies_list = frequencies(1):frequencies(2):frequencies(3);
num_frequencies = length(frequencies_list);
step_frequencies = frequencies(2);

%Phases
phases = data_desc.phases;  %List with relative phases
num_phases = size(phases,1);

%Shutters
shutters = data_desc.shutters; %List with shutters in microseconds
num_shutters = length(shutters);

%%Now load the actual image
img_loc = sprintf('%s/%s_datacube.mat', img_folder, img_name);
measurements = load( img_loc );
measurements = measurements.data_cube;  

%Select shutter (HDR not implemented here)
sel_shutter = num_shutters - 1;
fprintf('############## Used shutter : %d --> %3.3f', sel_shutter, shutters(sel_shutter) )

measurements = reshape( measurements(:,:,:,sel_shutter), imagedims(1), imagedims(2), num_frequencies * num_phases);
measurements = permute(measurements,[2 1 3]); %y, x, n

%Compute in double
measurements = double(measurements);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now load the calibration matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix_loc = sprintf('%s/%s.mat', cal_folder, cal_name);
M = load( matrix_loc );
M = M.calibration_matrix;

%Compute in double
M = double(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now crop the calibration matrix optionally (invalid measurements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rem_last = 2;
M = M(:, 1:end - rem_last);

%Debug
figure;
imagesc(M), colorbar, title('Measurement matrix');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save dims

n = size(M,2); %Timesteps
xs = imagedims(1); %Image Width
ys = imagedims(2); %Image Height

%Sanity check
if( size(M,1) ~= num_frequencies * num_phases ||  size(M,2) ~= n )
    error('Calibration matrix and dataset are not calibrated identically')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Solve I-subproblem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%% Debug
fprintf('\n\nProcessing i-subproblem full\n\n')

%Reconstruction

%Start timing
tic;

%% Reconstruct using primal-dual solver
lambda_residual = 1;
lambda_temporal = 5.0;
lambda_spatial = 0.3;
thresh_huber = 0.05;

I_reconstructed = pd_solve_full_svd_coherence(measurements, M, lambda_residual, lambda_temporal, lambda_spatial, thresh_huber, 50, 1e-5, verbose);

tt = toc;
disp(['  First subproblem took ' num2str(tt) ' secs']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve U-subproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare for subproblem processing (simply parallelize along the rows here)
num_total_pixels = size(I_reconstructed,1);
num_job_pixels = ceil( num_total_pixels / num_jobs );
start_y = (job - 1)*num_job_pixels + 1;
end_y = min( start_y + num_job_pixels - 1, size(I_reconstructed,1) );

%Crop
measurements = measurements(start_y:end_y,:,:);
I_reconstructed = I_reconstructed(start_y:end_y,:,:);

%% Debug
fprintf('\n\nProcessing u-subproblem job [%d / %d]\n\n', job, num_jobs)
fprintf('Processing rows [%d - %d]\n', start_y, end_y)

%% Store old solution
I_reconstructed_prev = I_reconstructed;

%% Find the local extrema for the constraints

%Cell array with nnz
nnz_positions = cell(size(I_reconstructed,1), size(I_reconstructed,2));
ymax_all = cell(size(I_reconstructed,1), size(I_reconstructed,2)); 

for x = 1:size(I_reconstructed,2)
    for y = 1:size(I_reconstructed,1)

        %fprintf('Finding zeros %d,%d', x, y)

        %Reshape
        ref_rec = reshape(I_reconstructed(y,x,:),1,[]);
        %meas_curr = reshape(measurements(y,x,:),[],1);

        %Print out maxes
        [ymax,imax,ymin,imin] = extrema(ref_rec);
        [imax,IX] = sort(imax);
        ymax = ymax(IX);

        %Thresholding
        imax = imax(ymax > 0.1);
        ymax = ymax(ymax > 0.1);

        ymax = ymax(imax > 30 & imax < n - 30);
        imax = imax(imax > 30 & imax < n - 30);

        %Save
        nnz_positions{y, x} = imax';
        ymax_all{y,x} = ymax';
    end
end

%% Reconstruct exponentials with spatially regularized solver

%Regularization parameters
beta_pos = 0.00001;
sigma_gauss = 3;
rho = 50.0; 

%Debug
fprintf('\n\n Gaussian/Exponential reconstruction with RHO=%3.3f of %s/%s: \n\n', rho, img_folder, img_name)

%Start timing
tic;

[I_reconstructed, U_reconstructed] =  ...
    non_linear_gaussian_exponential_independent(measurements, M, nnz_positions, ymax_all, sigma_gauss, beta_pos, rho, I_reconstructed_prev, verbose);

tt = toc;
disp([' Gaussian/exponential reconstruction took ' num2str(tt) ' secs']);

%Save the result for the current step
% save( sprintf('%s/reconstruction_I_job%03d.mat', output_folder, job), 'I_reconstructed' );
% save( sprintf('%s/reconstruction_U_job%03d.mat', output_folder, job), 'U_reconstructed' );

save( sprintf('%s/reconstruction_I_job%03d.mat', output_folder, job), 'I_reconstructed','tt' );
save( sprintf('%s/reconstruction_U_job%03d.mat', output_folder, job), 'U_reconstructed','tt' );

%Debug
disp('Done the reconstruction for current job.');

return;
