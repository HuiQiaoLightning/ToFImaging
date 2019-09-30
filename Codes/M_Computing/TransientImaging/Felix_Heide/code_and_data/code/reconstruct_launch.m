% ***Transient image reconstruction for one i-step***
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%
 
% This function is a launcher function for a single i-step. 
% This gives already a fairly good solution and is reasonably fast to compute.
% For a very efficient implementation using the GPU look at our project
% webpage (soon): http://www.cs.ubc.ca/labs/imager/tr/2013/TransientPMD/

% Copyright (C) 2013. Felix Heide
% Email: fheide@cs.ubc.ca

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DATASET and configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Debug parameters
verbose = 'all'; %Options: 'all', 'brief', 'none'

%Output
% output_folder = '.';                          % 注释 by whm  
output_folder='../../../Outputs';               % ../表示相对路径的层次。 % 文件输出 添加 by whm   


%Dataset
% img_folder = '../data/Images/Mirrors/2013-01-14.Mirrors-extracted';
% img_folder = '../data/Images/Wall/2013-01-12.Room-take1-extracted';
% img_folder = '../data/Images/DiscoBall/2013-01-12.Discoball-extracted';
% img_folder = '../data/Images/Corner/2013-01-13.Corner-extracted';
img_folder = '../data/Images/Bottles/2013-01-14.Bottles-extracted';


img_name = 'capture_take_avg';
img_desc = 'capture_take_desc';

%Calibration
cal_folder = '../data/Calibration/2013-01-12.Calibration6-results';
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
%%%% Show the result and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Show transient image
scaling = 1.0; %Amplitude scale
I_rec = max(I_reconstructed,0);
I_rec = I_rec./max(I_rec(:));
I_rec = min(scaling * I_rec,1);
I_rec = (I_rec * 1.0) .^ (1/2.2); %Gamma
scalefactor = 3;

%Show transient image
figure();
for t = 1:size(I_rec,3) - 100

    %Image
    I_show = imresize( I_rec(:,:,t), scalefactor,'bilinear');

    %Flip
    I_show = I_show(end:-1:1, :);

    %Show image
    imshow( I_show ), title(sprintf('Transient image frame %d', t));
    pause(0.1)

end

%Save the result
save( sprintf('%s/reconstruction_I.mat', output_folder), 'I_reconstructed' );

%Debug
disp('Done the reconstruction for current job.');
