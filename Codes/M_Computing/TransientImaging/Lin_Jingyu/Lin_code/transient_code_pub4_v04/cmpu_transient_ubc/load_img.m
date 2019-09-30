% load image data. input: sel_data,sel_shutter. output: measurements.

%% select data set
if ~exist('sel_data','var')
    sel_data = 1;
end
%Dataset
switch sel_data
    case 1
        img_path = 'Images/Mirrors/2013-01-14.Mirrors-extracted';
    case 2
        img_path = 'Images/Wall/2013-01-12.Room-take1-extracted';
    case 3
        img_path = 'Images/DiscoBall/2013-01-12.Discoball-extracted';
    case 4
        img_path = 'Images/Corner/2013-01-13.Corner-extracted';
    otherwise
        img_path = 'Images/Bottles/2013-01-14.Bottles-extracted';
end % switch
% ubc_path = '../data_ubc';             % ×¢ÊÍ by whm
ubc_path = '/data';               % Ìí¼Ó by whm
img_folder = sprintf('%s/%s',ubc_path,img_path);

img_name = 'capture_take_avg';
img_desc = 'capture_take_desc';

%% Load description

img_folder = ['../../../../Heide/code_and_data' img_folder];        % Ìí¼Ó by whm

data_desc = load( sprintf( '%s/%s.mat', img_folder, img_desc ) );
data_desc = data_desc.data_desc;
fprintf('\ndataset: %s\n', data_desc.dataset);

%Extract info
imagedims = data_desc.imagedims; %width, height
fprintf('image size: %d, %d\n', imagedims(1), imagedims(2));

%Freq
frequencies = data_desc.frequencies; % first step last
frequencies_list = frequencies(1):frequencies(2):frequencies(3);
num_frequencies = length(frequencies_list);
step_frequencies = frequencies(2);
fprintf('frequency: %d - %d, step: %f\n', frequencies(1), frequencies(3), frequencies(2));

%Phases
phases = data_desc.phases;  %List with relative phases
num_phases = size(phases,1);
fprintf('phases: %d, %d.\n', phases(1), phases(2));

%Shutters
shutters = data_desc.shutters; %List with shutters in microseconds
num_shutters = length(shutters);

%% Load image
img_loc = sprintf('%s/%s_datacube.mat', img_folder, img_name);
measurements = load( img_loc );
measurements = measurements.data_cube; 

%Select shutter (HDR not implemented here)
% if ~exist('sel_shutter','var')|| sel_shutter>num_shutters-1 || sel_shutter<1
%     sel_shutter = num_shutters - 1;
% end
sel_shutter = num_shutters - 1;

fprintf('Used shutter %d: %d us\n',sel_shutter,shutters(sel_shutter) )

%image
measurements = measurements(:,:,:,sel_shutter);
measurements = reshape(measurements, imagedims(1), imagedims(2), num_frequencies, num_phases);
fprintf('\n')


