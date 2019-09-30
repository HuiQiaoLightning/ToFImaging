% 最终测试代码时，应该直接调用该函数。在界面下，数据会占用内存。
function mM = LoadTransientDataWithUI
% 加载数据的数据形式：真实数据、仿真数据
%% Load the description
currpath = pwd;     % 获取当前路径
[img_desc, img_folder]=uigetfile({'*.*'},'Load the description');
if img_desc == 0
    mM = [];
    disp('没选择采集的多频率瞬态测量数据');
    return;
else
    data_desc = load( sprintf( '%s%s', img_folder, img_desc ) );
    data_desc = data_desc.data_desc;    
end

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

%% Load the actual image
cd(img_folder);     % 进入采集的数据的路径
[img_name, img_folder] = uigetfile({'*.*'},'Load the actual image');
if img_name == 0
    return;
else
    img_loc = sprintf('%s%s', img_folder, img_name);
    measurements = load( img_loc );
end
measurements = measurements.data_cube;  

%Select shutter (HDR not implemented here)
sel_shutter = num_shutters - 1;
fprintf('Used shutter : %d --> %3.3f\n', sel_shutter, shutters(sel_shutter) )

measurements = reshape( measurements(:,:,:,sel_shutter), imagedims(1), imagedims(2), num_frequencies * num_phases);
measurements = permute(measurements,[2 1 3]); %y, x, n

%Compute in double
measurements = double(measurements(end:-1:1,:,:));   % Heide采集的数据集是上下颠倒的，这里给反转过来

%% Load the calibration matrix
[cal_name, cal_folder] = uigetfile({'*.*'},'Load the calibration matrix');
if cal_name == 0
    return;
else
    matrix_loc = sprintf('%s%s', cal_folder, cal_name);
    M = load( matrix_loc );
end
M = load( matrix_loc );
M = M.calibration_matrix;

%Compute in double
M = double(M);

%% Crop the calibration matrix optionally (invalid measurements)
rem_last = 2;
M = M(:, 1:end - rem_last);

%Debug
% figure;
% imagesc(C), colorbar, title('Measurement matrix');

%Save dims
n = size(M,2); %Timesteps
xs = imagedims(1); %Image Width
ys = imagedims(2); %Image Height

%Sanity check
if( size(M,1) ~= num_frequencies * num_phases ||  size(M,2) ~= n )
    error('Calibration matrix and dataset are not calibrated identically')
end

filename = strfind(img_folder, '\');
filename = [filename(end-2) filename(end-1)];
filename = img_folder(filename(1)+1:filename(2)-1);
mM.filename = filename;
mM.measurements = measurements;
mM.M = M;
mM.imagedims = imagedims;
mM.frequencies = frequencies;
mM.num_frequencies = num_frequencies;
mM.phases = phases;
mM.num_phases = num_phases;
mM.shutters = shutters;

cd(currpath);       % 返回程序运行的当前路径      
disp(['读取采集数据 ' filename ' 结束。'])
return
