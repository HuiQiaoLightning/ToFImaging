%{ 
MatlabComputing_GUI：
1. Select a data set from the pop-up menu, 
2. Select experiment from the pop-up menu, 一个实验内容对应一个窗体

命名规则：函数名小写，对应被调用函数名各单词首字母大写
n,k：frame的计数，n表示计数，k为总的帧数
xs，ys：像素的行数和列数，x，y：像素的行计数和列计数
%}
function Matlab_Computing
% 初始化
init;

% 主窗体
MComputing = figure('MenuBar','None','Name','M Computing','NumberTitle','off');

% File菜单及其子项
mfile = uimenu(MComputing,'Label','File');

uimenu(mfile,'Label','Open_img','Callback',@open_img_Callback);
uimenu(mfile,'Label','Acquired_M_ToF(多频率信号)','Callback',@acquired_M_ToF_Callback, 'Separator','on');
uimenu(mfile,'Label','Acquired_S_ToF(单频率信号)','Callback',@acquired_S_ToF_Callback);
uimenu(mfile,'Label','Simu_data(时域信号)','Callback',@simu_data_Callback);

uimenu(mfile,'Label','Append_path','Callback',@append_path_Callback, 'Separator','on');
uimenu(mfile,'Label','Close_others','Callback',@close_others_Callback, 'Separator','on');
uimenu(mfile,'Label','Close_all','Callback',@close_all_Callback);
uimenu(mfile,'Label','Quit','Callback',@exit_Callback,... 
    'Separator','on','Accelerator','Q');

%{
其它可能用的菜单项    
uimenu(mfile,'Label','Save','Callback','disp(''save'')');

mh = uimenu(f,'Label','Find'); 
frh = uimenu(mh,'Label','Find and Replace ...',...
            'Callback','disp(''goto'')');
frh = uimenu(mh,'Label','Variable');                 
uimenu(frh,'Label','Name...', ...
          'Callback','disp(''variable'')');

      
uimenu(frh,'Label','Value...', ...
          'Callback','disp(''value'')');
%}    
%%
% Create the UICONTEXTMENU  建立上下文菜单，窗口中的右键菜单
cmenu = uicontextmenu;

% Create the parent menu  父菜单
closemenu = uimenu(cmenu,'label','Close');

% Create the submenus 子菜单
close_others = uimenu(closemenu,'label','Close_others',...
    'Callback',@close_others_Callback);
close_all = uimenu(closemenu,'label',...
    'Close_all','Callback',@close_all_Callback);
exit = uimenu(closemenu,'label','Quit',...
    'Callback',@exit_Callback);           
MComputing.UIContextMenu = cmenu;


% Lab菜单及其子项
mlab = uimenu(MComputing,'Label','Lab');
uimenu(mlab,'Label','Transient_imaging','Callback',@transient_imaging_Callback);
uimenu(mlab,'Label','Conventional_imaging','Callback',@conventional_imaging_Callback);
end

function init
%% 用数组的方式添加路径
needpath = genpath( './TransientImaging/Experiments' );
needpath = ["./TransientImaging"; needpath];
for n = 1:length(needpath)
    addpath (needpath(n));
end
% addpath ( './TransientImaging' );                            % 添加瞬态成像代码所在的路径
% addpath (genpath( './TransientImaging/Experiments' ));       % 添加瞬态成像实验代码所在的路径及其子路径

%% 设置数据结构
% mdata为应用程序的数据 使用 getappdata 和 setappdata 传递数据
% mdata.savefile = 'lastpath.mat'; 
mdata.filename = [];      % 图像(数据)文件名
mdata.images = [];        % 图像数据

mdata.acquired_M = [];    % 通过GUI把采集的多频率数据（如Heide的数据集）读入到matlabworkspace，实验时，不用每次都重新输入
mdata.acquired_S = [];    % 通过GUI把采集的单频率数据（如Kadambi的数据集）读入到matlabworkspace，实验时，不用每次都重新输入
mdata.simu = [];        % 通过GUI把仿真的数据读入到matlabworkspace，实验时，不用每次都重新输入
mdata.selected = [];    % 用于保存选择的采集数据或仿真数据
mdata.path = needpath;  % 程序退出时删除相应的路径

setappdata(0,'data',mdata);     % 数据读入到matlabworkspace

%% 画图设置
global lw fs                   % 定义画图时的线宽和字体
lw = 4;          %  绘图的线宽           3; 2; 1;          
fs = 19;         %  绘图标注的字体大小   19;14;11;  

%% 临时数据输出位置
global output_folder    
output_folder = '../TempData';

end

%% File
function open_img_Callback(hObject, eventdata)
mdata = getappdata(0,'data');
%%  打开传统图像
[img_desc, img_folder] = uigetfile({'*.*'},'Open an image');
if ~img_desc == 0
    images = imread( sprintf( '%s%s', img_folder, img_desc ) );
    mdata.images = images;
    setappdata(0,'data',mdata);
else
    disp('没选择图像呐...');
end
end

function acquired_M_ToF_Callback(hObject, eventdata)
mM = Load_Transient_Data_With_UI;           % 获取采集的数据
if ~isempty(mM)
    mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
    mdata.acquired_M = mM;
    setappdata(0,'data',mdata);     % 数据写入到matlabworkspace
end
end

function acquired_S_ToF_Callback(hObject, eventdata)
[img_desc, img_folder] = uigetfile({'*.*'},'Load cross-correlation function from camera');
if img_desc == 0
    disp('没选择采集的单频率瞬态测量数据');
    return;
else
    rawMean = load( sprintf( '%s%s', img_folder, img_desc ) );
    data_xcorr.rawMean = rawMean.rawMean;
    img_desc = 'DATA_kernel.mat';
    kernel = load( sprintf( '%s%s', img_folder, img_desc ) );
    data_xcorr.kernel = kernel.kernel;
    img_desc = 'DATA_ground_truth.mat';
    ground_truth = load( sprintf( '%s%s', img_folder, img_desc ) );
    data_xcorr.x_ground_truth = ground_truth.x_ground_truth;
    data_xcorr.y_measured = ground_truth.y_measured;
end

mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
mdata.acquired_S = data_xcorr;
setappdata(0,'data',mdata);     % 数据写入到matlabworkspace
end

function simu_data_Callback(hObject, eventdata)
mdata = getappdata(0,'data');
if ~isempty(mdata.acquired_M)
    M = mdata.acquired_M.M;
elseif ~isempty(mdata.simu)
    M = mdata.simu.M;
else
    [M, pathname] = uigetfile({'*.*'},'选择标定矩阵');
    M = load([pathname '\' M]);
    M = M.calibration_matrix;        
    M = double(M);  % Compute in double
    rem_last = 2;      % 用的是Heide的数据集，其它数据集需要根据具体情况进行修改
    M = M(:, 1:end - rem_last);  % TOF变换矩阵    
end
    

%% 打开视频文件,获取仿真的元数据

% 测试文件是 fruit_long1original.mat
[filename,pathname] = uigetfile({'*.*'},'选择用于仿真的视频文件或图像序列文件'); 
suffix = strfind(filename,'.');
suffix = filename(suffix(end):end);

switch lower(suffix)
    case {'.mp4','.avi','.mov'}                       % 输入是“.mp4”格式的视频文件 
        vidObj = VideoReader([pathname filename]);
        vidHeight = (vidObj.Height);
        vidWidth = (vidObj.Width);
        
        % 显示读入的视频            
        figure;
        currAxes = axes;
        while hasFrame(vidObj)
            vidFrame = readFrame(vidObj);
            image(vidFrame, 'Parent', currAxes);
            currAxes.Visible = 'off';
            pause(1/vidObj.FrameRate);
        end    
        close

        vid_data = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),'colormap',[]);

        % 读入特定时间段内的视频
        vidObj.CurrentTime = 00;   % 从第42s开始    25             default: 07
        k = 0;
        while vidObj.CurrentTime < vidObj.Duration    %17        % 读到第55秒     36 
            k = k+1;
            vid_data(k).cdata = readFrame(vidObj);                    
        end      

        parfor n=1:k
            true_data(:,:,n) = rgb2gray(vid_data(n).cdata);
        end
        
    case {'.mat'}     % 输入是“.mat”格式的视频文件     
        true_data = load([pathname filename]);
        true_data = true_data.true_data;

    case {'.png'}        % 输入是“.png”格式的图像文件,读取全部图像，然后连续选取199张，也可以像.mat文件那样随机选取；Heide用的仿真数据（443张图片） 
        img_files = dir(strcat(pathname,'*.png'));
        k = length(img_files);
        parfor n = 1:k
            true_data(:,:,n) = imread(strcat(pathname,img_files(n).name));
        end
    otherwise
        disp('输入数据格式不在考虑范围');
        return
end

%% 根据仿真数据构造标定矩阵
%{ 
if size( M,2 ) < size( true_data,3 )
    K = fix(size( true_data,3 ) / size( M,2 ));
    multi_M = M;
    for k = 1:K
        multi_M = [ multi_M M ];
    end
    multi_M = multi_M( :,1:size( true_data,3 ) );
elseif size( M,2 ) > size( true_data,3 )
    multi_M = M( :,1:size( true_data,3 ) );
end
simu_data = true_data;
M = multi_M;
%}

%% 构造仿真数据
 
if size(true_data,3) < size(M,2)    % 仿真数据不足199帧  动态视频 nature07980-s3.mov 有这个问题
    n = size(M,2) - size(true_data,3);
    k = zeros(size(true_data));
    k = k(:,:,1:n);
    tmp = zeros(size(true_data,1),size(true_data,2),size(M,2));
    m = ceil(n/2);
    tmp(:,:,1:m) = k(:,:,1:m);
    tmp(:,:,m+1:m+size(true_data,3)) = true_data;
    tmp(:,:,m+1+size(true_data,3):end) = k(:,:,m+1:end);
    simu_data = tmp;
else
    rand_select = 0;    % 1: 随机选取 
    switch rand_select
        case 0                   % 从第n张开始连续选取
            n = 96;            % 96：从Heide数据集选取的第一张图片
            simu_data = true_data(:,:,n:n+size(M,2)-1);
        case 1                   % 随机选取
            n = size(true_data,3);
            k = randperm(n);k = k(1:size(M,2));k = sort(k);    % 199是用ToF相机的采样的时间参数；k是随机选的，在实际中，“k随机选”意味着“拍照时，可能是不连续的”
            simu_data = true_data(:,:,k);                % max(simu_data(:)) = 255      max(beta_ideal(:)) = 1.4550e+03    由于，这两行的值供参考
        case 2             % 把时间轴缩小到1:199
            simu_data = reshape(true_data, [size(true_data,1)*size(true_data,2),size(true_data,3)]);
            simu_data = imresize(simu_data,[size(true_data,1)*size(true_data,2),size(M,2)]);
            simu_data = reshape(simu_data,[size(true_data,1),size(true_data,2),size(M,2)]);
            %{
            Show_Results(simu_data);
            sp = squeeze(simu_data(40,100,:));        % 双峰(96,86,:) (100,140,:) 低频：(120,140,:)  单峰：(100,40,:) (40,100,:)
            tp = squeeze(true_data(40,100,:));
            figure;plot(sp);figure;plot(tp);
            %}
    end
end


%{
 % 缩小尺寸，节省计算时间。
tmp = [];
parfor n = 1:size(simu_data,3)         
    tmp(:,:,n) = imresize(simu_data(:,:,n),[29,29]);   %  [40,55] [120,165]
end
simu_data = tmp;  
%}

mM.simu = double(simu_data);        % 时域信号
mM.M = M;

mM.imagedims(1) = size(simu_data,2);
mM.imagedims(2) = size(simu_data,1);
mM.num_frequencies = size(M,1)/2;
mM.num_phases = 2;;

% mM.filename = input('Filename is : ','s');
mM.filename = 'Tomato';


mdata.simu = mM;
setappdata(0,'data',mdata);     

%% debug   显示仿真输入
% tmp = mat2gray(simu_data);
% Show_Results(tmp);
disp('仿真数据读取结束。');
end

function append_path_Callback(hObject, eventdata)
% 添加需要的其它路径    
% Transient_Imaging_Path; 
end

function close_others_Callback(hObject, eventdata)
h = findobj('Type','Figure');
for n = 1:length(h)
    if isempty(h(n).Name)
        close (h(n)); % 除主界面和被调用界面外，其余窗口都被关闭。 length(h)>2
    end
end
end

function close_all_Callback(hObject, eventdata)
h = findobj('Type','Figure');  
for n = 1:length(h)
    if ~isempty(h(n).Name) && isequal(h(n).Name,'M Computing')
        break;
    end
end

% 除主界面外，关闭所有窗口。
if n == 1
    close(h(2:end)); 
else
    close(h(1:n-1));
    close(h(n+1:end));
end
end

function exit_Callback(hObject, eventdata)
h_main = findobj('Name','M Computing');
str1 = 'M Computing';
str2 = 'all';
selection = questdlg(['Quit ' str1 '?'],...
                     ['Close ' str2 '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'Yes')
    delete(h_main);
    close all;
    mdata = getappdata(0,'data');
    needpath = mdata.path;
    for k = 1:length(needpath)
        rmpath(needpath(k));             % 删除添加的路径；在运行他人的代码之前，要删除已添加的路径，以免函数重名，调用时出问题。
    end
end
end
%% Lab
function transient_imaging_Callback(hObject, eventdata)
Transient_Imaging;        % 调用瞬态成像窗体
end

function conventional_imaging_Callback(hObject, eventdata)
Conventional_Imaging;        % 调用传统成像窗体
end