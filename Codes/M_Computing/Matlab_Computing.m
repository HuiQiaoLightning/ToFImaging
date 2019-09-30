%{ 
MatlabComputing_GUI��
1. Select a data set from the pop-up menu, 
2. Select experiment from the pop-up menu, һ��ʵ�����ݶ�Ӧһ������

�������򣺺�����Сд����Ӧ�����ú���������������ĸ��д
n,k��frame�ļ�����n��ʾ������kΪ�ܵ�֡��
xs��ys�����ص�������������x��y�����ص��м������м���
%}
function Matlab_Computing
% ��ʼ��
init;

% ������
MComputing = figure('MenuBar','None','Name','M Computing','NumberTitle','off');

% File�˵���������
mfile = uimenu(MComputing,'Label','File');

uimenu(mfile,'Label','Open_img','Callback',@open_img_Callback);
uimenu(mfile,'Label','Acquired_M_ToF(��Ƶ���ź�)','Callback',@acquired_M_ToF_Callback, 'Separator','on');
uimenu(mfile,'Label','Acquired_S_ToF(��Ƶ���ź�)','Callback',@acquired_S_ToF_Callback);
uimenu(mfile,'Label','Simu_data(ʱ���ź�)','Callback',@simu_data_Callback);

uimenu(mfile,'Label','Append_path','Callback',@append_path_Callback, 'Separator','on');
uimenu(mfile,'Label','Close_others','Callback',@close_others_Callback, 'Separator','on');
uimenu(mfile,'Label','Close_all','Callback',@close_all_Callback);
uimenu(mfile,'Label','Quit','Callback',@exit_Callback,... 
    'Separator','on','Accelerator','Q');

%{
���������õĲ˵���    
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
% Create the UICONTEXTMENU  ���������Ĳ˵��������е��Ҽ��˵�
cmenu = uicontextmenu;

% Create the parent menu  ���˵�
closemenu = uimenu(cmenu,'label','Close');

% Create the submenus �Ӳ˵�
close_others = uimenu(closemenu,'label','Close_others',...
    'Callback',@close_others_Callback);
close_all = uimenu(closemenu,'label',...
    'Close_all','Callback',@close_all_Callback);
exit = uimenu(closemenu,'label','Quit',...
    'Callback',@exit_Callback);           
MComputing.UIContextMenu = cmenu;


% Lab�˵���������
mlab = uimenu(MComputing,'Label','Lab');
uimenu(mlab,'Label','Transient_imaging','Callback',@transient_imaging_Callback);
uimenu(mlab,'Label','Conventional_imaging','Callback',@conventional_imaging_Callback);
end

function init
%% ������ķ�ʽ���·��
needpath = genpath( './TransientImaging/Experiments' );
needpath = ["./TransientImaging"; needpath];
for n = 1:length(needpath)
    addpath (needpath(n));
end
% addpath ( './TransientImaging' );                            % ���˲̬����������ڵ�·��
% addpath (genpath( './TransientImaging/Experiments' ));       % ���˲̬����ʵ��������ڵ�·��������·��

%% �������ݽṹ
% mdataΪӦ�ó�������� ʹ�� getappdata �� setappdata ��������
% mdata.savefile = 'lastpath.mat'; 
mdata.filename = [];      % ͼ��(����)�ļ���
mdata.images = [];        % ͼ������

mdata.acquired_M = [];    % ͨ��GUI�Ѳɼ��Ķ�Ƶ�����ݣ���Heide�����ݼ������뵽matlabworkspace��ʵ��ʱ������ÿ�ζ���������
mdata.acquired_S = [];    % ͨ��GUI�Ѳɼ��ĵ�Ƶ�����ݣ���Kadambi�����ݼ������뵽matlabworkspace��ʵ��ʱ������ÿ�ζ���������
mdata.simu = [];        % ͨ��GUI�ѷ�������ݶ��뵽matlabworkspace��ʵ��ʱ������ÿ�ζ���������
mdata.selected = [];    % ���ڱ���ѡ��Ĳɼ����ݻ��������
mdata.path = needpath;  % �����˳�ʱɾ����Ӧ��·��

setappdata(0,'data',mdata);     % ���ݶ��뵽matlabworkspace

%% ��ͼ����
global lw fs                   % ���廭ͼʱ���߿������
lw = 4;          %  ��ͼ���߿�           3; 2; 1;          
fs = 19;         %  ��ͼ��ע�������С   19;14;11;  

%% ��ʱ�������λ��
global output_folder    
output_folder = '../TempData';

end

%% File
function open_img_Callback(hObject, eventdata)
mdata = getappdata(0,'data');
%%  �򿪴�ͳͼ��
[img_desc, img_folder] = uigetfile({'*.*'},'Open an image');
if ~img_desc == 0
    images = imread( sprintf( '%s%s', img_folder, img_desc ) );
    mdata.images = images;
    setappdata(0,'data',mdata);
else
    disp('ûѡ��ͼ����...');
end
end

function acquired_M_ToF_Callback(hObject, eventdata)
mM = Load_Transient_Data_With_UI;           % ��ȡ�ɼ�������
if ~isempty(mM)
    mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
    mdata.acquired_M = mM;
    setappdata(0,'data',mdata);     % ����д�뵽matlabworkspace
end
end

function acquired_S_ToF_Callback(hObject, eventdata)
[img_desc, img_folder] = uigetfile({'*.*'},'Load cross-correlation function from camera');
if img_desc == 0
    disp('ûѡ��ɼ��ĵ�Ƶ��˲̬��������');
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

mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
mdata.acquired_S = data_xcorr;
setappdata(0,'data',mdata);     % ����д�뵽matlabworkspace
end

function simu_data_Callback(hObject, eventdata)
mdata = getappdata(0,'data');
if ~isempty(mdata.acquired_M)
    M = mdata.acquired_M.M;
elseif ~isempty(mdata.simu)
    M = mdata.simu.M;
else
    [M, pathname] = uigetfile({'*.*'},'ѡ��궨����');
    M = load([pathname '\' M]);
    M = M.calibration_matrix;        
    M = double(M);  % Compute in double
    rem_last = 2;      % �õ���Heide�����ݼ����������ݼ���Ҫ���ݾ�����������޸�
    M = M(:, 1:end - rem_last);  % TOF�任����    
end
    

%% ����Ƶ�ļ�,��ȡ�����Ԫ����

% �����ļ��� fruit_long1original.mat
[filename,pathname] = uigetfile({'*.*'},'ѡ�����ڷ������Ƶ�ļ���ͼ�������ļ�'); 
suffix = strfind(filename,'.');
suffix = filename(suffix(end):end);

switch lower(suffix)
    case {'.mp4','.avi','.mov'}                       % �����ǡ�.mp4����ʽ����Ƶ�ļ� 
        vidObj = VideoReader([pathname filename]);
        vidHeight = (vidObj.Height);
        vidWidth = (vidObj.Width);
        
        % ��ʾ�������Ƶ            
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

        % �����ض�ʱ����ڵ���Ƶ
        vidObj.CurrentTime = 00;   % �ӵ�42s��ʼ    25             default: 07
        k = 0;
        while vidObj.CurrentTime < vidObj.Duration    %17        % ������55��     36 
            k = k+1;
            vid_data(k).cdata = readFrame(vidObj);                    
        end      

        parfor n=1:k
            true_data(:,:,n) = rgb2gray(vid_data(n).cdata);
        end
        
    case {'.mat'}     % �����ǡ�.mat����ʽ����Ƶ�ļ�     
        true_data = load([pathname filename]);
        true_data = true_data.true_data;

    case {'.png'}        % �����ǡ�.png����ʽ��ͼ���ļ�,��ȡȫ��ͼ��Ȼ������ѡȡ199�ţ�Ҳ������.mat�ļ��������ѡȡ��Heide�õķ������ݣ�443��ͼƬ�� 
        img_files = dir(strcat(pathname,'*.png'));
        k = length(img_files);
        parfor n = 1:k
            true_data(:,:,n) = imread(strcat(pathname,img_files(n).name));
        end
    otherwise
        disp('�������ݸ�ʽ���ڿ��Ƿ�Χ');
        return
end

%% ���ݷ������ݹ���궨����
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

%% �����������
 
if size(true_data,3) < size(M,2)    % �������ݲ���199֡  ��̬��Ƶ nature07980-s3.mov ���������
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
    rand_select = 0;    % 1: ���ѡȡ 
    switch rand_select
        case 0                   % �ӵ�n�ſ�ʼ����ѡȡ
            n = 96;            % 96����Heide���ݼ�ѡȡ�ĵ�һ��ͼƬ
            simu_data = true_data(:,:,n:n+size(M,2)-1);
        case 1                   % ���ѡȡ
            n = size(true_data,3);
            k = randperm(n);k = k(1:size(M,2));k = sort(k);    % 199����ToF����Ĳ�����ʱ�������k�����ѡ�ģ���ʵ���У���k���ѡ����ζ�š�����ʱ�������ǲ������ġ�
            simu_data = true_data(:,:,k);                % max(simu_data(:)) = 255      max(beta_ideal(:)) = 1.4550e+03    ���ڣ������е�ֵ���ο�
        case 2             % ��ʱ������С��1:199
            simu_data = reshape(true_data, [size(true_data,1)*size(true_data,2),size(true_data,3)]);
            simu_data = imresize(simu_data,[size(true_data,1)*size(true_data,2),size(M,2)]);
            simu_data = reshape(simu_data,[size(true_data,1),size(true_data,2),size(M,2)]);
            %{
            Show_Results(simu_data);
            sp = squeeze(simu_data(40,100,:));        % ˫��(96,86,:) (100,140,:) ��Ƶ��(120,140,:)  ���壺(100,40,:) (40,100,:)
            tp = squeeze(true_data(40,100,:));
            figure;plot(sp);figure;plot(tp);
            %}
    end
end


%{
 % ��С�ߴ磬��ʡ����ʱ�䡣
tmp = [];
parfor n = 1:size(simu_data,3)         
    tmp(:,:,n) = imresize(simu_data(:,:,n),[29,29]);   %  [40,55] [120,165]
end
simu_data = tmp;  
%}

mM.simu = double(simu_data);        % ʱ���ź�
mM.M = M;

mM.imagedims(1) = size(simu_data,2);
mM.imagedims(2) = size(simu_data,1);
mM.num_frequencies = size(M,1)/2;
mM.num_phases = 2;;

% mM.filename = input('Filename is : ','s');
mM.filename = 'Tomato';


mdata.simu = mM;
setappdata(0,'data',mdata);     

%% debug   ��ʾ��������
% tmp = mat2gray(simu_data);
% Show_Results(tmp);
disp('�������ݶ�ȡ������');
end

function append_path_Callback(hObject, eventdata)
% �����Ҫ������·��    
% Transient_Imaging_Path; 
end

function close_others_Callback(hObject, eventdata)
h = findobj('Type','Figure');
for n = 1:length(h)
    if isempty(h(n).Name)
        close (h(n)); % ��������ͱ����ý����⣬���ര�ڶ����رա� length(h)>2
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

% ���������⣬�ر����д��ڡ�
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
        rmpath(needpath(k));             % ɾ����ӵ�·�������������˵Ĵ���֮ǰ��Ҫɾ������ӵ�·�������⺯������������ʱ�����⡣
    end
end
end
%% Lab
function transient_imaging_Callback(hObject, eventdata)
Transient_Imaging;        % ����˲̬������
end

function conventional_imaging_Callback(hObject, eventdata)
Conventional_Imaging;        % ���ô�ͳ������
end