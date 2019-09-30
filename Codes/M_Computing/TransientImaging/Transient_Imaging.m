function Transient_Imaging
% 主窗体
Tran_Imgs = figure('MenuBar','None','Name','Transient   Imaging','NumberTitle','off');

% Experiments菜单及其子项
mexperiments = uimenu(Tran_Imgs,'Label','Experiments');

uimenu(mexperiments,'Label','Select_data','Callback',@select_data_Callback);
uimenu(mexperiments,'Label','Principle_simulation_analy','Callback',@principle_simulation_analy_Callback,...
    'Separator','on');
uimenu(mexperiments,'Label','Acquired_data_analy','Callback',@acquired_data_analy_Callback);

sub_experiments = uimenu(mexperiments,'Label','Experiment','Separator','on');
uimenu(sub_experiments,'Label','Experiment_1','Callback',@experiment_1_Callback);
uimenu(sub_experiments,'Label','Experiment_2','Callback',@experiment_2_Callback);
uimenu(sub_experiments,'Label','Experiment_3','Callback',@experiment_3_Callback);
uimenu(sub_experiments,'Label','Experiment_4','Callback',@experiment_4_Callback);
uimenu(sub_experiments,'Label','Experiment_5','Callback',@experiment_5_Callback);
uimenu(sub_experiments,'Label','Experiment_6','Callback',@experiment_6_Callback);
uimenu(sub_experiments,'Label','Experiment_7','Callback',@experiment_7_Callback);
uimenu(sub_experiments,'Label','Experiment_8','Callback',@experiment_8_Callback,'Separator','on');
uimenu(sub_experiments,'Label','Experiment_9','Callback',@experiment_9_Callback);
uimenu(sub_experiments,'Label','Experiment_10','Callback',@experiment_10_Callback);
uimenu(sub_experiments,'Label','Experiment_11','Callback',@experiment_11_Callback);
uimenu(sub_experiments,'Label','Experiment_12','Callback',@experiment_12_Callback);

sub_projects = uimenu(mexperiments,'Label','Project','Separator','on');
uimenu(sub_projects,'Label','Project_1','Callback',@project_1_Callback);
uimenu(sub_projects,'Label','Project_1_simu','Callback',@project_1_simu_Callback);
uimenu(sub_projects,'Label','Project_2','Callback',@project_2_Callback,'Separator','on');
uimenu(sub_projects,'Label','Project_3','Callback',@project_3_Callback);
uimenu(sub_projects,'Label','Project_4','Callback',@project_4_Callback);
uimenu(sub_projects,'Label','Project_5','Callback',@project_5_Callback);
uimenu(sub_projects,'Label','Project_6','Callback',@project_6_Callback,'Separator','on');
uimenu(sub_projects,'Label','Project_7','Callback',@project_7_Callback);


% uimenu(mexperiments,'Label','Acquired_data(频域信号)','Callback',@acquired_data_Callback,...
%     'Separator','on');


% Heide菜单及其子项
mheide = uimenu(Tran_Imgs,'Label','Heide');
uimenu(mheide,'Label','i_step','Callback',@i_step_Callback);
uimenu(mheide,'Label','u_step','Callback',@u_step_Callback);
uimenu(mheide,'Label','i_step_UI','Callback',@i_step_UI_Callback,...
    'Separator','on');
uimenu(mheide,'Label','i_step_UI_noise','Callback',@i_step_UI_noise_Callback);
uimenu(mheide,'Label','iu_step_UI','Callback',@iu_step_UI_Callback);


% Lin菜单及其子项
mlin = uimenu(Tran_Imgs,'Label','Lin');
uimenu(mlin,'Label','v03','Callback',@lin_v03_Callback);
uimenu(mlin,'Label','v04','Callback',@lin_v04_Callback);
uimenu(mlin,'Label','v05','Callback',@lin_v05_Callback);
uimenu(mlin,'Label','v03_v04_v05_UI','Callback',@lin_v03_v04_v05_UI_Callback,...
    'Separator','on');

% Kadambi菜单及其子项
mkadambi = uimenu(Tran_Imgs,'Label','Kadambi');
uimenu(mkadambi,'Label','xcorr','Callback',@kadambi_xcorr_Callback);
uimenu(mkadambi,'Label','cs_xcorr','Callback',@exp_cs_xcorr_Callback);

% Papers菜单及其子项
mpapers = uimenu(Tran_Imgs,'Label','Papers');
uimenu(mpapers,'Label','Model_Study_TI','Callback',@papers_Model_Study_TI_Callback);

end

%% Heide的算法
function i_step_Callback(hObject, eventdata)
oldFolder = pwd;                    % 获取当前路径
cd ./TransientImaging/Felix_Heide/code_and_data/code               % 使用相对路径
reconstruct_launch;
cd(oldFolder)                % 回到当前路径
end

function u_step_Callback(hObject, eventdata)
oldFolder = pwd;
cd ./TransientImaging/Felix_Heide/code_and_data/code/gaussian_exponential_model                % 使用相对路径
reconstruct_cluster_launch(1);     % <=128   代码虽然能运行，但得到的那些U究竟是什么呢？为什么是数组形式？而且计算速度也很快，不像传说中的那样慢？
cd(oldFolder)
end

function i_step_UI_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
if ~isempty(mdata.selected) 
    mM = mdata.selected;
    i_step_UI(mM);     
else
    disp('没有选择数据集...');
end
end

function iu_step_UI_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
if ~isempty(mdata.selected) 
    mM = mdata.selected;
    iu_step_UI(mM);     
else
    disp('没有选择数据集...');
end
end

function i_step_UI_noise_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
mM = mdata.selected;
if ~isempty(mdata.acquired)
    mM = mdata.acquired;
    i_step_UI_noise(mM);
else
    disp('没有采集数据...');
end
end
%% Lin的算法
function lin_v03_Callback(hObject, eventdata)
oldFolder = pwd;
cd ./TransientImaging/Lin_Jingyu/Lin_code/transient_pub3_v03/transient_ubc
rec_ubcdata;
cd(oldFolder)
end

function lin_v04_Callback(hObject, eventdata)
oldFolder = pwd;
cd ./TransientImaging/Lin_Jingyu/Lin_code/transient_code_pub4_v04/cmpu_transient_ubc
rec_ubcdata;
cd(oldFolder)
end

function lin_v05_Callback(hObject, eventdata)
oldFolder = pwd;
cd ./TransientImaging/Lin_Jingyu/Lin_code/sup_harmony_v05/sup_harmony/code
disp('start v05');               % 添加 by whm
% 从下面4个中选择一个测试，前2个是拐角重建的图像；后2个是相机建模 
% rec1_mirror;   
% rec2_mirror;
show_calib_amp;
show_calib_raw;
disp('end v05');               % 添加 by whm
cd(oldFolder)
end

function lin_v03_v04_v05_UI_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
if ~isempty(mdata.selected)
    mM = mdata.selected;
    Lin_v03_v04_v05_UI(mM);
else
    disp('没有采集数据...');
end
end

%% Experiments
function select_data_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace

if isempty(mdata.simu) && (~isempty(mdata.acquired_M))
    result = 2;
elseif (~isempty(mdata.simu)) && isempty(mdata.acquired_M)
    result = 1;
else
    result = input('选择实验数据来源...   \n 1:simu; \n 2:acquired \n');
end

switch result
    case 1
        if isempty(mdata.simu)
            display('请先输入仿真数据...');
            return
        end
        mM = mdata.simu;
        simu = mM.simu;
        imagedims = size(simu);
        simu = reshape(simu, [imagedims(1)*imagedims(2),imagedims(3)]); 
        M = mM.M;
%         noise_type = 0; %2;
        prompt = '选择仿真数据添加噪声的类型: ';
        display('0.  标定矩阵不加噪声');
        display('1.  标定矩阵加固定的Poisson噪声，在重建算法中也使用该噪声');
        display('2.  标定矩阵加泊松噪声');   
        display('3.  标定矩阵加均匀噪声');   
        display('6.  标定矩阵加均匀噪声，根据Heide的数据，得到的计算结果较好时采用的均匀噪声矩阵');
        display('4.  标定矩阵加 泊松噪声+均匀噪声，没实现');    
        display('5.  标定矩阵加高斯噪声1');
        display('7.  标定矩阵加高斯噪声2');
        noise_type = input(prompt);
        true_measurements = M * simu'; 
        switch noise_type
            case 0             % 不加噪声   
                noise_measurements = 0;
                display('仿真数据不加噪声');
            case 2              % Poisson噪声
                lambda = 2;
                noise = poissrnd(lambda,size(M));                 
%                 save( sprintf('%s/poisson_noise.mat', './Outputs'), 'noise' );               
                noise_measurements = (M + noise*0.02) * simu';    % 0.05   
                display('仿真数据加随机泊松噪声');
            case 1            % 加固定的Poisson噪声，实验16中，重建时也使用该噪声
                noise = load('poisson_noise.mat');  % 论文中用的测试数据，随机生成的，保留下来，在“E:\HM_org\papers\data”中同名备份。
                noise = noise.noise;                     
                noise_measurements = (M + noise*0.05) * simu';     % M + noise*0.05
                display('仿真数据加固定泊松噪声');
            case 3          % 均匀噪声
                noise = rand(size(M)); 
                noise_measurements = (M + noise*0.05) * simu';   % 0.05
                display('仿真数据加随机均匀噪声');
            case 6
                noise = load('./TransientImaging/Experiments/noise_Heide.mat');               % 根据Heide的数据，得到的计算结果较好时采用的均匀噪声矩阵。
                noise = noise.noise;
                noise_measurements = (M + noise*0.05) * simu';
                display('仿真数据加固定均匀噪声');
            case 4           % Poisson噪声 + 均匀噪声    感觉应该是这种噪声
            case 5          % 高斯噪声 1  主要用的还是 randn
                noise = normrnd(0,1,size(M));
                noise_measurements = (M + noise*0.05) * simu';
                display('仿真数据加随机高斯噪声');
            case 7          % 高斯噪声 2
                noise = randn(size(M)); 
                noise_measurements = (M + noise*0.03) * simu';      % 0.05           
                display('仿真数据加随机高斯噪声');
        end        
        true_measurements = reshape(true_measurements',[imagedims(1),imagedims(2),size(M,1)]);
        if noise_type ~= 0
            noise_measurements = reshape(noise_measurements',[imagedims(1),imagedims(2),size(M,1)]);
        end
        mM.true_measurements = true_measurements;
        mM.noise_measurements = noise_measurements;
        %{ 
        % 显示
        [y,x] = deal(120,140);             % 分析具体的像素点：(60,70);(30,40);(149,88);(83,136);(91,150);(131,67);(86,31);(95,148);(57,92);(45,53); (75,164)  
        true_m = squeeze(true_measurements(x,y,:));
        noise_m = squeeze(noise_measurements(x,y,:));     
        % 它们俩个的水平线不在一起，后面比较的时候，对每个像素要去掉均值，在同一个水平线上比较。但整个测量值去均值，每一个像素的测量值貌似也不在同一个水平线上。在比较时，还需要后处理。
        
        figure;plot( [ true_m noise_m ] );
        legend('True','Noise'); title([mM.filename ' measurements, at (', num2str(y) ',' num2str(x),')']);     
        
        tmp = mat2gray(true_measurements);
        Show_Results(tmp);
        %}
    case 2
        if isempty(mdata.acquired_M)
            display('请先输入观测数据...');
            return
        end        
        mM = mdata.acquired_M;
        display(['现使用Heide采集的数据：' mM.filename]);
end

mdata.selected = mM;
setappdata(0,'data',mdata);     % 数据写入到matlabworkspace
end

function principle_simulation_analy_Callback(hObject, eventdata)
% 使用仿真数据进行瞬态成像的原理分析
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
if ~isempty(mdata.simu)
    Principle_Simulation_Analy(mdata.simu);
else
    disp('没有加载仿真数据...');
end    
end

function acquired_data_analy_Callback(hObject, eventdata)
% 采集数据的分析
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
if ~isempty(mdata.acquired)
    Acquired_Data_Analy(mdata.acquired);
else
    disp('没有加载采集数据...');
end    
end

function experiment_1_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_1(mdata.selected);
end

function experiment_2_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_2(mdata.selected);
end

function experiment_3_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_3(mdata.selected);
end

function experiment_4_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_4(mdata.selected);
end

function experiment_5_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_5(mdata.selected);
end

function experiment_6_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_6_ver1(mdata.selected);
end

function experiment_7_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_7(mdata.selected);
end

function experiment_8_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_8(mdata.selected);
end

function experiment_9_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
Experiment_9(mdata.selected);
end

function experiment_10_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
% dbstop in Experiment_10 at 12  % 行号
% dbstop in Experiment_10 at Adaptive_Transient_Signal_n
Experiment_10(mdata.selected);
end

function experiment_11_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
% dbstop in Experiment_10 at 12  % 行号
% dbstop in Experiment_11 at Adaptive_Transient_Signal_n
Experiment_11(mdata.selected);
end

function experiment_12_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % 数据读入从matlabworkspace
% Experiment_12(mdata.selected);
Experiment_13(mdata.selected);
end

function kadambi_xcorr_Callback(hObject, eventdata)
oldFolder = pwd;
cd ./TransientImaging/Nanophotography_Demo_Code/demo_deconvolution
disp('start deconvolution');
demo_deconvolution;
disp('end deconvolution')
cd( oldFolder )
end

function exp_cs_xcorr_Callback(hObject, eventdata)
mdata = getappdata(0,'data');
exp_cs_xcorr(mdata.acquired_S);
end

function project_1_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_1(mdata.selected);
% Project_1_1(mdata.selected);
end

function project_2_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_2(mdata.selected);
end

function project_3_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_3(mdata.selected);  %  Project_3_0_1(mdata.selected); 同Project_3_0(mdata.selected);只是先把beta准备好，占内存
% Project_3_1(mdata.selected);
end

function project_6_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_6(mdata.selected);  
end

function project_5_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_5(mdata.selected);  
end

function project_7_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_7(mdata.selected);  
end

%{
function project_1_simu_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_1_simu(mdata.selected);
end

function project_2_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
% Project_2_0(mdata.selected);
Project_2_1(mdata.selected);
end

function project_4_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
Project_4(mdata.selected);
end

function papers_Model_Study_TI_Callback(hobject, eventdata)
mdata = getappdata(0,'data');
papers_Model_Study_TI(mdata.selected);
end
%}