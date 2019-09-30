function Transient_Imaging
% ������
Tran_Imgs = figure('MenuBar','None','Name','Transient   Imaging','NumberTitle','off');

% Experiments�˵���������
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


% uimenu(mexperiments,'Label','Acquired_data(Ƶ���ź�)','Callback',@acquired_data_Callback,...
%     'Separator','on');


% Heide�˵���������
mheide = uimenu(Tran_Imgs,'Label','Heide');
uimenu(mheide,'Label','i_step','Callback',@i_step_Callback);
uimenu(mheide,'Label','u_step','Callback',@u_step_Callback);
uimenu(mheide,'Label','i_step_UI','Callback',@i_step_UI_Callback,...
    'Separator','on');
uimenu(mheide,'Label','i_step_UI_noise','Callback',@i_step_UI_noise_Callback);
uimenu(mheide,'Label','iu_step_UI','Callback',@iu_step_UI_Callback);


% Lin�˵���������
mlin = uimenu(Tran_Imgs,'Label','Lin');
uimenu(mlin,'Label','v03','Callback',@lin_v03_Callback);
uimenu(mlin,'Label','v04','Callback',@lin_v04_Callback);
uimenu(mlin,'Label','v05','Callback',@lin_v05_Callback);
uimenu(mlin,'Label','v03_v04_v05_UI','Callback',@lin_v03_v04_v05_UI_Callback,...
    'Separator','on');

% Kadambi�˵���������
mkadambi = uimenu(Tran_Imgs,'Label','Kadambi');
uimenu(mkadambi,'Label','xcorr','Callback',@kadambi_xcorr_Callback);
uimenu(mkadambi,'Label','cs_xcorr','Callback',@exp_cs_xcorr_Callback);

% Papers�˵���������
mpapers = uimenu(Tran_Imgs,'Label','Papers');
uimenu(mpapers,'Label','Model_Study_TI','Callback',@papers_Model_Study_TI_Callback);

end

%% Heide���㷨
function i_step_Callback(hObject, eventdata)
oldFolder = pwd;                    % ��ȡ��ǰ·��
cd ./TransientImaging/Felix_Heide/code_and_data/code               % ʹ�����·��
reconstruct_launch;
cd(oldFolder)                % �ص���ǰ·��
end

function u_step_Callback(hObject, eventdata)
oldFolder = pwd;
cd ./TransientImaging/Felix_Heide/code_and_data/code/gaussian_exponential_model                % ʹ�����·��
reconstruct_cluster_launch(1);     % <=128   ������Ȼ�����У����õ�����ЩU������ʲô�أ�Ϊʲô��������ʽ�����Ҽ����ٶ�Ҳ�ܿ죬����˵�е���������
cd(oldFolder)
end

function i_step_UI_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
if ~isempty(mdata.selected) 
    mM = mdata.selected;
    i_step_UI(mM);     
else
    disp('û��ѡ�����ݼ�...');
end
end

function iu_step_UI_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
if ~isempty(mdata.selected) 
    mM = mdata.selected;
    iu_step_UI(mM);     
else
    disp('û��ѡ�����ݼ�...');
end
end

function i_step_UI_noise_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
mM = mdata.selected;
if ~isempty(mdata.acquired)
    mM = mdata.acquired;
    i_step_UI_noise(mM);
else
    disp('û�вɼ�����...');
end
end
%% Lin���㷨
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
disp('start v05');               % ��� by whm
% ������4����ѡ��һ�����ԣ�ǰ2���ǹս��ؽ���ͼ�񣻺�2���������ģ 
% rec1_mirror;   
% rec2_mirror;
show_calib_amp;
show_calib_raw;
disp('end v05');               % ��� by whm
cd(oldFolder)
end

function lin_v03_v04_v05_UI_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
if ~isempty(mdata.selected)
    mM = mdata.selected;
    Lin_v03_v04_v05_UI(mM);
else
    disp('û�вɼ�����...');
end
end

%% Experiments
function select_data_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace

if isempty(mdata.simu) && (~isempty(mdata.acquired_M))
    result = 2;
elseif (~isempty(mdata.simu)) && isempty(mdata.acquired_M)
    result = 1;
else
    result = input('ѡ��ʵ��������Դ...   \n 1:simu; \n 2:acquired \n');
end

switch result
    case 1
        if isempty(mdata.simu)
            display('���������������...');
            return
        end
        mM = mdata.simu;
        simu = mM.simu;
        imagedims = size(simu);
        simu = reshape(simu, [imagedims(1)*imagedims(2),imagedims(3)]); 
        M = mM.M;
%         noise_type = 0; %2;
        prompt = 'ѡ����������������������: ';
        display('0.  �궨���󲻼�����');
        display('1.  �궨����ӹ̶���Poisson���������ؽ��㷨��Ҳʹ�ø�����');
        display('2.  �궨����Ӳ�������');   
        display('3.  �궨����Ӿ�������');   
        display('6.  �궨����Ӿ�������������Heide�����ݣ��õ��ļ������Ϻ�ʱ���õľ�����������');
        display('4.  �궨����� ��������+����������ûʵ��');    
        display('5.  �궨����Ӹ�˹����1');
        display('7.  �궨����Ӹ�˹����2');
        noise_type = input(prompt);
        true_measurements = M * simu'; 
        switch noise_type
            case 0             % ��������   
                noise_measurements = 0;
                display('�������ݲ�������');
            case 2              % Poisson����
                lambda = 2;
                noise = poissrnd(lambda,size(M));                 
%                 save( sprintf('%s/poisson_noise.mat', './Outputs'), 'noise' );               
                noise_measurements = (M + noise*0.02) * simu';    % 0.05   
                display('�������ݼ������������');
            case 1            % �ӹ̶���Poisson������ʵ��16�У��ؽ�ʱҲʹ�ø�����
                noise = load('poisson_noise.mat');  % �������õĲ������ݣ�������ɵģ������������ڡ�E:\HM_org\papers\data����ͬ�����ݡ�
                noise = noise.noise;                     
                noise_measurements = (M + noise*0.05) * simu';     % M + noise*0.05
                display('�������ݼӹ̶���������');
            case 3          % ��������
                noise = rand(size(M)); 
                noise_measurements = (M + noise*0.05) * simu';   % 0.05
                display('�������ݼ������������');
            case 6
                noise = load('./TransientImaging/Experiments/noise_Heide.mat');               % ����Heide�����ݣ��õ��ļ������Ϻ�ʱ���õľ�����������
                noise = noise.noise;
                noise_measurements = (M + noise*0.05) * simu';
                display('�������ݼӹ̶���������');
            case 4           % Poisson���� + ��������    �о�Ӧ������������
            case 5          % ��˹���� 1  ��Ҫ�õĻ��� randn
                noise = normrnd(0,1,size(M));
                noise_measurements = (M + noise*0.05) * simu';
                display('�������ݼ������˹����');
            case 7          % ��˹���� 2
                noise = randn(size(M)); 
                noise_measurements = (M + noise*0.03) * simu';      % 0.05           
                display('�������ݼ������˹����');
        end        
        true_measurements = reshape(true_measurements',[imagedims(1),imagedims(2),size(M,1)]);
        if noise_type ~= 0
            noise_measurements = reshape(noise_measurements',[imagedims(1),imagedims(2),size(M,1)]);
        end
        mM.true_measurements = true_measurements;
        mM.noise_measurements = noise_measurements;
        %{ 
        % ��ʾ
        [y,x] = deal(120,140);             % ������������ص㣺(60,70);(30,40);(149,88);(83,136);(91,150);(131,67);(86,31);(95,148);(57,92);(45,53); (75,164)  
        true_m = squeeze(true_measurements(x,y,:));
        noise_m = squeeze(noise_measurements(x,y,:));     
        % ����������ˮƽ�߲���һ�𣬺���Ƚϵ�ʱ�򣬶�ÿ������Ҫȥ����ֵ����ͬһ��ˮƽ���ϱȽϡ�����������ֵȥ��ֵ��ÿһ�����صĲ���ֵò��Ҳ����ͬһ��ˮƽ���ϡ��ڱȽ�ʱ������Ҫ����
        
        figure;plot( [ true_m noise_m ] );
        legend('True','Noise'); title([mM.filename ' measurements, at (', num2str(y) ',' num2str(x),')']);     
        
        tmp = mat2gray(true_measurements);
        Show_Results(tmp);
        %}
    case 2
        if isempty(mdata.acquired_M)
            display('��������۲�����...');
            return
        end        
        mM = mdata.acquired_M;
        display(['��ʹ��Heide�ɼ������ݣ�' mM.filename]);
end

mdata.selected = mM;
setappdata(0,'data',mdata);     % ����д�뵽matlabworkspace
end

function principle_simulation_analy_Callback(hObject, eventdata)
% ʹ�÷������ݽ���˲̬�����ԭ�����
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
if ~isempty(mdata.simu)
    Principle_Simulation_Analy(mdata.simu);
else
    disp('û�м��ط�������...');
end    
end

function acquired_data_analy_Callback(hObject, eventdata)
% �ɼ����ݵķ���
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
if ~isempty(mdata.acquired)
    Acquired_Data_Analy(mdata.acquired);
else
    disp('û�м��زɼ�����...');
end    
end

function experiment_1_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_1(mdata.selected);
end

function experiment_2_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_2(mdata.selected);
end

function experiment_3_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_3(mdata.selected);
end

function experiment_4_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_4(mdata.selected);
end

function experiment_5_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_5(mdata.selected);
end

function experiment_6_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_6_ver1(mdata.selected);
end

function experiment_7_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_7(mdata.selected);
end

function experiment_8_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_8(mdata.selected);
end

function experiment_9_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
Experiment_9(mdata.selected);
end

function experiment_10_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
% dbstop in Experiment_10 at 12  % �к�
% dbstop in Experiment_10 at Adaptive_Transient_Signal_n
Experiment_10(mdata.selected);
end

function experiment_11_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
% dbstop in Experiment_10 at 12  % �к�
% dbstop in Experiment_11 at Adaptive_Transient_Signal_n
Experiment_11(mdata.selected);
end

function experiment_12_Callback(hObject, eventdata)
mdata = getappdata(0,'data');   % ���ݶ����matlabworkspace
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
Project_3(mdata.selected);  %  Project_3_0_1(mdata.selected); ͬProject_3_0(mdata.selected);ֻ���Ȱ�beta׼���ã�ռ�ڴ�
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