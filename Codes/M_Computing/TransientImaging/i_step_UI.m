function [I_rec, mdist] = i_step_UI(mM)
verbose = 'none'; %Options: 'all', 'brief', 'none'
lambda_residual = 1;
lambda_temporal = 5.0;
lambda_spatial = 0.3;
thresh_huber = 0.05;

if isfield(mM,'rec_measurements')
    measurements = mM.rec_measurements;
    M = mM.rec_M;    
elseif isfield(mM, 'simu')
    if mM.noise_measurements == 0
        measurements = mM.true_measurements;
    else
        measurements = mM.noise_measurements;   % 加噪声的数据
    end
    M = mM.M;
else
    measurements = mM.measurements;  % 采集的数据
    M = mM.M;
end

needpath = genpath('./TransientImaging/Felix_Heide/code_and_data/code'); 
addpath (needpath);         % 添加路径

tic;
%% Reconstruct using primal-dual solver
I_reconstructed = pd_solve_full_svd_coherence(measurements, M, lambda_residual, lambda_temporal, lambda_spatial, thresh_huber, 50, 1e-5, verbose);
telapsed = toc;
disp(['  First subproblem took ' num2str(telapsed) ' secs']);

rmpath(needpath);             % 删除添加的路径；在运行他人的代码之前，要删除已添加的路径，以免函数重名，调用时出问题。   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show the result and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Show transient image
scaling = 1.0; %Amplitude scale
I_rec = max(I_reconstructed,0);                    % I_rec 这个是真正重构的信号；后续的都是经过处理用于图像显示的

% mdist = TI_Postprocess(mM,I_rec,'_Heide',telapsed);   %这个函数实现的是下面代码的功能，但删掉了个别语句
mdist = 0;

global output_folder 
if isfield(mM,'rec_measurements')
    aux = '_CS';
elseif isfield(mM,'simu')
    aux = '_Poisson';
elseif 442 == size(mM.M,1)
    aux = '';
end
algo = '_Heide';

if isfield(mM,'simu')
    save( sprintf(['%s/' mM.filename algo aux '_QH.mat'], output_folder), 'I_rec','telapsed');         % 仿真数据
else
    save( sprintf(['%s/' mM.filename algo aux '_QH.mat'], output_folder), 'I_rec','telapsed');
end

%{
% 计算观测值的距离
if isfield(mM,'simu')
    PX2 = (reshape(I_rec,[mM.imagedims(1)*mM.imagedims(2),size(M,2)]))';
    simu = (reshape(mM.simu,[mM.imagedims(1)*mM.imagedims(2),size(M,2)]))';
    mdist = norm(simu - PX2)/norm(simu);
    %{
    simu = mM.simu;
    [y,x] = deal(100,40);             % 分析具体的像素点：低频：(136,83);(149,88);(150,91);(148,95);(140,120);单峰：(70,60);(31,86);(92,57);单峰2个return：(40,30);(53,45);双峰：(131,67);
    rp = squeeze(I_rec(x,y,:));
    sp = squeeze(simu(x,y,:));
    figure;plot([sp rp]);legend('True','Heide');title([mM.filename ', (', num2str(y) ',' num2str(x),')']);    
    xlabel('Time Series');ylabel( 'Amplitude' );    

    figure;plot(time, [sp rp],'LineWidth',lw);    % 回信中的图, time等的定义参见 function papers_Model_Study_TI(mM)
    title(['Signal at (', num2str(y) ',' num2str(x),')']);
    xlim([0 70]);ylabel( 'Amplitude' );xlabel('\tau (ns)');
    legend('True','Heide''s i-step');
    set(gca,'FontSize',fs);
    %}
else
    PX2 = (reshape(I_rec,[mM.imagedims(1)*mM.imagedims(2),size(M,2)]))';
    rec_m = M * PX2;
    true_m = (reshape(measurements,[mM.imagedims(1)*mM.imagedims(2),size(M,1)]))';
    mdist = norm(true_m - rec_m)/norm(true_m);
    %{
    rec_m = reshape(rec_m',[mM.imagedims(2),mM.imagedims(1),size(M,1)]);
    true_m = mM.measurements;
  
    [y,x] = deal(150,91);             % 分析具体的像素点：(60,70);(30,40);(149,88);(83,136);(150,91);(131,67);(86,31);(95,148);(57,92);(45,53); (75,164)
    rp = squeeze(I_rec(x,y,:));
    figure;plot(rp);title([mM.filename ', signal at (', num2str(y) ',' num2str(x),')']);
    
    rm = squeeze(rec_m(x,y,:));
    tm = squeeze(true_m(x,y,:));
    figure;plot([tm rm]);
    legend('True','Rec');
    title([mM.filename ', measurements at (', num2str(y) ',' num2str(x),')']);    
    %}
end
%
PX2 = I_rec./max(I_rec(:));
PX2 = min(scaling * PX2,1);
if isfield(mM,'simu')
    simu = mat2gray(mM.simu);
elseif isequal(mM.filename, 'DiscoBall') 
    PX2 = (PX2 * 1.0) .^ (1/5.2); %Gamma
else
    PX2 = (PX2 * 1.0) .^ (1/2.2); %Gamma
end

%Show transient image
figure;
if isfield(mM,'simu')
    for t = 1:size(PX2,3)    
        subplot(1,2,1);imshow(PX2(:,:,t));title(sprintf('Transient image frame %d', t));
        subplot(1,2,2);imshow(simu(:,:,t)); title(sprintf('True image frame %d', t));       
        pause(0.1)
    end
else
    scalefactor = 1;
    for t = 1:size(PX2,3)  
        % 放大后显示
        I_show = imresize( PX2(:,:,t), scalefactor,'bilinear');
        imshow( I_show )
        title(sprintf('Transient image frame %d', t));
        pause(0.1)
    end  
end


% Save the result  存成.mat格式

global output_folder 

aux = '_CS';
aux = '';
aux = '_poisson';
aux = '_poisson_CS';
if isfield(mM,'simu')
    save( sprintf(['%s/' mM.filename '_Heide' aux '.mat'], output_folder), 'I_rec','mdist','telapsed');         % 仿真数据
else
    save( sprintf(['%s/' mM.filename '_Heide' aux '.mat'], output_folder), 'I_rec','mdist','telapsed','true_m','rec_m');
end

% 以前保存的内容
% save( sprintf('%s/reconstruction_I.mat', output_folder), 'I_reconstructed','I_rec');


Show_Results(PX2,[mM.filename '_Heide' aux],'mp4');    % 生成视频文件


% 存成图片的格式
%{
currpath = pwd;     % 获取当前路径
global output_folder 
cd(output_folder );

for t = 1:size(PX2,3)  
    I_show = PX2(:,:,t);
    imwrite(I_show,[mM.filename '_Heide' aux '_' num2str(t) '.bmp']);      
end  
cd(currpath);       % 返回程序运行的当前路径
%}
disp('Done the reconstruction for current job.');
%}
  
end

%% 用Heide的tomato计算结果保存和显示
% save('Heide_tomato.mat', 'I_reconstructed','mM', 'I_rec' );
% 
% I_reconstructed = mat2gray(I_reconstructed);
% [n,m] = deal(96,86);             % 分析具体的像素点：(60,70);(30,40);(149,88);(136,83);(150,91);(131,67);(86,31);(95,148);(57,92);(45,53);
% rp = I_reconstructed(n,m,:);rp = rp(:);
% sp = simu(n,m,:);sp = sp(:);
% figure;plot([rp sp]);legend('Heide - reconstructed','original');title(['(' num2str(n) ',' num2str(m) ')']);  
% xlabel('Time Series');ylabel( 'Amplitude' );