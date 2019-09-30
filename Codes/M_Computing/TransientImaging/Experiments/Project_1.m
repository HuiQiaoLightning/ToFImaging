function Project_1(mM)
% 第一个paper中的Heide采集数据的实验，

%% Compressed sampling and recovery. Testing datasets are Heide's datasets.

[mM, err] = Compressed_Sampling_Recovery(mM);    % Algo1  压缩采集

i_step_UI(mM);        % Algo1 + Heide_i_step   这个放在这里是用来跑压缩测量的重构误差的。  可以测试Tomato。

Lin_v03_v04_v05_UI(mM);                  % Algo1 + Lin 这个放在这里是用来跑压缩测量的重构误差的。 可以测试Tomato。    

Project_6(mM);                           % Algo1 + Algo2.1

%{ 
Wavelet_based transient image reconstruction  注意：Heide的算法是归一化后，调整Gamma的。  
采用经验噪声矩阵，仅取第3层概貌信号最大的那个峰值重构瞬态信号。
Algo2.2 在指定尺度下（scale = 3）重构瞬态图像，Heide的数据集，噪声矩阵是基于最初做的时候的经验确定的。不能测试Tomato。   
%}

% [I_recW,mdistW] = Pseudoinverse_Wavelet_based_TI_rec(mM);     

return

%%
function [mM, err] = Compressed_Sampling_Recovery(mM)
%% 数据准备；
if isfield(mM, 'simu')
    if mM.noise_measurements == 0           % 如果是采样仿真数据，加的噪声如果太大，失真严重
        measurements = mM.true_measurements;
    else
        measurements = mM.noise_measurements;   % 加噪声的数据
    end
else
    measurements = mM.measurements;  % 采集的数据
end

xs = mM.imagedims(1);
ys = mM.imagedims(2);
    
Nf = size(mM.M,1)/2;Nt = size(mM.M,2);   % Nf：频率采样点;Nt：时域采样点。不同的实验环境，可能不一样，这是针对Heide的数据集的。
num_frequencies = Nt;

measurements(:,:,Nf+Nt+1:end) = [];
measurements(:,:,Nt+1:Nf) = [];

measurements = (reshape(measurements,[xs*ys,num_frequencies*2]))'; % 数据准备；

%% 获得压缩采样的参数
% sparsity basis
Psi{1} = mM.M(1:Nt,:);
Psi{2} = mM.M(Nf+1:Nf+Nt,:);


%% 压缩采样        重建方法需要修改，直接用这些采样点重构，而不是先重构输入信号，这个要怎么做？
% (1)
Phi = 2;            % 前3个是成熟的方法
switch Phi            
    case 1
        Phi = 'TranL20';      % 准备Phi，使用一个定制的矩阵，意味着只采样20个不同频率的值  采样点低频多
    case 2
        Phi = 'Random';       % 使用55个随机采样点，根据理论推导得到的
end

tstart = tic;
% (2)  获得压缩采样的参数
sparse_K = 28;  % input('稀疏度 sparse_k ：  ');     
[Phi, fre_pos, sampls_M] = Draw_L_fre(num_frequencies, Phi, sparse_K);  

%m = measurements(1,:);
%[hat_m(1,:),Phi,sparse_K,sampls_M,fre_pos] = Compressed_Sampling(m(1,1:num_frequencies),Phi);   % 把xpos保存下来，在做对比分析时，画stem图比较方便。
   
% 重构时迭代次数 
fre_sparse = 1;    % 为0时，按时域稀疏计算,这里按频域稀疏考虑
if sampls_M > 30
    K_loop = 9;        %  算法迭代次数(m>=K)，设x是K-sparse的,有时候，m太大，反而不收敛
else
    K_loop = 2*sparse_K;
end

[hat_m, rec_measurements]   = sub_Sampling_Recovery(xs,ys,measurements,sampls_M,Phi,Psi,num_frequencies,fre_sparse,K_loop);

telapsed = toc(tstart);

% 如果是计算每个像素点的重构误差，这种方法有时不太客观，因为有些点本身测量噪声大，而重构的结果通常都是光滑的。

err = (norm(rec_measurements(:)-measurements(:)))/norm(measurements(:));       % 计算重构的h和测量的h之间的欧氏距离/误差

%{
tm = reshape(measurements',[ys,xs,num_frequencies*2]);
tr = reshape(rec_measurements',[ys,xs,num_frequencies*2]);
th = reshape(hat_m',[ys,xs,sampls_M*2]);
[y,x] = deal(94,64);                  % noisy:(94,64),(70,60),(70,30),  noiseless: (150,30),  对于Heide的数据集，同样的像素点位置，不同数据集噪声不同。为什么会这样呢？
p1 = squeeze(tm(x,y,:));
p2 = squeeze(tr(x,y,:));
p3 = squeeze(th(x,y,:));
figure;
plot([p1 p2]);
hold on;stem([fre_pos fre_pos + 199], p3(:));
legend('Full measurements','Recovered measurements','Selected measurements');title([mM.filename ', (', num2str(y) ',' num2str(x),')']);
%}

mM.rec_measurements = reshape(rec_measurements',[ys,xs,num_frequencies*2]);
mM.measurements = reshape(measurements',[ys,xs,num_frequencies*2]);
mM.num_frequencies = num_frequencies;
mM.rec_M = [Psi{1};Psi{2}];

global output_folder 
aux = '_CS';
algo = '_Project_1';

save( sprintf(['%s/' mM.filename algo aux '_' num2str(sampls_M) '_QH.mat'], output_folder), 'err', 'telapsed','sampls_M'); 

return

function [Phi, fre_pos, sampls_M] = Draw_L_fre(num_frequencies, Phi, sparse_k)
N = num_frequencies;
switch Phi  
    case 'Random'     % K是根据信号的FFT计算的，观测信号随机给        
        sampls_M = ceil(sparse_k*log(N/sparse_k));     %  测量数(sampls_M>=sparse_K*log(N/sparse_K),至少40,但有出错的概率)，这个是理论值，采多了，效果也好不哪去。20个基本可以 
        tmp = randperm(N);                       % randi([1,N],sampls_M,1)这种方法生成的随机数有可能是重复的。 
        fre_pos = sort(tmp(1:sampls_M));   
        if fre_pos(1) ~= 1                          % 把起点和终点取上，避免端点发散
            fre_pos(1) = 1;
        end
        if fre_pos(sampls_M) ~= N  
            fre_pos(sampls_M) = N; 
        end
        Phi = zeros(sampls_M,N);
        for k = 1:size(Phi,1)        % 定制生成Phi
            Phi(k,fre_pos(k)) = 1;
        end  
end
return

%{
function [s,Phi,sparse_K,sampls_M,fre_pos] = Compressed_Sampling(x,Phi)
N = length(x);
switch Phi
    case 'Tran3'     % 3个观测信号
        sparse_K = 1;                  %  稀疏度(做FFT可以看出来)
        fre_pos = [29 75 158];     %  采样3个值   [58 75 93] [64 93 120]
        sampls_M = length(fre_pos);       %  测量数(sampls_M>=sparse_K*log(N/sparse_K),至少40,但有出错的概率) 
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)       % 定制生成Phi
            Phi(i,fre_pos(i)) = 1;
        end
    %         Phi = Phi+ (0.001 + 0.005.*randn(sampls_M,N));   % 给Phi添加均值为0.001，方差为0.005的随机数，避免方程不可解      

    case 'TranL20'     % 20个观测信号，使用低频多
        sparse_K = 10;     
        sampls_M = 20;     %  测量数(sampls_M>=sparse_K*log(N/sparse_K),至少40,但有出错的概率)，这个是理论值，采多了，效果也好不哪去。20个结果可以，30个没必要 
        fre_pos = [1 2 6 20 35 49 61 74 89 102 117 131 149 160 170 184 194 199 28 108];   % 199(205是采用Fourier basis用的，这里最大不超过199)后面的两个频率是在低频段选取的，担心高频部分测量误差大，没取。

        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)       % 定制生成Phi
            Phi(i,fre_pos(i)) = 1;
        end
    %         Phi = Phi+ (0.001 + 0.005.*randn(sampls_M,N));   % 给Phi添加均值为0.001，方差为0.005的随机数，避免方程不可解       
    case 'TranH20'     % 20个观测信号，使用高频多
        sparse_K = 10;     
        sampls_M = 20;     %  测量数(sampls_M>=sparse_K*log(N/sparse_K),至少40,但有出错的概率)，这个是理论值，采多了，效果也好不哪去。20个结果可以，30个没必要 
        fre_pos = [1 2 6 20 35 49 61 74 89 102 117 131 149 160 170 184 194 205 210 221];   % 在205后面，取了两个高频。        
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)        % 定制生成Phi
            Phi(i,fre_pos(i)) = 1;
        end
    %         Phi = Phi+ (0.001 + 0.055.*randn(sampls_M,N));        % 给Phi添加均值为0.001，方差为0.005的随机数，避免方程不可解 
    case 'Random'     % K是根据信号的FFT计算的，观测信号随机给
        sparse_K = input('稀疏度 sparse_K ：  ');     
        sampls_M = ceil(sparse_K*log(N/sparse_K));     %  测量数(sampls_M>=sparse_K*log(N/sparse_K),至少40,但有出错的概率)，这个是理论值，采多了，效果也好不哪去。20个结果可以，30个没必要 
        fre_pos = randi([1,N],sampls_M,1);   
        fre_pos = sort(fre_pos);
        fre_pos(1) = 1;fre_pos(sampls_M) = N;  % 把起点和终点取上，避免端点发散
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)        % 定制生成Phi
            Phi(i,fre_pos(i)) = 1;
        end
    case '7'
        sparse_K = 2;     
        sampls_M = 7;     %  测量数(sampls_M>=sparse_K*log(N/sparse_K),至少40,但有出错的概率)，这个是理论值，采多了，效果也好不哪去。20个结果可以，30个没必要 
        fre_pos = 3+[6 33 64 90 121 148 173];
        fre_pos(1) = 1;fre_pos(sampls_M) = N;  % 把起点和终点取上，避免端点发散
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)        % 定制生成Phi
            Phi(i,fre_pos(i)) = 1;
        end        
end
s = Phi*x.';             % 对信号压缩传感
return
%}

% 函数 CS_OMP_Psi 是根据 function [hat_x, M, tol] = Tran_Simu_CS_Principle_Analy_17(x,Phi,Psi,N,K)
% 中的 function hat_x = CS_OMP_1(K,fre_sparse,Phi,s,N,M,Psi)改写而来

function hat_x = CS_OMP_Psi(m,fre_sparse,Phi,s,N,M,Psi)         
%m = 2*K;                                            %  算法迭代次数(m>=K)，设x是K-sparse的

T = Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵) 

hat_y = zeros(1,N);                                 %  待重构的谱域(变换域)向量
Aug_t = [];                                         %  增量矩阵(初始值为空矩阵)
r_n = s;                                            %  残差值  

for times = 1:m                                    %  迭代次数(有噪声的情况下,该迭代次数为K)  
    for col = 1:N                                  %  恢复矩阵的所有列向量 
        product(col) = abs(T(:,col)'*r_n);         %  恢复矩阵的列向量和残差的投影系数(内积值)
    end
    [val,pos] = max(product);                       %  最大投影系数对应的位置，即找到一个其标记看上去与收集到的数据相关的小波 
    Aug_t = [Aug_t, T(:,pos)];                      %  矩阵扩充 
    T(:,pos) = zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零），在数据中去除这个标记的所有印迹  
    aug_y = (Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小
    r_n = s-Aug_t*aug_y;                            %  残差
    pos_array(times) = pos;                         %  纪录最大投影系数的位置    
end
hat_y(pos_array) = aug_y;                           %  重构的谱域向量 
if fre_sparse
    hat_x = real(Psi'*hat_y.');                     %  做逆傅里叶变换重构得到时域信号  
else
    hat_x = hat_y';                                 %  修正一下矩阵维数
end

return

function  [hat_measurements, rec_measurements] = sub_Sampling_Recovery(xs,ys,measurements,sampls_M,Phi,Psi,num_frequencies,fre_sparse,K_loop) 
parfor xy = 1:xs*ys      % 压缩采样
    m = measurements(:,xy);
    hat_m = zeros(2*sampls_M,1);
    hat_m(1:sampls_M) = Phi*m(1:num_frequencies);                    
    hat_m(sampls_M+1:2*sampls_M) = Phi*m(num_frequencies+1:end); 
    hat_measurements(:,xy) = hat_m;
end 

rec_measurements = zeros(size(measurements));
% tstart = tic;
parfor xy = 1:xs*ys      % 重构输入信号
    s = hat_measurements(:,xy);
    if any(s)             % 如果采样点全为0，不必计算
        hat_m = zeros(2*num_frequencies,1);
        hat_m(1:num_frequencies) = CS_OMP_Psi(K_loop,fre_sparse,Phi,s(1:sampls_M),num_frequencies,sampls_M,Psi{1});
        hat_m(num_frequencies+1:num_frequencies*2) = CS_OMP_Psi(K_loop,fre_sparse,Phi,s(sampls_M+1:end),num_frequencies,sampls_M,Psi{2}); 
        rec_measurements(:,xy) = hat_m; 
    end
end
% telapsed = toc(tstart);
return

function [I_rec, mdist] = Pseudoinverse_Wavelet_based_TI_rec(mM)
if isfield(mM,'rec_measurements')
    measurements = mM.rec_measurements;
    M = mM.rec_M;
else
    measurements = mM.measurements;
    M = mM.M;
end


xs = mM.imagedims(1);
ys = mM.imagedims(2);
measurements = (reshape(measurements,[xs*ys,size(M,1)]))'; % 数据准备；


weight = 1; 
load('noise_Heide.mat');

noise = noise/weight; 
if ~(size(noise,1) == size(M,1))
    noise = noise(1:size(M,1),:);
end
beta = pinv(M+noise) * measurements;

waveletname = 'db5';    %'sym5'
N = 5; 
PX = zeros(size(beta));

tic;
coef = 1;  % 2.2是根据结果反向估计的值，但实际上没啥用。因为整体都放到2.2倍，归一化后跟没放大是一样的。  
parfor xy = 1:size(beta,2)              % 小波分解与重构
    bp = beta(:,xy);
    [Cb,Lb] = wavedec(bp,N,waveletname);  % 用小波对信号进行N尺度小波分解
    ap = wrcoef('a',Cb,Lb,waveletname,3);  % cA+cD5+cD4=cA3   cA4：采样值为3时
    PX(:,xy) = coef * ap;   % 2.2是根据结果反向估计的值    
end

I_rec = zeros(size(beta));
parfor xy = 1:size(PX,2)
    ap = PX(:,xy);
    ap = Find_Main_Peak_without_Thresh(ap);  % 单纯的取一个大于0的主峰，与阈值无关。
    I_rec(:,xy) = ap;   % 2.2是根据结果反向估计的值
end
telapsed = toc;

I_rec = reshape(I_rec',[ys,xs,size(M,2)]);

mdist = TI_Postprocess(mM,I_rec,'_wavelet_Project_1',telapsed); 

%{
if isfield(mM,'simu')
    simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
    mdist = norm(simu - I_rec)/norm(simu);
    
    simu = mM.simu;
    I_rec = reshape(I_rec',[ys,xs,size(M,2)]);
    [y,x] = deal(86,96);             % 分析具体的像素点：低频：(136,83);(149,88);(150,91);(148,95);(120,140);单峰：(100,40) (70,60);(31,86);(92,57);单峰2个return：(40,30);(53,45);双峰：(86,96);
    rp = squeeze(I_rec(x,y,:));
    sp = squeeze(simu(x,y,:));
    figure;plot([sp rp]);legend('True','Wavelet');title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);    
    xlabel('Time Series');ylabel( 'Amplitude' );      
else    
    % 计算观测值的距离
    rec_m = M * I_rec;        
    mdist = norm(measurements - rec_m) / norm(measurements);   % 计算欧式距离，受观测误差值影响较大
    %{
    rec_m = reshape(rec_m',[ys,xs,size(M,1)]);
    true_m = mM.measurements;

    I_rec = reshape(I_rec',[ys,xs,size(M,2)]);
    
    [y,x] = deal(150,91);             % 分析具体的像素点：(60,70);(30,40);(149,88);(83,136);(150,91);(131,67);(86,31);(95,148);(57,92);(45,53); (75,164)
    rp = squeeze(I_rec(x,y,:));
    figure;plot(rp);title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);
    
    rm = squeeze(rec_m(x,y,:));
    tm = squeeze(true_m(x,y,:));
    figure;plot([tm rm]);
    legend('True','Rec');
    title([mM.filename ' measurements at (', num2str(y) ',' num2str(x),')']);  
    %}  
end
%
PX2 = I_rec ./ max(I_rec(:));  % 归一化
if isfield(mM,'simu')
    simu = mat2gray(mM.simu);
elseif isequal(mM.filename, 'DiscoBall') 
    PX2 = (PX2 * 1.0) .^ (1/5.2); %Gamma
else
    PX2 = (PX2 * 1.0) .^ (1/2.2); %Gamma
end

%Show transient image
if isfield(mM,'simu')
    figure;
    for t = 1:size(PX2,3)    
        subplot(1,2,1);imshow(PX2(:,:,t));title(sprintf('Transient image frame %d', t));
        subplot(1,2,2);imshow(simu(:,:,t)); title(sprintf('True image frame %d', t));       
        pause(0.1)
    end
else
    Show_Results(PX2);
end

Show_Results(PX2,[mM.filename '_wavelet_Project_1_0'],'mp4');    % 生成视频文件

output_folder = './TransientImaging/Outputs'; 
if isfield(mM,'simu')
    save( sprintf(['%s/' mM.filename '_wavelet_Project_1_0.mat'], output_folder), 'I_rec', 'mdist');
else
    save( sprintf(['%s/' mM.filename '_wavelet_Project_1_0.mat'], output_folder), 'I_rec','mdist','true_m','rec_m');    
end
% system('shutdown -s');
%}

%{ 
% 保存为图片
for t = 1:size(PX2,3)  
    I_show = PX2(:,:,t);
    imwrite(I_show,[mM.filename '_' num2str(t) '_Project_1_0.bmp']);      
end 
%}
return

