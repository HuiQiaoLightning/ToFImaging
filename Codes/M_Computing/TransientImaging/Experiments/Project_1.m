function Project_1(mM)
% ��һ��paper�е�Heide�ɼ����ݵ�ʵ�飬

%% Compressed sampling and recovery. Testing datasets are Heide's datasets.

[mM, err] = Compressed_Sampling_Recovery(mM);    % Algo1  ѹ���ɼ�

i_step_UI(mM);        % Algo1 + Heide_i_step   �������������������ѹ���������ع����ġ�  ���Բ���Tomato��

Lin_v03_v04_v05_UI(mM);                  % Algo1 + Lin �������������������ѹ���������ع����ġ� ���Բ���Tomato��    

Project_6(mM);                           % Algo1 + Algo2.1

%{ 
Wavelet_based transient image reconstruction  ע�⣺Heide���㷨�ǹ�һ���󣬵���Gamma�ġ�  
���þ����������󣬽�ȡ��3���ò�ź������Ǹ���ֵ�ع�˲̬�źš�
Algo2.2 ��ָ���߶��£�scale = 3���ع�˲̬ͼ��Heide�����ݼ������������ǻ����������ʱ��ľ���ȷ���ġ����ܲ���Tomato��   
%}

% [I_recW,mdistW] = Pseudoinverse_Wavelet_based_TI_rec(mM);     

return

%%
function [mM, err] = Compressed_Sampling_Recovery(mM)
%% ����׼����
if isfield(mM, 'simu')
    if mM.noise_measurements == 0           % ����ǲ����������ݣ��ӵ��������̫��ʧ������
        measurements = mM.true_measurements;
    else
        measurements = mM.noise_measurements;   % ������������
    end
else
    measurements = mM.measurements;  % �ɼ�������
end

xs = mM.imagedims(1);
ys = mM.imagedims(2);
    
Nf = size(mM.M,1)/2;Nt = size(mM.M,2);   % Nf��Ƶ�ʲ�����;Nt��ʱ������㡣��ͬ��ʵ�黷�������ܲ�һ�����������Heide�����ݼ��ġ�
num_frequencies = Nt;

measurements(:,:,Nf+Nt+1:end) = [];
measurements(:,:,Nt+1:Nf) = [];

measurements = (reshape(measurements,[xs*ys,num_frequencies*2]))'; % ����׼����

%% ���ѹ�������Ĳ���
% sparsity basis
Psi{1} = mM.M(1:Nt,:);
Psi{2} = mM.M(Nf+1:Nf+Nt,:);


%% ѹ������        �ؽ�������Ҫ�޸ģ�ֱ������Щ�������ع������������ع������źţ����Ҫ��ô����
% (1)
Phi = 2;            % ǰ3���ǳ���ķ���
switch Phi            
    case 1
        Phi = 'TranL20';      % ׼��Phi��ʹ��һ�����Ƶľ�����ζ��ֻ����20����ͬƵ�ʵ�ֵ  �������Ƶ��
    case 2
        Phi = 'Random';       % ʹ��55����������㣬���������Ƶ��õ���
end

tstart = tic;
% (2)  ���ѹ�������Ĳ���
sparse_K = 28;  % input('ϡ��� sparse_k ��  ');     
[Phi, fre_pos, sampls_M] = Draw_L_fre(num_frequencies, Phi, sparse_K);  

%m = measurements(1,:);
%[hat_m(1,:),Phi,sparse_K,sampls_M,fre_pos] = Compressed_Sampling(m(1,1:num_frequencies),Phi);   % ��xpos���������������Աȷ���ʱ����stemͼ�ȽϷ��㡣
   
% �ع�ʱ�������� 
fre_sparse = 1;    % Ϊ0ʱ����ʱ��ϡ�����,���ﰴƵ��ϡ�迼��
if sampls_M > 30
    K_loop = 9;        %  �㷨��������(m>=K)����x��K-sparse��,��ʱ��m̫�󣬷���������
else
    K_loop = 2*sparse_K;
end

[hat_m, rec_measurements]   = sub_Sampling_Recovery(xs,ys,measurements,sampls_M,Phi,Psi,num_frequencies,fre_sparse,K_loop);

telapsed = toc(tstart);

% ����Ǽ���ÿ�����ص���ع������ַ�����ʱ��̫�͹ۣ���Ϊ��Щ�㱾����������󣬶��ع��Ľ��ͨ�����ǹ⻬�ġ�

err = (norm(rec_measurements(:)-measurements(:)))/norm(measurements(:));       % �����ع���h�Ͳ�����h֮���ŷ�Ͼ���/���

%{
tm = reshape(measurements',[ys,xs,num_frequencies*2]);
tr = reshape(rec_measurements',[ys,xs,num_frequencies*2]);
th = reshape(hat_m',[ys,xs,sampls_M*2]);
[y,x] = deal(94,64);                  % noisy:(94,64),(70,60),(70,30),  noiseless: (150,30),  ����Heide�����ݼ���ͬ�������ص�λ�ã���ͬ���ݼ�������ͬ��Ϊʲô�������أ�
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
    case 'Random'     % K�Ǹ����źŵ�FFT����ģ��۲��ź������        
        sampls_M = ceil(sparse_k*log(N/sparse_k));     %  ������(sampls_M>=sparse_K*log(N/sparse_K),����40,���г���ĸ���)�����������ֵ���ɶ��ˣ�Ч��Ҳ�ò���ȥ��20���������� 
        tmp = randperm(N);                       % randi([1,N],sampls_M,1)���ַ������ɵ�������п������ظ��ġ� 
        fre_pos = sort(tmp(1:sampls_M));   
        if fre_pos(1) ~= 1                          % �������յ�ȡ�ϣ�����˵㷢ɢ
            fre_pos(1) = 1;
        end
        if fre_pos(sampls_M) ~= N  
            fre_pos(sampls_M) = N; 
        end
        Phi = zeros(sampls_M,N);
        for k = 1:size(Phi,1)        % ��������Phi
            Phi(k,fre_pos(k)) = 1;
        end  
end
return

%{
function [s,Phi,sparse_K,sampls_M,fre_pos] = Compressed_Sampling(x,Phi)
N = length(x);
switch Phi
    case 'Tran3'     % 3���۲��ź�
        sparse_K = 1;                  %  ϡ���(��FFT���Կ�����)
        fre_pos = [29 75 158];     %  ����3��ֵ   [58 75 93] [64 93 120]
        sampls_M = length(fre_pos);       %  ������(sampls_M>=sparse_K*log(N/sparse_K),����40,���г���ĸ���) 
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)       % ��������Phi
            Phi(i,fre_pos(i)) = 1;
        end
    %         Phi = Phi+ (0.001 + 0.005.*randn(sampls_M,N));   % ��Phi��Ӿ�ֵΪ0.001������Ϊ0.005������������ⷽ�̲��ɽ�      

    case 'TranL20'     % 20���۲��źţ�ʹ�õ�Ƶ��
        sparse_K = 10;     
        sampls_M = 20;     %  ������(sampls_M>=sparse_K*log(N/sparse_K),����40,���г���ĸ���)�����������ֵ���ɶ��ˣ�Ч��Ҳ�ò���ȥ��20��������ԣ�30��û��Ҫ 
        fre_pos = [1 2 6 20 35 49 61 74 89 102 117 131 149 160 170 184 194 199 28 108];   % 199(205�ǲ���Fourier basis�õģ�������󲻳���199)���������Ƶ�����ڵ�Ƶ��ѡȡ�ģ����ĸ�Ƶ���ֲ�������ûȡ��

        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)       % ��������Phi
            Phi(i,fre_pos(i)) = 1;
        end
    %         Phi = Phi+ (0.001 + 0.005.*randn(sampls_M,N));   % ��Phi��Ӿ�ֵΪ0.001������Ϊ0.005������������ⷽ�̲��ɽ�       
    case 'TranH20'     % 20���۲��źţ�ʹ�ø�Ƶ��
        sparse_K = 10;     
        sampls_M = 20;     %  ������(sampls_M>=sparse_K*log(N/sparse_K),����40,���г���ĸ���)�����������ֵ���ɶ��ˣ�Ч��Ҳ�ò���ȥ��20��������ԣ�30��û��Ҫ 
        fre_pos = [1 2 6 20 35 49 61 74 89 102 117 131 149 160 170 184 194 205 210 221];   % ��205���棬ȡ��������Ƶ��        
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)        % ��������Phi
            Phi(i,fre_pos(i)) = 1;
        end
    %         Phi = Phi+ (0.001 + 0.055.*randn(sampls_M,N));        % ��Phi��Ӿ�ֵΪ0.001������Ϊ0.005������������ⷽ�̲��ɽ� 
    case 'Random'     % K�Ǹ����źŵ�FFT����ģ��۲��ź������
        sparse_K = input('ϡ��� sparse_K ��  ');     
        sampls_M = ceil(sparse_K*log(N/sparse_K));     %  ������(sampls_M>=sparse_K*log(N/sparse_K),����40,���г���ĸ���)�����������ֵ���ɶ��ˣ�Ч��Ҳ�ò���ȥ��20��������ԣ�30��û��Ҫ 
        fre_pos = randi([1,N],sampls_M,1);   
        fre_pos = sort(fre_pos);
        fre_pos(1) = 1;fre_pos(sampls_M) = N;  % �������յ�ȡ�ϣ�����˵㷢ɢ
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)        % ��������Phi
            Phi(i,fre_pos(i)) = 1;
        end
    case '7'
        sparse_K = 2;     
        sampls_M = 7;     %  ������(sampls_M>=sparse_K*log(N/sparse_K),����40,���г���ĸ���)�����������ֵ���ɶ��ˣ�Ч��Ҳ�ò���ȥ��20��������ԣ�30��û��Ҫ 
        fre_pos = 3+[6 33 64 90 121 148 173];
        fre_pos(1) = 1;fre_pos(sampls_M) = N;  % �������յ�ȡ�ϣ�����˵㷢ɢ
        Phi = zeros(sampls_M,N);
        for i = 1:size(Phi,1)        % ��������Phi
            Phi(i,fre_pos(i)) = 1;
        end        
end
s = Phi*x.';             % ���ź�ѹ������
return
%}

% ���� CS_OMP_Psi �Ǹ��� function [hat_x, M, tol] = Tran_Simu_CS_Principle_Analy_17(x,Phi,Psi,N,K)
% �е� function hat_x = CS_OMP_1(K,fre_sparse,Phi,s,N,M,Psi)��д����

function hat_x = CS_OMP_Psi(m,fre_sparse,Phi,s,N,M,Psi)         
%m = 2*K;                                            %  �㷨��������(m>=K)����x��K-sparse��

T = Phi*Psi';                                       %  �ָ�����(��������*�������任����) 

hat_y = zeros(1,N);                                 %  ���ع�������(�任��)����
Aug_t = [];                                         %  ��������(��ʼֵΪ�վ���)
r_n = s;                                            %  �в�ֵ  

for times = 1:m                                    %  ��������(�������������,�õ�������ΪK)  
    for col = 1:N                                  %  �ָ���������������� 
        product(col) = abs(T(:,col)'*r_n);         %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ)
    end
    [val,pos] = max(product);                       %  ���ͶӰϵ����Ӧ��λ�ã����ҵ�һ�����ǿ���ȥ���ռ�����������ص�С�� 
    Aug_t = [Aug_t, T(:,pos)];                      %  �������� 
    T(:,pos) = zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩����������ȥ�������ǵ�����ӡ��  
    aug_y = (Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  ��С����,ʹ�в���С
    r_n = s-Aug_t*aug_y;                            %  �в�
    pos_array(times) = pos;                         %  ��¼���ͶӰϵ����λ��    
end
hat_y(pos_array) = aug_y;                           %  �ع����������� 
if fre_sparse
    hat_x = real(Psi'*hat_y.');                     %  ���渵��Ҷ�任�ع��õ�ʱ���ź�  
else
    hat_x = hat_y';                                 %  ����һ�¾���ά��
end

return

function  [hat_measurements, rec_measurements] = sub_Sampling_Recovery(xs,ys,measurements,sampls_M,Phi,Psi,num_frequencies,fre_sparse,K_loop) 
parfor xy = 1:xs*ys      % ѹ������
    m = measurements(:,xy);
    hat_m = zeros(2*sampls_M,1);
    hat_m(1:sampls_M) = Phi*m(1:num_frequencies);                    
    hat_m(sampls_M+1:2*sampls_M) = Phi*m(num_frequencies+1:end); 
    hat_measurements(:,xy) = hat_m;
end 

rec_measurements = zeros(size(measurements));
% tstart = tic;
parfor xy = 1:xs*ys      % �ع������ź�
    s = hat_measurements(:,xy);
    if any(s)             % ���������ȫΪ0�����ؼ���
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
measurements = (reshape(measurements,[xs*ys,size(M,1)]))'; % ����׼����


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
coef = 1;  % 2.2�Ǹ��ݽ��������Ƶ�ֵ����ʵ����ûɶ�á���Ϊ���嶼�ŵ�2.2������һ�����û�Ŵ���һ���ġ�  
parfor xy = 1:size(beta,2)              % С���ֽ����ع�
    bp = beta(:,xy);
    [Cb,Lb] = wavedec(bp,N,waveletname);  % ��С�����źŽ���N�߶�С���ֽ�
    ap = wrcoef('a',Cb,Lb,waveletname,3);  % cA+cD5+cD4=cA3   cA4������ֵΪ3ʱ
    PX(:,xy) = coef * ap;   % 2.2�Ǹ��ݽ��������Ƶ�ֵ    
end

I_rec = zeros(size(beta));
parfor xy = 1:size(PX,2)
    ap = PX(:,xy);
    ap = Find_Main_Peak_without_Thresh(ap);  % ������ȡһ������0�����壬����ֵ�޹ء�
    I_rec(:,xy) = ap;   % 2.2�Ǹ��ݽ��������Ƶ�ֵ
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
    [y,x] = deal(86,96);             % ������������ص㣺��Ƶ��(136,83);(149,88);(150,91);(148,95);(120,140);���壺(100,40) (70,60);(31,86);(92,57);����2��return��(40,30);(53,45);˫�壺(86,96);
    rp = squeeze(I_rec(x,y,:));
    sp = squeeze(simu(x,y,:));
    figure;plot([sp rp]);legend('True','Wavelet');title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);    
    xlabel('Time Series');ylabel( 'Amplitude' );      
else    
    % ����۲�ֵ�ľ���
    rec_m = M * I_rec;        
    mdist = norm(measurements - rec_m) / norm(measurements);   % ����ŷʽ���룬�ܹ۲����ֵӰ��ϴ�
    %{
    rec_m = reshape(rec_m',[ys,xs,size(M,1)]);
    true_m = mM.measurements;

    I_rec = reshape(I_rec',[ys,xs,size(M,2)]);
    
    [y,x] = deal(150,91);             % ������������ص㣺(60,70);(30,40);(149,88);(83,136);(150,91);(131,67);(86,31);(95,148);(57,92);(45,53); (75,164)
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
PX2 = I_rec ./ max(I_rec(:));  % ��һ��
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

Show_Results(PX2,[mM.filename '_wavelet_Project_1_0'],'mp4');    % ������Ƶ�ļ�

output_folder = './TransientImaging/Outputs'; 
if isfield(mM,'simu')
    save( sprintf(['%s/' mM.filename '_wavelet_Project_1_0.mat'], output_folder), 'I_rec', 'mdist');
else
    save( sprintf(['%s/' mM.filename '_wavelet_Project_1_0.mat'], output_folder), 'I_rec','mdist','true_m','rec_m');    
end
% system('shutdown -s');
%}

%{ 
% ����ΪͼƬ
for t = 1:size(PX2,3)  
    I_show = PX2(:,:,t);
    imwrite(I_show,[mM.filename '_' num2str(t) '_Project_1_0.bmp']);      
end 
%}
return

