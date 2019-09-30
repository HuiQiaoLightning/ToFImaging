function Project_6(mM)
% 第一个paper中的Heide仿真和采集数据的实验：tomato，对整幅图像用一个噪声矩阵进行重构，取平均值，第3、4、5层概貌信号都计算了。
% 可以保留noise矩阵；Project_3的升级版
% Algo 2.1
xs = mM.imagedims(1);
ys = mM.imagedims(2);

if isfield(mM,'rec_measurements')
    measurements = mM.rec_measurements;
    M = mM.rec_M; 
    weight = 1;
elseif isfield(mM,'simu')
    M = mM.M;
    if mM.noise_measurements == 0              % 这种情况没测，因为实际都是有噪声的
        measurements = mM.true_measurements;
        weight = 1;
        disp('这种情况没测，因为实际都是有噪声的');
        return;
    else
        measurements = mM.noise_measurements;
%     simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
        weight = 3;       % 这个权重是当仿真数据采用加固定的Poisson噪声时采用的，是一个需要根据具体数据进行估计的值
    end    
else
    M = mM.M;
    measurements = mM.measurements;
    weight = 1;  % 1;2;   
end

measurements = (reshape(measurements,[xs*ys,size(M,1)]))'; % 数据准备；
waveletname = 'db5';    %'sym5'
N = 5; 
num_noise = 100;  % 17; 100;         % 50次，15分钟  tomato:7次；

tstart = tic;

% 用3-D数组实现
rec_ti = zeros(3,size(M,2),size(measurements,2));
for nn = 1:num_noise     
%     noise{nn} = rand(size(M))/weight;                 % 保留noise矩阵
%     beta = pinv(M+noise{nn}) * measurements;   
    noise = rand(size(M))/weight;  
    beta = pinv(M+noise) * measurements;       
    parfor xy = 1:size(beta,2)   
        bp = beta(:,xy);
        [Cb,Lb] = wavedec(bp,N,waveletname); 
        ap3(:,xy) = wrcoef('a',Cb,Lb,waveletname,3);
        ap4(:,xy) = wrcoef('a',Cb,Lb,waveletname,4);
        ap5(:,xy) = wrcoef('a',Cb,Lb,waveletname,5);
%         figure;plot([ap3(:,pixel) ap4(:,pixel) ap5(:,pixel)]);grid on;  legend('ap3','ap4','ap5');      
    end
    rec_ti(1,:,:) = squeeze(rec_ti(1,:,:)) + ap3;
    rec_ti(2,:,:) = squeeze(rec_ti(2,:,:)) + ap4;
    rec_ti(3,:,:) = squeeze(rec_ti(3,:,:)) + ap5;
end

%{
% 用cell结构实现
for nn = 1:num_noise     
    noise = rand(size(M))/weight;  
    beta = pinv(M+noise) * measurements';    
    parfor pixel = 1:size(measurements,1)   
        rp = beta(:,pixel);
        [Cb,Lb] = wavedec(rp,N,waveletname); 
        ap3(:,pixel) = wrcoef('a',Cb,Lb,waveletname,3);
        ap4(:,pixel) = wrcoef('a',Cb,Lb,waveletname,4);
        ap5(:,pixel) = wrcoef('a',Cb,Lb,waveletname,5);
%         figure;plot([ap3(:,pixel) ap4(:,pixel) ap5(:,pixel)]);grid on;  legend('ap3','ap4','ap5');      
    end
    ap{nn,1} = ap3;
    ap{nn,2} = ap4;
    ap{nn,3} = ap5;
end

mean_ap3 = 0;
mean_ap4 = 0;
mean_ap5 = 0;
parfor nn = 1:num_noise
    mean_ap3 = mean_ap3 + ap{nn,1};
    mean_ap4 = mean_ap4 + ap{nn,2};
    mean_ap5 = mean_ap5 + ap{nn,3};
end

mean_ap3 = mean_ap3 ./ num_noise;
mean_ap4 = mean_ap4 ./ num_noise;
mean_ap5 = mean_ap5 ./ num_noise;
%}
ap3 = [];
ap4 = [];
ap5 = [];
PX = max(rec_ti,0);     % rec_ti:乘系数之前的重构信号   % 负数在下面做开方运算时，会出现复数，所以取大于零的部分
PX = PX / num_noise;
% 计算线性系数
parfor xy = 1:size(beta,2)
    single_measurement = measurements(:,xy);
    tmp = single_measurement - mean(single_measurement);
    
    tmp_ap = (squeeze(PX(1,:,xy)))';
    tmp1 = M * tmp_ap;
    tmp1 = tmp1 - mean(tmp1);
    coef = max(abs(tmp)) / max(abs(tmp1));
%     tmp1 = coef * tmp1;
%     figure;plot([tmp tmp1]);
    ap3(:,xy) = coef * tmp_ap;
%     figure;plot([tmp_ap ap3(:,xy)]);
    
    tmp_ap = (squeeze(PX(2,:,xy)))';
    tmp1 = M * tmp_ap;
    tmp1 = tmp1 - mean(tmp1);
    coef = max(abs(tmp)) / max(abs(tmp1));
%     tmp1 = coef * tmp1;
%     figure;plot([tmp tmp1]);
    ap4(:,xy) = coef * tmp_ap;
%     figure;plot([tmp_ap ap4(:,xy)]);
    
    tmp_ap = (squeeze(PX(3,:,xy)))';
    tmp1 = M * tmp_ap;
    tmp1 = tmp1 - mean(tmp1);
    coef = max(abs(tmp)) / max(abs(tmp1));
%     tmp1 = coef * tmp1;
%     figure;plot([tmp tmp1]);
    ap5(:,xy) = coef * tmp_ap;
%     figure;plot([tmp_ap ap5(:,xy)]);
end

telapsed = toc(tstart);

rec_ti_coef(1,:,:) = ap3;
rec_ti_coef(2,:,:) = ap4;
rec_ti_coef(3,:,:) = ap5;

global output_folder 
if isfield(mM,'rec_measurements')
    aux = '_CS';
elseif isfield(mM,'simu')
    aux = '_Poisson';
elseif 442 == size(mM.M,1)
    aux = '';
end
algo = '_Wavelet_Project_6';

%{
for level = 1:3
    I_rec = squeeze(rec_ti_coef(level,:,:));
    
    if isfield(mM,'simu')
        simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
        mdist(level) = norm(simu - I_rec)/norm(simu);
        
        simu = mM.simu;
        I_rec = reshape(I_rec',[ys,xs,size(M,2)]);
        [y,x] = deal(100,40);             % 分析具体的像素点：低频：(136,83);(149,88);(150,91);(148,95);(120,140);单峰：(100,40) (70,60);(31,86);(92,57);单峰2个return：(40,30);(53,45);双峰：(86,96);
        rp = squeeze(I_rec(x,y,:));
        sp = squeeze(simu(x,y,:));
        figure;plot([sp rp]);legend('True','Wavelet');title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);    
        xlabel('Time Series');ylabel( 'Amplitude' );            
    else    
        % 计算观测值的距离
        rec_m = M * I_rec;        % 把噪声矩阵保存后，这个可以再计算
        mdist(level) = norm(measurements - rec_m) / norm(measurements);   % 计算欧式距离，受观测误差值影响较大
   
        rec_m = reshape(rec_m',[ys,xs,size(M,1)]);
        true_m = mM.measurements;

        I_rec = reshape(I_rec',[ys,xs,size(M,2)]);

        [y,x] = deal(42,70);             % 分析具体的像素点：(60,70);(30,40);(149,88);(83,136);(150,91);(131,67);(86,31);(95,148);(57,92);(45,53); (75,164)
        rp = squeeze(I_rec(x,y,:));
        figure;plot(rp);title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);

        rm = squeeze(rec_m(x,y,:));
        tm = squeeze(true_m(x,y,:));
        figure;plot([tm rm]);
        legend('True','Rec');
        title([mM.filename ' measurements at (', num2str(y) ',' num2str(x),')']);          
    end
    
    %
    PX2 = I_rec./max(I_rec(:));
    if isfield(mM,'simu')
        simu = mat2gray(mM.simu);
    elseif isequal(mM.filename, 'DiscoBall') 
        PX2 = (PX2 * 1.0) .^ (1/5.2); %Gamma   开方运算
    else
        PX2 = (PX2 * 1.0) .^ (1/2.2); %Gamma   开方运算
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

    Show_Results(PX2,[mM.filename algo aux '_' num2str(num_noise) '_' num2str(level+2)],'mp4');    % 生成视频文件
    
    % 只保留某一层重构信号
    if isfield(mM,'simu')
        save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_' num2str(level+2) '.mat'], output_folder), 'I_rec', 'mdist', 'num_noise', 'telapsed');
    else
        save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_' num2str(level+2) '.mat'], output_folder), 'I_rec','mdist', 'num_noise', 'telapsed','true_m','rec_m');    
    end

%     Show_Results(PX2(:,:,1:100),['mean_wavelet_' mM.filename  '_' num2str(num_noise) '_' num2str(level)  '_Project_3_0'],'mp4');    % 生成视频文件
    
    %{
    % 存成图片的格式
    currpath = pwd;     % 获取当前路径
    global output_folder 
    cd(output_folder );

    for t = 1:size(PX2,3)  
        I_show = PX2(:,:,t);
        imwrite(I_show,[mM.filename algo aux '_' num2str(num_noise) '_' num2str(level+2) '_' num2str(t) '.bmp']);      
    end  
    cd(currpath);       % 返回程序运行的当前路径
    %}    
end
%}
save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_QH.mat'], output_folder), 'rec_ti', 'rec_ti_coef', 'telapsed', 'num_noise');   % 加载该数据，上面显示的内容可重新生成。
return

