function Project_5(mM)
% 第一个paper中的Heide仿真和采集数据的实验：tomato，对每个像素分别用一个噪声矩阵（20个）进行重构，从3、4、5层概貌信号中选取。这个计算结果有噪点，但有可能找到结果较好的那个信号。

xs = mM.imagedims(1);
ys = mM.imagedims(2);
M = mM.M;

if isfield(mM,'simu')   
    if mM.noise_measurements == 0
        disp('仿真的观测值需要加噪声');
        return;
    else
        measurements = mM.noise_measurements;
        simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
        weight = 3;       % 这个权重是当仿真数据采用加固定的Poisson噪声时采用的，是一个需要根据具体数据进行估计的值
    end
else
    measurements = mM.measurements;
    weight = 1;
end

measurements = (reshape(measurements,[xs*ys,size(M,1)]))'; % 数据准备；
 
waveletname = 'db5';    %'sym5'
N = 5; 
num_noise = 100;   % 70

sp = 0;  % 调试用
tstart = tic;
parfor xy = 1:size(measurements,2)
    tmp = measurements(:,xy);    

%     sp = simu(:,xy);     % 仿真用    
    
    tmp = rec_pixel_2(M, weight, waveletname, N, num_noise, tmp, sp);    
    
    I_rec(:,xy) = tmp;
end
telapsed = toc(tstart);

%{
% 计算观测值的距离
if isfield(mM,'simu')
    simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
    mdist = norm(simu - I_rec)/norm(simu);
    
    simu = mM.simu;
    I_rec = reshape(I_rec',[ys,xs,size(M,2)]);
    [y,x] = deal(86,96);             % 分析具体的像素点：低频：(136,83);(149,88);(150,91);(148,95);(120,140);单峰：(100,40) (70,60);(31,86);(92,57);单峰2个return：(40,30);(53,45);双峰：(86,96);
    rp = squeeze(I_rec(x,y,:));
    sp = squeeze(simu(x,y,:));
    figure;plot([sp rp]);legend('True','Rand wavelet');title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);    
    xlabel('Time Series');ylabel( 'Amplitude' );      
else    
    % 计算观测值的距离
    rec_m = M * I_rec;        % 把     噪声矩阵保存后，这个可以再计算
    mdist = norm(measurements - rec_m) / norm(measurements);   % 计算欧式距离，受观测误差值影响较大

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
end

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
%}

aux = '_Poisson';

algo = '_wavelet_Project_5';

% Show_Results(PX2,[mM.filename algo aux '_' num2str(num_noise)],'mp4');    % 生成视频文件


global output_folder 

if isfield(mM,'simu')
    save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_QH.mat'], output_folder), 'I_rec', 'num_noise', 'telapsed');
else
    save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_QH.mat'], output_folder), 'I_rec', 'num_noise', 'telapsed');    
end

%{
% 保存为图片
currpath = pwd;     % 获取当前路径
cd(output_folder );
for t = 1:size(PX2,3)  
    I_show = PX2(:,:,t);
    imwrite(I_show,[mM.filename algo aux '_' num2str(num_noise) '_' num2str(t) '.bmp']);      
end 

cd(currpath);       % 返回程序运行的当前路径

% system('shutdown -s');
%}
return

% 逐个像素计算逆矩阵，版本2，也能实现版本1的功能。目前是只考虑第4、5层概貌信号，由小波分解引起的小的振荡可通过仅取大于0的部分来去掉。
function rec_ap = rec_pixel_2(M, weight, waveletname, N, num_noise, single_measurement, sp)
nn = 1;
while nn <= num_noise     
    noise = rand(size(M))/weight;  
    rp = pinv(M+noise) * single_measurement; 
    
    %{ 
    figure;plot([rp sp]);
    figure;plot(rp);
    %}
    
    [Cb,Lb] = wavedec(rp,N,waveletname); 
    level_start = 3;     % 从第3层概貌信号开始，如果是从第4层概貌信号开始的话，就是 4
    level_step = level_start - 1;
    for level = level_start:N        
        tmp = wrcoef('a',Cb,Lb,waveletname,level);   % 乘以系数2并不是最佳的，可根据细节的最大值出现的概率确定可能的峰值点。
        [tmpv,tmpi] = max(tmp);
        if tmpi <= 3 || tmpi >= numel(tmp)-2     % 如果极值点靠近两个端点，则认为是由小波振荡引起的，不考虑。       
            nn = nn - 1;        % for 循环后有nn+1，这里，nn不应该递增
            break;
        end
        % tmp,  max( tmp,0 ) 观察一下取0和不取0的区别，按取平均值的情况来看，不取0应该也是可以的。其实是不应该有取大于0的操作的。但取大于0后，能减少因振荡引起的噪声
        if level == 3
            ap{ nn,level - level_step } = Find_Main_Peak_without_Thresh(tmp);  % 对于第3层概貌信号，单纯的取一个主峰，与阈值无关。 因为第3层因振荡引起的波峰较多。             
        else
            ap{ nn,level - level_step } = max(tmp,0);
        end
        %{ 
        figure;plot([tmp ap{ nn,level - level_step } sp]);grid on;title(['Level = ' num2str(level)]);
        figure;plot([tmp ap{ nn,level - level_step }]);grid on;title(['Level = ' num2str(level)]);
        %}
    end
    nn = nn + 1;
end

tmp = single_measurement - mean(single_measurement);

for nn = 1:num_noise
    for level = 1:size( ap,2 )  
        tmp1 =  M * ap{nn,level};
        tmp1 = tmp1 - mean(tmp1);        % 去均值，否则，它们不在一个水平线上
%         figure;plot([tmp tmp1]);
        similarity( nn,level ) = norm(tmp - tmp1);
%         similarity( nn,level ) = std( tmp - tmp1 );      % 此处可以考虑使用 modwtxcorr       
    end
    
%     figure;plot( ap{ nn,end } );grid on;
    if size( ap,2 ) == 3
        %{ 
        if length( find( ap{ nn,3 }>0 ) ) > 150         % 对于低频信号，在第4、5层概貌信号中选；阈值是150
            similarity( nn,1 ) = Inf;
        end   
        %}
        %{
        if length( find( ap{ nn,3 }>0 ) ) > 150          % 对于低频信号，在第5层概貌信号中选
            similarity( nn,2 ) = Inf;
            similarity( nn,1 ) = Inf;
        elseif length( find( ap{ nn,3 }>0 ) ) > 120     % 对于有多个峰值的低频信号，在4、5层概貌信号中选
            similarity( nn,1 ) = Inf;
        end
        %}
         
        if length( find( ap{ nn,3 }>0 ) ) > 120         % % 对于低频信号，在第4、5层概貌信号中选；阈值是120；这个目前需要根据信号进行调整，应该是自适应
            similarity( nn,1 ) = Inf;
        end   
                
    elseif size( ap,2 ) == 2                             % 只用第4、5层概貌信号 
        if length( find( ap{ nn,2 }>0 ) ) > 150          % 对于低频信号，在第5层概貌信号中选
            similarity( nn,1 ) = Inf;
        end
    end
end

[ tmp,tmp1 ] = min( similarity, [  ], 2 );
[ tmp,ind ] = min( tmp );

rec_ap = ap{ ind,tmp1( ind ) };
%{
figure;plot([rec_ap sp]);grid on;
figure;plot(rec_ap);grid on;
%}
return
