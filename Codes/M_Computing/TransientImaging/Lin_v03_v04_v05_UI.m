function Lin_v03_v04_v05_UI(mM)
%{
03和04的函数名一样，内容没核实是否也一样。所以，两个路径分别指定。05的版本没在这里。
此处显示详细说明  重构结果：ver03的比ver04的延迟10多个帧，采用ver04。
但印象中，第一次下载的那个版本计算结果最好。
%}

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

num_frequencies = mM.num_frequencies;
num_phases = mM.num_phases;

%{
measurements = mM.measurements;
if isfield(mM,'simu')
    imagedims = mM.imagedims;
    num_phases = 2;    
    num_frequencies = size(measurements,3)/num_phases;
else
    measurements = double(measurements(end:-1:1,:,:,:));
    imagedims = mM.imagedims;
    num_frequencies = mM.num_frequencies;
    num_phases = mM.num_phases;
end
measurements = reshape(measurements, imagedims(2), imagedims(1), num_frequencies, num_phases);
%}
measurements = reshape(measurements, mM.imagedims(2), mM.imagedims(1), num_frequencies, num_phases);


ver = 'paper';
switch ver
    case '03' 
        needpath = genpath('./TransientImaging/Lin_Jingyu/Lin_code/transient_pub3_v03');    
        addpath (needpath);         % 添加路径
        load('./TransientImaging/Lin_Jingyu/Lin_code/transient_pub3_v03/transient_ubc/CM_2.mat')  % CM, tau0, An, phin, fre0, fre_step, tau_step, Kw
        k0 = 2*pi/1000;

        % freq_len = size(CM,1);
        freq_len = 210;
        fre_len0 = floor(fre0/fre_step);
        % nnum = size(An,2);  % number of harmonic components
        freq1 = fre0;  % MHz
        freqstep = fre_step; 
        freq2 = fre0+(freq_len-1)*fre_step;

        measurements = measurements(:,:,1:freq_len,:);
        cdata = measurements(:,:,:,1) + 1i*measurements(:,:,:,2); % complex image
        sz = size(cdata);

        tau0 = 22.44;    % 源代码中是 18 对应的tau_step=0.1; 22.44对应的tau_step=0.33; 
        tau_step = 0.33;   % 源代码中是 0.1;        time step. Heide et al. set tau_step = 0.33.
        tau_len = 1000;
        % tau0 = 18;
        % tau_step = 0.01;
        % tau_len = 3000;

        %% reconstruct
        fprintf('\nStart Running at %s.\n',datestr(now))
        tic; % atart timer

        % data rectification
        Hc1 = cdata;
        ep = exp(1i*(phin(:,1)));
        for i=1:freq_len
            Hc1(:,:,i) = cdata(:,:,i)./An(i,1)*ep(i);
        end

        % windowing
        wHc1 = Hc1;
        beta = 6;
        wind = kaiser(floor(freq2/freqstep)*2+1,beta);
        wind = wind(end-freq_len+1:end);
        for i=1:freq_len
            wHc1(:,:,i) = Hc1(:,:,i)*wind(i);
        end

        % trans_img = TransientwithLow(wHc1,freq1,freqstep,tau0,tau_step,tau_len);
        trans_img = ComputeTransient2(wHc1,freq1,freqstep,tau0,tau_step,tau_len);

        tElapsed = toc;
        fprintf('\nRunning Time: %f\n', tElapsed)

        % save(sprintf('trans_%s',sel_data),'trans_img','scenename',...
        %     'tau_step','tau_len','freq1','freq2','freqstep','freq_len');

        %% post processing
        trans_iron = IronTransient3(trans_img,1,2);
        % trans_iron = PulsizeTransient(trans_img,freq2,tau_step,0.1,0.2);
        % norm_img = NormalizeImage(trans_iron,2);
        norm_img = CompressHDR(trans_iron);
        % norm_img = CompressHDR(trans_img);

        %% plots
        y = 100;
        x = 40;

        a = permute(cdata(y,x,:),[3 2 1]);
        b = permute(Hc1(y,x,:),[3 2 1])*max(abs(An(:,1)));
        c = permute(wHc1(y,x,:),[3 2 1]);
        figure(1); plot([real(a) real(b) real(c)])

        tr1 = permute(trans_img(y,x,:),[3 2 1]);
        tr2 = permute(trans_iron(y,x,:),[3 2 1]);

        figure(2); plot([tr1/max(tr1) tr2/max(tr2)])
        % plot([tr1 tr2])

        %% show video
        % i0 = 600;
        i0 = 0;
        i_len = 600;
        for i=1:1:i_len
            figure(10)
        %     image(trans_iron(:,:,i0+i)*2*256); colormap(gray(256))
        %     image(trans_img(:,:,i0+i)*256); colormap(gray(256))
        %     image(ourFrame(:,:,i0+i)*256); colormap(gray(256))
            image(norm_img(:,:,i0+i)*256); colormap(gray(256))
            ttl = sprintf('scene frame: %.3d',i); title(ttl)
            pause(0.1)
        end
        
        disp('end v03');
        rmpath(needpath);             % 删除添加的路径；在运行他人的代码之前，要删除已添加的路径，以免函数重名，调用时出问题。    

        
        [y,x] = deal(20,40);  
        tr1 = permute(trans_img(y,x,:),[3 2 1]);
        tr2 = permute(trans_iron(y,x,:),[3 2 1]);
        tr3 = permute(norm_img(y,x,:),[3 2 1]);        
        figure;plot([tr1/max(tr1) tr2/max(tr2) tr3/max(tr3)]);
        legend('trans-img','trans-iron','norm-img');
        title(['(' num2str(y) ',' num2str(x) ')']);

    case '04'
        needpath = genpath('./TransientImaging/Lin_Jingyu/Lin_code/transient_code_pub4_v04');    
        addpath (needpath);         % 添加路径
        
        % time step. Heide et al. set tau_step = 0.33.
        tau0 = 22.44;  % ns   % 20 对应的tau_step=0.1; 22.44对应的tau_step=0.33;
        tau_step = 0.33;   % 0.1;   源代码中是0.1
        tau_len = 1000;
        % tau0 = 20;
        % tau_step = 0.01;
        % tau_len = 3600;
        
        load('./TransientImaging/Lin_Jingyu/Lin_code/transient_code_pub4_v04/calib/cfs/PmCF_ubc.mat','An','phin','fre0','fre_step')  
        % CM, tau0, An, phin, fre0, fre_step, tau_step, Kw

        % freq_len = size(CM,1);
        freq_len = 210;       % 210是作者全采样使用的频率长度； 压缩感知中使用标定矩阵时，为199；
        fre_len0 = floor(fre0/fre_step);
        % nnum = size(An,2);  % number of harmonic components
        freq1 = fre0;  % MHz
        freqstep = fre_step; 
        freq2 = fre0+(freq_len-1)*fre_step;

        measurements = measurements(:,:,1:freq_len,:);
        cdata = measurements(:,:,:,1) + 1i*measurements(:,:,:,2); % complex image
        sz = size(cdata);

        %% reconstruct
        fprintf('\nStart Running at %s.\n',datestr(now))
        tic; % atart timer

        % data rectification
        Hc1 = cdata;
        ep = exp(1i*(phin(:,1)));
        for i=1:freq_len
            Hc1(:,:,i) = cdata(:,:,i)./An(i,1)*ep(i);
        end

        % windowing
        wHc1 = Hc1;
        beta = 6;
        wind = kaiser(floor(freq2/freqstep)*2+1,beta);
        wind = wind(end-freq_len+1:end);
        for i=1:freq_len
            wHc1(:,:,i) = Hc1(:,:,i)*wind(i);
        end

%         trans_img = TransientwithLow(wHc1,freq1,freqstep,tau0,tau_step,tau_len);
        trans_img = ComputeTransient2(wHc1,freq1,freqstep,tau0,tau_step,tau_len);

        tElapsed = toc;
        fprintf('\nRunning Time: %f\n', tElapsed)

        % save(sprintf('trans_%s',sel_data),'trans_img','scenename',...
        %     'tau_step','tau_len','freq1','freq2','freqstep','freq_len');

        %% post processing
        trans_iron = IronTransient3(trans_img,1,2);
        % trans_iron = PulsizeTransient(trans_img,freq2,tau_step,0.1,0.2);
        % norm_img = NormalizeImage(trans_iron,2);
        norm_img = CompressHDR(trans_iron);
        % norm_img = CompressHDR(trans_img);

        %% show video
        i0 = floor(0/tau_step);
        i_len = floor(36/tau_step);
        for i=1:1:i_len
            figure(10)
            image(norm_img(:,:,i0+i)*256); colormap(gray(256))
            ttl = sprintf('scene frame: %.3d',i); title(ttl)
            pause(0.1*tau_step)
        end
        rmpath(needpath);             % 删除添加的路径；在运行他人的代码之前，要删除已添加的路径，以免函数重名，调用时出问题。    
        
        [y,x] = deal(20,40);  
        tr1 = permute(trans_img(y,x,:),[3 2 1]);
        tr2 = permute(trans_iron(y,x,:),[3 2 1]);
        tr3 = permute(norm_img(y,x,:),[3 2 1]);        
        figure;plot([tr1/max(tr1) tr2/max(tr2) tr3/max(tr3)]);
        legend('trans-img','trans-iron','norm-img');
        title(['(' num2str(y) ',' num2str(x) ')']);     
  case 'paper'
    needpath = genpath('./TransientImaging/Lin_Jingyu/Lin_code/transient_code_pub4_v04');    
    addpath (needpath);         % 添加路径
    
    % time step. Heide et al. set tau_step = 0.33.
    tau0 = 22.44;  % ns   % 20 对应的tau_step=0.1; 22.44对应的tau_step=0.33;
    tau_step = 0.33;   % 0.1;   源代码中是0.1
    tau_len = 1000;

    
    load('./TransientImaging/Lin_Jingyu/Lin_code/transient_code_pub4_v04/calib/cfs/PmCF_ubc.mat','An','phin','fre0','fre_step')  
    % CM, tau0, An, phin, fre0, fre_step, tau_step, Kw

    if size(M,1) < 442
        freq_len = num_frequencies;  % 压缩采样重构的频率个数用的是199个
    else
        freq_len = 210;       % 210是作者全采样使用的频率长度； 压缩感知中使用标定矩阵时，为199；
    end
    fre_len0 = floor(fre0/fre_step);
    % nnum = size(An,2);  % number of harmonic components
    freq1 = fre0;  % MHz
    freqstep = fre_step; 
    freq2 = fre0+(freq_len-1)*fre_step;

    measurements = measurements(:,:,1:freq_len,:);
    cdata = measurements(:,:,:,1) + 1i*measurements(:,:,:,2); % complex image
    sz = size(cdata);

    %% reconstruct
    fprintf('\nStart Running at %s.\n',datestr(now))
    tic; % atart timer

    % data rectification
    Hc1 = cdata;
    ep = exp(1i*(phin(:,1)));
    parfor i=1:freq_len
        Hc1(:,:,i) = cdata(:,:,i)./An(i,1)*ep(i);
    end

    % windowing
    wHc1 = Hc1;
    beta = 6;
    wind = kaiser(floor(freq2/freqstep)*2+1,beta);
    wind = wind(end-freq_len+1:end);
    parfor i=1:freq_len
        wHc1(:,:,i) = Hc1(:,:,i)*wind(i);
    end

    trans_img = ComputeTransient2(wHc1,freq1,freqstep,tau0,tau_step,tau_len);

    telapsed = toc;
    fprintf('\nRunning Time: %f\n', telapsed)

    rmpath(needpath);             % 删除添加的路径；在运行他人的代码之前，要删除已添加的路径，以免函数重名，调用时出问题。    

    % 计算观测值的距离

    I_rec = trans_img(:,:,1:size(M,2)); I_rec = max(I_rec,0);

%    TI_Postprocess(mM,I_rec,'_Lin',tElapsed);  
   
   global output_folder 
    if isfield(mM,'rec_measurements')
        aux = '_CS';
    elseif isfield(mM,'simu')
        aux = '_Poisson';
    elseif 442 == size(mM.M,1)
        aux = '';
    end
    algo = '_Lin';

    if isfield(mM,'simu')
        save( sprintf(['%s/' mM.filename algo aux '_QH.mat'], output_folder), 'I_rec','telapsed');         % 仿真数据
    else
        save( sprintf(['%s/' mM.filename algo aux '_QH.mat'], output_folder), 'I_rec','telapsed');
    end

end

%{
% My debug
if isfield(mM,'simu')
    norm_img = mat2gray(norm_img);
    simu = mat2gray(mM.simu);
    figure;
    for t = 1:size(simu,3)    
        subplot(1,2,1);imshow(norm_img(:,:,t));title(sprintf('Transient image frame %d', t));
        subplot(1,2,2);imshow(simu(:,:,t)); title(sprintf('simu image frame %d', t));       
        pause(0.1)
    end    
    [y,x] = deal(45,53);  
    tr1 = norm_img(y,x,:);tr1 = tr1(:);
    tr2 = simu(y,x,:);
    figure;plot([tr1(1:199)/6 tr2(:)]);legend(['Lin ',ver,' - rec'],'simu');title(['(' num2str(y) ',' num2str(x) ')']);   
    xlabel('Time Series');ylabel( 'Amplitude' );
    set(gca,'FontSize',14);
else
    [y,x] = deal(129,71);           
    tm1 = squeeze(measurements(x,y,:,1));
    tm2 = squeeze(measurements(x,y,:,2));
    figure;plot(tm1);hold on;
    figure;plot(tm2);hold on;
    tp1 = squeeze(trans_img(x,y,:));
    tp2 = squeeze(trans_iron(x,y,:));
    figure;plot([tp1 tp2]);grid on;

    time =[0:0.1:20];time=time(1:end-2);      % 时间采样点，去掉最后2个
    time = time/3*10;   
    figure;plot(time,tp1(1:length(time)),time,tp2(1:length(time)));grid on;  
end
%}

return

