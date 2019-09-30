function Project_5(mM)
% ��һ��paper�е�Heide����Ͳɼ����ݵ�ʵ�飺tomato����ÿ�����طֱ���һ����������20���������ع�����3��4��5���ò�ź���ѡȡ���������������㣬���п����ҵ�����Ϻõ��Ǹ��źš�

xs = mM.imagedims(1);
ys = mM.imagedims(2);
M = mM.M;

if isfield(mM,'simu')   
    if mM.noise_measurements == 0
        disp('����Ĺ۲�ֵ��Ҫ������');
        return;
    else
        measurements = mM.noise_measurements;
        simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
        weight = 3;       % ���Ȩ���ǵ��������ݲ��üӹ̶���Poisson����ʱ���õģ���һ����Ҫ���ݾ������ݽ��й��Ƶ�ֵ
    end
else
    measurements = mM.measurements;
    weight = 1;
end

measurements = (reshape(measurements,[xs*ys,size(M,1)]))'; % ����׼����
 
waveletname = 'db5';    %'sym5'
N = 5; 
num_noise = 100;   % 70

sp = 0;  % ������
tstart = tic;
parfor xy = 1:size(measurements,2)
    tmp = measurements(:,xy);    

%     sp = simu(:,xy);     % ������    
    
    tmp = rec_pixel_2(M, weight, waveletname, N, num_noise, tmp, sp);    
    
    I_rec(:,xy) = tmp;
end
telapsed = toc(tstart);

%{
% ����۲�ֵ�ľ���
if isfield(mM,'simu')
    simu = (reshape(mM.simu,[xs*ys,size(M,2)]))';
    mdist = norm(simu - I_rec)/norm(simu);
    
    simu = mM.simu;
    I_rec = reshape(I_rec',[ys,xs,size(M,2)]);
    [y,x] = deal(86,96);             % ������������ص㣺��Ƶ��(136,83);(149,88);(150,91);(148,95);(120,140);���壺(100,40) (70,60);(31,86);(92,57);����2��return��(40,30);(53,45);˫�壺(86,96);
    rp = squeeze(I_rec(x,y,:));
    sp = squeeze(simu(x,y,:));
    figure;plot([sp rp]);legend('True','Rand wavelet');title([mM.filename ' signal at (', num2str(y) ',' num2str(x),')']);    
    xlabel('Time Series');ylabel( 'Amplitude' );      
else    
    % ����۲�ֵ�ľ���
    rec_m = M * I_rec;        % ��     �������󱣴����������ټ���
    mdist = norm(measurements - rec_m) / norm(measurements);   % ����ŷʽ���룬�ܹ۲����ֵӰ��ϴ�

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
end

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
%}

aux = '_Poisson';

algo = '_wavelet_Project_5';

% Show_Results(PX2,[mM.filename algo aux '_' num2str(num_noise)],'mp4');    % ������Ƶ�ļ�


global output_folder 

if isfield(mM,'simu')
    save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_QH.mat'], output_folder), 'I_rec', 'num_noise', 'telapsed');
else
    save( sprintf(['%s/' mM.filename algo aux '_' num2str(num_noise) '_QH.mat'], output_folder), 'I_rec', 'num_noise', 'telapsed');    
end

%{
% ����ΪͼƬ
currpath = pwd;     % ��ȡ��ǰ·��
cd(output_folder );
for t = 1:size(PX2,3)  
    I_show = PX2(:,:,t);
    imwrite(I_show,[mM.filename algo aux '_' num2str(num_noise) '_' num2str(t) '.bmp']);      
end 

cd(currpath);       % ���س������еĵ�ǰ·��

% system('shutdown -s');
%}
return

% ������ؼ�������󣬰汾2��Ҳ��ʵ�ְ汾1�Ĺ��ܡ�Ŀǰ��ֻ���ǵ�4��5���ò�źţ���С���ֽ������С���񵴿�ͨ����ȡ����0�Ĳ�����ȥ����
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
    level_start = 3;     % �ӵ�3���ò�źſ�ʼ������Ǵӵ�4���ò�źſ�ʼ�Ļ������� 4
    level_step = level_start - 1;
    for level = level_start:N        
        tmp = wrcoef('a',Cb,Lb,waveletname,level);   % ����ϵ��2��������ѵģ��ɸ���ϸ�ڵ����ֵ���ֵĸ���ȷ�����ܵķ�ֵ�㡣
        [tmpv,tmpi] = max(tmp);
        if tmpi <= 3 || tmpi >= numel(tmp)-2     % �����ֵ�㿿�������˵㣬����Ϊ����С��������ģ������ǡ�       
            nn = nn - 1;        % for ѭ������nn+1�����nn��Ӧ�õ���
            break;
        end
        % tmp,  max( tmp,0 ) �۲�һ��ȡ0�Ͳ�ȡ0�����𣬰�ȡƽ��ֵ�������������ȡ0Ӧ��Ҳ�ǿ��Եġ���ʵ�ǲ�Ӧ����ȡ����0�Ĳ����ġ���ȡ����0���ܼ����������������
        if level == 3
            ap{ nn,level - level_step } = Find_Main_Peak_without_Thresh(tmp);  % ���ڵ�3���ò�źţ�������ȡһ�����壬����ֵ�޹ء� ��Ϊ��3����������Ĳ���϶ࡣ             
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
        tmp1 = tmp1 - mean(tmp1);        % ȥ��ֵ���������ǲ���һ��ˮƽ����
%         figure;plot([tmp tmp1]);
        similarity( nn,level ) = norm(tmp - tmp1);
%         similarity( nn,level ) = std( tmp - tmp1 );      % �˴����Կ���ʹ�� modwtxcorr       
    end
    
%     figure;plot( ap{ nn,end } );grid on;
    if size( ap,2 ) == 3
        %{ 
        if length( find( ap{ nn,3 }>0 ) ) > 150         % ���ڵ�Ƶ�źţ��ڵ�4��5���ò�ź���ѡ����ֵ��150
            similarity( nn,1 ) = Inf;
        end   
        %}
        %{
        if length( find( ap{ nn,3 }>0 ) ) > 150          % ���ڵ�Ƶ�źţ��ڵ�5���ò�ź���ѡ
            similarity( nn,2 ) = Inf;
            similarity( nn,1 ) = Inf;
        elseif length( find( ap{ nn,3 }>0 ) ) > 120     % �����ж����ֵ�ĵ�Ƶ�źţ���4��5���ò�ź���ѡ
            similarity( nn,1 ) = Inf;
        end
        %}
         
        if length( find( ap{ nn,3 }>0 ) ) > 120         % % ���ڵ�Ƶ�źţ��ڵ�4��5���ò�ź���ѡ����ֵ��120�����Ŀǰ��Ҫ�����źŽ��е�����Ӧ��������Ӧ
            similarity( nn,1 ) = Inf;
        end   
                
    elseif size( ap,2 ) == 2                             % ֻ�õ�4��5���ò�ź� 
        if length( find( ap{ nn,2 }>0 ) ) > 150          % ���ڵ�Ƶ�źţ��ڵ�5���ò�ź���ѡ
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
