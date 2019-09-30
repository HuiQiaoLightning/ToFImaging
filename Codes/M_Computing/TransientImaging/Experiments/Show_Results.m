%  显示成像结果
function Show_Results(X,varargin)

if X == 0
    for t = 1:size(I_rec,3)  
        I_show = I_rec(:,:,t);
        imwrite(I_show,['DiscoBall' num2str(t) '.bmp']);      
    end  
    
elseif isempty(varargin)    % 显示采集数据重构的图像
    scalefactor = 1;           % 3
    figure;
    for t = 1:size(X,3)  
        I_show = imresize( X(:,:,t), scalefactor,'bilinear');
%         I_show = medfilt2(imresize( X(:,:,t), scalefactor,'bilinear'));   % 对图像中值滤波
        imshow(I_show);
        title(sprintf('image frame %d', t));
        %title(sprintf('Transient image frame %d', t));
        pause(0.01)
%         pause;
    end   
    
elseif ischar(varargin{1}) && isequal(varargin{2}, 'avi')
    currpath = pwd;     % 获取当前路径
    cd('./TransientImaging/Outputs');
    v = VideoWriter(varargin{1},'Motion JPEG AVI');           % 使用 Motion JPEG 编码的 AVI 文件  
    v.Quality = 95;
    open(v)
    for n = 1:size(X,3)   
        writeVideo(v,X(:,:,n))
    end
    close(v)
    cd(currpath);       % 返回程序运行的当前路径
    
elseif ischar(varargin{1}) && isequal(varargin{2}, 'mp4')
    currpath = pwd;     % 获取当前路径
    global output_folder 
    cd(output_folder );
    v = VideoWriter(varargin{1},'MPEG-4');           % 使用 Motion JPEG 编码的mp4 文件      
    v.Quality = 95;
    open(v)
    for n = 1:size(X,3)   
        writeVideo(v,X(:,:,n))
    end
    close(v)
    cd(currpath);       % 返回程序运行的当前路径
end

return

%{
    case 2
        figure;             % Show transient image
        for t = 1:size(X,3)    
            subplot(1,2,1);imshow(X(:,:,t));title(sprintf('Transient image frame %d', t));
            subplot(1,2,2);imshow(simu(:,:,t)); title(sprintf('simu image frame %d', t));       
            pause(0.1);
        end
        % 未归一化时的显示
        [n,m] = deal(30,40);             % 分析具体的像素点：(60,70);    (30,40)       (149,88)                  
        rp = X(n,m,:);rp = rp(:);    
        sp = simu(n,m,:);sp = sp(:);
        figure;
        yyaxis left
        plot(rp);
        yyaxis right
        plot(sp);
        legend('reconstructed','original');title(['(' num2str(n) ',' num2str(m) ')']);    
        
        % 都归一化时的显示
        [n,m] = deal(30,40);             % 分析具体的像素点：(60,70);(149,88);(30,40)    
        rp = X(n,m,:);rp = rp(:);
        sp = simu(n,m,:);sp = sp(:);
        figure;plot([rp sp]);legend('reconstructed','original');title(['(' num2str(n) ',' num2str(m) ')']);
end


switch nargin
    case 1        
       %% 采用不同的显示方式，显示效果也不一样。
        X=varargin{1};
        prompt = '选择显示的方式: ';
        display('1.  用imshow显示图像，没有去掉小于零的部分，相当于全部都是交流分量');
        display('2.  用imshow显示图像，去掉小于零的部分，相当于人为构造了直流分量');
        display('3.  用image显示图像，没有去掉小于零的部分，相当于全部都是交流分量');
        display('4.  用image显示图像，去掉小于零的部分，相当于人为构造了直流分量');
        display('5.  直接显示，彩色图像；可以直观地观察亮度的分布');
        display('6.  显示某些像素');

        result = input(prompt,'s');

        switch result
            case '1'        
                % 用imshow显示图像，没有去掉小于零的部分，相当于全部都是交流分量，这个应该是对的，因为采集的时候，也没有交流分量。但是...
                scaling = 1.0; %Amplitude scale
                beta = GammaCorrection(X,scaling,2.2);            % 没有去掉小于零的部分

                scalefactor = 3;
                figure;
                for time = 1:size(beta,3)                
                    I_show = imresize( beta(end:-1:1,:,time), scalefactor,'bilinear');
                    imshow( I_show );
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end
            case '2'
                % 用imshow显示图像，去掉小于零的部分，相当于人为构造了直流分量，这个应该是不对的，因为采集的时候，没有直流分量。但是...
                scaling = 1.0; %Amplitude scale
                beta = TakePositiveAndGammaCorrection(X,scaling,2.2);           % 去掉小于零的部分

                scalefactor = 3;
                figure;
                for time = 1:size(beta,3)                
                    I_show = imresize( beta(end:-1:1,:,time), scalefactor,'bilinear');
                    imshow( I_show );
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end
            case '3'
                % 用image显示图像，没有去掉小于零的部分，相当于全部都是交流分量，这个应该是对的，因为采集的时候，也没有交流分量。但是...
                scaling = 1.0; %Amplitude scale
                beta = GammaCorrection(X,scaling,2.2)*255;  % 0~255
                figure;
                for time = 1:size(beta,3)        
                    image(beta(end:-1:1,:,time)); colormap(gray(256))
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end  
            case '4'
                % 用image显示图像，去掉小于零的部分，相当于人为构造了直流分量，这个应该是不对的，因为采集的时候，没有直流分量。但是...
                scaling = 1.0; %Amplitude scale
                beta = TakePositiveAndGammaCorrection(X,scaling,2.2)*255;  % 0~255
                figure;
                for time = 1:size(beta,3)        
                    image(beta(end:-1:1,:,time)); colormap(gray(256))
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end 
            case '5'
                % 直接显示，彩色图像；根据X的结果，可以直观地观察亮度的分布
                figure;
                for time = 1:size(X,3)        
                    image(X(end:-1:1,:,time),'CDataMapping','scaled'); colorbar     
                    title(['Frame ' num2str(time) '   Ex?']);
                    pause(0.1);
                end       	   
            case '6'
                % 显示某些像素
                tp1 = X(40,40,:); tp1 = tp1(:);        
                tp2 = X(60,60,:); tp2 = tp2(:);        
                tp3 = X(120,40,:); tp3 = tp3(:);           
                tp4 = X(120,60,:); tp4 = tp4(:);           
                figure;
                subplot(4,1,1); plot(tp1); title('X(40,40,:)');
                subplot(4,1,2); plot(tp2); title('X(60,60,:)');
                subplot(4,1,3); plot(tp3); title('X(120,40,:)');
                subplot(4,1,4); plot(tp4); title('X(120,60,:)');

        end  
    case 2        

        else

        end
end

% output_folder='../Outputs';      % 输出文件保存的位置   

if isequal(varargin{1},'w')         % 先写.bmp图形文件到硬盘，再做成视频文件；再删除.bmp文件 
    beta=varargin{2}/256;           % varargin{2}已被映射到[0,256]，imwrite中的范围是[0,1] 
    % 写成磁盘文件
    filename='Mirrors_Tran3_P';
    for i=1:size(beta,3)
        I_show = imresize(beta(end:-1:1,:,i), 3,'bilinear');
        imwrite(I_show,[filename num2str(i) '.bmp'],'bmp');
    end

    % 保存成avi格式视频
    vidObj = VideoWriter([filename '.avi']);    % 创建VideoWrite对象
    open(vidObj);                               % 打开该对象
    %创建Matlab电影
%     figure;
    for i=1:199
        img=imread([filename num2str(i) '.bmp']);
        imshow(img);
        currFrame = getframe;            % 复制当前图形
        writeVideo(vidObj,currFrame);    % 将当前图形写到vidObj对象中
    end
    close(vidObj);  
    !del *.bmp                          % 删除所有的.bmp文件
elseif nargin==3                        % 算法比较
    FX=varargin{1};LX=varargin{2};PX=varargin{3};
    figure;
    for time=1:size(FX,3)
        subplot(2,2,1);image(FX(end:-1:1,:,time)); colormap(gray(256));title(['Felix: ' num2str(time)]);
        subplot(2,2,2);image(LX(end:-1:1,:,time)); colormap(gray(256));title(['Lin: ' num2str(time)]);
        subplot(2,2,3);image(PX(end:-1:1,:,time)); colormap(gray(256));title(['Wavelet: ' num2str(time)]);
        pause(0.1);
    end
end
%}

    %% 第一种显示方式
%     beta=mat2gray(varargin{1}); 
%     beta=beta.^ (1/2.2);  %Gamma    % 对DiscoBall需要取3.2，视觉效果好
%     % 显示
%     figure;
%     for time=1:size(beta,3)
%         I_show = imresize(beta(end:-1:1,:,time), 3,'bilinear');
%         imshow(I_show);title(['Frame ' num2str(time)]);                
%         pause(1/25);
% %         pause;
%     end 

            i=100;j=30;
%             i=50;j=60;
%             i=70;j=80;
            fx=FX(i,j,:);
            lx=LX(i,j,:);
            px=PX(i,j,:);
            figure;plot([fx(:) lx(:) px(:)]);title(['At (' num2str(i) ',' num2str(j) ')']);
            legend('Felix','Lin','CS');
 
positionVector1 = [0.1, 0.2, 0.3, 0.3];
subplot('Position',positionVector1);
imshow(t(:,:,50));


