%  ��ʾ������
function Show_Results(X,varargin)

if X == 0
    for t = 1:size(I_rec,3)  
        I_show = I_rec(:,:,t);
        imwrite(I_show,['DiscoBall' num2str(t) '.bmp']);      
    end  
    
elseif isempty(varargin)    % ��ʾ�ɼ������ع���ͼ��
    scalefactor = 1;           % 3
    figure;
    for t = 1:size(X,3)  
        I_show = imresize( X(:,:,t), scalefactor,'bilinear');
%         I_show = medfilt2(imresize( X(:,:,t), scalefactor,'bilinear'));   % ��ͼ����ֵ�˲�
        imshow(I_show);
        title(sprintf('image frame %d', t));
        %title(sprintf('Transient image frame %d', t));
        pause(0.01)
%         pause;
    end   
    
elseif ischar(varargin{1}) && isequal(varargin{2}, 'avi')
    currpath = pwd;     % ��ȡ��ǰ·��
    cd('./TransientImaging/Outputs');
    v = VideoWriter(varargin{1},'Motion JPEG AVI');           % ʹ�� Motion JPEG ����� AVI �ļ�  
    v.Quality = 95;
    open(v)
    for n = 1:size(X,3)   
        writeVideo(v,X(:,:,n))
    end
    close(v)
    cd(currpath);       % ���س������еĵ�ǰ·��
    
elseif ischar(varargin{1}) && isequal(varargin{2}, 'mp4')
    currpath = pwd;     % ��ȡ��ǰ·��
    global output_folder 
    cd(output_folder );
    v = VideoWriter(varargin{1},'MPEG-4');           % ʹ�� Motion JPEG �����mp4 �ļ�      
    v.Quality = 95;
    open(v)
    for n = 1:size(X,3)   
        writeVideo(v,X(:,:,n))
    end
    close(v)
    cd(currpath);       % ���س������еĵ�ǰ·��
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
        % δ��һ��ʱ����ʾ
        [n,m] = deal(30,40);             % ������������ص㣺(60,70);    (30,40)       (149,88)                  
        rp = X(n,m,:);rp = rp(:);    
        sp = simu(n,m,:);sp = sp(:);
        figure;
        yyaxis left
        plot(rp);
        yyaxis right
        plot(sp);
        legend('reconstructed','original');title(['(' num2str(n) ',' num2str(m) ')']);    
        
        % ����һ��ʱ����ʾ
        [n,m] = deal(30,40);             % ������������ص㣺(60,70);(149,88);(30,40)    
        rp = X(n,m,:);rp = rp(:);
        sp = simu(n,m,:);sp = sp(:);
        figure;plot([rp sp]);legend('reconstructed','original');title(['(' num2str(n) ',' num2str(m) ')']);
end


switch nargin
    case 1        
       %% ���ò�ͬ����ʾ��ʽ����ʾЧ��Ҳ��һ����
        X=varargin{1};
        prompt = 'ѡ����ʾ�ķ�ʽ: ';
        display('1.  ��imshow��ʾͼ��û��ȥ��С����Ĳ��֣��൱��ȫ�����ǽ�������');
        display('2.  ��imshow��ʾͼ��ȥ��С����Ĳ��֣��൱����Ϊ������ֱ������');
        display('3.  ��image��ʾͼ��û��ȥ��С����Ĳ��֣��൱��ȫ�����ǽ�������');
        display('4.  ��image��ʾͼ��ȥ��С����Ĳ��֣��൱����Ϊ������ֱ������');
        display('5.  ֱ����ʾ����ɫͼ�񣻿���ֱ�۵ع۲����ȵķֲ�');
        display('6.  ��ʾĳЩ����');

        result = input(prompt,'s');

        switch result
            case '1'        
                % ��imshow��ʾͼ��û��ȥ��С����Ĳ��֣��൱��ȫ�����ǽ������������Ӧ���ǶԵģ���Ϊ�ɼ���ʱ��Ҳû�н�������������...
                scaling = 1.0; %Amplitude scale
                beta = GammaCorrection(X,scaling,2.2);            % û��ȥ��С����Ĳ���

                scalefactor = 3;
                figure;
                for time = 1:size(beta,3)                
                    I_show = imresize( beta(end:-1:1,:,time), scalefactor,'bilinear');
                    imshow( I_show );
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end
            case '2'
                % ��imshow��ʾͼ��ȥ��С����Ĳ��֣��൱����Ϊ������ֱ�����������Ӧ���ǲ��Եģ���Ϊ�ɼ���ʱ��û��ֱ������������...
                scaling = 1.0; %Amplitude scale
                beta = TakePositiveAndGammaCorrection(X,scaling,2.2);           % ȥ��С����Ĳ���

                scalefactor = 3;
                figure;
                for time = 1:size(beta,3)                
                    I_show = imresize( beta(end:-1:1,:,time), scalefactor,'bilinear');
                    imshow( I_show );
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end
            case '3'
                % ��image��ʾͼ��û��ȥ��С����Ĳ��֣��൱��ȫ�����ǽ������������Ӧ���ǶԵģ���Ϊ�ɼ���ʱ��Ҳû�н�������������...
                scaling = 1.0; %Amplitude scale
                beta = GammaCorrection(X,scaling,2.2)*255;  % 0~255
                figure;
                for time = 1:size(beta,3)        
                    image(beta(end:-1:1,:,time)); colormap(gray(256))
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end  
            case '4'
                % ��image��ʾͼ��ȥ��С����Ĳ��֣��൱����Ϊ������ֱ�����������Ӧ���ǲ��Եģ���Ϊ�ɼ���ʱ��û��ֱ������������...
                scaling = 1.0; %Amplitude scale
                beta = TakePositiveAndGammaCorrection(X,scaling,2.2)*255;  % 0~255
                figure;
                for time = 1:size(beta,3)        
                    image(beta(end:-1:1,:,time)); colormap(gray(256))
                    title(['Frame ' num2str(time)]);
                    pause(0.1);
                end 
            case '5'
                % ֱ����ʾ����ɫͼ�񣻸���X�Ľ��������ֱ�۵ع۲����ȵķֲ�
                figure;
                for time = 1:size(X,3)        
                    image(X(end:-1:1,:,time),'CDataMapping','scaled'); colorbar     
                    title(['Frame ' num2str(time) '   Ex?']);
                    pause(0.1);
                end       	   
            case '6'
                % ��ʾĳЩ����
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

% output_folder='../Outputs';      % ����ļ������λ��   

if isequal(varargin{1},'w')         % ��д.bmpͼ���ļ���Ӳ�̣���������Ƶ�ļ�����ɾ��.bmp�ļ� 
    beta=varargin{2}/256;           % varargin{2}�ѱ�ӳ�䵽[0,256]��imwrite�еķ�Χ��[0,1] 
    % д�ɴ����ļ�
    filename='Mirrors_Tran3_P';
    for i=1:size(beta,3)
        I_show = imresize(beta(end:-1:1,:,i), 3,'bilinear');
        imwrite(I_show,[filename num2str(i) '.bmp'],'bmp');
    end

    % �����avi��ʽ��Ƶ
    vidObj = VideoWriter([filename '.avi']);    % ����VideoWrite����
    open(vidObj);                               % �򿪸ö���
    %����Matlab��Ӱ
%     figure;
    for i=1:199
        img=imread([filename num2str(i) '.bmp']);
        imshow(img);
        currFrame = getframe;            % ���Ƶ�ǰͼ��
        writeVideo(vidObj,currFrame);    % ����ǰͼ��д��vidObj������
    end
    close(vidObj);  
    !del *.bmp                          % ɾ�����е�.bmp�ļ�
elseif nargin==3                        % �㷨�Ƚ�
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

    %% ��һ����ʾ��ʽ
%     beta=mat2gray(varargin{1}); 
%     beta=beta.^ (1/2.2);  %Gamma    % ��DiscoBall��Ҫȡ3.2���Ӿ�Ч����
%     % ��ʾ
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


