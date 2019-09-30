% reconstruct transient

addpath('..\functions')

sel_data = 'bottle16'; % data name
cf_name = 'PmCF1_pix'; % correlation function name
rawpath = '..\data'; % raw data path

% time step (ns)
tau0 = 33.3;
tau_step = 0.1;
tau_len = 1000;
% tau0 = 33.3;
% tau_step = 0.01;
% tau_len = 3000;

%% acquaring parameters
freq1 = 5;  % MHz
freq2 = 165;
freqstep = 1; 

phs1 = 0;
phs2 = 1; % period
phs_len = 10;

w = 165;
h = 120;
shuttersz = 3;

freq = freq1:freqstep:freq2;
freq_len = length(freq);
phsstep = (phs2-phs1)/phs_len;
phs = phs1 + (0:phs_len-1)*phsstep;

sel_sh = shuttersz;

takes = 1; % number of taking

%% load correlation function
if  ~exist('phi1_pix','var')
    load(sprintf('..\\calib\\cfs\\%s',cf_name),'A1_pix','phi1_pix');
end
%   A1, phi1: Amplitudes and phases of fundamental component
%             dimensions: frequency,0/90deg,w,h,shutter.
%   shutter£º shutter time.
%   freq1, freq2, freqstep, freq_len, freq: frequency setting.

A1 = A1_pix(1:freq_len,:,:,:,sel_sh);
phi1 = phi1_pix(1:freq_len,:,:,:,sel_sh);

%% read data
dat = zeros(w,h,2,freq_len,phs_len);
outputpath = sprintf('%s\\rawdata\\scene_%s',rawpath,sel_data);
for i=1:takes
    dataname = sprintf('meas_data_take%.3d.dat',i-1);
    fn = sprintf('%s\\%s',outputpath,dataname);
    fprintf('reading %s ...\n',dataname);
    dat = dat + ReadDump_Shutter(fn,phs_len,freq_len,shuttersz,h,w,sel_sh);
end

% check
showdata = false; %%%% switch %%%%
if showdata
    % show captured images
    shw = permute(dat,[2 1 4 5 3]);
    shw = shw(end:-1:1,:,:,:,:)/max(shw(:));
    for i=1:freq_len
        figure(1); 
        image([shw(:,:,i,1,1) shw(:,:,i,1,2)]*256*100); 
        colormap(gray(256))
        ttl = sprintf('frequency: %d',i); title(ttl) 
        pause(0.1)
    end
end % if showdata

%% reconstruct
fprintf('\nStart Running at %s.\n',datestr(now))
tic; % atart timer

dat1 = RemoveHarmonicbyMP(dat(:,:,:,1:freq_len,:),1);

Hc1 =  ComplexImageRectify(dat1,A1,phi1);
wHc1 = Hc1;

beta = 6;
wind = kaiser(floor(freq2/freqstep)*2+1,beta);
wind = wind(end-freq_len+1:end);
for i=1:freq_len
    wHc1(:,:,i) = Hc1(:,:,i)*wind(i);
end

trans_img = TransientwithLow(wHc1,freq1,freqstep,tau0,tau_step,tau_len);
% trans_img = ComputeTransient2(Hc1,freq1,freqstep,tau0,tau_step,tau_len);

tElapsed = toc;
fprintf('\nRunning Time: %f\n', tElapsed)

% save(sprintf('results\\trans_%s',sel_data),'trans_img',...
%     'tau_step','tau_len','freq','freqstep','freq_len');

%% plots
x = 100;
y = 80;
% x = 80;
% y =60;

a = permute(dat(x,y,:,1:freq_len,1),[4 3 2 1 5]);
b = permute(dat1(x,y,:,:),[4 3 2 1]);
c = permute(Hc1(y,x,:),[3 2 1]);
tr1 = permute(trans_img(y,x,:),[3 2 1]);

figure(1);
plot([a b real(c) imag(c)])
legend('dat_r','dat_i','dat1_r','dat1_i','Hc_{1r}','Hc_{1i}')

figure(2);
plot(tr1)

%% show video
trans_iron = IronTransient3(trans_img,1,10);
% norm_img = NormalizeImage(trans_iron,2);
norm_img = CompressHDR(trans_iron);

i0 = 0;
i_len = floor(6/tau_step);
for i=1:i_len
    figure(10)
%     image(trans_iron(:,:,i)*256); colormap(gray(256))
%     image(trans_img(:,:,i)*256); colormap(gray(256))
%     image(ourFrame(:,:,i)*256); colormap(gray(256))
    image(norm_img(:,:,i0+i)*256); colormap(gray(256))
    ttl = sprintf('scene %s: %.3d',sel_data,i); title(ttl)
    pause(0.5*tau_step)
end
