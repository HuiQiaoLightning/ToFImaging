% reconstruct transient images from UBC data
% 
% Download the data package from the website of this paper
% F. Heide, M. Hullin, J. Gregson, W. Heidrich. Low-budget Transient Imaging using Photonic Mixer Devices. SIGGRAPH 2013.
% http://www.cs.ubc.ca/labs/imager/tr/2013/TransientPMD/
% 
% Then copy the content in the 'data/images' in the package to 'data_ubc'.
%
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

disp('start v04');               % Ìí¼Ó by whm

addpath('..\functions')           % ×¢ÊÍ by whm

% time step. Heide et al. set tau_step = 0.33.
tau0 = 20;  % ns
tau_step = 0.1;
tau_len = 1000;
% tau0 = 20;
% tau_step = 0.01;
% tau_len = 3600;

%% load data
sel_data = 1; % 1.Mirrors, 2.Wall, 3.Discoball, 4.Corner, 5.Bottles.
sel_shutter = -1; % use default valuemax(trans_img,[],3)
load_img % load data

% mx = max(abs(measurements(:)));
% measurements = measurements/mx; % normalizing
measurements = permute(measurements,[2 1 3 4]);
measurements = double(measurements(end:-1:1,:,:,:));

%% initializing
cf_name = 'PmCF_ubc';

load(sprintf('..\\calib\\cfs\\%s',cf_name),'An','phin','fre0','fre_step')     

% CM, tau0, An, phin, fre0, fre_step, tau_step, Kw

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

%% show video
i0 = floor(0/tau_step);
i_len = floor(36/tau_step);
for i=1:1:i_len
    figure(10)
    image(norm_img(:,:,i0+i)*256); colormap(gray(256))
    ttl = sprintf('scene frame: %.3d',i); title(ttl)
    pause(0.1*tau_step)
end

disp('end v04');               % Ìí¼Ó by whm