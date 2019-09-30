% depth estimation

addpath('..\functions')

sel_data = 'mirror1'; % data name
cf_name = 'PmCF5_stripe'; % correlation function name
rawpath = '..\data'; % raw data path

% time step (ns)
tau0 = 30;
tau_step = 0.1;
tau_len = 1000;
% tau0 = 28;
% tau_step = 0.00667;
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
shuttersz = 1;

freq = freq1:freqstep:freq2;
freq_len = length(freq);
phsstep = (phs2-phs1)/phs_len;
phs = phs1 + (0:phs_len-1)*phsstep;

sel_sh = shuttersz;

takes = 1; % number of taking

%% load correlation function
% if  ~exist('phi1_pix','var')
%     load(sprintf('..\\calib\\cfs\\%s',cf_name),'A1_pix','phi1_pix');
if  ~exist('phi1_str','var')
    load(sprintf('..\\calib\\cfs\\%s',cf_name),'A1_str','phi1_str');
end
%   A1, phi1: Amplitudes and phases of fundamental component
%             dimensions: frequency,0/90deg,w,h,shutter.
%   shutter£º shutter time.
%   freq1, freq2, freqstep, freq_len, freq: frequency setting.

if shuttersz==3 % same as the calibration data
    s = sel_sh;
else
    s = 3;
end
A1 = A1_str(freq1:freq2,:,:,:,s);
phi1 = phi1_str(freq1:freq2,:,:,:,s);

%% read data
outputpath = sprintf('%s\\rawdata\\scene_%s',rawpath,sel_data);
dataname = 'meas_data_take000.dat';
fn = sprintf('%s\\%s',outputpath,dataname);
fprintf('reading data ...\n');
dat = ReadDump_Shutter(fn,phs_len,freq_len,shuttersz,h,w,sel_sh);

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

% trans_img = TransientwithLow(wHc1,freq1,freqstep,tau0,tau_step,tau_len);
% trans_img = ComputeTransient2(wHc1,freq1,freqstep,tau0,tau_step,tau_len);
[depthmap,grayvalue] = DepthnGray(wHc1,freq1,freqstep,tau0,tau_step,tau_len);
depthmap = depthmap(:,1:163);
grayvalue = grayvalue(:,1:163);

% trans_iron = IronTransient3(trans_img,1,6);

tElapsed = toc;
fprintf('\nRunning Time: %f\n', tElapsed)

%% plots
x = 80;
y = 60;

a = permute(dat(x,y,:,1:freq_len,1),[4 3 2 1 5]);
b = permute(dat1(x,y,:,:),[4 3 2 1]);
c = permute(Hc1(y,x,:),[3 2 1]);

figure(1);
plot([a b real(c) imag(c)])
legend('dat_r','dat_i','dat1_r','dat1_i','Hc_{1r}','Hc_{1i}')

% return

%% depthmap
% [grayvalue,depthmap] = max(trans_img(:,1:163,:),[],3);
% [grayvalue,depthmap] = max(trans_img(:,1:163,200:300),[],3);
% grayvalue = grayvalue/max(grayvalue(:));
figure(1); image(grayvalue/max(grayvalue(:))*256); 
colormap(gray(256));colorbar

dep1 = (tau0+0)/tau_step; dep2 = dep1+10/tau_step;
objdepth = depthmap;
objdepth = min(objdepth,dep2);
objdepth = max(objdepth,dep1);

thr = max(grayvalue(:))*0.004;
mask = grayvalue>thr;
objdepth = objdepth.*mask + max(objdepth(:))*(~mask);
figure(2); imagesc(objdepth); 
colormap(jet(256));colorbar

fprintf('depthmap range %d ~ %d.\n',...
    min(depthmap(:)),max(depthmap(:)));

% SPEEDOFLIGHT = 299792458e-9; % mm/ps
% dep = objdepth*SPEEDOFLIGHT;
% figure(3); imagesc(dep); 
% colormap(jet(256));colorbar

% return

%% save
fn = sprintf('results\\depth_%s',sel_data);
save(fn,'depthmap','objdepth','grayvalue','tau_step');
