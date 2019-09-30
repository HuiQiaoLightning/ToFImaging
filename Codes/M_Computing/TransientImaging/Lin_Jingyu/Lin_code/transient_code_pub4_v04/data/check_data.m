% check data
addpath('..\functions')

sel_data = 'bottle16';
takeno = 1;
phs_len = 10;
freq_len = 161;
shuttersz = 3;

%% read data
rawpath = '..\data';
outputpath = sprintf('%s\\rawdata\\scene_%s',rawpath,sel_data);
% dataname = 'meas_data_take000.dat';
dataname = sprintf('meas_data_take%.3d.dat',takeno-1);
fn = sprintf('%s\\%s',outputpath,dataname);
sel_sh = shuttersz;
h = 120; w = 165;
fprintf('reading %s ...\n',outputpath);
dat = ReadDump_Shutter(fn,phs_len,freq_len,shuttersz,h,w,sel_sh);
% w, h, 2, freq_len, phs_len

%% plots
x = 80; y = 100;
d = permute( dat(x,y,1,50,:), [5 4 3 2 1]); % along phase
% d = permute( dat(x,y,1,:,1), [4 5 3 2 1]); % along frequency
figure(1); plot(d); legend('1','2','3','4','5','6')
