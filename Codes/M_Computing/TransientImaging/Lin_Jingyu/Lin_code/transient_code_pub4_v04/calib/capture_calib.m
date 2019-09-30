% invoke capturetool to capture images
% multiple frequencies, multiple fixed phase points

%% capturetool
cmdpath = '..\TOF_tools';
captool = 'capturetool_v2.exe';

%% acquaring parameters (important!)
freq1 = 1;  % MHz
freq2 = 180;
freqstep = 1; 

phase1 = 0; 
phasepoint = 160;
phase2 = 2;  % period

shutter = '480 2 1920';  % shutter time (us)
otherstring = 'com3 8'; % com port and number of taking
outputpath = 'calib_5';

%% default settings
phasestep = (phase2-phase1)/phasepoint;
phase2 = phase2 - phasestep/2;
% SPEEDOFLIGHT = 299792458; % m/s
% delaystep = phasestep/freq1*1e-6*SPEEDOFLIGHT; % m
% delay1 = 0; 
% delay2 = delay1+delaystep*phasepoint;

%% system command
% syscmd = sprintf('%s\\%s %s %f %f %f %f %f %f %s %s %s', ...
%     cmdpath,captool,datapath,...
%     freq1,freqstep,freq2,phase1,phasestep,phase2,...
%     shutter, outputpath, otherstring);
% system(syscmd);
for f = freq1:freqstep:freq2
    datapath = sprintf('calib_freq_%05.1f',f);
    syscmd = sprintf('%s\\%s %s %f %f %f %f %f %f %s %s %s', ...
        cmdpath,captool,datapath,...
        f,freqstep,f,phase1,phasestep,phase2,...
        shutter, outputpath, otherstring);
    system(syscmd);
end
