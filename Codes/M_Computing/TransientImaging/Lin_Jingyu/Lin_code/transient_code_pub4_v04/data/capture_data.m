% invoke capturetool to capture images
% multiple frequencies, multiple fixed phase points

sel_data = 'mirror';

%% capturetool
cmdpath = '..\TOF_tools';
captool = 'capturetool_v2.exe';

%% acquaring parameters (important!)
freq1 = 1;  % MHz
freq2 = 180;
freqstep = 1; 

phase1 = 0; 
phasepoint = 10;
phase2 = 1;  % period

shutter = '480 2 1920'; % shutter time (us)
otherstring = 'com5 1'; % com port and number of taking

%% default settings
phasestep = (phase2-phase1)/phasepoint;
phase2 = phase2 - phasestep/2;
% SPEEDOFLIGHT = 299792458; % m/s
% delaystep = phasestep/freq1*1e-6*SPEEDOFLIGHT; % m
% delay1 = 0; 
% delay2 = delay1+delaystep*phasepoint;

outputpath = sprintf('rawdata\\scene_%s',sel_data);
datapath = '.\';

%% system commandy
syscmd = sprintf('%s\\%s %s %f %f %f %f %f %f %s %s %s', ...
    cmdpath,captool,datapath,...
    freq1,freqstep,freq2,phase1,phasestep,phase2,...
    shutter, outputpath, otherstring);
system(syscmd);
