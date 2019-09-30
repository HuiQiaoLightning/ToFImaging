% invoke capturetool to capture images
% testing

%% capturetool
cmdpath = '..\TOF_tools';
captool = 'capturetool_v2.exe';

%% important parameters
datapath = 'testdata';
freq1 = 100; 
freqstep=1; 
freq2 = freq1;
phase1 = 0; 
phasepoint = 100;
phase2 = 1;  % period
phasestep = (phase2-phase1)/phasepoint;
phase2 = phase2 - phasestep/2;

otherstring = '120 2 1920 .\ com5 1'; % shutter and com port

%% system command
syscmd = sprintf('%s\\%s %s %f %f %f %f %f %f %s &', ...
    cmdpath,captool,datapath,...
    freq1,freqstep,freq2,phase1,phasestep,phase2,...
    otherstring);
system(syscmd)
