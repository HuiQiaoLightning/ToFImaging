function [dat] = ReadDump_all(fn,delaysz,freqsz,shuttersz,h,w)
% ReadDump_all  Read all dumpdata, output 6D metrix 'dat'.
% Source dimensions: w,h,shutter,phase(0,90),freq,delay
% Output dimensions: w,h,shutter,phase(0,90),freq,delay
% Parameters:
%   fn - dump data name
%   delaysz - delay dimension
%   freqsz - frequency dimension
%   shuttersz - shutter dimension (120 ~ 1920)
%   h - image height (120)
%   w - image width (165)
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

[fid, msg] = fopen(fn,'r');
if fid==-1
    fprintf(msg);
    fprintf('\n');
    dat = 0;
    return;
end

count = delaysz*freqsz*2*shuttersz*h*w;
fdat = fread(fid, count, 'uint16');
fclose(fid);

dat = reshape(fdat,[w,h,shuttersz,2,freqsz,delaysz]);
dat = dat - 32768;  % ushort to short
