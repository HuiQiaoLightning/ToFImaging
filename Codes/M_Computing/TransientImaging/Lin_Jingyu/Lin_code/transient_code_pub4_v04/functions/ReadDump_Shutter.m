function [dat] = ReadDump_Shutter(fn,delaysz,freqsz,shuttersz,h,w,shutno)
% ReadDump_Shutter  Read dumpdata at a shutter no, output 5D metrix 'dat'.
% Source dimensions: w,h,shutter,phase(0,90),freq,delay
% Output dimensions: w,h,phase(0,90),freq,delay
% Parameters:
%   fn - dump data name
%   delaysz - delay dimension
%   freqsz - frequency dimension
%   shuttersz - shutter dimension (120 ~ 1920)
%   h - image height (120)
%   w - image width (165)
%   shutno - selected shutter no.
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

% count = delaysz*freqsz*2*shuttersz*h*w;
count = shuttersz*h*w;
if shutno>shuttersz
    shutno = shuttersz;
end
dat = zeros([w,h,2,freqsz,delaysz]);
for i=1:delaysz
    for j=1:freqsz
        for k=1:2 % phases
            fdat = fread(fid, count, 'uint16');
            fdat = reshape(fdat,[w h shuttersz]);
            dat(:,:,k,j,i) = fdat(:,:,shutno);
        end
    end
end
fclose(fid);
dat = dat - 32768;  % ushort to short
