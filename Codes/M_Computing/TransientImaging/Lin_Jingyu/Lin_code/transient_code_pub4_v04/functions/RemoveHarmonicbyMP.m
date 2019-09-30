function [dat1] = RemoveHarmonicbyMP(dat,periods)
% RemoveHarmonicbyMP  Remove harmonic components by multiphase imaging.
% Source dimensions: w,h,phase(0,90),freq,delay
% Output dimensions: w,h,phase(0,90),freq
% Parameters:
%   dat - multiphase dumpdata (5D).
%         delay times must cover exact integer periods.
%   periods - period number covered by delay times.
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

% initializing
sz = size(dat);
if length(sz)<5 % not multiphase
    dat1 = dat;
    return;
end
w = sz(1); 
h = sz(2);
freq_len = sz(4);
K = sz(5); % phase number
if ~exist('periods','var')
    periods = 1;
end

% integration
cs = cos(2*pi*periods*(0:K-1)/K).';
dat1 = reshape(dat,[w*h*freq_len*2,K])*cs;
dat1 = reshape(dat1,[w,h,2,freq_len])*2/K;

