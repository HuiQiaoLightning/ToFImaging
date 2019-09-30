function [Hc1] = RemoveHarmonicbyMP(Hc,periods)
% RemoveHarmonicbyMP  Remove harmonic components by multiphase imaging.
% Source dimensions: h,w,freq,delay. 
% Output dimensions: h,w,freq
% Parameters:
%   Hc - multiphase complex images (4D)
%        delay times must cover exact integer periods.
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

%% initializing
h = size(Hc,1);
w = size(Hc,2);
fre_len = size(Hc,3); % phase number
K = size(Hc,4); % phase number
if ~exist('periods','var')
    periods = 1;
end

%% integration
cs = cos(2*pi*periods*(0:K-1)/K);
dat = permute(Hc,[4 1 2 3]);
Hc1 = cs*reshape(dat,[K,h*w*fre_len]);
Hc1 = reshape(Hc1,[h,w,fre_len])*2/K;