function [Hc] = ComplexImage(dat)
% ComplexImageCollection  Convert dumpdata to Hc.
% Source dimensions: w,h,phase(0,90),freq,delay
% Output dimensions: h,w,freq,delay
% Parameters:
%   dat - 5D dumpdata
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

dat = permute(dat,[2 1 4 5 3]);
Hc = dat(:,:,:,:,1) + 1i*dat(:,:,:,:,2);
Hc = Hc(end:-1:1,:,:,:);
