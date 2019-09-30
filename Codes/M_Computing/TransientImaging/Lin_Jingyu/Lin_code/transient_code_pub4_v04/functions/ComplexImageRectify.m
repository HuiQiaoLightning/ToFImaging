function [Hc] = ComplexImageRectify(dat,A1,phi1)
% ComplexImageRectify  Convert dumpdata to Hc and rectify.
% Source dimensions: w,h,phase(0,90),freq
% Output dimensions: h,w,freq
% Parameters:
%   dat - 4D dumpdata
%   A1 - amplitude of fundamental component
%   phi1 - phase of fundamental component
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

% initializing
sz = size(dat);
w = sz(1); 
h = sz(2);
freq_len = sz(4);
% w = size(dat,1); 
% h = size(dat,2);
% freq_len = size(dat,4);
% phs_len = size(dat,5);
Hc = zeros([h,w,freq_len]);

% converting
if length(size(A1))==2 % 2D A1,phi1
    % dat = permute(dat(:,end:-1:1,:,:,:),[2 1 4 5 3]);
    ep = exp(-1i*(phi1(:,1)+phi1(:,2))/2);
    for i=1:freq_len
        Hc_re = dat(:,end:-1:1,1,i).'/A1(i,1);
        Hc_im = dat(:,end:-1:1,2,i).'/A1(i,2);
        Hc(:,:,i) = (Hc_re + 1i*Hc_im)*ep(i);
    %     Hc(:,:,i) = (Hc_re + 1i*Hc_im);
    end
else  % 4D A1,phi1
    A1 = permute(A1,[3 4 2 1]); % w,h,phase,freq
    dat = dat./A1;
    ep = exp(-1i*(phi1(:,1,:,:)+phi1(:,2,:,:))/2);
    ep = permute(ep,[3 4 2 1]);
    for i=1:freq_len
        a = (dat(:,:,1,i) + 1i*dat(:,:,2,i)).*ep(:,:,1,i);
        Hc(:,:,i) = a(:,end:-1:1).';
    end
end
