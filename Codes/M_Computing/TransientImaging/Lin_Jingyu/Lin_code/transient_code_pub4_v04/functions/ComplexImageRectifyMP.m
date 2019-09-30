function [Hc] = ComplexImageRectifyMP(dat,A1,phi1)
% ComplexImageRectifyMP  Convert multiphase dumpdata to Hc and rectify.
% Source dimensions: w,h,phase(0,90),freq,delay
% Output dimensions: h,w,freq,delay
% Parameters:
%   dat - 5D dumpdata
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
phs_len = sz(5);
% w = size(dat,1); 
% h = size(dat,2);
% freq_len = size(dat,4);
% phs_len = size(dat,5);
A0 = sum(dat,5)/phs_len; % DC component
Hc = zeros([h,w,freq_len,phs_len],'single');

% converting
if length(size(A1))==2 % 2D A1,phi1
    % dat = permute(dat(:,end:-1:1,:,:,:),[2 1 4 5 3]);
    ep = exp(-1i*(phi1(:,1)+phi1(:,2))/2);
    for i=1:freq_len
        for j=1:phs_len
            a = dat(:,:,:,i,j) - A0(:,:,:,i);
            Hc_re = a(:,end:-1:1,1).'/A1(i,1);
            Hc_im = a(:,end:-1:1,2).'/A1(i,2);
            Hc(:,:,i,j) = (Hc_re + 1i*Hc_im)*ep(i);
    %         Hc(:,:,i,j) = (Hc_re + 1i*Hc_im);
        end
    end
else  % 4D A1,phi1
    A1 = permute(A1,[3 4 2 1]); % w,h,phase,freq
    ep = exp(-1i*(phi1(:,1,:,:)+phi1(:,2,:,:))/2);
    ep = permute(ep,[3 4 2 1]);
    for i=1:freq_len
        for j=1:phs_len
            a = dat(:,:,:,i,j) - A0(:,:,:,i);
            Hc_re = a(:,:,1)./A1(:,:,1,i);
            Hc_im = a(:,:,2)./A1(:,:,2,i);
            a = (Hc_re + 1i*Hc_im).*ep(:,:,1,i);
            Hc(:,:,i,j) = a(:,end:-1:1).';
        end
    end
end
