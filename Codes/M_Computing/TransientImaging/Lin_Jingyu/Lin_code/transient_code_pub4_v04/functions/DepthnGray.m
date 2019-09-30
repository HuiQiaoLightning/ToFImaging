function [depthmap,grayvalue] = DepthnGray(cap_img,fre0,fre_step,tau0,tau_step,tau_len)
% DepthnGray - Compute depth image and gray image from transient image by iFFT. 
% Parameters:
%   cap_img - complex images captured with increasing frequencies
%   fre0,fre_step - frequency parameters (MHz)
%   tau0,tau_step,tau_len - time domain settings (ns)
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

% initializing
sz = size(cap_img);
N = floor(1000/fre_step/tau_step); % number of FFT points
tau_n0 = floor(tau0/tau_step); % number of points before tau0,
fre_n0 = floor(fre0/fre_step); % number of points before fre0,

% computing depthmap
depthmap = zeros(sz(1),sz(2));
grayvalue = zeros(sz(1),sz(2));
lowZ = zeros(fre_n0,1); % zeros padding in low spectrum
for i=1:sz(1)
    for j=1:sz(2)
        dat = cap_img(i,j,:);
        A = [lowZ;dat(:)];
        a = ifft(A,N);
        a = 2*real(a(tau_n0+(1:tau_len)));
        [Y,I]=max(a);
        grayvalue(i,j) = Y/tau_step; % normalize on tau_step
        depthmap(i,j) = I+tau_n0;
    end
end
