function [trans_img] = ComputeTransient2(Hc,fre0,fre_step,tau0,tau_step,tau_len)
% ComputeTransient2 - Compute transient image by iFFT. 
% Parameters:
%   Hc - complex images captured with increasing frequencies
%   fre0,fre_step - frequency parameters (MHz)
%   tau0,tau_step,tau_len - time domain settings (ns)
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

% initializing
sz = size(Hc);
N = floor(1000/fre_step/tau_step); % number of FFT points
% k0 = 2*pi/1000;
tau_n0 = floor(tau0/tau_step); % number of points before tau0,
% tau = tau0+(0:tau_len-1)*tau_step; % ns
fre_n0 = floor(fre0/fre_step); % number of points before fre0,
% fre_len = sz(3);
% fre = fre0+(0:fre_len-1)*fre_step; % MHz

% computing beta
beta = zeros(sz(1),sz(2),tau_len);
lowZ = zeros(fre_n0,1); % zeros padding in low spectrum
for i=1:sz(1)
    for j=1:sz(2)
        dat = Hc(i,j,:);
        A = [lowZ;dat(:)];
        a = ifft(A,N);
        a = 2*real(a(tau_n0+(1:tau_len)));
        beta(i,j,:) = permute(a,[3 2 1]);
    end
end

% trans_img = beta;
trans_img = beta/tau_step; % normalize on tau_step
