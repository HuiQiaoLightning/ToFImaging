function [trans_img] = TransientwithLow_90(Hc,fre0,fre_step,tau0,tau_step,tau_len)
% TransientwithLow - Compute transient image and recover low spectrum.
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
N1 = floor(1000/fre0/tau_step); % number of FFT points
N2 = floor(1000/10/tau_step); % number of FFT points
tau_n0 = floor(tau0/tau_step); % number of points before tau0,
% tau = tau0+(0:tau_len-1)*tau_step; % ns
fre_n0 = floor(fre0/fre_step); % number of points before fre0,
% fre_len = sz(3);
% fre = fre0+(0:fre_len-1)*fre_step; % MHz

% computing transient image
trans_img = zeros(sz(1),sz(2),tau_len);
for i=1:sz(1)
    for j=1:sz(2)
        % iFFT with frequency step fre0
        pix = Hc(i,j,:);
        pix = pix(:);
        A = [abs(pix(1));pix(1:fre_n0:90)]; % pad a zero and downsample
        a = ifft(A,N1);
        a = 2*real(a); % length N1
        B = [0;pix((10-fre0+1):10:end)]; % pad a zero and downsample
        b = ifft(B,N2);
        b = 2*real(b); % length N2
        % FFT with frequency step fre_step
        Hc_low = fft(a-min(a),N);
        Hc_low = Hc_low(1:fre_n0);
        Hc_90 = fft(b-min(b),N);
        % iFFT with frequency step fre_step
        A = [Hc_low;pix]; % padding recovered low sprectrum
        A(91:100) = Hc_90(91:100);
        a = ifft(A,N);
        a = 2*real(a(tau_n0+(1:tau_len)));
        trans_img(i,j,:) = permute(a,[3 2 1]);
    end
end
trans_img = trans_img/tau_step; % normalize on tau_step
