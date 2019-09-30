function [beta1] = ExtractFundamental(beta,fre0,fre_step,tau0,tau_step,An,phin)
% ExtractFundamental  Extract fundamental component / remove dilated components 
% Parameters:
%   beta - complex images captured with increasing frequencies
%   fre0,fre_step - frequency parameters (MHz)
%   tau0,tau_step - time domain settings (ns)
%   An, phin - correlation matrix parameters
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

% initializing
N = floor(1000/fre_step/tau_step); % number of FFT points
sz = size(beta);
tau_len = sz(3);
tau_n0 = floor(tau0/tau_step); % number of points before tau0
fre_n0 = floor(fre0/fre_step); % number of points before fre0
fre_len = size(An,1);

nnum = size(An,2); % number of harmonic components
Bn = An.*exp(-1i*phin);
for n=2:nnum
    Bn(:,n) = Bn(:,n)./Bn(:,1)/n;
end

% removing dilated components
iter = ceil(log2(tau_len/tau_n0+1))-1;
for i=1:sz(1)
    for j=1:sz(2)
        timprof = permute(beta(i,j,:),[3 2 1]); % time profile
        m = 1;
        for it = 1:iter % iterations
            a = fft(timprof,N);
%             freprof = a(1:(fre_n0+fre_len)); % frequency profile
            freprof = a(1:(fre_n0+fre_len*2)); % frequency profile
            for n = 2:nnum
                dil_tau_n0 = tau_n0*n*m;
                dil_taulen = min(dil_tau_n0,tau_n0+tau_len-dil_tau_n0);
                if dil_taulen<=0 % dil_tau_n0 exceed time profile
                    break;
                end
                dil_fre = downsample(freprof,n)*n; % shrink spectrum 
                dil_frelen = length(dil_fre)-fre_n0;
                dil_fre(1:fre_n0) = 0;
                dil_fre(fre_n0+1:end) = dil_fre(fre_n0+1:end).*Bn(1:dil_frelen,n);
                a = ifft(dil_fre,N);
                dil_tau = 2*real(a(dil_tau_n0+(1:dil_taulen)));
                timprof(dil_tau_n0-tau_n0+(1:dil_taulen)) = ...
                    timprof(dil_tau_n0-tau_n0+(1:dil_taulen)) - dil_tau;
            end
            m = m*2;
        end % end of iterations
        beta(i,j,:) = permute(timprof,[3 2 1]);
    end
end 
    
beta1 = beta;
