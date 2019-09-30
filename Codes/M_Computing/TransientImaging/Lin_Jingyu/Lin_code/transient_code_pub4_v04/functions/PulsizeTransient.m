function [pulse] = PulsizeTransient(trans_img,f_H,tau_step,psc,kthr)
% PulsizeTransient - Find and compress sparse pulses from transient images.
% Parameters:
%   trans_img - images for processing.
%   f_H - the highest frequency (MHz)
%   tau_step - time step (ns)
%   psc - scale for pulses.
%   kthr - relative threshold for supressing (default 0.1).
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

if ~exist('kthr','var')
    kthr = 0.1;
end
if ~exist('psc','var')
    psc = 1; % no scale
end
sz = size(trans_img);
% min_trans = min(trans_img,[],3);
% for i=1:sz(3);
%     trans_img(:,:,i) = trans_img(:,:,i) - min_trans;
% end
% thr = max(trans_img(:))*kthr;
min_trans = min(trans_img,[],3);
max_trans = max(trans_img,[],3);
thr = (max_trans-min_trans)*kthr;
pw = 1000/f_H/tau_step; % width of pulses
win = floor(pw/6)*2; % searching window
pulse = trans_img;
x = 1:sz(3);

for i=1:sz(1)
    for j=1:sz(2)
        pix = zeros(1,sz(3));
        dat = trans_img(i,j,:)-min_trans(i,j);
        dat = dat(:);
        k = 0;
        while k<=sz(3)-win
            [y,idx] = max(dat(k+(1:win)));
            if y>thr(i,j) && idx<=win*3/4
                % find a peak
                k = k+idx;
                pix = pix + y*gauss(x,k,pw*psc/1.414);
                while k<=sz(3)-win
                    [~,idx] = min(dat(k+(1:win)));
                    if idx<=win*3/4
                        % find a valley
                        k = k+idx;
                        break;                        
                    else
                        k = k + win/2;
                    end
                end
            else
                k = k + win/2;
            end
        end
        pulse(i,j,:) = reshape(pix,[1 1 sz(3)]);
    end % sz(2)
end % sz(1)
pulse = pulse/max(pulse(:)); % normalizing
return

function g = gauss(x,x0,sigma)
% gaussian kernel
x = (x-x0)/sigma;
g = exp(-x.*x);
return