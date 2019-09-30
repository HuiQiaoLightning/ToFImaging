function [trans_smooth] = SmoothTransient1(trans_img,f_len,sigma,sigma_pk)
% SmoothTransient  Smoothing transient image by bilateral filter
% Parameters:
%   trans_img - images for processing
%   f_len - half filter length
%   sigma - for edge reservation, high to suppress great edge
%   sigma_pk - for peak reservation, high to suppress great peak
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

% sz = size(trans_img);
% trans_smooth = zeros(sz);
% for i=1:sz(3)
%     W = zeros(sz(1:2)); % weights
%     for j=max(i-f_len,1):min(i+f_len,sz(3))
%         d = (trans_img(:,:,i) - trans_img(:,:,j))/sigma;
%         gau = exp(-d.*d);
%         W = W + gau;
%         trans_smooth(:,:,i) = trans_smooth(:,:,i) + trans_img(:,:,j).*gau;
%     end
%     trans_smooth(:,:,i) = trans_smooth(:,:,i)./W;
% end

% peak-researved bilateral filter 
sz = size(trans_img);
trans_smooth = trans_img;
% sigma2 = sigma/10;
for i=1:sz(3)
    W = ones(sz(1:2)); % weights
    dat = trans_img(:,:,i);
    pk = max(dat/sigma_pk,0); % preventing peak from filtered
    for j=max(i-f_len,1):min(i+f_len,sz(3))
        if j==i 
            continue;
        end
        d = (dat - trans_img(:,:,j))/sigma;
        gau = exp(-d.*d-pk); % gaussian weight
        W = W + gau;
        trans_smooth(:,:,i) = trans_smooth(:,:,i) + trans_img(:,:,j).*gau;
    end
    trans_smooth(:,:,i) = trans_smooth(:,:,i)./W;
end
