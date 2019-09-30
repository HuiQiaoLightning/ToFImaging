function [trans_smooth] = SmoothTransient2(trans_img,f_len,ksigma,ksigma_pk)
% SmoothTransient  Smoothing transient image by bilateral filter, adaptive sigma
% Parameters:
%   trans_img - images for processing
%   f_len - half filter length
%   ksigma - ratio for edge reserving (<1), high to suppress great edge
%   ksigma_pk - ratio for peak reserving (<1), high to suppress great peak

sz = size(trans_img);
trans_smooth = trans_img;
sigma = (max(trans_img,[],3)-min(trans_img,[],3))*ksigma;
sigma_pk = (max(trans_img,[],3)-min(trans_img,[],3))*ksigma_pk;
for i=1:sz(3)
    W = ones(sz(1:2)); % weights
    dat = trans_img(:,:,i);
    pk = max(dat./sigma_pk,0); % preventing peak from filtered
    for j=max(i-f_len,1):min(i+f_len,sz(3))
        if j==i 
            continue;
        end
        d = (dat - trans_img(:,:,j))./sigma;
        gau = exp(-d.*d-pk); % gaussian weight
        W = W + gau;
        trans_smooth(:,:,i) = trans_smooth(:,:,i) + trans_img(:,:,j).*gau;
    end
    trans_smooth(:,:,i) = trans_smooth(:,:,i)./W;
end
