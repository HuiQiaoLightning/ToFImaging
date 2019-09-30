function [trans_boost] = BoostTransient1(trans_img,llow,pk_sigma)
% BoostTransient  Compensating low frequency components, fast and good
% Method: simply detect peak above zero and boost in time.
% Parameters:
%   trans_img - images for processing
%   llow - sin(k0*fre0*t)./sin(k0*fre_step/2*t); odd length
%   pk_sigma - for peak detection
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

sz = size(trans_img);
l_len = length(llow);
l_ctr = floor(l_len/2)+1;
bm = min(trans_img,[],3);
trans_boost = zeros(sz);
for i=1:sz(3)
    dat = trans_img(:,:,i);
    cond = max(dat/pk_sigma,0);
    pk = (dat-bm).*(1-exp(-cond.*cond)); % peak
    posc = l_ctr - i;
    p1 = max(1-posc,1);
    p2 = min(l_len-posc,sz(3));
    for j=p1:p2 % convolution
        trans_boost(:,:,j) = trans_boost(:,:,j) + pk*llow(posc+j) ;
    end
end
% bst = (-bm)./max(trans_boost,[],3);
% for i=1:sz(3)
%     trans_boost(:,:,i) = trans_boost(:,:,i).*bst + trans_img(:,:,i);
% end
trans_boost = trans_boost + trans_img;
