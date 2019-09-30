function [trans_boost] = BoostTransient6s(trans_img,llow,pk_sigma,edg,iter)
% BoostTransient  Compensating low frequency components, stable
% Method: boost peak iteratively, SOFT threshold, edge omitting
% Parameters:
%   trans_img - images for processing
%   llow - sin(k0*fre0*t)./sin(k0*fre_step/2*t); low envelop filter, odd length
%   thr - threshold for peak detection
%   edg - not processed edge 
%   iter - iteration time
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

sz = size(trans_img);
l_len = length(llow);
l_ctr = floor(l_len/2)+1; % filter center
trans_boost = trans_img;
% thr = max(trans_img,[],3)/6;
% edg = 50;

% iterations
for n = 1:iter
    % find peak and convolute
    bm = min(trans_boost,[],3);
    low_env = zeros(sz);
    for i=edg+1:sz(3)-edg
        dat = trans_boost(:,:,i);
        cond = max(dat./pk_sigma,0);
        pk = (dat-bm).*(1-exp(-cond.*cond)); % peak
        posc = l_ctr - i;
        p1 = max(1-posc,1);
        p2 = min(l_len-posc,sz(3));
        for j=p1:p2 % convolution
            low_env(:,:,j) = low_env(:,:,j) + pk*llow(posc+j) ;
        end
    end
    % boost transient images
    trans_boost = low_env + trans_img;
end
