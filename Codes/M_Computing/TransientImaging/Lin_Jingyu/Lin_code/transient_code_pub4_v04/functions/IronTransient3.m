function [trans_iron] = IronTransient3(trans_img,kthr,sigma)
% IronTransient - Reshape transient image.
% Method: softly suppress value below threshold.
% Parameters:
%   trans_img - images for processing.
%   kthr - relative threshold for supressing (default 0.1).
%   sigma - high value for heavier suppressing (default 1).
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

if ~exist('kthr','var')
    kthr = 0.1;
end
if ~exist('sigma','var')
    sigma = 1;
end
sz = size(trans_img);
min_trans = min(trans_img,[],3);
max_trans = max(trans_img,[],3);
thr = (max_trans-min_trans)*kthr;
sigma2 = sigma*sigma;
trans_iron = trans_img;
for i=1:sz(3)
    dat = trans_iron(:,:,i)-min_trans;
    cond = dat>thr;
    d = dat./thr;
    K = (exp(d.*d*sigma2)-1)/(exp(sigma2)-1);
%     K = max(min(K,1),cond);
    K = K.*~cond + cond;
    trans_iron(:,:,i) = dat.*K;
end
