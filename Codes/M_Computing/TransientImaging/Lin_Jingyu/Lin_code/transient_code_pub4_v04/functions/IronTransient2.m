function [trans_iron] = IronTransient2(trans_img,ksigma)
% IronTransient  Iron transient images from both ends.
% Method: softly suppress leading slope.
% Parameters:
%   trans_img - images for processing.
%   ksigma - ratio for supressing(<1), high for great leading slope
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

sz = size(trans_img);
trans_iron = trans_img;
diff_img = diff(trans_img,1,3);

% forward
thr = max(trans_img,[],3)*ksigma;
cond = zeros(sz(1),sz(2));
for i=1:sz(3)
    dat = trans_iron(:,:,i);
    cond = dat>thr | cond;
    d = dat./thr;
    K = (exp(d.*d)-1)/(exp(1)-1);
    K = max(min(K,1),cond);
    trans_iron(:,:,i) = dat.*K;
end

% backward
dat = zeros(sz(1),sz(2));
trans_iron(:,:,sz(3)) = dat;
for i=sz(3)-1:-1:2
    dat = max(dat-diff_img(:,:,i),0);
    trans_iron(:,:,i) = min(dat,trans_iron(:,:,i));
end
