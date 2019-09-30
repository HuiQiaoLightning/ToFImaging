function [trans_iron] = IronTransient(trans_img)
% IronTransient  Iron transient images from both ends 
% Method: Both ends start from zeros, keep positive vlaues.
% Parameters:
%   trans_img - images for processing.
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

sz = size(trans_img);
trans_iron = zeros(sz);
diff_img = diff(trans_img,1,3);
% dat1 = zeros(sz(1),sz(2));
% dat2 = zeros(sz(1),sz(2));
% hftime = ceil(sz(3)/2);
% for i=2:hftime
%     dat1 = max(dat1+diff_img(:,:,i-1),0);
%     dat2 = max(dat2-diff_img(:,:,end-i+2),0);
%     trans_iron(:,:,i) = dat1;
%     trans_iron(:,:,end-i+1) = dat2;
% end
dat = zeros(sz(1),sz(2));
for i=2:sz(3)-1
    dat = max(dat+diff_img(:,:,i-1),0);
    trans_iron(:,:,i) = dat;
end
dat = zeros(sz(1),sz(2));
for i=sz(3)-1:-1:2
    dat = max(dat-diff_img(:,:,i),0);
    trans_iron(:,:,i) = min(dat,trans_iron(:,:,i));
end
