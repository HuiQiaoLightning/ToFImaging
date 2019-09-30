function [norm_img] = CompressHDR(img,thr)
% CompressHDR - Compress HDR transient image to 256 levels and normalize.
% Method: histogram equalization
% Parameters:
%   img - HDR transient images for processing
%   thr - normalized threshold for gray level 0.
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

if ~exist('thr','var')
    thr = 0;
end
img = max(img,0);
steady_img = max(img,[],3);
thr = thr*max(steady_img(:));

% sort gray values
s = sort(steady_img(:));
N = length(s);
for i=1:N  % remove zero pixels
    if s(i)>thr
        s = s(i:end);
        break;
    end
end
zeromask = (steady_img<=thr);  % gray level 0

% compute mapping ratio
% N = floor(length(s)/255);
N = length(s)/255;  % gray level step
inv_img = 1./(steady_img + zeromask);
wgh = zeros(size(steady_img));
mask0 = wgh;
for i=254:-1:1  % gray level 2~255
    n = ceil(i*N);
    mask = (steady_img>s(n));
    wgh = wgh + inv_img.*(mask&~mask0)*(i+1);
    mask0 = mask; % record last mask
end
wgh = wgh + inv_img.*~mask;  % gray level 1
wgh = wgh.*~zeromask/255;

% compressing range
tau_len = size(img,3);
norm_img = img;
for i=1:tau_len
    norm_img(:,:,i) = img(:,:,i).*wgh;
end
