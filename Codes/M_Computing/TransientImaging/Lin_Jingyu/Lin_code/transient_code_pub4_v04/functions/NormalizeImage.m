function [norm_img] = NormalizeImage(img,sigma)
% NormalizeImage - Normalize transient image with [1-exp(-t/sigma)]/N.
% Parameters:
%   img - transient images for processing
%   sigma - scaling factor, + for dark improving, - for dark suppressing
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

if ~exist('sigma','var')
    sigma = 0.5;
end
N = 1-exp(-1/sigma); % normalizing factor

% computing weight
img = max(img,0);
steady_img = max(img,[],3);
M = max(steady_img(:));
norm_val = (1-exp(-steady_img/M/sigma))/N;
mask = (steady_img>0);
wgh = norm_val./(steady_img+~mask).*mask;

% normalizing
tau_len = size(img,3);
norm_img = img;
for i=1:tau_len
    norm_img(:,:,i) = img(:,:,i).*wgh;
end
