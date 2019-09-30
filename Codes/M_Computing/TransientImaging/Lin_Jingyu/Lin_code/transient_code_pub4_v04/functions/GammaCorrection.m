function [new_img] = GammaCorrection(img,scaling)
% GammaCorrection  Normalize images and implement gamma correction
% Parameters:
%   img - images for processing
%   scaling - scaling factor
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

% scaling = 1.0; %Amplitude scale
img = max(img,0);
img = img/max(img(:));
img = min(scaling*img,1);
new_img = img.^(1/2.2); %Gamma
