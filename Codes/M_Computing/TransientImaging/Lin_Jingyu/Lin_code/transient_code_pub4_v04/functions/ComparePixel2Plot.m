function ComparePixel2Plot(imgs1, imgs2, y, x)
% ComparePixel2Plot  Compare plots of a pixel at two image sequences
% Parameters:
%   imgs1, imgs2 - image sequences for comparing
%   y, x - pixel coordinates at the image sequences
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

pixel = [imgs1(y,x,:) imgs2(y,x,:)];
pl = permute(pixel,[3 2 1]);
pl = pl(:,:);
plot(pl); legend('1','2')