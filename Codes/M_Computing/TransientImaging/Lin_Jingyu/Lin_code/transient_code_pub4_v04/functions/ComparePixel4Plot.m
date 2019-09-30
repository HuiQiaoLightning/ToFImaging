function ComparePixel4Plot(imgs1, imgs2, imgs3, imgs4, y, x)
% ComparePixel6Plot  Compare plots of a pixel at 4 image sequences
% Parameters:
%   imgs - image sequences for comparing
%   y, x - pixel coordinates at the image sequences
%
% Copyright (C) 2013 Jinyu Lin, linjy02@hotmail.com
%

pixel = [imgs1(y,x,:), imgs2(y,x,:), imgs3(y,x,:), imgs4(y,x,:)];
pl = permute(pixel,[3 2 1]);
pl = pl(:,:);
plot(pl); legend('1', '2', '3', '4')