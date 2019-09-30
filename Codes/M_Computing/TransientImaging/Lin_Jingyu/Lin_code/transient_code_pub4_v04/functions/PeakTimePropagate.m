function [prog_map, c_map] = PeakTimePropagate(peak_map,dark,n_clr,n_sf,wp)
% PeakTimeShow - convert peak time map to propagation.
% Parameters:
%   peak_map - input peak time map 
%   dark - dark region
%   n_clr - number of color
%   n_sf - for smooth filter
%   wp - width of a propagation wave
% Output:
%   prog_map - progation map 
%   c_map - color map
% Command:
%   figure(1); imagesc(prog_map); colormap(c_map); colorbar
%   imwrite(prog_map,colormap(c_map),'prog_map.jpg','jpg')
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%
if ~exist('r','var')
    wp = 0.16; % width of a propagation wave
end
if ~exist('n_sf','var')
    n_sf = 4; % for smooth filter
end
sf = ones(n_sf,n_sf)/n_sf/n_sf; % smooth filter
if ~exist('n_clr','var')
    n_clr = 16; % number of color
end
c_map = jet(n_clr); c_map(1,:) = [0 0 0]; % colormap

a = conv2(ones(size(peak_map)),sf,'same');
peak_map = conv2(peak_map,sf,'same')./a;
peak_map = peak_map - min(peak_map(:)); 
peak_map = peak_map.*dark;
p_max = max(peak_map(:));
c_step = p_max/n_clr; % color step
a = peak_map/c_step;
prog_map = round(a).*(abs(a-round(a))<wp);
prog_map = prog_map + 1;
% figure(3); imagesc(prog_map); colormap(c_map); colorbar
