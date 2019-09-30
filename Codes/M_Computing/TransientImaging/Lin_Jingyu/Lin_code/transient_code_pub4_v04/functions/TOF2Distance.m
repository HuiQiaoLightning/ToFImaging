function distances = TOF2Distance(objdepth,tau_step,d0)
% TOF2Distance - convert TOF depth to real distance 
% Parameters:
%   objdepth - depthmap in TOF units
%   tau_step - TOF unit
%   d0 - synchronizing signal delay offset
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

if ~exist('d0','var')
    d0 = 3660; 
end

SPEEDOFLIGHT = 299792458e-6; % mm/ns
k = SPEEDOFLIGHT*tau_step/2; % mm/depth
distances = objdepth*k-d0;
