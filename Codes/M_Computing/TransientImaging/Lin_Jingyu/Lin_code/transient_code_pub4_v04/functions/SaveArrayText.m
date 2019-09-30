function SaveArrayText(fn,arr,fmt)
% SaveArrayText - Save a 2D array to a txt file. 
% Parameters:
%   fn - file name
%   arr - 2D array to be saved
%   fmt - fprintf format string for each element
%
% Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com
%

if ~exist('fmt','var')
    fmt = '%f';
end

dat_row = size(arr,1);
dat_col = length(arr(:))/dat_row; 
f_dat = reshape(arr,[dat_row dat_col]);
[fid, msg] = fopen(fn,'w');
if fid==-1
    fprintf(msg);
    fprintf('\n');
    return;
end
for i=1:dat_row
    for j=1:dat_col
        fprintf(fid,fmt,f_dat(i,j));
        fprintf(fid,'  ');        
    end
    fprintf(fid,'\r\n');        
end
fclose(fid);
