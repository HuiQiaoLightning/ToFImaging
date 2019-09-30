function X = reshapeinput(Xin, X)
% RESHAPEINPUT reshape X to match the shape of Xin 

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:27:05 $


[m,n] = size(Xin);
% Scalar input is always col major
if m == n && m == 1 
%    X = X(:); % Will make it a row major
    return; % Retain shape
end

% Single point evaluation
if isvector(X) % X is a vector so use shape information
    Xin(:) = X;
    X = Xin;
    return;
elseif isvector(Xin)
    if  m > n && n == 1  % col major input
        return;
    elseif n > m && m == 1 % row major input
        X = X';
    end
else % Xin is a matrix; not a documented feature
    p = size(Xin,2);
    if p > 1          % Matrix Xin with vectorized 'on' 
        X = reshape(X,m,n,p);
    else
        X = reshape(X,m,n);
    end
end
