function NullSpace = eqconstrnullspace(A,nDims)
%EQCONSTRNULLSPACE is private to PATTERNSEARCH
% Function that returns the null space of the input matrix A .
% nDims is the dimensionality of the optimization space.
%
% Private to PATTERNSEARCH
%
% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/11/05 16:59:43 $
%

% One can also do this by finding the non trivial orthonormal basis from QR itself
% but it currently causes scaling issues in some of the test problems. 
% The NULL function uses SVD which is less efficient  and more appropriate if A 
% might not be full row rank, but rank(A) = size(A,1) is guaranteed by pre processing.

numEqConstr = size(A,1);
if ~isempty(A)
    % Find Null Space of equality constraints
    [Q,~] = qr(A',0);
    if (numEqConstr <= nDims)
		NullSpace = eye(nDims) - Q*Q';
    else
        error('globaloptim:eqconstrnullspace:Infeasible', ...
		'Cannot proceed further, too many independent equality constraints.\n No feasible points found.'); 
    end
else
    NullSpace = eye(nDims);
end

