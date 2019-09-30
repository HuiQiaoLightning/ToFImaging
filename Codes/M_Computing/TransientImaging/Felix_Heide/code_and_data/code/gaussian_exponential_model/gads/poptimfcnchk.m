function [Iterate,count] = poptimfcnchk(FUN,nonlcon,Xin,Iterate,Vectorized,objFcnArg,conFcnArg)
%POPTIMFCNCHK Calls objective function 'FUN' and nonlinear constraint.
%   function 'nonlcon' for the first time at the start point 'Iterate.x'. If
%   'Vectorized' option in 'on' the objective function is called with two
%   points.

%   Private to PATTERNSEARCH.


%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:36 $

y = NaN;
count = 0;
X = Iterate.x;
% Check the objective function. 
[FUN,msg] = fcnchk(FUN);
if ~isempty(msg)
    error('globaloptim:poptimfcnchk:objFcnArgCheck',msg);
end
if strcmpi(Vectorized,'off')  % The function is not vectorized. 
    try
        [y,count] = funevaluate(FUN,Xin,X,'init',[],[],objFcnArg{:});
    catch userFcn_ME
        gads_ME = MException('globaloptim:poptimfcnchk:objfunCheck', ...
            'Failure in initial user-supplied objective function evaluation. PATTERNSEARCH cannot continue.');
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
   end
    if numel(y) ~=1
        error('globaloptim:poptimfcnchk:objfunCheck','%s\n', ...
            'Your objective function must return a scalar value.');
    end
elseif strcmpi(Vectorized,'on') % If vectorized is 'on', Completepoll MUST be 'on' too
    X2 = [X, X];
    try
        [f,count] = funevaluate(FUN,Xin,X2,'init',[],[],objFcnArg{:});
    catch userFcn_ME
        gads_ME = MException('globaloptim:poptimfcnchk:objfunCheck', ...
            ['Failure in initial user-supplied objective function evaluation ' ...
            '(hint: ''vectorized'' option is ''on''). PATTERNSEARCH cannot continue.']);
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
    end
    if 2 ~=numel(f)
        error('globaloptim:poptimfcnchk:objfunCheck','%s\n', ...
            ['When ''Vectorized'' is ''on'', your objective function must ' ...
             'return a vector of length equal to number of input points.']);
    end
    y = f(1);
end
% We want the function value to be real
if isnan(y)
    error('globaloptim:poptimfcnchk:objfunNaN','Objective function must be real and not NaN at the starting point.');
end
Iterate.f = y;

% Evaluate nonlinear constraints for the first time
if  ~isempty(nonlcon)
    try
        [cineq,ceq] = feval(nonlcon,reshapeinput(Xin,X),conFcnArg{:});
        Iterate.cineq = zeros(numel(cineq),1);
        Iterate.ceq = zeros(numel(ceq),1);
        Iterate.cineq(:) = cineq;
        Iterate.ceq(:) = ceq;
    catch userFcn_ME
        gads_ME = MException('globaloptim:poptimfcnchk:confunCheck', ...
            'Failure in initial user-supplied nonlinear constraint function evaluation. PATTERNSEARCH cannot continue.');
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
    end
    c = [Iterate.cineq;Iterate.ceq];
    if ~all(isreal(c) & isfinite(c))
        error('globaloptim:poptimfcnchk:confunNotReal','Constraint function must return real value.');
    end
end
