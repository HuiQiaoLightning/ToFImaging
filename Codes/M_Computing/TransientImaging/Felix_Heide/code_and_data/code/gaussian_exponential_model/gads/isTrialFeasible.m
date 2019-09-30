function feasible = isTrialFeasible(X,Aineq,bineq,Aeq,beq,lb,ub,tol)
%isTrialFeasible Checks if X is feasible w.r.t. linear constraints. 
% 	

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:26:24 $

feasibleIneq = true;
feasibleEq = true;

% Inequality constraints
if ~isempty(Aineq)
    feasibleIneq = max(Aineq*X-bineq) <= tol;
end
% Upper bounds
argub = ~eq(ub,inf);
if any(argub) && ~isempty(argub)
    feasibleIneq = feasibleIneq && all(X(argub) <= ub(argub));
end
% Lower bounds
arglb = ~eq(lb,-inf);
if any(arglb) && ~isempty(arglb)
    feasibleIneq = feasibleIneq && all(X(arglb) >= lb(arglb));
end
% Equality constraints
if ~isempty(Aeq)
    constrViolation = (Aeq*X-beq);
    feasibleEq   = all(abs(constrViolation(~isinf(constrViolation))) <= tol);
end
feasible = feasibleIneq && feasibleEq;