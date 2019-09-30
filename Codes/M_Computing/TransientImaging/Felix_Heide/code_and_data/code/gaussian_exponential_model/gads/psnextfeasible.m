function [msg,nextIterate,direction,FunEval] = psnextfeasible(ObjFunc,Xin,sites,Iterate, ...
    Aineq,bineq,Aeq,beq,lb,ub,tol,constr,objFcnArg,FunEval,options)
%PSNEXTFEASIBLE finds next feasible iterate among 'sites' in a pattern.
%   Input argument 'objfunc' is the objective function on which POLL step
%   is implemented, 'sites' is a set of points generated by POLL technique.
%   We would like to test the function at these sites hoping that one of them
%   have less function value than the current iterate.
%
%   Input argument 'Iterate' is the incumbent around which polling is done.
%   Iterate is a structure which have fields 'X', the current point and 'f',
%   the function value 'f' at the current point. The tolerance 'tol' is used
%   to determine if a point is feasible w.r.t. linear constraints.
%
%   Output argument 'msg' is a binary flag indicating whether a better iterate
%   is found or not, 'nextIterate' is a successful iterate after polling is done.
%   If POLL is unsuccessful 'nextIterate' is same as the input argument 'Iterate'.
%
%   The successful search direction vector is returned as 'direction' and number of
%   function evaluation is stored in 'FunEval'.
%
% 	Example:
% 	If there are 4 points in 2 dimension space then
%    X is     [2  1  9 -2
%              0  1 -4  5]

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:26:58 $


direction= [];
msg = false;
nextIterate = Iterate;
% Remove infeasible points
feasible = true(length(sites),1);
for k = 1:length(sites)
    %Initialize the default function value to Inf.
    sites(k).f = inf;
    if constr
        feasible(k) = isTrialFeasible(sites(k).x,Aineq,bineq,Aeq,beq,lb,ub,tol);
    end
end
sites(~feasible) = [];
if (isempty(sites))    % No design sites to evaluate
    return
end
% maxeval will make sure that we respect this limit on function evaluation.
maxeval = min(length(sites),options.MaxFunEvals - FunEval);
% Local variables to use in parfor loop
Cache = options.Cache;
CacheTol = options.CacheTol;
CacheSize = options.CacheSize;

if options.NotVectorizedPoll
    if strcmpi(options.CompletePoll,'on')
        trialX = [sites(:).x]; trialF = [sites.f];
        if options.SerialUserFcn
            for k = 1:maxeval
                trialF(k) = funevaluate(ObjFunc,Xin,trialX(:,k),Cache, ...
                    CacheTol,CacheSize,objFcnArg{:});
            end
        else
            parfor (k = 1:maxeval)
                % Do not use Cache when running parfor
                trialF(k) = funevaluate(ObjFunc,Xin,trialX(:,k),'off', ...
                    [],[],objFcnArg{:});
            end
        end
        FunEval = FunEval+maxeval;
        % This if condition will complete the parfor loop (get the best fval)
        [unused,index] = min(trialF);
        if nextIterate.f > trialF(index)
            nextIterate.x = trialX(:,index);
            nextIterate.f = trialF(index);
            direction = index;
            msg = true;
        end
    else % completePoll is 'off'
        for k = 1:maxeval
            [sites(k).f,count] = funevaluate(ObjFunc,Xin,sites(k).x,Cache, ...
                CacheTol,CacheSize,objFcnArg{:});
            FunEval = FunEval+count;
            if nextIterate.f > sites(k).f
                nextIterate.x = sites(k).x;
                nextIterate.f = sites(k).f;
                direction =k;
                msg = true;
                return;
            end
        end
    end
    % vectorized is 'on' (guaranteed from the previous if-condition),
    % if CompletePoll is 'on' then evaluate all the points
elseif strcmpi(options.CompletePoll,'on')
    Xfeas = [sites(:).x];
    feasible = true(size(Xfeas,2),1);
    % we must check all the points for feasibility before evaluating them
    if constr
        [feasible,Xfeas] = getFeasiblePoints([sites(1:maxeval).x],Aineq,bineq,Aeq,beq,lb,ub,tol);
    end
    [f,count] = funevaluate(ObjFunc,Xin,Xfeas,options.Cache,options.CacheTol, ...
        options.CacheSize,objFcnArg{:});
    maxeval = size(Xfeas,2);
    FunEval = FunEval+count;
    for i =1:maxeval
        sites(i).x = Xfeas(:,i);
        sites(i).f = f(i);
        if (nextIterate.f > sites(i).f)
            nextIterate.x  = sites(i).x;
            nextIterate.f  = sites(i).f;
            direction = i;
            msg = true;
        end
    end
    % offset direction to match the right order
    direction = direction + (length(feasible) - nnz(feasible));
end
