function options  = checkoptions_gads(options,defaultopt,numberOfVariables)
%CHECKOPTIONS validates all PATTERNSEARCH options before they are used by
%   solver
%
%   private to pfminlcon, pfminbnd, and pfminunc.

%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:01 $

% Sanity check for the options structure
options = psoptimset_gads(options);

options.Display = psoptimget_gads(options,'Display',defaultopt,'fast');

% Define verbosity here (Later we can use options structure)
switch options.Display  
    case {'off','none'}
        options.Verbosity = 0;
    case 'final'
        options.Verbosity = 1;    
    case 'iter'
        options.Verbosity = 2;
    case 'diagnose'
        options.Verbosity = 3;
    otherwise
        options.Verbosity = 1;
end

% Retrieve options using PSOPTIMGET
options.MeshExpansion     = psoptimget_gads(options,'MeshExpansion',defaultopt,'fast'); 
options.MeshContraction   = psoptimget_gads(options,'MeshContraction',defaultopt,'fast'); 
options.CompleteSearch    = psoptimget_gads(options,'CompleteSearch',defaultopt,'fast');
options.MeshAccelerator   = psoptimget_gads(options,'MeshAccelerator',defaultopt,'fast');
options.TolMesh           = psoptimget_gads(options,'TolMesh',defaultopt,'fast');
options.TolCon            = psoptimget_gads(options,'TolCon',defaultopt,'fast');
options.MaxMeshSize       = psoptimget_gads(options,'MaxMeshSize',defaultopt,'fast');
options.MaxIter           = psoptimget_gads(options,'MaxIter',defaultopt,'fast');
options.MaxFunEvals       = psoptimget_gads(options,'MaxFunEvals',defaultopt,'fast');
options.TimeLimit         = psoptimget_gads(options,'TimeLimit',defaultopt,'fast');
options.TolBind           = psoptimget_gads(options,'TolBind',defaultopt,'fast');
options.TolFun            = psoptimget_gads(options,'TolFun',defaultopt,'fast');
options.TolX              = psoptimget_gads(options,'TolX',defaultopt,'fast');
options.InitialMeshSize   = psoptimget_gads(options,'InitialMeshSize',defaultopt,'fast');
options.PollMethod        = psoptimget_gads(options,'PollMethod',defaultopt,'fast');
options.PollingOrder      = psoptimget_gads(options,'PollingOrder',defaultopt,'fast');
options.CompletePoll      = psoptimget_gads(options,'CompletePoll',defaultopt,'fast');
options.PlotInterval      = psoptimget_gads(options,'PlotInterval',defaultopt,'fast');
options.Vectorized        = psoptimget_gads(options,'Vectorized',defaultopt,'fast');
options.Cache             = psoptimget_gads(options,'Cache',defaultopt,'fast');
options.CacheTol          = psoptimget_gads(options,'CacheTol',defaultopt,'fast');
options.CacheSize         = psoptimget_gads(options,'CacheSize',defaultopt,'fast');
options.ScaleMesh         = psoptimget_gads(options,'ScaleMesh',defaultopt,'fast');
options.MeshRotate        = psoptimget_gads(options,'MeshRotate',defaultopt,'fast');
options.InitialPenalty    = psoptimget_gads(options,'InitialPenalty',defaultopt,'fast');
options.PenaltyFactor     = psoptimget_gads(options,'PenaltyFactor',defaultopt,'fast');
options.UseParallel       = psoptimget_gads(options,'UseParallel',defaultopt,'fast');
% These options will be stuffed in the structure later (after some
% processing)
outputFcns        = psoptimget_gads(options,'OutputFcns',defaultopt,'fast');
plotFcns          = psoptimget_gads(options,'PlotFcns',defaultopt,'fast');
searchFcn        = psoptimget_gads(options,'SearchMethod',defaultopt,'fast');
           
% Modify some fields if they are not yet assigned
if ischar(options.MaxFunEvals) && isequal(lower(options.MaxFunEvals),'2000*numberofvariables')
        options.MaxFunEvals = 2000*numberOfVariables;
end
if ischar(options.MaxIter) && isequal(lower(options.MaxIter),'100*numberofvariables')
        options.MaxIter = 100*numberOfVariables;
end

options.MaxFunEvals  = floor(options.MaxFunEvals);
options.MaxIter = floor(options.MaxIter);

% If searchFcn is a cell array with additional arguments, handle them
if iscell(searchFcn)
    searchFcnArg = searchFcn(2:end);
    searchFcn = searchFcn{1};
else
    searchFcnArg = {};
end
% Search technique could be [], char, or function_handle
if isempty(searchFcn)
    searchFcn = [];
elseif isa(searchFcn,'function_handle')
    [searchFcn,msg] = fcnchk(searchFcn);
    if ~isempty(msg)
        error('globaloptim:checkoptions:searchFcnArgCheck',msg);
    end
    searchFcnString = func2str(searchFcn);
elseif ischar(searchFcn) 
    [searchFcn,msg] = fcnchk(searchFcn);
    if ~isempty(msg)
        error('globaloptim:checkoptions:searchFcnArgCheck',msg);
    end
    searchFcnString = func2str(searchFcn);
else
    error('globaloptim:checkoptions:invalidSearchMethod','Invalid choice of Search method: See psoptimset for SearchMethod.\n');    
end

% Make sure that search method is different from poll method and not 'none'
if ~isempty(searchFcn) && any(strcmpi(searchFcnString,{options.PollMethod,'none'}))
   searchFcn = [];
end

% Only some choices can be strings (special case)
if isa(searchFcn,'function_handle') && any(strcmpi(searchFcnString,{'positivebasisnp1', 'positivebasis2n', ...
            'gpspositivebasisnp1','gpspositivebasis2n','madspositivebasisnp1', 'madspositivebasis2n', ...
            'gsspositivebasisnp1','gsspositivebasis2n'}))
        searchFcn = searchFcnString; % Convert to a string because these are not functions
end


options.SearchMethod = searchFcn;
options.SearchMethodArg = searchFcnArg;

% If options.MaxMeshSize is less than options.Meshsize (This should not happen)
if options.MaxMeshSize < options.InitialMeshSize
    warning('globaloptim:checkoptions:maxMeshSize','MaxMeshSize should be greater than InitialMeshSize;changing initial mesh size to MaxMeshSize.\n');
    options.InitialMeshSize = options.MaxMeshSize;
end

% It is NOT vectorized in these conditions
options.NotVectorizedPoll   = (strcmpi(options.Vectorized,'off') || ...
    (strcmpi(options.Vectorized, 'on') && strcmpi(options.CompletePoll,'off')));
options.NotVectorizedSearch = (strcmpi(options.Vectorized,'off') || ...
    (strcmpi(options.Vectorized, 'on') && strcmpi(options.CompleteSearch,'off')));

% If using 2N basis or MADS RotatePattern has no effect.
if any(strcmpi(options.PollMethod,{'positivebasis2n','gpspositivebasis2n', ...
    'madspositivebasisnp1','madspositivebasis2n','gsspositivebasis2n'}))
    options.MeshRotate = 'off';
end
% Error checking on InitialPenalty and PenaltyFactor
if options.InitialPenalty < 1
    warning('globaloptim:checkoptions:smallPenalty','InitialPenalty must be greater than or equal to one; using default.\n');
    options.InitialPenalty = defaultopt.InitialPenalty;
end
% Penalty factor to increase penalty
if options.PenaltyFactor <= 1
    warning('globaloptim:checkoptions:smallPenaltyFactor','PenaltyFactor must be greater than one; using default.\n');
    options.PenaltyFactor = defaultopt.PenaltyFactor;
end
% If outputFcns is a cell array with additional arguments, handle them
[options.OutputFcns,options.OutputFcnsArg] = functionHandleOrCellArray('OutputFcns',outputFcns);

 if isempty(options.OutputFcns)
     options.OutputTrue = false;
 else
     options.OutputTrue = true;
 end
 % If plotFcns is a cell array with additional arguments, handle them
[options.PlotFcns,options.PlotFcnsArgs] = functionHandleOrCellArray('PlotFcns',plotFcns);

 if isempty(options.PlotFcns)
     options.PlotTrue = false;
 else
     options.PlotTrue = true;
 end
 
% Test for valid strings
if ~isempty(options.UseParallel)
    stringSet('UseParallel',options.UseParallel,{'never','always'});
    options.SerialUserFcn = strcmpi(options.UseParallel,'never');
else
    options.SerialUserFcn = true;
end

