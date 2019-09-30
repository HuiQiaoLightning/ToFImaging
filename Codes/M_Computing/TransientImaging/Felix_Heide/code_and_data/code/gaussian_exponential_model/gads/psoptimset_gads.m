function options = psoptimset_gads(varargin)
%PSOPTIMSET Create/alter PATTERNSEARCH OPTIONS structure.
%   OPTIONS = PSOPTIMSET('PARAM1',VALUE1,'PARAM2',VALUE2, ...) creates an
%   optimization options structure OPTIONS in which the named parameters
%   have the specified values. Any unspecified parameters are set to []
%   (parameters with value [] indicate to use the default value for that
%   parameter when  OPTIONS is passed to the optimization function). It is
%   sufficient to type  only the leading characters that uniquely identify
%   the parameter. Case is  ignored for parameter names. NOTE: For
%   values that are strings, correct case and the complete string are
%   required.
%   
%   OPTIONS = PSOPTIMSET(OLDOPTS,'PARAM1',VALUE1, ...) creates a copy of
%   OLDOPTS with the named parameters altered with the specified values.
%   
%   OPTIONS = PSOPTIMSET(OLDOPTS,NEWOPTS) combines an existing options
%   structure OLDOPTS with a new options structure NEWOPTS. Any parameters
%   in NEWOPTS with non-empty values overwrite the corresponding old
%   parameters in OLDOPTS. 
%   
%   PSOPTIMSET with no input arguments and no output arguments displays all
%   parameter names and their possible values, with defaults shown in {}.
%
%   OPTIONS = PSOPTIMSET (with no input arguments) creates an options
%   structure OPTIONS where all the fields are set to defaults for
%   PATTERNSEARCH.
%
%PSOPTIMSET PARAMETERS
%   TolMesh              - Termination tolerance on mesh size 
%                          [ positive scalar | {1e-6} ]
%   TolCon               - Termination tolerance on constraints 
%                          [ positive scalar | {1e-6} ]
%   TolFun               - Termination tolerance on function value
%                          [ positive scalar | {1e-6} ] 
%   TolX                 - Termination tolerance on X 
%                          [ positive scalar | {1e-6} ] 
%   TolBind              - Binding Tolerance
%                          [ positive scalar | {1e-3} ]
%   MaxIter              - Maximum number of iterations allowed 
%                          [ positive scalar | {100*numberOfVariables} ]
%   MaxFunEvals          - Maximum number of function (objective)
% 	                        evaluations allowed 
%                          [ positive scalar | {2000*numberOfVariables} ]
%   TimeLimit            - Total time (in seconds) allowed for optimization
%                          [ positive scalar | {Inf} ]
% 	
%   MeshContraction      - Mesh refining factor used when an iteration is not
%                          successful 
%                          [ positive scalar | {0.5} ]
%   MeshExpansion        - Mesh coarsening factor used when an iteration is
%                          successful 
%                          [ positive scalar | {2.0} ] 
%   MeshAccelerator      - Accelerate convergence when near a minima (may
%                          lose some accuracy)
%                          [ 'on',{'off'} ]
%   MeshRotate           - Rotate the pattern before declaring a point to
%                          be optimum
%                          [ 'off',{'on'} ]
%   InitialMeshSize      - Initial mesh size at start
%                          [ positive scalar | {1.0} ]
%   ScaleMesh            - Mesh is scaled if 'on'
%                          [ 'off', {'on'} ]
%   MaxMeshSize          - Maximum mesh size used in a POLL/SEARCH step
%                          [ positive scalar | {Inf} ]
%
%   InitialPenalty       - Initial value of penalty parameter
%                          [ positive scalar | {10} ]
%   PenaltyFactor        - Penalty update parameter
%                          [ positive scalar | {100} ]
%
%   PollMethod           - Polling method
%                          [ 'GSSPositiveBasisNp1' | 'GSSPositiveBasis2N' |
%                            'GPSPositiveBasisNp1' | {'GPSPositiveBasis2N'} |
%                            'MADSPositiveBasisNp1' | 'MADSPositiveBasis2N' ]
%   CompletePoll         - Complete polling around the current iterate
%                          [ 'on',{'off'} ]
%   PollingOrder         - Ordering of the POLL directions
%                          [ 'Random' |'Success'| {'Consecutive'} ]
% 	
%   SearchMethod         - Search method
%                          [ @GSSPositiveBasisNp1 | @GSSPositiveBasis2N |
%                            @GPSPositiveBasisNp1 | @GPSPositiveBasis2N  |
%                            @MADSPositiveBasisNp1 | @MADSPositiveBasis2N |
%                            @searchlhs        | @searchneldermead | 
%                            @searchga         | {[]} ]
%   CompleteSearch       - Complete search around the current iterate
%                          [ 'on',{'off'} ]
%
%   Display              - Level of display 
%                          [ 'off' | 'iter' | 'diagnose' | {'final'} ]
%   OutputFcns           - A set of functions called in every iteration 
%                          [ @psoutputhistory | {[]} ]
%   PlotFcns             - A set of specialized plot functions called in
%                          every iteration 
%                          [ @psplotbestf    |  @psplotmeshsize | 
%                            @psplotfuncount |  @psplotbestx    | {[]} ]
%   PlotInterval         - Plot functions will be called every interval
%                          [ {1} ]
%
%   Cache                - Use a CACHE of the points evaluated. This
% 	                        options could be expensive in terms of memory 
%                          and speed . If the objective function is
%                          stochastic, it is advised not to use this option.
%                          [ 'on', {'off'} ]
%   CacheSize            - Limit the CACHE size. A typical choice of 1e4 is
%                          suggested but it depends on computer speed and 
%                          memory. Used if 'Cache' option is 'on'.
%                          [ positive scalar | {1e4} ]
%   CacheTol             - Tolerance used to determine if two points are
%                          close enough to be declared same. 
%                          Used if 'Cache' option is 'on'. 
%                          [ positive scalar | {eps} ]
% 	
%   Vectorized           - Objective function is vectorized and it can
%                          evaluate more than one point in one call 
%                          [ 'on' | {'off'} ]
% 	 	
%   UseParallel          - Use PARFOR to evaluate objective and nonlinear 
%                          constraint functions.
%                          [ 'always' | {'never'} ]
% 	 See also PSOPTIMGET.


% Print out possible values of properties, when called with no input/output arguments
if (nargin == 0) && (nargout == 0)
    fprintf('                 TolMesh: [ positive scalar | {1e-6} ]\n');
    fprintf('                  TolCon: [ positive scalar | {1e-6} ]\n');
    fprintf('                    TolX: [ positive scalar | {1e-6} ]\n');
    fprintf('                  TolFun: [ positive scalar | {1e-6} ]\n');
    fprintf('                 TolBind: [ positive scalar | {1e-3} ]\n');
    fprintf('                 MaxIter: [ positive scalar | {100*numberOfVariables} ]\n');
    fprintf('             MaxFunEvals: [ positive scalar | {2000*numberOfVariables} ]\n');
    fprintf('               TimeLimit: [ positive scalar | {Inf} ]\n\n');
    
    fprintf('         MeshContraction: [ positive scalar | {0.5} ]\n');
    fprintf('           MeshExpansion: [ positive scalar | {2.0} ]\n');
    fprintf('         MeshAccelerator: [ on  | {off} ]\n');
    fprintf('              MeshRotate: [ off | {on} ]\n');
    fprintf('         InitialMeshSize: [ positive scalar | {1.0} ]\n');
    fprintf('               ScaleMesh: [ off | {on} ]\n');
    fprintf('             MaxMeshSize: [ positive scalar | {Inf} ]\n\n');

    fprintf('          InitialPenalty: [ positive scalar | {10} ]\n');
    fprintf('           PenaltyFactor: [ positive scalar | {100} ]\n\n');

    fprintf('              PollMethod: [ MADSPositiveBasisNp1 | MADSPositiveBasis2N |\n');
    fprintf('                            GPSPositiveBasisNp1 | {GPSPositiveBasis2N} ]\n');
    fprintf('            CompletePoll: [ on  | {off} ]\n');
    fprintf('            PollingOrder: [ Random | Success | {Consecutive} ]\n\n');
    
    fprintf('            SearchMethod: [ function_handle  | @MADSPositiveBasisNp1 |\n');
    fprintf('                            @MADSPositiveBasis2N | @GPSPositiveBasisNp1 |\n'); 
    fprintf('                            @GPSPositiveBasis2N | @GSSPositiveBasisNp1 |\n');
    fprintf('                            @GSSPositiveBasis2N | @searchga |\n'); 
    fprintf('                            @searchlhs       | @searchneldermead | {[]} ]\n');
    fprintf('          CompleteSearch: [ on  | {off} ]\n\n');
    
    fprintf('                 Display: [ off | iter | diagnose | {final} ]\n');
    fprintf('              OutputFcns: [ function_handle | @psoutputhistory | {[]} ]\n\n');
    fprintf('                PlotFcns: [ function_handle | @psplotbestf |\n');
    fprintf('                            @psplotmeshsize | @psplotfuncount |\n');
    fprintf('                            @psplotbestx    | {[]} ]\n');
    fprintf('            PlotInterval: [ positive scalar | {1} ]\n\n');
    
    fprintf('                   Cache: [ on  | {off} ]\n');
    fprintf('               CacheSize: [ positive scalar | {1e4} ]\n');
    fprintf('                CacheTol: [ positive scalar | {eps} ]\n\n');
    
    fprintf('              Vectorized: [ on | {off} ]\n\n');
    
    fprintf('             UseParallel: [ ''always'' | {''never''} ]\n');
    return; 
end     

numberargs = nargin; 

%Return options with default values and return it when called with one
%output argument
options= struct('TolMesh', [], ...
                'TolCon', [], ...
                'TolX', [] , ...
                'TolFun',[] , ...
                'TolBind',[], ...
                'MaxIter', [], ...
                'MaxFunEvals', [], ...
                'TimeLimit', [], ...
                'MeshContraction', [], ... 
                'MeshExpansion', [], ...
                'MeshAccelerator',[], ... 
                'MeshRotate',[], ...
                'InitialMeshSize', [], ...
                'ScaleMesh', [], ...
                'MaxMeshSize', [], ...
                'InitialPenalty',[], ...
                'PenaltyFactor',[], ...
                'PollMethod', [], ...
                'CompletePoll',[], ...
                'PollingOrder', [], ...
                'SearchMethod', [], ...
                'CompleteSearch',[], ...
                'Display', [], ...
                'OutputFcns', [], ... 
                'PlotFcns',[]', ...
                'PlotInterval', [], ...
                'Cache', [], ...
                'CacheSize',[], ... 
                'CacheTol',[], ...
                'Vectorized',[], ...
                'UseParallel', [] ...
               ); 

% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname)
            msg = sprintf( ...
                'No default options available: the function ''%s'' does not exist on the path.',funcname);
            error('globaloptim:PSOPTIMSET:functionNotFound',msg)
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
        optionsfcn = feval(varargin{1},'defaults');
    catch
        msg = sprintf( ...
            'No default options available for the function ''%s''.',funcname);
        error('globaloptim:PSOPTIMSET:noDefaultOptions',msg)
    end
    % To get output, run the rest of psoptimset as if called with psoptimset(options, optionsfcn)
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = fieldnames(options);
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error('globaloptim:PSOPTIMSET:invalidArgument',['Expected argument %d to be a string parameter name ' ...
                    'or an options structure\ncreated with PSOPTIMSET.'], i);
        end
        argFieldnames = fieldnames(arg);
        if strmatch('MaxIteration',argFieldnames,'exact')
            val = arg.MaxIteration;
            rmfield(arg,'MaxIteration');
            arg.MaxIter = val;
            warning('globaloptim:PSOPTIMSET:obsoleteProperty','Options ''MaxIteration'' is obsolete; use ''MaxIter'' instead.');
        end
        
        for j = 1:m
            if any(strcmp(argFieldnames,Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                [valid, errmsg] = checkfield(Names{j,:},val);
                if valid
                    options.(Names{j,:}) = val;
                else
                    error('globaloptim:PSOPTIMSET:invalidOptionField',errmsg);
                end
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('globaloptim:PSOPTIMSET:invalidArgPair','Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error('globaloptim:PSOPTIMSET:invalidArgFormat','Expected argument %d to be a string parameter name.', i);
        end
        
        lowArg = lower(arg);
        if strcmp(lowArg,'maxiteration')
            arg = 'maxiter';
            lowArg = lower(arg);
            warning('globaloptim:PSOPTIMSET:obsoleteProperty','Options ''MaxIteration'' is obsolete; use ''MaxIter'' instead.');
        end
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error('globaloptim:PSOPTIMSET:invalidParamName','Unrecognized parameter name ''%s''.', arg);
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ', ' Names{k,:}];
                end
                msg = sprintf('%s).', msg);
                error('globaloptim:PSOPTIMSET:AmbiguousParamName',msg);
            end
        end
        expectval = 1;                      % we expect a value next
        
    else           
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        [valid, errmsg] = checkfield(Names{j,:},arg);
        if valid
            options.(Names{j,:}) = arg;
        else
            error('globaloptim:PSOPTIMSET:invalidParamVal',errmsg);
        end
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error('globaloptim:PSOPTIMSET:invalidParamVal','Expected value for parameter ''%s''.', arg);
end


%-------------------------------------------------
function [valid, errmsg] = checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   [VALID, MSG] = CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%

valid = 1;
errmsg = '';
% empty matrix is always valid
if isempty(value)
    return
end

switch field
    case {'SearchMethod'}
        if iscell(value)  || isa(value,'function_handle')
        elseif isa(value,'char') && any(strcmp(value,{'gpspositivebasisnp1','gpspositivebasis2n', ...
                'positivebasisnp1','positivebasis2n','madspositivebasisnp1','madspositivebasis2n', ...
                'gsspositivebasisnp1','gsspositivebasis2n','searchlhs','searchga','searchneldermead'}))
        else
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''@GPSPositiveBasisNp1'',''@GPSPositiveBasis2N'',''@MADSPositiveBasisNp1'',''@MADSPositiveBasis2N'',''@searchlhs'',\n''@GSSPositiveBasisNp1'',''@GSSPositiveBasis2N'',''@searchga'',''@searchneldermead'' or a function_handle.',field); 
        end
        
    case {'PollMethod'}  
        if ~isa(value,'char') || ~any(strcmp(value,{'gpspositivebasisnp1', 'gpspositivebasis2n' ...
                'positivebasisnp1', 'positivebasis2n','madspositivebasisnp1', 'madspositivebasis2n', ...
                'gsspositivebasisnp1', 'gsspositivebasis2n'}))
            valid=0;
            if length(strmatch(value,{'gpspositivebasisnp1','gpspositivebasis2n', ...
                    'positivebasisnp1','positivebasis2n','madspositivebasisnp1','madspositivebasis2n', ...
                    'gsspositivebasisnp1', 'gsspositivebasis2n'})) >1
                valid =0;
                errmsg = sprintf('Ambiguous value for OPTIONS parameter %s: No match found''.',field); 
            else
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''GSSPositiveBasis2N'',''GSSPositiveBasisNp1'',''GPSPositiveBasis2N'',''GPSPositiveBasisNp1'',''MADSPositiveBasis2N'', or ,''MADSPositiveBasisNp1''.',field);
            end
        end
        
    case {'PollingOrder'}
        if ~isa(value,'char') || ~any(strcmp(value,{'random','success','consecutive'}))
            valid=0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''Random'',''Success'', or ,''Consecutive''.',field);
        end
        
    case {'Display'}
        if ~isa(value,'char') || ~any(strcmpi(value,{'off','iter','none', ...
                    'diagnose','final'}))
            valid=0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'',''iter'', ''diagnose'',''final'', or ,''none''.',field);
        end
        
    case {'OutputFcns','PlotFcns'}
        if iscell(value) ||  isa(value,'function_handle')
            valid = 1;
        else
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a function handle or cell array of function handles.',field);
        end
        
    case {'InitialMeshSize','MeshContraction','MeshExpansion','MaxMeshSize', ...
                'CacheTol','CacheSize','TimeLimit' ...
                'InitialPenalty','PenaltyFactor'} 
        if ~(isa(value,'double') && value > 0) 
            valid = 0;
            if ischar(value)
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar (not a string).',field);
            else
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar.',field);
            end
        end
    case {'TolMesh','TolBind','TolX','TolFun','TolCon', }
        if ~(isa(value,'double') && value >= 0) 
            valid = 0;
            if ischar(value)
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar (not a string).',field);
            else
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
            end
        end
        
    case {'MaxIter'} % integer including inf or default string
        if ~(isa(value,'double') && value >= 0) ...
                && ~isequal(value, '100*numberofvariables') ...
                && ~isequal(value, '200*numberofvariables')
            valid = 0;
            if ischar(value)
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive numeric (not a string).',field);
            else
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive numeric.',field);
            end
        end
        
    case {'MaxFunEvals'} % integer including inf or default string
        if ~(isa(value,'double') && value >= 0) ...
                && ~isequal(value,'1000*numberofvariables') ...
                && ~isequal(value, '2000*numberofvariables') 
            valid = 0;
            if ischar(value)
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive numeric (not a string).',field);
            else
                errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive numeric.',field);
            end
        end 
        
    case {'CompletePoll','CompleteSearch','MeshAccelerator','MeshRotate','Vectorized','Cache'}  %off, on
        if ~isa(value,'char') || ~any(strcmp(value,{'on','off'}))
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
        end
    case 'ScaleMesh'  %off, on, dynamic
        if ~isa(value,'char') || ~any(strcmp(value,{'on','off','dynamic'}))
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
        end
    case 'PlotInterval'
        valid =  isnumeric(value) && (1 == length(value)) && (value > 0) && (value == floor(value));
        if(~valid)
            errmsg = sprintf('The field ''%s'' must be a positive integer.',field);
        end
     
    case 'UseParallel'
        if ischar(value)
            valid = 1;
        else
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''never'' or ''always''.',field);
        end
    
    otherwise
        error('globaloptim:PSOPTIMSET:unknownOptionsField','Unknown field name for Options structure.')
end    


