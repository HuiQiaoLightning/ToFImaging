function     [MeshSize,Iterate,X,optimState] = updateparam(successPoll,successSearch, ...
        MeshSize,nextIterate,Iterate,X,optimState,options)
%UPDATEPARAM Updates the parameters of pattern search for next iteration
% 	
% 	SUCCESSPOLL,SUCCESSSEARCH: Flag indicating if poll or search was
% 	successful in last iteration
% 	
% 	MeshAccelerator: Used for fast mesh contraction when close to minima
% 	(might loose accuracy)
% 	
% 	MAXMESHSIZE,MINMESH: Maximum and minimum mesh size which can be used
% 	
% 	MeshExpansion,MeshContraction: These factors(scalar) are used to coarsen
% 	or refining the mesh
% 	
% 	MESHSIZE: Current mesh size used.
% 	
% 	NEXTITERATE,ITERATE: Iterates in last iteration and current iteration
                 
%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:27:28 $


factor = 1e-2;
MeshExpansion = options.MeshExpansion;
Contract = optimState.MeshCont;
MaxMeshSize = options.MaxMeshSize;
AdaptiveMesh = any(strcmpi(options.PollMethod,{'madspositivebasisnp1','madspositivebasis2n'}));
% Override these options for MADS
if AdaptiveMesh
    MeshExpansion = 4;
    Contract = 0.25;
    MaxMeshSize = 1;
end

if successPoll
    MeshSize = min(MaxMeshSize,MeshExpansion*MeshSize);
else 
    if  (nextIterate.f < 0 && ~isfinite(nextIterate.f))
        optimState.infMessage ='Function has reached -Inf value';
    elseif ~successSearch
        optimState.how = 'Refine Mesh';
        MeshSize = Contract*MeshSize;
    end
end

% MeshAccelerator step; Speed up convergence when near the optimizer. The
% new value of MeshContraction will be used in the next successive
% unsuccessful iteration
if strcmpi(options.MeshAccelerator,'on')
    if abs(MeshSize) < min(factor,options.TolMesh*(1/factor))
        % MeshSize is less than factor, we can use a fast convergence
        optimState.MeshCont  =  max(factor,optimState.MeshCont/2);
    else
         % MeshSize is greater than factor, reset the contraction factor
         optimState.MeshCont = options.MeshContraction;
    end
end

%How did last iteration go
if successSearch
    optimState.how = 'Successful Search';
elseif successPoll
    optimState.how = 'Successful Poll';
end

if strcmpi(options.MeshRotate,'on') && MeshSize < options.TolMesh*1/factor && ~successPoll
    optimState.scale = -optimState.scale ;
    optimState.how = [optimState.how,'\Rotate'];
end

%Other termination criteria
optimState.deltaX    = norm(nextIterate.x-Iterate.x);
optimState.deltaF   = abs(nextIterate.f-Iterate.f);
Iterate  = nextIterate;
optimState.Iter     = optimState.Iter + 1;
X(:)     = Iterate.x;

