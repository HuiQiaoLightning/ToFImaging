function [ u, p_u, p_u_g, val ] = ls_nonlinear_pos( M, b, nnz_pos, nnz_val, beta_pos, sigma, rho, i_prev )
   
% ***LS_NONLINEAR Non-linear gaussian/exponential fitting***
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%

% Copyright (C) 2013. Felix Heide
% Email: fheide@cs.ubc.ca

    %Three parameters in each position (gaussian amplitude, exponential ampitude, decay time)
    
    %%%%Initial point
    u_0 = zeros( size(nnz_pos,1), 4);
    u_0(:,3) = 0.05;
    u_0(:,4) = nnz_pos; %Initialize with previous pos
        
    %%%%Bounds
    lb = zeros(size(reshape(u_0, [], 4)));
    ub = ones(size(reshape(u_0, [], 4))) * 100; %Not bounded
    
    %Amplitude
    for i = 1:size(nnz_pos,1)
        u_0(i,1:2) = min(max(nnz_val(i,1), 0),1);
    end
    
    %Vectorize
    u_0 = u_0(:);
    
    %Decay
    %lb(:,3) = 0.01;
    %ub(:,3) = 0.1;
    
    lb(:,3) = 0.001;
    ub(:,3) = 0.5;
    
    %Position
    lb(:,4) = 1;
    ub(:,4) = size(M,2);
 
    %%% HYBRID SEARCH
    tic;
    funObj = @(x) minimzation_obj_handle(x, M, b, nnz_pos, beta_pos, sigma, rho, i_prev );
    pattern_opts = psoptimset_gads( 'SearchMethod', {@MADSPositiveBasis2N}, ...
                                'MaxFunEvals', 60000, 'MaxIter', 1000, 'PollingOrder', 'Random', ... %'CompletePoll', 'on'...
                                'Display', 'off');
    u_t = patternsearch_gads(funObj, u_0, [],[],[],[],lb(:),ub(:),[],pattern_opts);
    
    options = [];
    options.MaxIter = 500;
    options.verbose = 0;
    funProj = @(x)boundProject(x,lb(:),ub(:));
    funObj = @(x) minimzation_obj_handle(x, M, b, nnz_pos, beta_pos, sigma, rho, i_prev );
    u = minConf_SPG(funObj,u_t,funProj,options);

    f_u_hybrid= funObj(u);
    
    hybrid_time = toc;
    
    %Function value
    val= funObj(u);
    
    %Reshape
    u = reshape(u, [], 4);
    
    
    %Optionally reconstruct
    if nargout > 1

        %Dimensions
        t_dim = size(M,2);

        %Generate current reconstruction vector
        p_u = zeros(1, t_dim);
        p_u_g = zeros(1, t_dim);

        %Iterate over positions
        for i = 1:size(nnz_pos,1)
            %Spatial position
            pos = u(i,4);
            pos = max(min(pos, t_dim),1); %Clamp

            %Gaussians
            p_u = p_u + u(i, 1) * exp(-0.5 * (([1:t_dim] - pos)/sigma).^2);
            p_u_g = p_u_g + u(i, 1) * exp(-0.5 * (([1:t_dim] - pos)/sigma).^2); %Just the gaussians

            %Exponentials
            exponential = u(i, 2) * exp( -u(i, 3) * ([1:t_dim] - pos) );
            exponential(1:floor(pos)) = 0;
            p_u = p_u + exponential;
        end
    end
    
return;

function [f, df, H] = minimzation_obj_handle(u, M, b, nnz_pos, beta_pos, sigma, rho, i_prev )

    % This function returns the function value, partial derivatives.

    f = obj_fun( u, M, b, nnz_pos, beta_pos, sigma, rho, i_prev );

    if nargout > 1
      df = grad_fun( u, M, b, nnz_pos, beta_pos, sigma, rho, i_prev );
    end
    
    if nargout > 2
      H = [];
    end

return;

function [ f_u ] = obj_fun( u, M, b, nnz_pos, beta_pos, sigma, rho, i_prev )

    %Vector reshape
    u = reshape(u, [], 4);

    %Dimensions
    t_dim = size(M,2);

    %Generate current reconstruction vector
    p_u = zeros(1, t_dim);
    
    %Iterate over positions
    for i = 1:size(nnz_pos,1)
        
        %Spatial position
        pos = u(i,4);
        pos = max(min(pos, t_dim),1); %Clamp
        
        %Gaussians
        p_u = p_u + u(i, 1) * exp(-0.5 * (([1:t_dim] - pos)/sigma).^2);
        
        %Exponentials
        exponential = u(i, 2) * exp( -u(i, 3) * ([1:t_dim] - pos) );
        exponential(1:floor(pos)) = 0;
        p_u = p_u + exponential;

    end
    
    f_u = norm( M * p_u' - b, 2).^2 + beta_pos * norm( u(:,4) - nnz_pos, 2 ).^2 + rho * norm( p_u' - i_prev, 2).^2 ;
    
return;

function [ g_u ] = grad_fun( u, M, b, nnz_pos, beta_pos, sigma, rho, i_prev )

    %Vector reshape
    u = reshape(u, [], 4);

    %Dimensions
    t_dim = size(M,2);

    %Generate current reconstruction vector
    p_u = zeros(1, t_dim);
    
    %Iterate over positions
    for i = 1:size(nnz_pos,1)
        
        %Spatial position
        pos = u(i,4);
        pos = max(min(pos, t_dim),1); %Clamp
        
        %Gaussians
        p_u = p_u + u(i, 1) * exp(-0.5 * (([1:t_dim] - pos)/sigma).^2);
        
        %Exponentials
        exponential = u(i, 2) * exp( -u(i, 3) * ([1:t_dim] - pos) );
        exponential(1:floor(pos)) = 0;
        p_u = p_u + exponential;

    end
    
    %Right hand term
    rh_term = M' * (M * p_u' - b);
    rh_prev_term = (p_u' - i_prev);
    
    %Compute derivatives
    g_u = zeros( size(u) );
    
    %Iterate over positions
    for i = 1:size(nnz_pos,1)
        
        %Spatial position
        pos = u(i,4);
        pos = max(min(pos, t_dim),1); %Clamp
                
        %Gaussian amplitude
        gauss_amp_d = exp(-0.5 * (([1:t_dim] - pos)/sigma).^2);
        g_u(i,1) = gauss_amp_d * rh_term + rho * gauss_amp_d * rh_prev_term;
        
        %Exponential amplitude
        exponential_amp_d = exp( -u(i, 3) * ([1:t_dim] - pos) );
        exponential_amp_d(1:floor(pos)) = 0;
        g_u(i,2) = exponential_amp_d * rh_term + rho *exponential_amp_d * rh_prev_term;
        
        %Exponential decay
        exponential_decay_d = u(i, 2) * exp( -u(i, 3) * ([1:t_dim] - pos) ) .* (-([1:t_dim] - pos));
        exponential_decay_d(1:floor(pos)) = 0;
        g_u(i,3) = exponential_decay_d * rh_term + rho * exponential_decay_d * rh_prev_term;
        
        %Gaussian and exponential position
        gaussian_pos_d = u(i, 1) * exp(-0.5 * (([1:t_dim] - pos)/sigma).^2) .* ( ([1:t_dim] - pos)/(sigma.^2) );
        
        exponential_pos_d = u(i, 2) * exp( -u(i, 3) * ([1:t_dim] - pos) ) .* u(i, 3);
        exponential_pos_d(1:floor(pos)) = 0;
        
        pos_d = gaussian_pos_d + exponential_pos_d;
        g_u(i,4) = pos_d * rh_term + rho * pos_d * rh_prev_term;

    end
    
    %Muliplicator from squared norm
    g_u = 2 * g_u;
    
    %Compute derivative component for positions
    g_u(:,4) = g_u(:,4) + 2 * beta_pos * ( u(:,4) -  nnz_pos );
    
    %Vectorize
    g_u = g_u(:);
    
return;