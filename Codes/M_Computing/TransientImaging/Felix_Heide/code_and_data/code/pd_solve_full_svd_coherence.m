function [ res ] = pd_solve_full_svd_coherence(B, M, lambda_residual, lambda, lambda_spatial, thresh_huber, max_it, tol, verbose)

% ***Forward-backward splitting algorithm for solving the i-step***
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%
% This algorithm is based on:
% CHAMBOLLE , A., AND POCK , T. 2011. A first-order primal-dual algorithm for
% convex problems with applications to imaging. J. Math. Imaging Vis. 40 , 120–14
%
% which is essentially a (preconditioned) ADMM method. Stephen Boyd has
% published a nice manuscript on ADMM: Distributed Optimization and Statistical 
% Learning via the Alternating Direction Method of Multipliers
% S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein
% Foundations and Trends in Machine Learning, 3(1):1–122, 2011. 
 
% Copyright (C) 2013. Felix Heide
% Email: fheide@cs.ubc.ca

    
    %Shortcut for the l1 norm.
    Amplitude = @(u)sqrt(u.^2);
    F = @(u)sum(sum(Amplitude(u)));
    
    %Convert dims
    n_meas = size(B, 3);
    B_vec = reshape( shiftdim(B,2), n_meas, [] );
    
    %Precompute constants
    MtB = M' * B_vec;
    MtM = M' * M;
    
    %Minimization weights
    %PRECOMPUTE THE VALUES OF L IN PRACTICE FOR FAST CONVERGENCE (No lookup table given here for completeness)
    if strcmp(verbose, 'brief') || strcmp( verbose, 'all')
            fprintf('Computing linear operator norm... ')
    end
    L = compute_operator_norm(@(x)Kmult(x, lambda, lambda_spatial), ...
                              @(x)KSmult(x, lambda, lambda_spatial), [size(B,1), size(B,2), size(M,2)]);
    if strcmp(verbose, 'brief') || strcmp( verbose, 'all')
            fprintf(' done! Starting iterations.\n')
    end
    
    sigma = 10.0;
    tau = .7/(sigma*L^2);
    theta = 1.0;
    
    %Precompute prox inversion
    A = tau * lambda_residual * MtM + eye(size(MtM));
    Apinv = pinv(A);
    
    %Prox operators
    %L1
    %ProxFS = @(u,s) u./ max(1, Amplitude(u));
    
    %Huber
    ProxFS = @(u,s) (u/(1 + thresh_huber * s)) ./ max(1, Amplitude(u/(1 + thresh_huber * s)));
    ProxG = @(f,tau,l) solve_SVD(MtB, MtM, Apinv, f, tau, l );
    %ProxG = @(f,tau,l) solve_quadprog(MtB, MtM, f, tau, l );
    
    
    %Set initial iterate
    f = zeros(size(B,1), size(B,2), size(M,2)); %Simply start with 0-vector
    g = Kmult(f, lambda, lambda_spatial);
    f1 = f;
    
    %Example of one iterations.
    for i = 1:max_it
        
        fold = f;
        g = ProxFS( g + sigma * Kmult(f1, lambda, lambda_spatial), sigma);
        f = ProxG( f - tau * KSmult(g, lambda, lambda_spatial), tau, lambda_residual);
        f1 = f + theta * (f-fold);

        diff = f - fold;
        diff_norm = norm(diff(:), 2);
        if norm(f(:), 2) > eps()
            diff_norm = diff_norm / norm(f(:), 2);
        end
        if strcmp(verbose, 'brief') || strcmp( verbose, 'all')
            fprintf('iter %d, diff %5.5g\n', i, diff_norm)
        end
        if diff_norm < tol
            break;
        end
    end
    
    res = f1;
 
return;

function x = solve_SVD(MtB, MtM, Apinv, f, tau, l )
    %Solves Ax = b with
    % A = (tau*lambda* M'* M + eye ) and b = tau * lambda * M' * B + f
    % Solve quadratic program here if you want constraints on x
    
    %Compute b vector
    n = size(f, 3);
    f_vec = reshape( shiftdim(f,2), n, [] );
    
    %Compute b
    b = tau * l * MtB + f_vec;

    %Solve
    x = Apinv*b;
    
    %Reshape
    x = shiftdim( reshape(x, [n, size(f,1), size(f,2)]), 1 );

return;

function Kmultf = Kmult(f, lambda, lambda_sp )

%1) Gradient in time
dtf=[-1 1];
dtf_t = fliplr(flipud(dtf));
d_t = dtf_t(1,2)* f(:,:,[1 1:end]) + dtf_t(1,1)*f(:,:,[1:end end]);
d_t = lambda * d_t(:,:, 2:end);

%2) Gradients
d_x = zeros( size(f) );
d_y = zeros( size(f) );

%Iterate over time
for t = 1:size(f,3)

    % derivative filters
    dxf=[-1 1];
    dyf=[-1;1];

    %Compute tv terms
    fx = imconv( f(:,:,t), fliplr(flipud(dxf)), 'full');
    d_x(:,:,t) = lambda_sp * fx(:, 2:end);

    fy = imconv( f(:,:,t), fliplr(flipud(dyf)), 'full');
    d_y(:,:,t) = lambda_sp * fy(2:end, :);

end

%3) Put all together (in vectorized form)
Kmultf = cat(4, d_t, d_x, d_y);

return;


function KSmultf = KSmult(f, lambda, lambda_sp )

%1) Gradient in time
dtf=[-1 1];
f_g = f(:,:,:,1);
d_t = dtf(1,2) * lambda * f_g(:,:,[1 1:end]) + dtf(1,1) * lambda * f_g(:,:,[1:end end]);
d_t = d_t(:,:, 1:end-1);

%2) Gradients
f_x = f(:,:,:,2);
f_y = f(:,:,:,3);

d_x = zeros( size(f_x) );
d_y = zeros( size(f_y) );

%Iterate over time
for t = 1:size(f,3)
    
    %1) Gradient 
    % derivative filters
    dxf=[-1 1];
    dyf=[-1;1];

    %Compute tv terms
    fx = imconv( lambda_sp * f_x(:,:,t), dxf, 'full');
    d_x(:,:,t) = fx(:, 1:end-1);

    fy = imconv( lambda_sp * f_y(:,:,t), dyf, 'full');
    d_y(:,:,t) = fy(1:end-1,:);

end

% gather result
KSmultf = d_t + d_y + d_x;

return;

function F_filt = imconv(F,K,output)

%Convolution
%General: F_filt = imfilter(F, K, 'full', 'conv', 'replicate');

%Speedup for small two entry kernels (full)
if size(K,1) == 1 && size(K,2) == 2 && strcmp('full', output)
    F_filt = K(1,2)* F(:,[1 1:end],:) + K(1,1)*F(:,[1:end end],:);
elseif size(K,1) == 2 && size(K,2) == 1 && strcmp('full', output)
    F_filt = K(2,1)* F([1 1:end],:,:) + K(1,1)*F([1:end end],:,:);
else
    %General model
    F_filt = imfilter(F, K, output, 'conv', 'replicate');
end

return;

