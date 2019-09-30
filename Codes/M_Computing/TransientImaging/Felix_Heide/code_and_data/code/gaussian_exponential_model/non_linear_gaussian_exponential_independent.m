function [ f_u, u_res ] = non_linear_gaussian_exponential_independent(B, M, nnz_positions, ymax_all, sigma, beta_pos, rho, I_prev, verbose)
   
% ***LS_NONLINEAR Non-linear gaussian/exponential fitting with variable positions***
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%

% Copyright (C) 2013. Felix Heide
% Email: fheide@cs.ubc.ca
    
    %Individual derivatives are only depending on single time pixel. So
    %simply iterate over all pixels and solve the prox
    
    %We have here u
    f_u = zeros(size(nnz_positions,1), size(nnz_positions,2), size(M,2) ); 
    
    u_res = cell(size(nnz_positions,1), size(nnz_positions,2)); 

    for y = 1:size(nnz_positions,1)
        for x = 1:size(nnz_positions,2)
            
            %Debug
            fprintf('Gauss/exp fit for %d, %d \n', x, y)
         
            %Get local params
            nnz_pos_curr = nnz_positions{y, x};
            ymax_curr = ymax_all{y, x};
            b_curr = reshape(B(y,x,:),[],1);
            i_prev_curr = reshape(I_prev(y,x,:),[],1);
            
            if isempty( nnz_pos_curr )
                continue;
            end
            
            %Compute 
            tic;

            [ u, img_rec, ~, val] = ls_nonlinear_pos( M, b_curr, nnz_pos_curr, ymax_curr, beta_pos, sigma, rho, i_prev_curr );
            
            u_res{y, x} =  u;
            
            rectime = toc;
            
            %Function value output
            fprintf('--> val %3.3f in time %3.2f sec\n', val, rectime)

            %Save
            f_u(y, x,:) = img_rec';            
        end
    end

return;