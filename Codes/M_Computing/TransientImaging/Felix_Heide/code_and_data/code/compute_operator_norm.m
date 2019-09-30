function L = compute_operator_norm(A, AS, vec_size)

% ***Operator norm***
% computes the operator norm for a linear operator AS on images with size sx, 
% which is the square root of the largest eigenvector of AS*A.
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%

    %Compute largest eigenvalue
    opts.tol = 1.0e-3;
    lambda_largest = eigs(@(x)ASAfun(x, A, AS, vec_size), prod(vec_size(:)), 1,'lm', opts);
    L = sqrt(lambda_largest);

return;

function ASAx = ASAfun(x, A, AS,vec_size)
    x_img = reshape(x,vec_size);
    ASAx = AS(A(x_img));
    ASAx = ASAx(:);
return;