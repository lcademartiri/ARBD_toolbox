function Draw = build_noise_library(dims, std) 
% BUILD_NOISE_LIBRARY Generates a whitened, zero-mean noise library.
%
% Inputs:
%   dims : Vector of target dimensions, e.g., [1000, 10, 50]
%          dims(1) is the number of samples (rows).
%          The product of dims(2:end) defines the columns to decorrelate.
%   std  : Target standard deviation for the output
%
% Output:
%   Draw : Array of specified 'dims' with fully decorrelated columns.

    % 1) Parse dimensions and flatten trailing dimensions
    M = dims(1);                     % Number of samples (rows)
    total_cols = prod(dims(2:end));  % e.g., 10 * 50 = 500 columns

    % 2) Generate random values directly in a 2D matrix for processing
    D = normrnd(0, std, M, total_cols);
    D = double(D);                   % Ensure double precision

    % 3) Exact zero-mean centering
    mu = mean(D, 1);                 % 1 x total_cols
    D0 = D - mu;                     % Center the data

    % 4) Compute empirical covariance matrix (total_cols x total_cols)
    C = (D0' * D0) / (M - 1);        

    % 5) Pick target isotropic variance = average of diagonal of C
    avgvar = trace(C) / total_cols;  

    % 6) Compute symmetric whitening/rescaling transform T
    % Use eigendecomposition (stable for a few hundred/thousand dimensions)
    [V, L] = eig((C + C') / 2);      % Force strict symmetry
    lambda = diag(L);
    
    % Protect against tiny or negative eigenvalues due to numerical noise
    lambda(lambda <= 0) = eps;
    scale = sqrt(avgvar ./ lambda);
    T = V * diag(scale) * V';         % total_cols x total_cols transform

    % 7) Apply transform 
    % Note: (T * D0')' simplifies mathematically to (D0 * T) because T is symmetric
    Draw_flat = D0 * T;              % M x total_cols

    % 8) Reshape the flat 2D matrix back to the requested N-D array dimensions
    Draw = reshape(Draw_flat, dims);
end