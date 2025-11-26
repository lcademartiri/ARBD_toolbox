function Draw = build_noise_library(std,nodisp) 
% S : simulation struct (will be updated with noise parameters)
% Draw : Mx3 matrix of precomputed displacement vectors (each row a sample)
% Result: stores S.noise.{mean,whiten,avgvar,M}
    Draw=normrnd(0,std,nodisp,3);
    D = double(Draw);          % ensure double precision
    M = size(D,1);
    if size(D,2) ~= 3
        error('Noise library must be M x 3');
    end

    % 1) exact zero-mean
    mu = mean(D,1);            % 1x3
    D0 = D - mu;               % center

    % 2) covariance -> compute empirical covariance (3x3)
    C = (D0' * D0) / (M - 1);  % same as cov(D0) but slightly faster

    % 3) pick target isotropic variance = average of diagonal of C
    avgvar = trace(C) / 3;     % scalar

    % 4) compute symmetric whitening/rescaling transform T such that:
    %    T * C * T' = avgvar * I
    % Use eigendecomposition (stable for small 3x3)
    [V, L] = eig((C + C')/2);  % ensure symmetric
    lambda = diag(L);
    % Protect against tiny negative eigenvalues due to numerical noise:
    lambda(lambda <= 0) = eps;
    scale = sqrt(avgvar ./ lambda);
    T = V * diag(scale) * V';  % 3x3 transform

    % 5) apply transform to produce final library
    Draw = (T * D0')';      % M x 3
end
