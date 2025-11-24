function disp_vecs = mic_all_pair_displacements(p, S)
% MIC_ALL_PAIR_DISPLACEMENTS Computes all pairwise MIC displacement vectors.
%
%   disp_vecs = mic_all_pair_displacements(p, S)
%
%   Input:
%       p - (N x 3) particle positions
%       S - System struct with S.fcc.A and S.fcc.invA initialized
%           (Works for both Cubic and FCC)
%
%   Output:
%       disp_vecs - (M x 3) matrix of displacement vectors, where M = N*(N-1).
%                   Contains p(i) - p(j) for all i != j.
%                   Format matches the reshaped PDFD output (long column list).

    N = size(p, 1);
    
    % 1. Generate Indices for all i != j
    % We use meshgrid to generate the full permutation N x N
    [J_grid, I_grid] = meshgrid(1:N, 1:N);
    
    % Create a mask to exclude the diagonal (self-interactions)
    mask = I_grid ~= J_grid;
    
    % Extract valid indices as column vectors
    I = I_grid(mask);
    J = J_grid(mask);
    
    % 2. Coordinate Transformation (Cartesian -> Fractional)
    % Works for both Cubic (invA is 1/L diag) and FCC (invA is lattice inverse)
    f = p * S.fcc.invA';
    
    % 3. Calculate Pairwise Difference in Fractional Space
    % df = f_i - f_j
    df = f(I, :) - f(J, :);
    
    % 4. Apply Minimum Image Convention (MIC)
    % In fractional space, the nearest image is always within [-0.5, 0.5].
    % We simply subtract the nearest integer.
    df = df - round(df);
    
    % 5. Transform Back to Cartesian Space
    % dr = df * A^T
    disp_vecs = df * S.fcc.A';

end