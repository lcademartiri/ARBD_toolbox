function p_wrapped = mic_wrap_positions(p, S)
% MIC_WRAP_POSITIONS Wraps particle coordinates into the primary simulation cell.
%
%   p_wrapped = mic_wrap_positions(p, S)
%
%   Input:
%       p - (N x 3) matrix of particle positions
%       S - Simulation structure (contains S.bc, S.br, S.fcc matrices)
%
%   Output:
%       p_wrapped - (N x 3) wrapped positions

    p_wrapped = p;

    if S.bc == 2  % --- Cubic PBC ---
        L = 2 * S.br;
        % Wrap to interval [-L/2, L/2]
        % This keeps the cluster centered at the origin
        p_wrapped = p - L .* round(p ./ L);
        
    elseif S.bc == 3  % --- FCC PBC ---
        % 1. Transform to Fractional Coordinates (Basis of lattice vectors)
        % f = p * invA'
        f = (S.fcc.invA * p.').';
        
        % 2. Wrap Fractional Coordinates
        % Use 'round' to wrap to [-0.5, 0.5] (Centered Unit Cell)
        % Use 'floor' to wrap to [0, 1] (Positive Quadrant Unit Cell)
        % Since your simulation seems centered at 0, 'round' is safer.
        f_wrapped = f - round(f);
        
        % 3. Transform back to Cartesian
        % p = f * A'
        p_wrapped = (S.fcc.A * f_wrapped.').';
        
    elseif S.bc == 4  % --- Big Box (Teleport Check) ---
        % BB doesn't wrap continuously, but we check for hard limits.
        % (Usually handled by bigboxTeleport, but strictly speaking:
        % particles outside BB should be removed or reflected, not wrapped).
        % This function does nothing for BB or SBC.
    end
end