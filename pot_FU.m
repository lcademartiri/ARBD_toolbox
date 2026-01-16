function [H, H_interpolant, U_interpolant_trunc, U_interpolant_full] = pot_FU(potentialtype, cutoff, table_limit, points, eqdistance, welldepth, kappa)
% pot_FU
% Generates force and TWO energy interpolants (Metric vs Physical).
%
% Inputs:
%   potentialtype: 1=LJ, 2=WCA, 3=Quadratic, 4=Hertzian, 5=Yukawa
%   cutoff:        Simulation cutoff (Forces and Metric U drop to 0 here)
%   table_limit:   Physical limit (e.g., 50*sigma) for the Full U interpolant
%   points:        Grid points
%   eqdistance:    sigma
%   welldepth:     epsilon
%   kappa:         (Optional) Yukawa screening
%
% Outputs:
%   H:                   Table [r, Force_trunc, Energy_trunc] (for inspection)
%   H_interpolant:       Force F(r) (Truncated at cutoff)
%   U_interpolant_trunc: Metric U(r) (Truncated at cutoff)
%   U_interpolant_full:  Physical U(r) (extends to table_limit)

    if nargin < 7
        kappa = 1 / eqdistance; 
    end

    % 1. Generate Grid to the "Physical" Limit (table_limit)
    x = linspace(0, table_limit, points)';
    sigma = eqdistance;
    xtilde = x / sigma;
    
    % Initialize arrays (Full Range)
    F_full = zeros(points, 1);
    U_full = zeros(points, 1);

    % --- 2. Calculate Full Potentials ---
    if potentialtype == 1 % Lennard-Jones
        F_full = (welldepth/sigma) .* 24 .* (2./xtilde.^13 - 1./xtilde.^7);
        U_full = 4 * welldepth .* (xtilde.^-12 - xtilde.^-6);
        
    elseif potentialtype == 2 % WCA
        r_min = 2^(1/6) * sigma;
        F_full = (welldepth/sigma) .* 24 .* (2./xtilde.^13 - 1./xtilde.^7);
        U_full = 4 * welldepth .* (xtilde.^-12 - xtilde.^-6) + welldepth;
        
        % WCA is naturally truncated, so Full = Truncated
        mask = x >= r_min;
        F_full(mask) = 0; U_full(mask) = 0;
        
    elseif potentialtype == 3 % Quadratic
        U_full = 0.5 * welldepth .* (1 - x./sigma).^2;
        dx = x(2)-x(1);
        F_full = -gradient(U_full, dx);
        mask = x >= sigma;
        F_full(mask) = 0; U_full(mask) = 0;
        
    elseif potentialtype == 4 % Hertzian
        U_full = welldepth .* (1 - x./sigma).^(5/2);
        dx = x(2)-x(1);
        F_full = -gradient(U_full, dx);
        mask = x >= sigma;
        F_full(mask) = 0; U_full(mask) = 0;

    elseif potentialtype == 5 % Yukawa
        r_safe = x; r_safe(r_safe==0) = eps;
        term_exp = exp(-kappa .* r_safe);
        term_inv = 1 ./ r_safe;
        U_full = welldepth .* sigma .* term_inv .* term_exp;
        F_full = welldepth .* sigma .* term_exp .* (term_inv.^2 + kappa .* term_inv);
    end

    % Handle singularity at r=0 (common for LJ/Yukawa)
    if potentialtype == 1 || potentialtype == 5
        % Just remove the first point if it's 0 or very close to singularity
        if x(1) < 1e-9
            x(1) = []; F_full(1) = []; U_full(1) = [];
        end
    end
    
    % --- 3. Create TRUNCATED versions (The "Internal Metric") ---
    % Everything beyond 'cutoff' becomes 0.
    F_trunc = F_full;
    U_trunc = U_full;
    
    mask_cutoff = x > cutoff;
    F_trunc(mask_cutoff) = 0;
    U_trunc(mask_cutoff) = 0;

    % --- 4. Build Output Table H ---
    % We typically store the TRUNCATED version in H for quick inspection
    H = [x, F_trunc, U_trunc];

    % --- 5. Generate Interpolants ---
    % H_interpolant: Force (Truncated)
    H_interpolant = griddedInterpolant(x, F_trunc, 'linear', 'none');
    
    % U_interpolant_trunc: The Stability Metric (Zero beyond cutoff)
    U_interpolant_trunc = griddedInterpolant(x, U_trunc, 'linear', 'none');
    
    % U_interpolant_full: The Physical Energy (Valid out to table_limit)
    U_interpolant_full = griddedInterpolant(x, U_full, 'linear', 'none');
end