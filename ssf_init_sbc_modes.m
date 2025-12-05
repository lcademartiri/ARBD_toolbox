function SSF = ssf_init_sbc_modes(SSF, S, opts)
% SSF = ssf_init_sbc_modes(SSF, S, opts)
%
% Computes the "Allowed" K-values for a Spherical Cavity (SBC).
% Solves j_l'(k*R) = 0.
%
% Inputs:
%   S.br : Sphere Radius
%   opts.l_max : Maximum angular momentum (e.g. 20)
%   opts.n_max : Number of radial roots per l (e.g. 50)
%   opts.store_rho_ts : (boolean) whether to store time series

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'l_max'), opts.l_max = 25; end
    if ~isfield(opts, 'n_max'), opts.n_max = 50; end
    if ~isfield(opts, 'store_rho_ts'), opts.store_rho_ts = false; end
    
    R = S.br;
    
    fprintf('SSF SBC: Computing spherical Bessel roots (l=0..%d, n=1..%d)... ', opts.l_max, opts.n_max);
    
    % Pre-allocate
    % We store: k_val, l, n
    count = (opts.l_max + 1) * opts.n_max;
    modes = zeros(count, 3);
    idx = 1;
    
    for l = 0:opts.l_max
        for n = 1:opts.n_max
            % Find the n-th root of j_l'(x) = 0
            % Guess: Roots are roughly spaced by pi.
            % Offset depends on l. approx x ~ (n + l/2)*pi
            guess = (n + l/2) * pi;
            
            % Root finder
            f = @(x) spherical_bessel_deriv(l, x);
            
            % Robust search
            try
                root_x = fzero(f, guess);
            catch
                % Fallback search bracket if guess is bad
                root_x = fzero(f, [guess-pi/2, guess+pi]);
            end
            
            k_val = root_x / R;
            
            modes(idx, :) = [k_val, l, n];
            idx = idx + 1;
        end
    end
    
    % Sort by k magnitude
    [~, sort_ord] = sort(modes(:,1));
    modes = modes(sort_ord, :);
    
    % Filter only those within our global thermal limit (if it exists)
    % This relies on SSF.kmax_thermal being set by the previous unified init
    if isfield(SSF, 'kmax_thermal')
        mask = modes(:,1) <= SSF.kmax_thermal;
        modes = modes(mask, :);
    end
    
    SSF.sbc_modes.data = modes; % [k, l, n]
    SSF.sbc_modes.kvals = modes(:,1);
    
    % Storage for accumulation
    % We accumulate the "Power Spectrum" C_ln = sum_m |Rho_lnm|^2
    SSF.sbc_modes.sum_power = zeros(size(modes,1), 1);
    
    % --- FIXED LINE BELOW (Uses 'opts' directly) ---
    if opts.store_rho_ts
        SSF.sbc_modes.power_ts = zeros(0, size(modes,1));
    else
        SSF.sbc_modes.power_ts = [];
    end
    
    fprintf('Done. Found %d modes.\n', size(modes,1));
end

% --- MATH HELPER ---
function val = spherical_bessel_deriv(l, x)
    % Derivative of j_l(x)
    % Relation: j_l'(x) = j_{l-1}(x) - ((l+1)/x) * j_l(x)
    
    % Compute j_l(x) and j_{l-1}(x)
    % Matlab's besselj is cylindrical J_nu.
    % j_l(x) = sqrt(pi/(2x)) * J_{l+0.5}(x)
    
    if x < 1e-9
        if l == 0
            % j0(x) = sin(x)/x -> deriv is 0 at 0? No, limit is 0.
            val = 0; return; 
        elseif l == 1
            % j1(x) start linear, deriv const
            val = 1/3; return;
        else
            val = 0; return;
        end
    end

    pre = sqrt(pi/(2*x));
    jl      = pre * besselj(l + 0.5, x);
    
    if l == 0
        % j0'(x) = -j1(x)
        val = - (pre * besselj(1.5, x));
    else
        jlm1    = pre * besselj(l - 0.5, x);
        val = jlm1 - ((l+1)/x) * jl;
    end
end