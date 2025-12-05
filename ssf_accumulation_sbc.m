function SSF = ssf_accumulation_sbc(p, SSF)
% Projects density onto Spherical Bessel Eigenmodes.
% This is the "Spectral Analysis" for a sphere.

    % Coordinate transform: Cartesian -> Spherical
    % p is N x 3
    x = p(:,1); y = p(:,2); z = p(:,3);
    [phi, elev, r] = cart2sph(x, y, z);
    theta = pi/2 - elev; % Convert elevation to zenith (0 to pi)

    N = size(p, 1);
    modes = SSF.sbc_modes.data; % [k, l, n]
    num_modes = size(modes, 1);
    
    % Prepare array for power in this snapshot
    current_power = zeros(1, num_modes);
    
    % We loop over L because Y_lm computation is expensive and shared
    unique_ls = unique(modes(:,2));
    
    for i = 1:numel(unique_ls)
        l = unique_ls(i);
        
        % 1. Compute Spherical Harmonics Y_lm for all particles at this l
        % (Requires a Ylm function - see below)
        % Y is N_part x (2l+1)
        Y = compute_ylm_vectorized(l, theta, phi); 
        
        % 2. Find all n-roots associated with this l
        idx_mask = (modes(:,2) == l);
        ks_for_l = modes(idx_mask, 1);
        indices  = find(idx_mask);
        
        % 3. Compute Radial Bessel Part for these k's
        % R_mat is (Number_of_n_roots) x N_part
        % We calculate jl(k * r_particle)
        R_mat = compute_bessel_matrix(l, ks_for_l, r);
        
        % 4. Perform the Projection (The "Hankel-Fourier" transform)
        % Rho_lm = Sum_over_particles( R(r) * Y(theta,phi) )
        % We want Sum_m |Rho_lm|^2
        
        % Matrix mult: (Num_n x N_part) * (N_part x Num_m) -> (Num_n x Num_m)
        % This gives us rho_{nm} for this specific l
        Rho_nm = R_mat * Y;
        
        % Sum square modulus over m to get rotationally invariant power
        Power_n = sum(abs(Rho_nm).^2, 2); % Sum over m columns
        
        current_power(indices) = Power_n;
    end
    
    % Accumulate
    SSF.sbc_modes.sum_power = SSF.sbc_modes.sum_power + current_power';
    
    if ~isempty(SSF.sbc_modes.power_ts)
        SSF.sbc_modes.power_ts(SSF.nsnap, :) = current_power;
    end
end

% --- HELPERS ---

function R_mat = compute_bessel_matrix(l, k_vals, r_vals)
    % Returns matrix J(i,j) = j_l( k_vals(i) * r_vals(j) )
    % k: M x 1, r: N x 1
    
    % Argument matrix
    Z = k_vals * r_vals'; % M x N
    
    % Compute spherical bessel
    % j_l(z) = sqrt(pi/2z) J_{l+0.5}(z)
    
    % Handle small z to avoid NaN
    safe_Z = Z;
    safe_Z(Z<1e-9) = 1e-9; 
    
    pre = sqrt(pi ./ (2*safe_Z));
    cyl = besselj(l + 0.5, safe_Z);
    R_mat = pre .* cyl;
    
    % Fix z=0 limit
    if l==0
        R_mat(Z<1e-9) = 1;
    else
        R_mat(Z<1e-9) = 0;
    end
end

function Y = compute_ylm_vectorized(l, theta, phi)
    % Wrapper for Legendre Polynomials to get Y_lm
    % Returns N x (2l+1) matrix
    % Columns correspond to m = -l, ..., 0, ..., +l
    
    N = numel(theta);
    % legendre(l, x) returns matrix of size (l+1) x N for m=0..l
    P = legendre(l, cos(theta)); % shape: (l+1) x N
    
    % We need to construct full Y_lm including negative m and phase
    % Y_lm = N_lm * P_lm(cos theta) * exp(i*m*phi)
    
    % Allocation (2l+1 cols)
    Y = zeros(N, 2*l+1);
    
    % m = 0
    m = 0;
    norm_factor = sqrt((2*l+1)/(4*pi));
    Y(:, l+1) = norm_factor * P(1, :)'; 
    
    % Loop m = 1 to l
    for m = 1:l
        % Normalization factor
        % N_lm = sqrt( (2l+1)/(4pi) * (l-m)! / (l+m)! )
        fac = factorial(l-m) / factorial(l+m);
        norm = sqrt( ((2*l+1)/(4*pi)) * fac );
        
        P_lm = P(m+1, :)'; % Fetch P_lm (m is row index m+1)
        
        % Positive m
        Y(:, l+1+m) = norm * P_lm .* exp(1i * m * phi) * (-1)^m; 
        
        % Negative m (symmetry relation)
        % Y_{l,-m} = (-1)^m Y_{l,m}^*
        Y(:, l+1-m) = conj(Y(:, l+1+m)) * (-1)^m;
    end
end