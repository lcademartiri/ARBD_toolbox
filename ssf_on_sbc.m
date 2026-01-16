function SSF = ssf_on_sbc(S, p, opts)
% SSF = ssf_offline_sbc_optimized(S, p, opts)
%
% Robust Offline Static Structure Factor for SBC (Cartesian + Bessel)
%
% INPUTS
%   S : struct with fields:
%       S.N, S.br (box radius), S.rp (particle radius)
%       (Optional) S.stdx
%   p : N x 3 x T array (Particle Trajectories)
%   opts : Configuration struct (optional)
fprintf('### Initializing: Static Structure Factor and Compressibility Calculation ###\n');

if nargin < 3, opts = struct(); end

% --- Defaults ---
if ~isfield(opts, 'num_shells'),     opts.num_shells = 1; end
if ~isfield(opts, 'dirs_per_shell'), opts.dirs_per_shell = 1; end
if ~isfield(opts, 'l_max'),          opts.l_max = 3; end
if ~isfield(opts, 'n_max'),          opts.n_max = 30; end
if ~isfield(opts, 'cacheSizeMB'),    opts.cacheSizeMB = S.cacheSizeMB; end
if ~isfield(opts, 'force_structural_k'), opts.force_structural_k = true; end

[N, dim, T] = size(p);
assert(dim == 3, 'Particle array must be N x 3 x T');

fprintf('========================================================\n');
fprintf('SSF OFFLINE SBC (Optimized)\n');
fprintf('Particles: %d | Snapshots: %d\n', N, T);

% --- 0. DATA PREP & CENTERING CHECK ---
% SBC Bessel modes assume the center of the container is (0,0,0).
% If the data is shifted (e.g., 0 to L), we must center it.
mean_pos = mean(mean(p, 3), 1); % Average over Time and Particles
if norm(mean_pos) > (S.br / 10)
    fprintf('WARNING: Data does not seem centered at origin.\n');
    fprintf('         Mean pos: [%.2f %.2f %.2f]. Auto-centering...\n', mean_pos);
    p = p - mean_pos;
end

% --- 1. INTELLIGENT K-LIMITS ---
% The thermal limit is often too high for structural visualization.
% We define a structural limit (~15-20 sphere diameters).
k_struct_max = 25 * (pi / S.rp); 

if isfield(S, 'stdx') && S.stdx > 1e-9
    k_thermal_max = pi / S.stdx;
else
    k_thermal_max = k_struct_max * 10;
end

kmin = pi / S.br;

if opts.force_structural_k
    % Stick to the structural range for a useful plot
    kmax = k_struct_max;
    fprintf('Mode: STRUCTURAL Display (k_max set to %.2e)\n', kmax);
else
    % Go all the way to thermal limit
    kmax = k_thermal_max;
    fprintf('Mode: THERMAL Display (k_max set to %.2e)\n', kmax);
end

% Diagnostic: Expected first peak location
k_peak_est = 2*pi / (2*S.rp);
fprintf('Diagnostic: First structural peak expected around k = %.2e\n', k_peak_est);


% --- 2. CARTESIAN PROJECTION (FIBONACCI) ---
fprintf('\n--- Cartesian Projection ---\n');
[K_cart, kmag_cart] = generate_fibonacci_k(kmin, kmax, opts.num_shells, opts.dirs_per_shell);
fprintf('Sampling %d k-vectors...\n', size(K_cart,1));

[cart_mean_rho, cart_mean_abs2] = project_cartesian_chunked(p, K_cart, opts.cacheSizeMB,S);

% S(k) calculation
S_cart_tot  = cart_mean_abs2 / N;           % < |rho|^2 > / N
S_cart_stat = abs(cart_mean_rho).^2 / N;    % | <rho> |^2 / N (Bragg/Wall)
S_cart_fluc = S_cart_tot - S_cart_stat;     % Pure structure

SSF.cartesian.k     = kmag_cart;
SSF.cartesian.Stot  = S_cart_tot;
SSF.cartesian.Sstat = S_cart_stat;
SSF.cartesian.Sfluc = S_cart_fluc;


% --- 3. BESSEL EIGENMODE GENERATION ---
fprintf('\n--- SBC Bessel Projection ---\n');
fprintf('Generating roots/modes (l_max=%d, n_max=%d)...\n', opts.l_max, opts.n_max);
[modes, W_ln] = generate_sbc_modes(S.br, kmax, opts.l_max, opts.n_max);
fprintf('Generated %d Bessel modes.\n', size(modes,1));

if isempty(modes)
    error('No Bessel modes found in range! Check units of S.br vs S.rp.');
end


% --- 4. OPTIMIZED BESSEL PROJECTION ---
fprintf('Projecting (Matrix Optimized)...\n');
S_bessel_raw = project_bessel_optimized(p, modes, W_ln, opts.cacheSizeMB,S);

SSF.bessel.k = modes(:,1);
SSF.bessel.l = modes(:,2);
SSF.bessel.n = modes(:,3);
SSF.bessel.S = S_bessel_raw;

% --- 5. CALCULATE COMPRESSIBILITY ---
idx0=modes(:,2)==0;
temp=[modes(idx0,1),S_bessel_raw(idx0)];
temp=sortrows(temp,1);
temp=temp(1:5,:);
temp(:,1)=temp(:,1).^2;
X = [ones(size(temp(:,1))), temp(:,1)];
% weights
w = 1 ./ temp(:,1);
% weighted least squares
X = [ones(size(temp(:,1))), temp(:,1)];
W = diag(w);
beta = (X' * W * X) \ (X' * W * y);
SSF.S0=beta(1);

fprintf('=== Completed: Static Structure Factor and Compressibility Calculation ===\n');

end


%% ========================================================================
%  HELPER FUNCTIONS
% ========================================================================

function [K, kmag] = generate_fibonacci_k(kmin, kmax, nshell, ndir)
    % Log-spaced shells to capture low-k structure and high-k tails
    k_rad = logspace(log10(kmin), log10(kmax), nshell)';
    
    % Golden spiral points on sphere
    i = (0:ndir-1)';
    phi = acos(1 - 2*(i+0.5)/ndir);
    theta = pi*(1+sqrt(5))*(i+0.5);
    dirs = [sin(phi).*cos(theta), sin(phi).*sin(theta), cos(phi)];
    
    % Kronecker product equivalent
    Nk = nshell * ndir;
    K  = zeros(Nk,3);
    kmag = zeros(Nk,1);
    
    idx = 1;
    for i = 1:nshell
        K(idx:idx+ndir-1, :) = k_rad(i) * dirs;
        kmag(idx:idx+ndir-1) = k_rad(i);
        idx = idx + ndir;
    end
end

function [mean_rho, mean_abs2] = project_cartesian_chunked(p, K, cacheMB,S)
	fprintf('Calculating Cartesian Projection...\n');
    [N,~,T] = size(p);
    Nk = size(K,1);
    
    mean_rho  = zeros(Nk,1);
    mean_abs2 = zeros(Nk,1);
    
    % Memory chunking
    bytes_per_num = 16; % Complex double
    target_bytes = cacheMB * 1024^2 * 0.4; % Use 40% of cache per buffer
    
    % We compute matrix [ChunkK x N].
    chunk_k = floor(target_bytes / (N * bytes_per_num));
    chunk_k = max(1, min(chunk_k, 5000));
    
    % Pre-permute p for faster access: [3 x N x T]
    p_perm = permute(p, [2 1 3]);
    tStart=tic;
    for s = 1:chunk_k:Nk
		counterstruct = struct('Stage','Cartesian Projection', 'Chunk', s, 'Total_Chunks', ceil(Nk/chunk_k));
        e = min(s+chunk_k-1, Nk);
        K_sub = K(s:e,:); % [Mk x 3]
        
        rho_sum_t = zeros(e-s+1, 1);
        abs2_sum_t = zeros(e-s+1, 1);
        
        % Iterate snapshots
        for t = 1:T
            % Vectorized phase: (Mk x 3) * (3 x N) -> (Mk x N)
            ptemp=squeeze(p_perm(:,:,t));
            ptemp(:,vecnorm(ptemp,2,1)>S.br)=[]; % exclusion of reals in the fuzz layer
            ptemp=ptemp-mean(ptemp,2); % COM subtraction
            phase = -1i * (K_sub * ptemp);
            rho_t = sum(exp(phase), 2); % Sum over particles -> (Mk x 1)
            
            rho_sum_t = rho_sum_t + rho_t;
            abs2_sum_t = abs2_sum_t + abs(rho_t).^2;
			progressUpdate(t, T, tStart, 100, counterstruct)
        end
        
        mean_rho(s:e) = rho_sum_t / T;
        mean_abs2(s:e) = abs2_sum_t / T;
    end
end

function [modes, W_ln] = generate_sbc_modes(R, kmax, lmax, nmax)
	fprintf('Generating Bessel Basis...\n');
    % Solves j_l'(k*R) = 0 for Neumann boundary
    modes_list = [];
    
    for l = 0:lmax
        roots = [];
        % Approximate roots for jl'
        for n = 1:nmax
            % Empirical guess for derivative roots
            guess = (n + l/2 - 0.5) * pi; 
            if l==0, guess = n*pi; end % j0' = -j1 -> roots of j1
            
            f = @(x) sph_bessel_deriv(l,x);
            
            try
                rt = fzero(f, guess);
                if rt > 1e-4
                    roots = [roots; rt];
                end
            catch
                % Silent catch for root finding
            end
        end
        roots = unique(roots);
        k_vals = roots / R;
        
        valid = k_vals <= kmax;
        if any(valid)
            k_valid = k_vals(valid);
            n_indices = (1:length(k_valid))';
            % Append: [k, l, n]
            modes_list = [modes_list; [k_valid, repmat(l, numel(k_valid),1), n_indices]];
        end
    end
    
    if isempty(modes_list)
        modes = []; W_ln = []; return;
    end
    
    % Sort by k
    [~, idx] = sort(modes_list(:,1));
    modes = modes_list(idx,:);
    
    % Normalization W_ln
    % For Neumann (jl'(kR)=0): W = (R^3/2) * jl(kR)^2 * (1 - l(l+1)/(kR)^2)
    V = (4/3)*pi*R^3;
    W_ln = zeros(size(modes,1),1);
    
    for i = 1:size(modes,1)
        k = modes(i,1);
        l = modes(i,2);
        x = k*R;
        
        val = sph_bessel_j(l, x);
        if l==0
            geo = 1; 
        else
            geo = (1 - (l*(l+1))/(x^2));
        end
        W_ln(i) = (R^3 / 2) * val^2 * geo / V; % Normalized by volume
    end
end

function S_final = project_bessel_optimized(p, modes, W_ln, cacheMB,S)
	fprintf('Calculating Bessel Basis Projection...\n');
    [N, ~, T] = size(p);
    
    unique_l = unique(modes(:,2));
    
    % Accumulator for Power Spectrum C_ln = < sum_m |rho_lnm|^2 >
    C_ln_acc = zeros(size(modes,1), 1);
    
    % Map modes to indices for fast lookup
    % We process one 'l' at a time to reuse Y_lm
    tStart=tic;
    for il = 1:length(unique_l)
        l = unique_l(il);
        counterstruct = struct('Stage','Bessel Projection', 'l', l, 'total_ls', length(unique_l));
        % 1. Slice the mode-specific vectors for this 'l'
        mask = (modes(:,2) == l);
        idxs_l = find(mask);
        ks_l = modes(mask, 1);
        num_k = length(ks_l);
        
        % SLICE THE NORMALIZATION VECTORS HERE
        W_ln_sub = W_ln(mask); 
        degen_sub = (2*l + 1); % Or degen(mask) if you pre-calculated it
        
        if num_k == 0, continue; end
        
        for t = 1:T
			% 1. Extract particles in window
			r_vec = p(:, :, t); 
			r_vec(vecnorm(r_vec,2,2) > S.br, :) = [];
			% 2. GET INSTANTANEOUS COUNT
			Nt = size(r_vec, 1);
			if Nt == 0, continue; end % Safety
            r_vec=r_vec-mean(r_vec,1); % COM subtraction
            [phi, elev, r] = cart2sph(r_vec(:,1), r_vec(:,2), r_vec(:,3));
            theta = pi/2 - elev;
            
            % 2. Spherical Harmonics Y_lm [N x (2l+1)]
            Y = compute_ylm(l, theta, phi);
            
            % 3. Radial Bessel Functions [NumK x N]
            % J(k_idx, particle_idx)
            % This is the most expensive part usually.
            kr = ks_l * r'; % Outer product: (NumK x 1) * (1 x N) -> NumK x N
            J = sph_bessel_j_matrix(l, kr); 
            
            % 3. Project and get Power [num_k x 2l+1] -> [num_k x 1]
            Rho_lnm = J * Y;
            Power_ln = sum(abs(Rho_lnm).^2, 2);
            
            % 4. NORMALIZE IMMEDIATELY (Using the sliced vectors)
            % Nt is size [1x1], Power_ln is [num_k x 1], W_ln_sub is [num_k x 1]
            S_snapshot = (Power_ln ./ W_ln_sub) ./ Nt ./ degen_sub;
            
            % 5. Accumulate into the main array at the correct indices
            C_ln_acc(idxs_l) = C_ln_acc(idxs_l) + S_snapshot;
			progressUpdate(t, T, tStart, 100, counterstruct)
        end
    end
    
    % Time Average
    C_ln_avg = C_ln_acc / T;
    
    % Final S(k) Normalization
    % S = (C_ln / W_ln) / N / (Degeneracy)
    % Degeneracy for RIPS is (2l+1)
    
    
    S_final = C_ln_acc / T;

end

% --- MATH HELPERS ---

function val = sph_bessel_j(l, x)
    % Scalar/Vector wrapper
    val = sqrt(pi./(2*x)) .* besselj(l+0.5, x);
    if l==0
        % Fix singularity at 0
        val(x<1e-6) = 1;
    else
        val(x<1e-6) = 0;
    end
end

function J = sph_bessel_j_matrix(l, Z)
    % Optimized for large matrix Z
    Z(Z<1e-9) = 1e-9; % Avoid div by zero
    J = sqrt(pi./(2*Z)) .* besselj(l+0.5, Z);
    if l==0
        J(Z<1e-6) = 1;
    end
end

function val = sph_bessel_deriv(l, x)
    % Derivative of j_l(x)
    if x < 1e-5
        if l==1, val=0.33; else, val=0; end % Rough approx for roots check
        return; 
    end
    % Recurrence: j_l'(x) = j_{l-1}(x) - (l+1)/x * j_l(x)
    % Note: j_{-1} = -y_0 (spherical neumann)? No, use explicit form:
    % j_l'(x) = ( l*j_{l-1} - (l+1)*j_{l+1} ) / (2l+1) ... wait, simpler:
    jl = sph_bessel_j(l,x);
    if l==0
        val = -sph_bessel_j(1,x);
    else
        jlm1 = sph_bessel_j(l-1,x);
        val = jlm1 - (l+1)/x * jl;
    end
end

