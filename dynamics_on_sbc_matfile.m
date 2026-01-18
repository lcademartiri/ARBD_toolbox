function [FSTD,GAMMASTD,DEFFSTD]=dynamics_on_sbc_matfile(S,filelist,opts)
% file_list: Cell array of the 10 'POS_partXX.mat' filenames
if nargin < 3, opts = struct(); end
% --- Defaults ---
if ~isfield(opts, 'vecs'),     				opts.vecs = 600; end
if ~isfield(opts, 'do_HReq'), 				opts.do_HReq = false; end
if ~isfield(opts, 'HReq_vecs'),          	opts.HReq_vecs = 180; end
if ~isfield(opts, 'lag_ratio'),          	opts.lag_ratio = 100; end
if ~isfield(opts, 'cacheSizeMB'),    		opts.cacheSizeMB = S.cacheSizeMB; end

	fprintf('### Initializing: Collective Dynamics Analysis ###\n');
	
	% K-WINDOW
    k_fundamental = 2*pi/(2*S.br);
    k_max=pi/S.rp;
    k_mags=(k_fundamental:k_fundamental:k_max)';
    k_mags=sort([k_mags;k_fundamental.*[sqrt(2);sqrt(3);pi]]);
    nK = length(k_mags);
    
    % Fibonacci vectors
    az=fibonacci_sphere(opts.vecs);
    [az,el,~]=cart2sph(az(:,1),az(:,2),az(:,3));
    azel=[az,el];
    azel(azel(:,1)<0,:)=[];

	if opts.do_HReq
		% add equator vectors
		equator=linspace(0,pi,opts.HReq_vecs+1)';
		equator(end,:)=[];
		equator(1,2)=0;
		azel_rad=[azel;equator];
	else
		azel_rad=azel;
	end
	
    azel_deg=rad2deg(azel_rad);
    nAzel=size(azel_rad,1);
    clear az el equator azel
	
	% --- 1. PRE-CALCULATE ALL 3D WAVE-VECTORS (Q) ---
    [K_grid, Az_grid, El_grid] = ndgrid(k_mags, azel_rad(:,1), azel_rad(:,2));
    % Flatten into an M x 3 matrix of [qx, qy, qz]
    qx = K_grid(:) .* (cos(Az_grid(:)) .* cos(El_grid(:)));
    qy = K_grid(:) .* (sin(Az_grid(:)) .* cos(El_grid(:)));
    qz = K_grid(:) .* sin(El_grid(:));
    Q = [qx, qy, qz]; 
    num_k_total = size(Q, 1);
    T_total = 100000; % Total snapshots
    
    % --- 2. ALLOCATE SPECTRAL MATRIX (4.8 GB) ---
    RHO = complex(zeros(num_k_total, T_total, 'single'));
    
    % --- 3. SINGLE-PASS FILE PROCESSING ---
    t_global = 1;
    for f = 1:numel(file_list)
        fprintf('Processing File %d/10...\n', f);
        m = matfile(file_list{f});
        % POS is [N x 3 x T_part]
        POS_part = m.POS; 
        T_part = size(POS_part, 3);
        
        for t_local = 1:T_part
            pos_t = POS_part(:,:,t_local); % [N x 3]
            
            % Truncate to observation window
            idx_in = sum(pos_t.^2, 2) <= S.br^2;
            pos_t = pos_t(idx_in, :);
            
            % VECTORIZED DENSITY CALCULATION
            phase = -1i * (Q * pos_t');
            RHO(:, t_global) = sum(exp(phase), 2);
            
            t_global = t_global + 1;
        end
        clear POS_part % Free RAM for next file
    end

    % --- 4. DYNAMICS ANALYSIS (On the Spectral Matrix) ---
    fprintf('Calculating ACFs...\n');
    max_lag = floor(T_total / opts.lag_ratio);
    % Pre-allocate outputs
    F_AVG_raw = zeros(num_k_total, max_lag+1, 'single');
    
    % Use a simple loop over k-vectors (fast because data is in RAM)
    for k_idx = 1:num_k_total
        rho_t = double(RHO(k_idx, :)); % Convert to double for xcorr
        rho_c = rho_t - mean(rho_t);
        acf = xcorr(rho_c, max_lag, 'biased');
        acf = acf(max_lag+1:end);
        if abs(acf(1)) > 1e-9
            F_AVG_raw(k_idx, :) = single(abs(acf) / abs(acf(1)));
        end
    end
	% --- 5. RESHAPE AND FIT ---
    fprintf('Fitting decay curves for D and Gamma...\n');
    
    % F_AVG_raw was [num_k_total x max_lag+1]
    % Q was built from ndgrid(k_mags, azel_rad(:,1), azel_rad(:,2))
    % So the order of indices in Q (and thus RHO) is:
    % (k=1,dir=1), (k=2,dir=1) ... (k=nK,dir=1), (k=1,dir=2) ...
    
    nK = numel(k_mags);
    nAzel = size(azel_rad, 1);
    
    % Pre-allocate the result maps
    DEFFSTD = zeros(nAzel, nK);
    GAMMASTD = zeros(nAzel, nK);
    FSTD = zeros(nAzel, nK, max_lag + 1, 'single');
    
    % Explicit Loop for Fitting
    % We iterate through the Q-vector list and map results back to the Azel/K grid
    for idir = 1:nAzel
        for k_idx = 1:nK
            % 1. Linear Index in the flattened Spectral Matrix
            % Note: This mapping MUST match how you generated the Q matrix
            q_idx = (idir - 1) * nK + k_idx;
            
            % 2. Extract the specific F(k, t) curve
            F_curve = double(F_AVG_raw(q_idx, :));
            k_val = k_mags(k_idx);
            
            % 3. Store in the 3D output array
            FSTD(idir, k_idx, :) = single(F_curve);
            
            % 4. Call the Fit Helper
            % S.timestep is the time between snapshots (e.g., 100 * dt_internal)
            [D_val, G_val] = fit_dynamics_internal(F_curve, S.timestep, k_val);
            
            % 5. Store in the scalar maps
            DEFFSTD(idir, k_idx) = D_val;
            GAMMASTD(idir, k_idx) = G_val;
        end
    end
    
    fprintf('=== Complete: Collective Dynamics Analysis ===\n');
end

%% --- UPDATED HELPER FUNCTIONS ---

function [D, Gamma] = fit_dynamics_internal(F_norm, dt, k_val)
    % Create time vector for the lags
    lags = (0:(length(F_norm)-1))';
    time = lags * dt;
    
    % LOG-LINEAR FITTING: F(k,t) = exp(-Gamma * t)
    % We fit in the 'Linear' region of the decay
    % 0.9: skip initial ballistic/numerical noise
    % 0.2: skip the noise tail where data is unreliable
    idx = (F_norm < 0.9) & (F_norm > 0.2);
    
    if sum(idx) < 5
        % Not enough points to fit a reliable line
        D = NaN; Gamma = NaN;
    else
        % p = polyfit(x, y, 1) -> y = p(1)x + p(2)
        % log(F) = -Gamma * t
        p = polyfit(time(idx), log(F_norm(idx)), 1);
        
        Gamma = -p(1); % Slope is -Gamma
        
        % Hydrodynamic limit: D = Gamma / k^2
        D = Gamma / (k_val^2);
        
        % Physical sanity check: D and Gamma cannot be negative
        if Gamma < 0
            D = NaN; Gamma = NaN;
        end
    end
end