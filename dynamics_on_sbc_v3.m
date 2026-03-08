function [FSTD,GAMMASTD,DEFFSTD,opts,azel_rad,k_mags]=dynamics_on_sbc_v3(S,POS,opts)

    if nargin < 3, opts = struct(); end
    % --- Defaults ---
    if ~isfield(opts, 'vecs'),          opts.vecs = 600; end
    if ~isfield(opts, 'do_HReq'),       opts.do_HReq = true; end
    if ~isfield(opts, 'HReq_vecs'),     opts.HReq_vecs = 180; end
    if ~isfield(opts, 'lag_ratio'),     opts.lag_ratio = 100; end
	% --- 1. DATA PREPARATION ---
    % Convert cell array back to a 3D matrix (T x N x 3) for vectorization
    % Use NaN for particles outside the boundary to keep matrix shape
    T = size(POS, 3);
    max_lag = T/opts.lag_ratio;
	
	% Pre-extract coordinates [N x T]
    pux = squeeze(POS(1:S.N,1,:)); 
    puy = squeeze(POS(1:S.N,2,:)); 
    puz = squeeze(POS(1:S.N,3,:));
    
    % Mask particles outside the observation window
    R_sq = pux.^2 + puy.^2 + puz.^2;
	mask = R_sq > S.br^2;
	pux(mask) = NaN; puy(mask) = NaN; puz(mask) = NaN;
	clear R_sq mask
    
    % --- 2. K-WINDOW & VECTORS ---
    k_fundamental = 2*pi/(2*S.br);
	% The thermal limit is often too high for structural visualization.
	% We define a structural limit (~15-20 sphere diameters).
	k_struct_max = 25 * (pi / S.rp); 
    
	if isfield(S, 'stdx') && S.stdx > 1e-9
		k_max = pi / S.stdx;
	else
		k_max = k_struct_max/10;
	end
    k_mags = sort([(k_fundamental:k_fundamental:k_max)'; k_fundamental.*[sqrt(2);sqrt(3);pi]]);
    
    az_vecs = fibonacci_sphere(opts.vecs);
    [az, el, ~] = cart2sph(az_vecs(:,1), az_vecs(:,2), az_vecs(:,3));
    azel = [az, el];
    azel(azel(:,1) < 0, :) = [];
    if opts.do_HReq
        equator = [linspace(0, pi, opts.HReq_vecs+1)', zeros(opts.HReq_vecs+1, 1)];
        azel_rad = [azel; equator(1:end-1, :)];
    else
        azel_rad = azel;
    end
	
    nAzel=size(azel_rad,1);
    nK = length(k_mags);
    [~,rotmat]=FCCrotate([1,0,0],[1,1,1]./norm([1,1,1]));
    rotmat=rotmat';
    
    % Scalar Maps (Keep Replicates for Error Bars)
        % Dimensions: Theta/Phi x K x Replicates
    DEFFSTD   = zeros(nAzel, nK);
    GAMMASTD  = zeros(nAzel, nK);
    FSTD = zeros(nAzel, nK, max_lag+1);
    
    if S.bc == 2 || S.bc==3
        DEFFMASK  = zeros(nAzel, nK);
        GAMMAMASK = zeros(nAzel, nK);
        FMASK= zeros(nAzel, nK, max_lag+1);
    end
    
    % SUPERLOOP
    tStart=tic;
    for idir=1:nAzel
        counterstruct = struct('Stage','Dynamics', 'k_vector', idir, 'Total_k_vectors', nAzel);	
        az = azel_rad(idir,1);
        sin_az = sin(az); 
        cos_az = cos(az);
        el = azel_rad(idir,2);
        cos_el = cos(el);
        sin_el = sin(el);
        nx = cos_az * cos_el; % Unit Direction Vector n_hat
        ny = sin_az * cos_el; % Unit Direction Vector n_hat
        nz = sin_el; % Unit Direction Vector n_hat
        if S.bc==3
            ntemp=rotmat*[nx;ny;nz];
            nx=ntemp(1);
            ny=ntemp(2);
            nz=ntemp(3);
        end            
        for k_idx = 1:nK
            k_val = k_mags(k_idx);
            qx = k_val * nx; % 3D q-vector
            qy = k_val * ny; % 3D q-vector
            qz = k_val * nz; % 3D q-vector
            
            % -- DYNAMICS Deff(k) --
            % Use unwrapped coords (put)
            
            phase_dyn = -(qx.*pux + qy.*puy + qz.*puz);
            E_dyn = exp(1i * phase_dyn);
			% Determine who is inside (not NaN)
            inside_mask = ~isnan(E_dyn);
            E_dyn(~inside_mask) = 0; % Set outside to zero
            
            % Calculate Rho normalized by instantaneous particle count
            % This removes the "Shot Noise" of entries/exits
            rho_dyn = sum(E_dyn, 1) ./ (sum(inside_mask, 1) + eps);
            [F_curve, ~] = compute_acf(rho_dyn', max_lag);
            [D_val, G_val] = fit_dynamics(F_curve, S.timestep, k_val);
            
            DEFFSTD(idir, k_idx) = D_val;
            GAMMASTD(idir, k_idx) = G_val;
            FSTD(idir, k_idx, :) = FSTD(idir, k_idx, :) + reshape((F_curve), 1, 1, []);
            
            if S.bc == 2 || S.bc==3
                rho_dyn_m = sum(mask{irep,1} .* E_dyn, 1);
                [F_curve_m, ~] = compute_acf(rho_dyn_m, max_lag);
                [D_m, G_m] = fit_dynamics(F_curve_m, S.timestep, k_val);
                DEFFMASK(idir, k_idx) = D_m;
                GAMMAMASK(idir, k_idx) = G_m;
                FMASK(idir, k_idx, :) = FMASK(idir, k_idx, :) + reshape(single(F_curve_m), 1, 1, []);
            end           
        end
        progressUpdate(idir, nAzel, tStart, 1, counterstruct);
    end
end

%% FAST HELPER FUNCTIONS

function [F_norm, lags_vec] = compute_acf(rho_t, max_lag)
    % Subtract mean to get fluctuations
    rho_c = rho_t - mean(rho_t);
    % Fast FFT-based autocorrelation
    n = length(rho_c);
    m = 2^nextpow2(2*n-1);
    Rho_f = fft(rho_c, m);
    acf = ifft(Rho_f .* conj(Rho_f));
    acf = real(acf(1:max_lag+1));
    
    % Biased normalization (matches xcorr 'biased')
    lags_vec = 0:max_lag;
    counts = n - lags_vec';
    acf = acf ./ counts; % Normalization per lag
    
    if acf(1) > 1e-12
        F_norm = acf / acf(1);
    else
        F_norm = zeros(max_lag+1, 1);
    end
end

function [D, Gamma] = fit_dynamics(F_norm, dt, k_val)
    time = (0:(length(F_norm)-1))' * dt;
    % Fit only the reliable decay range
    idx = F_norm < 0.9 & F_norm > 0.15;
    if sum(idx) < 5
        D = NaN; Gamma = NaN;
    else
        % Linear fit of log(F) vs time: log(F) = -Gamma * t
        p = polyfit(time(idx), log(F_norm(idx)), 1);
        Gamma = -p(1);
        D = Gamma / k_val^2;
    end
end