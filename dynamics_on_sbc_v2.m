function [FSTD,GAMMASTD,DEFFSTD,opts,azel_rad,k_mags]=dynamics_on_sbc_v2(S,POS,opts)

    if nargin < 3, opts = struct(); end
    % --- Defaults ---
    if ~isfield(opts, 'vecs'),          opts.vecs = 600; end
    if ~isfield(opts, 'do_HReq'),       opts.do_HReq = true; end
    if ~isfield(opts, 'HReq_vecs'),     opts.HReq_vecs = 180; end
    if ~isfield(opts, 'lag_ratio'),     opts.lag_ratio = 100; end

    fprintf('### Initializing: Collective Dynamics Analysis ###\n');

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
    k_max = pi/S.rp;
    k_mags = sort([(k_fundamental:k_fundamental:k_max)'; k_fundamental.*[sqrt(2);sqrt(3);pi]]);
    nK = length(k_mags);
    
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
    
    % Scalar Maps (Keep Replicates for Error Bars)
        % Dimensions: Theta/Phi x K x Replicates
    DEFFSTD   = zeros(nAzel, nK, size(POS,1));
    GAMMASTD  = zeros(nAzel, nK, size(POS,1));
    FSTD = zeros(nAzel, nK, max_lag+1);
    tStart = tic;
    for idir = 1:nAzel
        az = azel_rad(idir, 1); el = azel_rad(idir, 2);
        nx = cos(az)*cos(el); ny = sin(az)*cos(el); nz = sin(el);
        
        % Projection [N x T]
        U_proj = nx*pux + ny*puy + nz*puz; 

        for k_idx = 1:nK
            k_val = k_mags(k_idx);
            
            % 1. Compute Phase Matrix [N x T]
            % We keep the particles separate!
            E_dyn = exp(-1i * k_val * U_proj); 
            
            % 2. THE VECTORIZED ACF
            % Instead of one ACF of a sum, we do the sum of ACFs.
            % This is mathematically equivalent to the Self-Scattering Function.
            F_sum = 0;
            for i = 1:S.N
                % Correlation for single particle 'i'
                [Fi, ~] = compute_acf(E_dyn(i,:).', max_lag);
                F_sum = F_sum + Fi;
            end
            F_curve = F_sum / S.N; % Average F(k,t)

            % Fit and Store
            [D_val, G_val] = fit_dynamics(F_curve, S.timestep, k_val);
            DEFFSTD(idir, k_idx) = D_val;
            GAMMASTD(idir, k_idx) = G_val;
            FSTD(idir, k_idx, :) = F_curve;
        end
        progressUpdate(idir, nAzel, tStart, 1, struct('Stage','SBC Dynamics'));
    end
	
    fprintf('=== Complete: Collective Dynamics Analysis ===\n');
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