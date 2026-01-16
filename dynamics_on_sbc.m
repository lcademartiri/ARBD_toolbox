function [FSTD,GAMMASTD,DEFFSTD]=dynamics_on_sbc(S,POS,opts)

if nargin < 3, opts = struct(); end
% --- Defaults ---
if ~isfield(opts, 'vecs'),     				opts.vecs = 600; end
if ~isfield(opts, 'do_HReq'), 				opts.do_HReq = false; end
if ~isfield(opts, 'HReq_vecs'),          	opts.HReq_vecs = 180; end
if ~isfield(opts, 'lag_ratio'),          	opts.lag_ratio = 100; end
if ~isfield(opts, 'cacheSizeMB'),    		opts.cacheSizeMB = S.cacheSizeMB; end

	fprintf('### Initializing: Collective Dynamics Analysis ###\n');
	T_steps=size(POS{1,1},3);
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
    
    % DYNAMICS
    max_lag = T_steps/opts.lag_ratio;
    nAzel=size(azel_rad,1);
    nK = length(k_mags);
    
    % Scalar Maps (Keep Replicates for Error Bars)
        % Dimensions: Theta/Phi x K x Replicates
    DEFFSTD   = zeros(nAzel, nK, size(POS,1));
    GAMMASTD  = zeros(nAzel, nK, size(POS,1));
    
    % Full Time Correlation Function (Average Only)
    % Dimensions: Theta/Phi x K x Time
    % Using 'single' to save 50% RAM
    F_AVG_STD = zeros(nAzel, nK, max_lag+1);
    
    % SUPERLOOP
	
    for irep=1:size(POS,1)
        % initialize replicate storage
        seg_D = zeros(nAzel, nK);
        seg_G = zeros(nAzel, nK);
        for idir=1:nAzel
			counterstruct = struct('Stage','Dynamics', 'Slice', irep, 'k_vector', idir, 'Total_k_vectors', nAzel);
			tStart=tic;
            az = azel_rad(idir,1);
            sin_az = sin(az); 
            cos_az = cos(az);
            el = azel_rad(idir,2);
            cos_el = cos(el);
            sin_el = sin(el);
            nx = cos_az * cos_el; % Unit Direction Vector n_hat
            ny = sin_az * cos_el; % Unit Direction Vector n_hat
            nz = sin_el; % Unit Direction Vector n_hat          
            for k_idx = 1:nK
                k_val = k_mags(k_idx);
                qx = k_val * nx; % 3D q-vector
                qy = k_val * ny; % 3D q-vector
                qz = k_val * nz; % 3D q-vector
                
                % -- DYNAMICS Deff(k) --
				pux = squeeze(POS{irep,1}(:,1,:));
				puy = squeeze(POS{irep,1}(:,2,:));
				puz = squeeze(POS{irep,1}(:,3,:));
                
                phase_dyn = -(qx.*pux + qy.*puy + qz.*puz);
                E_dyn = exp(1i * phase_dyn);
                rho_dyn = sum(E_dyn, 1);
                [F_curve, ~] = compute_acf(rho_dyn, max_lag);
                [D_val, G_val] = fit_dynamics(F_curve, S.timestep, k_val);
                
                seg_D(idir, k_idx) = D_val;
                seg_G(idir, k_idx) = G_val;
                F_AVG_STD(idir, k_idx, :) = F_AVG_STD(idir, k_idx, :) + reshape((F_curve), 1, 1, []);
                
            end
			progressUpdate(idir, nAzel, tStart, 1, counterstruct)
        end
        % Store segment
        FSTD(:,:,:,irep)=F_AVG_STD;
        DEFFSTD(:,:,irep) = seg_D;
        GAMMASTD(:,:,irep) = seg_G;
        if S.bc == 2 || S.bc==3
            FMASK(:,:,:,irep)=F_AVG_MASK;
            DEFFMASK(:,:,irep) = seg_Dm;
            GAMMAMASK(:,:,irep) = seg_Gm;
        end
    end
	fprintf('=== Complete: Collective Dynamics Analysis ===\n');
end

%% HELPER FUNCTIONS

function [F_norm, lags_vec] = compute_acf(rho_t, max_lag)
    if any(isnan(rho_t))
        F_norm = nan(max_lag+1, 1); lags_vec = 0:max_lag; return; 
    end
    rho_c = rho_t - mean(rho_t);
    acf = xcorr(rho_c, max_lag, 'biased');
    acf = acf(max_lag+1:end);
    if abs(acf(1)) > 1e-9
        F_norm = abs(acf) / abs(acf(1));
    else
        F_norm = zeros(size(acf));
    end
    lags_vec = 0:max_lag;
end

function [D, Gamma] = fit_dynamics(F_norm, dt, k_val)
    time = (0:(length(F_norm)-1)) * dt;
    idx = F_norm < 0.9 & F_norm > 0.2;
    if sum(idx) < 5
        D = NaN; Gamma = NaN;
    else
        p = polyfit(time(idx), log(F_norm(idx)), 1);
        Gamma = -p(1);
        D = Gamma / k_val^2;
    end
end