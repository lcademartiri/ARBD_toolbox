function [FSTD,GAMMASTD,DEFFSTD]=dynamics_on_sbc(S,sput)
	
	T_steps=size(POS{1,1},3);
	% K-WINDOW
    k_fundamental = 2*pi/(2*S.br);
    k_max=pi/S.rp;
    k_mags=(k_fundamental:k_fundamental:k_max)';
    k_mags=sort([k_mags;k_fundamental.*[sqrt(2);sqrt(3);pi]]);
    nK = length(k_mags);
    
    % Fibonacci vectors
    az=fibonacci_sphere(600);
    [az,el,~]=cart2sph(az(:,1),az(:,2),az(:,3));
    azel=[az,el];
    azel(azel(:,1)<0,:)=[];

    % add equator vectors
    equator=linspace(0,pi,181)';
    equator(end,:)=[];
    equator(1,2)=0;
    azel_rad=[azel;equator];
    azel_deg=rad2deg(azel_rad);
    nAzel=size(azel_rad,1);
    clear az el equator azel
    
    % DYNAMICS
    max_lag = T_steps/100;
	
    nAzel=size(azel_rad,1);
    nK = length(k_mags);
    rotmat=rotmat';
    
    % Scalar Maps (Keep Replicates for Error Bars)
        % Dimensions: Theta/Phi x K x Replicates
    DEFFSTD   = zeros(nAzel, nK, size(sput,1));
    GAMMASTD  = zeros(nAzel, nK, size(sput,1));
    
    % Full Time Correlation Function (Average Only)
    % Dimensions: Theta/Phi x K x Time
    % Using 'single' to save 50% RAM
    F_AVG_STD = zeros(nAzel, nK, max_lag+1);
    
    % SUPERLOOP
    for irep=1:size(sput,1)
        % initialize replicate storage
        seg_D = zeros(nAzel, nK);
        seg_G = zeros(nAzel, nK);
        for idir=1:nAzel
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
				pux = squeeze(spt{irep,1}(:,1,:));
				puy = squeeze(spt{irep,1}(:,2,:));
				puz = squeeze(spt{irep,1}(:,3,:));
                
                phase_dyn = -(qx.*pux + qy.*puy + qz.*puz);
                E_dyn = exp(1i * phase_dyn);
                rho_dyn = sum(E_dyn, 1);
                [F_curve, ~] = compute_acf(rho_dyn, max_lag);
                [D_val, G_val] = fit_dynamics(F_curve, S.timestep, k_val);
                
                seg_D(idir, k_idx) = D_val;
                seg_G(idir, k_idx) = G_val;
                F_AVG_STD(idir, k_idx, :) = F_AVG_STD(idir, k_idx, :) + reshape((F_curve), 1, 1, []);
                disp([irep,idir, k_idx])
            end
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
end