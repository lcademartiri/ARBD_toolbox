function [FSTD,GAMMASTD,DEFFSTD,FMASK,GAMMAMASK,DEFFMASK]=dynamics(S,sput,azel_rad,k_mags,max_lag,mask,N_eff)

    nAzel=size(azel_rad,1);
    nK = length(k_mags);
    [~,rotmat]=FCCrotate([1,0,0],[1,1,1]./norm([1,1,1]));
    rotmat=rotmat';
    
    % Scalar Maps (Keep Replicates for Error Bars)
        % Dimensions: Theta/Phi x K x Replicates
    DEFFSTD   = zeros(nAzel, nK, size(sput,1));
    GAMMASTD  = zeros(nAzel, nK, size(sput,1));
    
    % Full Time Correlation Function (Average Only)
    % Dimensions: Theta/Phi x K x Time
    % Using 'single' to save 50% RAM
    F_AVG_STD = zeros(nAzel, nK, max_lag+1);
    
    if S.bc == 2 || S.bc==3
        DEFFMASK  = zeros(nAzel, nK, size(sput,1));
        GAMMAMASK = zeros(nAzel, nK, size(sput,1));
        F_AVG_MASK= zeros(nAzel, nK, max_lag+1);
    end
    
    % SUPERLOOP
    for irep=1:size(sput,1)
        % initialize replicate storage
        seg_D = zeros(nAzel, nK);
        seg_G = zeros(nAzel, nK);
        if S.bc==2 || S.bc==3
            seg_Dm = zeros(nAzel, nK);
            seg_Gm = zeros(nAzel, nK);
        end
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
                if S.bc~=1
                    pux = squeeze(sput{irep,1}(:,1,:));
                    puy = squeeze(sput{irep,1}(:,2,:));
                    puz = squeeze(sput{irep,1}(:,3,:));
                else
                    pux = squeeze(spt{irep,1}(:,1,:));
                    puy = squeeze(spt{irep,1}(:,2,:));
                    puz = squeeze(spt{irep,1}(:,3,:));
                end
                
                phase_dyn = -(qx.*pux + qy.*puy + qz.*puz);
                E_dyn = exp(1i * phase_dyn);
                rho_dyn = sum(E_dyn, 1);
                [F_curve, ~] = compute_acf(rho_dyn, max_lag);
                [D_val, G_val] = fit_dynamics(F_curve, S.timestep, k_val);
                
                seg_D(idir, k_idx) = D_val;
                seg_G(idir, k_idx) = G_val;
                F_AVG_STD(idir, k_idx, :) = F_AVG_STD(idir, k_idx, :) + reshape((F_curve), 1, 1, []);
                
                if S.bc == 2 || S.bc==3
                    rho_dyn_m = sum(mask{irep,1} .* E_dyn, 1);
                    [F_curve_m, ~] = compute_acf(rho_dyn_m, max_lag);
                    [D_m, G_m] = fit_dynamics(F_curve_m, S.timestep, k_val);
                    seg_Dm(idir, k_idx) = D_m;
                    seg_Gm(idir, k_idx) = G_m;
                    F_AVG_MASK(idir, k_idx, :) = F_AVG_MASK(idir, k_idx, :) + reshape(single(F_curve_m), 1, 1, []);
                end
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