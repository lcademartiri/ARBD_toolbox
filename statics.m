function [SSTD,SMASK]=statics(S,spt,azel_rad,k_mags,mask,N_eff)
    nAzel=size(azel_rad,1);
    nK = length(k_mags);
    
    % Scalar Maps (Keep Replicates for Error Bars)
        % Dimensions: Theta/Phi x K x Replicates
    SSTD      = zeros(nAzel, nK, size(spt,1));    
    if S.bc == 2 || S.bc==3
        SMASK     = zeros(nAzel, nK, size(spt,1));
    end
    [~,rotmat]=FCCrotate([1,0,0],[1,1,1]./norm([1,1,1]));
    rotmat=rotmat';
    % SUPERLOOP
    for irep=1:size(spt,1)
        % initialize replicate storage
        seg_S = zeros(nAzel, nK);
        if S.bc==2 || S.bc==3
            seg_Sm = zeros(nAzel, nK);
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
                
                % -- STATIC S(k) --
                % Use wrapped coords (pt)
                px = squeeze(spt{irep,1}(:,1,:)); 
                py = squeeze(spt{irep,1}(:,2,:)); 
                pz = squeeze(spt{irep,1}(:,3,:));
                
                % Dot product q*r
                phase = -(qx.*px + qy.*py + qz.*pz);
                E = exp(1i * phase);
                
                % Standard
                rho_stat = sum(E, 1);
                seg_S(idir, k_idx) = mean(abs(rho_stat).^2) / S.N;
                
                if S.bc == 2 || S.bc==3
                    rho_mask = sum(mask{irep,1} .* E, 1);
                    seg_Sm(idir, k_idx) = mean(abs(rho_mask).^2) / N_eff{irep,1};
                end
                disp([irep,idir, k_idx])
            end
        end
        % Store segment
        SSTD(:,:,irep) = seg_S;
        if S.bc == 2 || S.bc==3
            SMASK(:,:,irep) = seg_Sm;
        end
    end
end