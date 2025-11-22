function [p,pgp,Nreal,idxgp] = ghostParticles_all(S,p,GPMAT)

    if S.bc==1
        % create antipodal ghost particles for the initial positions
        pnorms=vecnorm(p(:,1:3),2,2);
        idxgp=pnorms>(S.br-S.gtrig); % find if any particle touches border
        pgp=p(:,1:3)-(2*S.br).*(p(:,1:3)./pnorms);               
    elseif S.bc==2
        % create ghost particles for the initial positions
        idxgp=abs(p(:,1:3))>(S.br-S.gtrig); % find if any particle touches border
        % identification and classification of escape mode
        row=find(sum(idxgp,2)>0); % identify which particle touch border
        idxgp=idxgp.*sign(p(:,1:3)); % label of how particles are leaving the border for complete list of particles
        escape=double(idxgp(row,:)); % label of how particles are leaving the border only for those who are leaving
        gpmatrixtouse=14*ones(length(row),1)+escape(:,1)*9+escape(:,2)*3+escape(:,3); % find out, for each type of escape, which matrix of gp to use
        % calculation of ghost particle positions
        gpmattemp=GPMAT(gpmatrixtouse,:,:).*2.*S.br; % calculate ghost displacements (3D matrix with particle;xyz;# of ghost particles to build)
        gpi=p(row,1:3)+gpmattemp; % initialize ghost particle position matrix
        pgp=reshape(permute(gpi,[2 3 1]),[3 size(gpi,1)*7])';
        % code verified
    elseif S.bc==3
        N = size(p,1);
        % Convert original particles to fractional coordinates
        frac_coords_orig = p * S.fcc.invA;
        % identify particles that are within a interaction cutoff from the
        % walls
        is_candidate = any((frac_coords_orig(:,1:3) <= S.fcc.frac_thresh) | (frac_coords_orig(:,1:3) >= 1 - S.fcc.frac_thresh), 2);        
        % list of the boundary particles
        frac_candidates = frac_coords_orig(is_candidate, :);
        % 26-fold replication of the boundary particles, reshaped for
        % summation
        frac_candidates_replicated = repmat(frac_candidates, 26, 1);
        % adding unit cell displacements to the various replicas of
        % boundary particles
        frac_shifts_replicated = kron(S.fcc.shift_coeffs, ones(size(frac_candidates, 1), 1));
        frac_replicated_subset = frac_candidates_replicated(:,1:3) + frac_shifts_replicated;
        final_selection_idx = (all(frac_replicated_subset(:,1:3) >= -S.fcc.frac_thresh & frac_replicated_subset(:,1:3) <= (1 + S.fcc.frac_thresh), 2));
        pgp=(frac_replicated_subset(final_selection_idx,1:3)) * S.fcc.A;
        idxgp=is_candidate;
    end
    Nreal=size(p,1)+size(pgp,1);
    % ----

end