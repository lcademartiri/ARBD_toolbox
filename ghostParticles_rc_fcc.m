function [p,Nreal] = ghostParticles_rc_fcc(S,p,GPMAT,displacements)
    if S.bc==1
        % create ghost particles for the initial positions
        idxgpi=vecnorm(p(:,1:3),2,2)>(S.br-S.gtrig); % find if any particle touches border
        gpi=[];
        if sum(idxgpi)>0 % if any particle is touching the boundary
            finaltemp=p(idxgpi,8:10);
            initialtemp=p(idxgpi,1:3);
            final=finaltemp-(finaltemp./vecnorm(finaltemp,2,2))*2*S.br; % calculate the post-displacement positions of the ghost particles
            initial=initialtemp-((initialtemp./vecnorm(initialtemp,2,2))*2*S.br); % calculate the pre-displacement positions of the ghost particles
            ids=p(idxgpi,4); % find the ids of the real particles touching the boundary
            dp=final-initial; % calculate the displacements of the ghost particles
            gpi=[initial,ids,dp,final]; % compiling the position matrix for the ghost particles
        end
        % create ghost particles for the final positions
        idxgpf=vecnorm(p(:,8:10),2,2)>(S.br-S.gtrig); % find if any particle touches border
        gpf=[];
        if sum(idxgpf)>0 % if any particle is touching the boundary
            finaltemp=p(idxgpf,8:10);
            initialtemp=p(idxgpf,1:3);
            final=finaltemp-(finaltemp./vecnorm(finaltemp,2,2))*2*S.br; % calculate the post-displacement positions of the ghost particles
            initial=initialtemp-((initialtemp./vecnorm(initialtemp,2,2))*2*S.br); % calculate the pre-displacement positions of the ghost particles
            ids=p(idxgpf,4); % find the ids of the real particles touching the boundary
            dp=final-initial; % calculate the displacements of the ghost particles
            gpf=[initial,ids,dp,final]; % compiling the position matrix for the ghost particles
        end
    elseif S.bc==2
        % create ghost particles for the initial positions
        idxgpi=abs(p(:,1:3))>(S.br-S.gtrig); % find if any particle touches border
        gpi=[];
        if sum(sum(idxgpi))>0
            % identification and classification of escape mode
            row=find(sum(idxgpi,2)>0); % identify which particle touch border
            idxgpi=idxgpi.*sign(p(:,1:3)); % label of how particles are leaving the border for complete list of particles
            escape=double(idxgpi(row,:)); % label of how particles are leaving the border only for those who are leaving
            gpmatrixtouse=14*ones(length(row),1)+escape(:,1)*9+escape(:,2)*3+escape(:,3); % find out, for each type of escape, which matrix of gp to use
            % calculation of ghost particle positions
            gpmattemp=GPMAT(gpmatrixtouse,:,:).*2.*S.br; % calculate ghost displacements (3D matrix with particle;xyz;# of ghost particles to build)
            gpi=p(row,1:3)+gpmattemp; % initialize ghost particle position matrix
            gpi=reshape(permute(gpi,[2 3 1]),[3 size(gpi,1)*7])';
            % calculation of ghost particle displacements
            idtemp=round(p(row,4)+((squeeze(gpmattemp(:,1,:))).*0+1)-1,0); % find ids of ghosted particles
            idtemp=reshape(idtemp',[1 size(idtemp,1)*size(idtemp,2)])'; % reshape idtemp matrix
            dptemp=zeros(size(idtemp,1),3); % initialize temporary displacement matrix
            dptemp(isnan(idtemp)==0,:)=displacements(idtemp(isnan(idtemp)==0),:); % compile displacements for all the ghost particles
            % cleanup of temp matrices and assembly
            idnan=isnan(idtemp)==1;
            gpi(idnan,:)=[];
            idtemp(idnan,:)=[];
            dptemp(idnan,:)=[];
            gpi=[gpi,idtemp,dptemp,gpi+dptemp];
            % code verified
        end
        % create ghost particles for the final positions
        idxgpf=abs(p(:,8:10))>(S.br-S.gtrig); % find if any particle touches border
        gpf=[];
        if sum(sum(idxgpf))>0
            % identification and classification of escape mode
            row=find(sum(idxgpf,2)>0); % identify which particle touch border
            idxgpf=idxgpf.*sign(p(:,8:10)); % label of how particles are leaving the border for complete list of particles
            escape=double(idxgpf(row,:)); % label of how particles are leaving the border only for those who are leaving
            gpmatrixtouse=14*ones(length(row),1)+escape(:,1)*9+escape(:,2)*3+escape(:,3); % find out, for each type of escape, which matrix of gp to use
            % calculation of ghost particle positions
            gpmattemp=GPMAT(gpmatrixtouse,:,:).*2.*S.br; % calculate ghost displacements (3D matrix with particle;xyz;# of ghost particles to build)
            gpf=p(row,8:10)+gpmattemp; % initialize ghost particle position matrix
            gpf=reshape(permute(gpf,[2 3 1]),[3 size(gpf,1)*7])';
            % calculation of ghost particle displacements
            idtemp=round(p(row,4)+((squeeze(gpmattemp(:,1,:))).*0+1)-1,0); % find ids of ghosted particles
            % idtemp=round(p(row,4)+((squeeze(gpmattemp(:,1,:))).*0+1)-1,0); % find ids of ghosted particles
            idtemp=reshape(idtemp',[1 size(idtemp,1)*size(idtemp,2)])'; % reshape idtemp matrix
            dptemp=zeros(size(idtemp,1),3); % initialize temporary displacement matrix
            dptemp(isnan(idtemp)==0,:)=displacements(idtemp(isnan(idtemp)==0),:); % compile displacements for all the ghost particles
            % cleanup of temp matrices and assembly
            idnan=isnan(idtemp)==1;
            gpf(idnan,:)=[];
            idtemp(idnan,:)=[];
            dptemp(idnan,:)=[];
            gpf=[gpf-dptemp,idtemp,dptemp,gpf];
            % code verified
        end
        idtemp=[];
        dptemp=[];
    % ----
    elseif S.bc==3

        N = size(p,1);
        
        % INITIAL ghost particles
        % Convert original particles to fractional coordinates
        frac_coords_orig = [p(:,1:3) * S.fcc.invA,p(:,4)];

        % --- Identify Candidate Particles and Required Shifts ---
        % Find original particles that are close to any face. These are our candidates.
        % A particle is near a "low" face if its fractional coordinate c < thresh.
        % A particle is near a "high" face if its fractional coordinate c > 1 - thresh.
        is_candidate = any((frac_coords_orig(:,1:3) <= S.fcc.frac_thresh) | (frac_coords_orig(:,1:3) >= 1 - S.fcc.frac_thresh), 2);        
        frac_candidates = frac_coords_orig(is_candidate, :);
        frac_candidates_replicated = repmat(frac_candidates, 26, 1);
        frac_shifts_replicated = kron(S.fcc.shift_coeffs, ones(size(frac_candidates, 1), 1));
        frac_replicated_subset = [frac_candidates_replicated(:,1:3) + frac_shifts_replicated,frac_candidates_replicated(:,4)];
        final_selection_idx = (all(frac_replicated_subset(:,1:3) >= -S.fcc.frac_thresh & frac_replicated_subset(:,1:3) <= (1 + S.fcc.frac_thresh), 2));
        gpi = zeros(sum(final_selection_idx),4);   % preallocate
        gpi=[(frac_replicated_subset(final_selection_idx,1:3)) * S.fcc.A,frac_replicated_subset(final_selection_idx,4)];
        gpi=[gpi,p(gpi(:,4),5:7),gpi(:,1:3)+p(gpi(:,4),5:7)];

        % FINAL ghost particles
        % Convert original particles to fractional coordinates
        frac_coords_orig = [p(:,8:10) * S.fcc.invA,p(:,4)];

        % --- Identify Candidate Particles and Required Shifts ---
        % Find original particles that are close to any face. These are our candidates.
        % A particle is near a "low" face if its fractional coordinate c < thresh.
        % A particle is near a "high" face if its fractional coordinate c > 1 - thresh.
        is_candidate = any((frac_coords_orig(:,1:3) <= S.fcc.frac_thresh) | (frac_coords_orig(:,1:3) >= 1 - S.fcc.frac_thresh), 2);
        frac_candidates = frac_coords_orig(is_candidate, :);
        frac_candidates_replicated = repmat(frac_candidates, 26, 1);
        frac_shifts_replicated = kron(S.fcc.shift_coeffs, ones(size(frac_candidates, 1), 1));
        frac_replicated_subset = [frac_candidates_replicated(:,1:3) + frac_shifts_replicated,frac_candidates_replicated(:,4)];
        final_selection_idx = (all(frac_replicated_subset(:,1:3) >= -S.fcc.frac_thresh & frac_replicated_subset(:,1:3) <= (1 + S.fcc.frac_thresh), 2));
        gpf = zeros(sum(final_selection_idx),4);   % preallocate
        gpf=[(frac_replicated_subset(final_selection_idx,1:3)) * S.fcc.A,frac_replicated_subset(final_selection_idx,4)];
        gpf=[gpf(:,1:3)-p(gpf(:,4),5:7),gpf(:,4),p(gpf(:,4),5:7),gpf(:,1:3)];
        p=[p;gpi;gpf];
        Nreal=size(p,1);
    end
    p=[p;gpi;gpf];
    % --- remove ghost particle duplicates
    if size(p,1)>S.N
        idxdupl=((p(S.N+1:end,1)-p(S.N+1:end,1)').^2+(p(S.N+1:end,2)-p(S.N+1:end,2)').^2+(p(S.N+1:end,3)-p(S.N+1:end,3)').^2)<(2*S.rp/100)^2;
        [row,col]=find(idxdupl);
        idxdupl=[row,col];
        idxdupl(idxdupl(:,1)>=idxdupl(:,2),:)=[];
        if isempty(idxdupl)==0
            p(idxdupl(:,2)+S.N,:)=[];
        end
    end
    Nreal=size(p,1);
end