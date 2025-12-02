clear all
%% cell list testing
S.timestep=1.1411e-8;
S.esdiff=2.4526e-11;
S.alpha=0.7045;
S.rp=1e-8;
S.stdx=6.2793e-10;
S.pot_sigma=2*S.rp+0.25e-9;
S.kbT=1.38e-23*298;
S.pot_epsilon=S.kbT;
S.rc=3*S.pot_sigma;
S.fdxconv=(S.esdiff*S.alpha/S.kbT)*S.timestep;
S.bc=1;
[~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /value');
tokens = regexp(cmdout, '\d+', 'match');
S.cacheSizeMB = (max(str2double(tokens))/1024)/feature('numCores');
S.opt_ppc=floor(0.5*sqrt((S.cacheSizeMB*1024^2)/(27*36)));

% parameters
phis=logspace(log10(1e-5),log10(1e-1),5)';
Ns=ceil(logspace(log10(1e2),log10(1e5),5)');

c=[];
q=1;
for i1=1:numel(phis)
    for i2=1:numel(Ns)
            phi=phis(i1);
            N=Ns(i2);
            pv=(4/3)*pi*S.rp^3;
            bv=(N*pv)/phi;
            br=((3*bv)/(4*pi))^(1/3);
            ndens=N/bv;
            cv=S.opt_ppc/ndens;
            cl=cv^(1/3);
            if cl<S.rc
                cl=S.rc;
            end
            c(q,:)=[S.rp,S.stdx,N,phi,br,bv,S.rc,cl];
            q=q+1;
    end
end
noconds=q-1;
clear phis Ns clms N phi clm br bv cl q i*

for ic=1:noconds
    S.N=c(ic,3);
    S.br=c(ic,5);
    S.cl=c(ic,8);
    S.icl=1/S.cl;
    S.clm=1;
    % displacement library
    DISP=build_noise_library(S.stdx,1e6);
    DISPcl=single(DISP./S.cl);
    % create cells padded to 2^n for morton indexing
    maxr=S.br+S.rc;
    if S.bc==1 || S.bc==4
        pcs=2*(ceil(maxr/S.cl)+1);
    else
        pcs=2*ceil(maxr/S.cl);
    end

    lcs=2^nextpow2(pcs);
    padding_cells=(lcs-pcs)/2;
    lce=(0:S.cl:0.5*lcs*S.cl)';
    lce=sort(unique([-lce;lce]));
    lcc=lce(1:end-1,:)+S.cl/2;
    needed_bits = nextpow2(lcs); 
    
    % morton indices calculated on the cellcenters^3 meshgrid which
    % includes padding 
    morton_grid_dim=[lcs,lcs,lcs];
    num_morton_cells = prod(morton_grid_dim);
    [X, Y, Z] = ndgrid(0:morton_grid_dim(1)-1, 0:morton_grid_dim(2)-1, 0:morton_grid_dim(3)-1);
    MI = zeros(size(X));
    for b = 0:(needed_bits-1)
        MI = bitor(MI, bitshift(bitget(X, b+1), 3*b + 0));
        MI = bitor(MI, bitshift(bitget(Y, b+1), 3*b + 1));
        MI = bitor(MI, bitshift(bitget(Z, b+1), 3*b + 2));
    end
    MI=MI+1;
    % calculate cell cartesian offsets within rc - include self
    [dX, dY, dZ] = meshgrid(-S.clm:S.clm, -S.clm:S.clm, -S.clm:S.clm);
    offsets = single([dX(:), dY(:), dZ(:)]);
    no=size(offsets,1);
    if max(MI,[],'all') > 2e9
        NL = zeros(num_morton_cells, no, 'uint32'); % Use UINT32
    else
        NL = zeros(num_morton_cells, no, 'int32'); 
    end
    % determine neighbor table only for physical meshgrid (without padding)
    idxphys=histcounts([-maxr,maxr],lce);
    vi=find(idxphys); 

    for ix=vi(1):vi(2)
        for iy=vi(1):vi(2)
            for iz=vi(1):vi(2)
                mi=MI(ix,iy,iz); % Morton index of center cell
                for io=1:no % loop over all cells in rc range
                    nx = ix + offsets(io, 1);
                    ny = iy + offsets(io, 2);
                    nz = iz + offsets(io, 3);
                    % append Morton index of the cells in rc range at row 
                    % associated with Morton index of the center cell
                    neigh_morton_idx = MI(nx, ny, nz);
                    NL(mi, io) = neigh_morton_idx; 
                end                
            end
        end
    end
    MIP=MI(vi(1):vi(2),vi(1):vi(2),vi(1):vi(2));
    XP=X(vi(1):vi(2),vi(1):vi(2),vi(1):vi(2))+1;
    XP=XP(:);
    YP=Y(vi(1):vi(2),vi(1):vi(2),vi(1):vi(2))+1;
    YP=YP(:);
    ZP=Z(vi(1):vi(2),vi(1):vi(2),vi(1):vi(2))+1;
    ZP=ZP(:);
    vXYZP=[XP,YP,ZP];
    vMIP=MIP(:);
    
    % cleanup and organization
    
    cc_cl=single(lcc/S.cl); % cell center coordinates in single reduced dims
    ce_cl=single(lce/S.cl); % cell center coordinates in single reduced dims
    cc=double(cc_cl*S.cl); % cell center coordinates in double real dims
    ce=double(ce_cl*S.cl); % cell center coordinates in double real dims
    ccp=cc(vi(1):vi(2),1); % physical cell center coordinates in double real dims
    ccp_cl=single(ccp.*S.icl);  % physical cell center coordinates in single reduced dims
    vMIP=sort(vMIP); % sorted list of physically reachable morton indices
    b=size(MIP,1); % number of physical cells
    NL=NL; % neighbor list in morton indices;
    vi=vi; % map of physical cells on logical cells

    % generate starting particle configuration
    p=fillSphere(S.N,S.br,S.rp,0);
    %%% DUMMY DEBUGGING P
        % p=[ceil(mean(vi)),ceil(mean(vi)),ceil(mean(vi))]+offsets;
        % p=ccp(p);
    %%%
    p_cl=single(p.*S.icl);
    pcli=floor(p_cl);
    pcll=p_cl-pcli;
    offset_to_physical = 1+b/2; % From your setup code
    ix = pcli(:,1) + offset_to_physical;
    iy = pcli(:,2) + offset_to_physical;
    iz = pcli(:,3) + offset_to_physical;
    lin_idx = ix + (iy-1)*b + (iz-1)*b^2; % correct
    p_morton_id = MIP(lin_idx);
    
    [sorted_p_mid, sort_perm] = sort(p_morton_id);
    sorted_p = p(sort_perm, :);
    sorted_p_cl = p_cl(sort_perm, :);
    sorted_pcli = pcli(sort_perm, :);
    sorted_pcll = pcll(sort_perm, :);
    sorted_p_morton_id=sort(p_morton_id);
    max_morton_id = size(NL, 1);
    for it=1:1e2
        tic
        counts = (histcounts(sorted_p_mid, 1:(max_morton_id+1)))';
        starts = [1; cumsum(counts(1:end-1)) + 1];
        sorted_forces = zeros(size(p,1), 3, 'single');
        for k=1:length(vMIP)
            i_cell = vMIP(k);
            % Skip if cell is empty
            count_i = counts(i_cell);
            if count_i == 0, continue; end
            % 1. Get Cell i Data
            start_i = starts(i_cell);
            range_i = start_i : (start_i + count_i - 1);
             
            % 2. Identify Valid Neighbors
            neigh_ids = NL(i_cell, :)';          % 1x27 vector of Morton IDs
            neigh_counts = counts(neigh_ids);   % 1x27 vector of particle counts
            
            % 3. Filter out empty neighbors immediately
            if sum(neigh_counts)==0, continue; end
            
            % 4. Harvest valid neighbors using reduced global coordinates
            valid_mask = neigh_counts > 0;
            valid_ids=neigh_ids(valid_mask);
            valid_starts=starts(valid_ids);
            valid_counts=counts(valid_ids);
            
            temp_map = repelem(1:length(valid_counts), valid_counts)';
            temp_seq = (1:sum(valid_counts))';
            temp_cum_counts = [0; cumsum(valid_counts)];
            temp_adj = temp_cum_counts(temp_map); 
            valid_indices = valid_starts(temp_map) + temp_seq - temp_adj - 1;
            valid_shifts_base = offsets(valid_mask, :);
            expanded_shifts = repelem(valid_shifts_base, valid_counts, 1);
    
            self_local = sorted_pcll(range_i, :);
            neighbors_local = sorted_pcll(valid_indices, :);
            neighbors_projected = neighbors_local + expanded_shifts;
            %%% VERIFIED AS WORKING IN SBC TILL HERE
            dr_red = reshape(neighbors_projected, [1, numel(valid_indices), 3]) - ...
            reshape(self_local, [count_i, 1, 3]);
            dr_phys = double(dr_red * S.cl);
            r2 = sum(dr_phys.^2, 3);
    
            % --- 7. PHYSICS CALCULATION ---
    
            mask = (r2 < S.rc^2) & (r2 > 1e-30);
            if ~any(mask(:)), continue; end
            r2 = r2(mask); % Vector of valid squared distances
            r2_reduced = r2 / (S.pot_sigma^2);
            min_safe_r2 = (0.5)^2; 
            r2_reduced = max(r2_reduced, min_safe_r2);
        
            % Example: Lennard-Jones
            inv_r2 = 1 ./ r2_reduced;
            inv_r6 = inv_r2.^3;
            
            % Force Scalar (Magnitude / r)
            f_scalar_val = 48 * S.pot_epsilon* (inv_r6.^2 - 0.5 * inv_r6) .* (inv_r2./S.pot_sigma^2);
            
            % --- 8. FORCE ACCUMULATION ---
            % We only write to 'sorted_forces(range_i)'. 
            % We do NOT write back to the neighbors.
            % This linear memory write is optimal for cache.
            
            for d = 1:3
                % Get the coordinate component (dx, dy, or dz)
                % d_comp is automatically squeezed to 2D: [Ni, Nj]
                d_comp = dr_phys(:,:,d);
                
                % --- FIX: Initialize using 2D size ---
                % We use size(d_comp) instead of size(r2) to ensure it is 2D
                f_comp_matrix = zeros(size(d_comp), 'single');
                
                % Map the calculated forces back to their positions
                % (Linear indexing works fine even if mask is 3D and matrix is 2D)
                f_comp_matrix(mask) = f_scalar_val .* d_comp(mask);
                
                % Sum across neighbors (columns)
                % Result is [Ni, 1] (2D), which matches sorted_forces dimensions
                sorted_forces(range_i, d) = sorted_forces(range_i, d) - sum(f_comp_matrix, 2);
            end
        end
        disppot=sorted_forces.*S.fdxconv;
        idxovershoot=vecnorm(disppot,2,2)>S.rp;
        disppot(idxovershoot,:)=S.rp.*(disppot(idxovershoot,:)./vecnorm(disppot(idxovershoot,:)));
        displacements_cl=disppot.*S.fdxconv*S.icl+DISPcl(randi(1e6,[size(p,1) 1]),:);
        sorted_p_cl = sorted_p_cl + displacements_cl;
        sorted_p = sorted_p_cl * S.cl;
        pnorms=vecnorm(sorted_p,2,2);
        idxwrap=pnorms>S.br;
        if any(idxwrap)
            sorted_p(idxwrap,:)=sorted_p(idxwrap,:)-2*S.br.*(sorted_p(idxwrap,:)./pnorms(idxwrap,:));
            sorted_p_cl(idxwrap,:)=sorted_p(idxwrap,:).*S.icl;
        end
        pcli = floor(sorted_p_cl); 
        pcll = sorted_p_cl - pcli;
        ix = pcli(:,1) + offset_to_physical;
        iy = pcli(:,2) + offset_to_physical;
        iz = pcli(:,3) + offset_to_physical;
        lin_idx = ix + (iy-1)*b + (iz-1)*b^2;
        p_morton_id = MIP(lin_idx);
        [sorted_p_mid, sort_perm] = sort(p_morton_id);
        sorted_p = p(sort_perm, :);
        sorted_p_cl = p_cl(sort_perm, :);
        sorted_pcli = pcli(sort_perm, :);
        sorted_pcll = pcll(sort_perm, :);
        if mod(it,10)==0
            disp([ic,log10(it)])
        end
        TIMES(ic,it)=toc;
    end
end