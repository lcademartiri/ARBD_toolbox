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

% parameters
phis=logspace(log10(1e-5),log10(1e-3),3)';
Ns=ceil(logspace(log10(1e3),log10(1e3),1)');
clms=(1)';

c=[];
q=1;
for i1=1:3
    for i2=1:1
        for i3=1:1
            phi=phis(i1);
            N=Ns(i2);
            clm=clms(i3);
            bv=(N*(4/3)*pi*S.rp^3)/phi;
            br=((3*bv)/(4*pi))^(1/3);
            cl=S.rc/clm+1e-16;
            c(q,:)=[S.rp,S.stdx,N,phi,br,bv,S.rc,cl,clm];
            q=q+1;
        end
    end
end
noconds=q-1;
clear phis Ns clms N phi clm br bv cl q i*

for ic=1:noconds
    S.N=c(ic,3);
    S.br=c(ic,5);
    S.cl=c(ic,8);
    S.icl=1/S.cl;
    S.clm=c(ic,9);
    DISP=build_noise_library(S.stdx,1e6);
    DISPcl=single(DISP./S.cl);
    % create cells padded to 2^n for morton indexing
    maxr=S.br+S.rc;
    cellcenters=(S.cl/2:S.cl:maxr)';
    if max(cellcenters-S.cl/2)<maxr
        cellcenters(end+1,:)=cellcenters(end,:)+S.cl;
    end
    padded_grid_size=2^nextpow2(numel(cellcenters)*2);
    padding_cells=(padded_grid_size/2)-numel(cellcenters);
    if padding_cells<1+S.clm
        padding_cells=1+S.clm;
        padded_grid_size=2*(numel(cellcenters)+padding_cells);
    end
    needed_bits = nextpow2(padded_grid_size); 
    padding=(max(cellcenters):S.cl:padding_cells*S.cl+max(cellcenters))';
    
    cellcenters=[cellcenters;padding(2:end,1)];
    cellcenters=sort(unique([-cellcenters;cellcenters]));
    celledges=[cellcenters(1)-S.cl/2;cellcenters+S.cl/2];
    % morton indices
    morton_grid_dim=[padded_grid_size,padded_grid_size,padded_grid_size];
    num_morton_cells = prod(morton_grid_dim);
    [X, Y, Z] = ndgrid(0:morton_grid_dim(1)-1, 0:morton_grid_dim(2)-1, 0:morton_grid_dim(3)-1);
    morton_idx = zeros(size(X));
    % FIX: Loop up to needed_bits-1 instead of hardcoded 9
    for b = 0:(needed_bits-1)
        morton_idx = bitor(morton_idx, bitshift(bitget(X, b+1), 3*b + 0));
        morton_idx = bitor(morton_idx, bitshift(bitget(Y, b+1), 3*b + 1));
        morton_idx = bitor(morton_idx, bitshift(bitget(Z, b+1), 3*b + 2));
    end
    % neighbor list
    [dX, dY, dZ] = meshgrid(-S.clm:S.clm, -S.clm:S.clm, -S.clm:S.clm);
    offsets = single([dX(:), dY(:), dZ(:)]);
    % offsets(ceil(size(offsets,1)/2),:)=[];
    no=size(offsets,1);
    if max(morton_idx,[],'all') > 2e9
        NeighborTable = zeros(num_morton_cells, no, 'uint32'); % Use UINT32
    else
        NeighborTable = zeros(num_morton_cells, no, 'int32'); 
    end
    for ix=padding_cells:padded_grid_size-padding_cells
        for iy=padding_cells:padded_grid_size-padding_cells
            for iz=padding_cells:padded_grid_size-padding_cells
                mi=morton_idx(ix,iy,iz);
                for io=1:no
                    nx = ix + offsets(io, 1);
                    ny = iy + offsets(io, 2);
                    nz = iz + offsets(io, 3);
                    neigh_morton_idx = morton_idx(nx, ny, nz);
                    NeighborTable(mi, io) = neigh_morton_idx;
                end                
            end
        end
    end
    vi=[find(abs(cellcenters-S.cl/2)<maxr,1)-1,find((cellcenters+S.cl/2)>maxr,1)]; % physical span of reachable cells
    morton_idx_phys=morton_idx(vi(1):vi(2),vi(1):vi(2),vi(1):vi(2));
    vmi=morton_idx_phys(:);
    % cleanup and organization
    
    cc_cl=single(cellcenters/S.cl); % cell center coordinates in single reduced dims
    ce_cl=single(celledges/S.cl); % cell center coordinates in single reduced dims
    vind=vi; % physical span of reachable cells
    vcc_cl=cc_cl(vind(1):vind(2),1); % cell center coordinates in single reduced dims
    vce_cl=ce_cl(vind(1):vind(2),1); % cell center coordinates in single reduced dims
    vmi=sort(vmi); % sorted list of physically reachable morton indices
    pmmi=morton_idx_phys; % map between coordinate indices and morton indices
    pmmig=size(pmmi,1);
    nl=NeighborTable; % neighbor list in morton indices;
    clear NeighborTable morton_idx_phys celledges cellcenters vi dX dY dZ X Y Z padding* io ix iy iz b neigh_morton_idx

    % generate starting particle configuration
    p=fillSphere(S.N,S.br,S.rp,0);
    p_cl=single(p.*S.icl);
    pcli=floor(p_cl);
    pcll=p_cl-pcli;
    offset_to_physical = 1+pmmig/2; % From your setup code
    ix = pcli(:,1) + offset_to_physical;
    iy = pcli(:,2) + offset_to_physical;
    iz = pcli(:,3) + offset_to_physical;
    lin_idx = ix + (iy-1)*pmmig + (iz-1)*pmmig^2;
    p_morton_id = pmmi(lin_idx);
    [sorted_p_mid, sort_perm] = sort(p_morton_id);
    sorted_p = p(sort_perm, :);
    sorted_p_cl = p_cl(sort_perm, :);
    sorted_pcli = pcli(sort_perm, :);
    sorted_pcll = pcll(sort_perm, :);
    max_morton_id = size(nl, 1);
    for it=1:1e2
        tic
        counts = (histcounts(sorted_p_mid, 1:(max_morton_id+1)))';
        starts = [1; cumsum(counts(1:end-1)) + 1];
        sorted_forces = zeros(S.N, 3, 'single');
        for k=1:length(vmi)
            i_cell = vmi(k);
            % Skip if cell is empty
            count_i = counts(i_cell);
            if count_i == 0, continue; end
            % 1. Get Cell i Data
            start_i = starts(i_cell);
            range_i = start_i : (start_i + count_i - 1);
             
            % 2. Identify Valid Neighbors
            neigh_ids = nl(i_cell, :);          % 1x27 vector of Morton IDs
            neigh_counts = counts(neigh_ids);   % 1x27 vector of particle counts
            
            % 3. Filter out empty neighbors immediately
            if sum(neigh_counts)==0, continue; end
            
            % 4. Harvest valid neighbors using reduced global coordinates
            valid_mask = neigh_counts > 0;
            valid_ids=neigh_ids(valid_mask);
            valid_starts=starts(valid_ids);
            
            valid_counts=counts(valid_ids);
            
            temp_map = repelem(1:length(valid_counts), valid_counts);
            temp_seq = (1:sum(valid_counts))';
            temp_cum_counts = [0; cumsum(valid_counts)];
            temp_adj = temp_cum_counts(temp_map); 
            valid_indices = valid_starts(temp_map) + temp_seq - temp_adj - 1;
            valid_shifts_base = offsets(valid_mask, :);
            expanded_shifts = repelem(valid_shifts_base, valid_counts, 1);
    
            self_local = sorted_pcll(range_i, :);
            neighbors_local = sorted_pcll(valid_indices, :);
            neighbors_projected = neighbors_local + expanded_shifts;
    
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
            f_scalar_val = 48 * S.pot_epsilon* (inv_r6.^2 - 0.5 * inv_r6) .* inv_r2;
            
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
                sorted_forces(range_i, d) = sorted_forces(range_i, d) + sum(f_comp_matrix, 2);
            end
        end
        disppot=sorted_forces.*S.fdxconv;
        idxovershoot=vecnorm(disppot,2,2)>S.rp;
        disppot(idxovershoot,:)=S.rp.*(disppot(idxovershoot,:)./vecnorm(disppot(idxovershoot,:)));
        displacements_cl=sorted_forces.*S.fdxconv*S.icl+DISPcl(randi(1e6,[S.N 1]),:);
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
        lin_idx = ix + (iy-1)*pmmig + (iz-1)*pmmig^2;
        p_morton_id = pmmi(lin_idx);
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