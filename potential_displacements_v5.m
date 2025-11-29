function [disppot, pairs_i, pairs_j, d_mic] = potential_displacements_v5(p, S, H, H_interpolant, ~)
   % Optimized specifically for CUBIC PBC (S.bc == 2).
% FIXED: Handles S.neighbor_linear as a cell array correctly.

% V4: HIGH-DENSITY OPTIMIZED CUBIC PBC.
% - Replaces sparse 'accumarray' with fast dense 'sum'.
% - Flattens neighbor cell-arrays to matrices for cache efficiency.
% - Handles N > 16,000 scaling issues.

if S.bc ~= 2
    error('This function only supports Cubic PBC (S.bc == 2)');
end

return_pairs = (nargout > 1);

N = size(p, 1);
L = 2*S.br;
invL = 1/L; 
rc = S.rc;
rc2 = rc^2;

% Potential constants
pot_r_min = H(1,1);
pot_r_max = H(end,1);
pot_F_min = H(1,2);

% -------------------------------------------------------------------------
% 1. NEIGHBOR OFFSETS
% -------------------------------------------------------------------------
ncell = S.ncell;
[offX, offY, offZ] = meshgrid(-1:1, -1:1, -1:1);
offsets = [offX(:), offY(:), offZ(:)]; % 27x3

% -------------------------------------------------------------------------
% 2. CELL ASSIGNMENT & SORTING
% -------------------------------------------------------------------------
p_shift = p + S.br; 
p_shift = p_shift - floor(p_shift .* invL) * L;

c_idx = floor(p_shift / S.cellsize) + 1;
c_idx(c_idx > ncell) = ncell;
c_idx(c_idx < 1) = 1;

cell_ids = c_idx(:,1) + ncell*(c_idx(:,2)-1) + ncell*ncell*(c_idx(:,3)-1);

% SORT (Structure of Arrays)
[sorted_cell_ids, sort_perm] = sort(cell_ids);
p_sorted = p(sort_perm, :);

cell_counts = accumarray(sorted_cell_ids, 1, [ncell^3, 1]);
cell_starts = [1; cumsum(cell_counts(1:end-1)) + 1]; 

% -------------------------------------------------------------------------
% 3. OUTPUT ARRAYS
% -------------------------------------------------------------------------
fx_sorted = zeros(N, 1);
fy_sorted = zeros(N, 1);
fz_sorted = zeros(N, 1);

if return_pairs
    est_pairs = N * 150; 
    pairs_i = zeros(est_pairs, 1);
    pairs_j = zeros(est_pairs, 1);
    d_mic_out = zeros(est_pairs, 3);
    pair_count = 0;
else
    pairs_i=[]; pairs_j=[]; d_mic=[];
end

% -------------------------------------------------------------------------
% 4. INTERACTION LOOP (Full Shell Gather)
% -------------------------------------------------------------------------
active_cells = find(cell_counts > 0);
n_active = length(active_cells);

% Convert grid indices
[cx_all, cy_all, cz_all] = ind2sub([ncell, ncell, ncell], active_cells);

for k = 1:n_active
    cid = active_cells(k);
    
    count_h = cell_counts(cid);
    start_h = cell_starts(cid);
    range_h = start_h : (start_h + count_h - 1);
    
    p_home = p_sorted(range_h, :);
    
    % Accumulator for this block
    f_home_x = zeros(count_h, 1);
    f_home_y = zeros(count_h, 1);
    f_home_z = zeros(count_h, 1);
    
    cx = cx_all(k); cy = cy_all(k); cz = cz_all(k);
    
    for n = 1:27
        % Neighbor ID with PBC
        nx = cx + offsets(n,1);
        ny = cy + offsets(n,2);
        nz = cz + offsets(n,3);
        
        if nx > ncell, nx = 1; elseif nx < 1, nx = ncell; end
        if ny > ncell, ny = 1; elseif ny < 1, ny = ncell; end
        if nz > ncell, nz = 1; elseif nz < 1, nz = ncell; end
        
        nid = nx + ncell*(ny-1) + ncell*ncell*(nz-1);
        
        if cell_counts(nid) == 0, continue; end
        
        count_n = cell_counts(nid);
        start_n = cell_starts(nid);
        range_n = start_n : (start_n + count_n - 1);
        p_neigh = p_sorted(range_n, :);
        
        % Compute Diffs
        dx = p_home(:,1) - p_neigh(:,1).';
        dy = p_home(:,2) - p_neigh(:,2).';
        dz = p_home(:,3) - p_neigh(:,3).';
        
        % MIC
        dx = dx - L .* round(dx .* invL);
        dy = dy - L .* round(dy .* invL);
        dz = dz - L .* round(dz .* invL);
        
        r2 = dx.^2 + dy.^2 + dz.^2;
        
        % Mask: r <= rc AND r > 0 (to avoid self i=i)
        mask = (r2 <= rc2) & (r2 > 0); 
        
        if ~any(mask(:)), continue; end
        
        [fx_mat, fy_mat, fz_mat] = get_force_mat(dx, dy, dz, r2, mask, H, H_interpolant, pot_r_min, pot_r_max, pot_F_min, rc);
        
        % Accumulate on HOME only
        f_home_x = f_home_x + sum(fx_mat, 2);
        f_home_y = f_home_y + sum(fy_mat, 2);
        f_home_z = f_home_z + sum(fz_mat, 2);
        
        if return_pairs
            % --- PAIR EXTRACTION (ROBUST) ---
            [loc_i, loc_j] = find(mask);
            
            % Map to global sorted indices
            % FORCE COLUMN VECTORS immediately for safety
            g_i = range_h(loc_i); g_i = g_i(:);
            g_j = range_n(loc_j); g_j = g_j(:);
            
            % Unique Pair Filter (Store i < j only)
            valid_p = g_i < g_j;
            
            if any(valid_p)
                % Filter indices
                loc_i = loc_i(valid_p); 
                loc_j = loc_j(valid_p);
                g_i = g_i(valid_p); 
                g_j = g_j(valid_p);
                
                n_add = length(loc_i);
                rng = pair_count+1 : pair_count+n_add;
                pairs_i(rng) = sort_perm(g_i);
                pairs_j(rng) = sort_perm(g_j);
                
                % Extract Distances
                lin_idx = sub2ind(size(dx), loc_i, loc_j);
                
                % CRITICAL FIX: Explicitly reshape extracted values to columns
                % This works even if dx(lin_idx) returns a row vector.
                vx = dx(lin_idx); vx = vx(:);
                vy = dy(lin_idx); vy = vy(:);
                vz = dz(lin_idx); vz = vz(:);
                
                d_mic_out(rng, :) = [vx, vy, vz];
                
                pair_count = pair_count + n_add;
            end
        end
    end
    
    % Write to Global Array
    fx_sorted(range_h) = f_home_x;
    fy_sorted(range_h) = f_home_y;
    fz_sorted(range_h) = f_home_z;
end

% -------------------------------------------------------------------------
% 5. FINALIZE
% -------------------------------------------------------------------------
fx_final = zeros(N, 1); fx_final(sort_perm) = fx_sorted;
fy_final = zeros(N, 1); fy_final(sort_perm) = fy_sorted;
fz_final = zeros(N, 1); fz_final(sort_perm) = fz_sorted;

potdisp = [fx_final, fy_final, fz_final] * (S.esdiff / S.kbT) * S.timestep;

maxstep = S.pot_clamp * sqrt(3) * S.stdx;
norms = sqrt(sum(potdisp.^2, 2));
overshoot = norms > maxstep;
if any(overshoot)
    potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
end
disppot = potdisp;

if return_pairs
    pairs_i = pairs_i(1:pair_count);
    pairs_j = pairs_j(1:pair_count);
    d_mic = d_mic_out(1:pair_count, :);
end

end

% -------------------------------------------------------------------------
% HELPER
% -------------------------------------------------------------------------
function [fx, fy, fz] = get_force_mat(dx, dy, dz, r2, mask, H, H_interpolant, p_min, p_max, f_min, rc)
    r = zeros(size(r2));
    r(mask) = sqrt(r2(mask));
    
    r_c = r;
    r_c(r_c < p_min) = p_min;
    r_c(r_c > p_max) = p_max;
    
    Fij = zeros(size(r));
    Fij(mask) = H_interpolant(r_c(mask));
    
    Fij(r >= rc | r >= p_max) = 0;
    
    bad_low = isnan(Fij) & (r < p_min) & mask;
    if any(bad_low(:)), Fij(bad_low) = f_min; end
    
    inv_r = zeros(size(r));
    inv_r(mask) = 1 ./ r(mask);
    
    f_scalar = Fij .* inv_r; 
    
    fx = f_scalar .* dx;
    fy = f_scalar .* dy;
    fz = f_scalar .* dz;
end