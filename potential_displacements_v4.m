function [disppot, pairs_i, pairs_j, d_mic_out] = potential_displacements_v4(p, S, H, H_interpolant, ghostghost)
% Optimized V4: Block-processing + Fast MIC
% drastically reduces memory bandwidth usage.

if nargin < 5, ghostghost = 0; end

N = size(p,1);
L = 2*S.br;
halfL = S.br; % L/2
rc = S.rc;
rc2 = rc^2;

% Prepare potential constants
pot_r_min = H(1,1);
pot_r_max = H(end,1);
pot_F_min = H(1,2);

% Pre-constants for force conversion
force_scale = (S.esdiff / S.kbT) * S.timestep;
maxstep = S.pot_clamp * sqrt(3) * S.stdx;

% Allocate outputs
fx_total = zeros(N,1);
fy_total = zeros(N,1);
fz_total = zeros(N,1);

% Optional: Accumulate pairs if requested (Note: keeping track of pairs slows things down)
return_pairs = (nargout > 1);
if return_pairs
    % Preallocate estimation (max 100 neighbors per particle)
    est_neighbors = N * 50; 
    pairs_i = zeros(est_neighbors, 1);
    pairs_j = zeros(est_neighbors, 1);
    d_mic_out = zeros(est_neighbors, 3);
    pair_count = 0;
end

% BLOCK SIZE: Tunable parameter. 
% 500 to 1000 is usually optimal for L2/L3 cache retention.
BLOCK_SIZE = 512; 

px = p(:,1); 
py = p(:,2); 
pz = p(:,3);

% Loop over blocks of i
for i_start = 1:BLOCK_SIZE:N
    i_end = min(i_start + BLOCK_SIZE - 1, N);
    
    % Indices for this block
    idx_i = i_start:i_end; 
    
    % Get coordinates for this block (Size: B x 1)
    bx = px(idx_i);
    by = py(idx_i);
    bz = pz(idx_i);
    
    % Compute differences against ALL particles (Size: B x N)
    % Using implicit expansion (MATLAB R2016b+)
    dx = bx - px.'; 
    dy = by - py.'; 
    dz = bz - pz.';
    
    % ---  MIC ---
    % Only works if particles have not moved > L/2 in one step (safe assumption)
    dx = dx - L .* round(dx ./ L);
    dy = dy - L .* round(dy ./ L);
    dz = dz - L .* round(dz ./ L);
    
    % --- Distance Squared ---
    r2 = dx.^2 + dy.^2 + dz.^2;
    
    % --- Filter ---
    % 1. Within cutoff
    % 2. Avoid self-interaction (i > j) to do Upper Triangle only
    %    We create a mask of indices.
    
    % Create a column index matrix for the block
    j_indices = 1:N;
    
    % We only want j > i to avoid double counting and self-interaction
    % This logic is tricky in blocks, so simpler approach:
    % Calculate full force, divide by 2? No, forces are directional.
    % Newton's 3rd law: Force on i by j is -Force on j by i.
    % To vectorize efficiently, we usually compute full N x N logic in blocks 
    % but only act where r2 < rc2.
    
    mask = (r2 < rc2) & (r2 > 0);
    
    % To implement Newton's 3rd law (Upper Triangle) efficiently in blocks:
    % It is often faster to compute ALL interactions (double work) and sum 
    % than to handle the complex indexing of the upper triangle in blocks, 
    % UNLESS we only need to write to the 'i' accumulator.
    
    % Let's stick to computing the force on 'i' from all 'j'.
    % This means we calculate ij and ji separately. 
    % This is O(N^2) work, but avoids 'accumarray' which is slow.
    % The previous code used upper triangle + accumarray. 
    % Here, we just sum directly into fx_total(idx_i).
    
    if ~any(mask(:)), continue; end
    
    % Extract valid pairs for this block
    % (Linear indexing within the block matrices)
    valid_idx = find(mask);
    
    if isempty(valid_idx), continue; end
    
    r_val = sqrt(r2(valid_idx));
    
    dx_val = dx(valid_idx);
    dy_val = dy(valid_idx);
    dz_val = dz(valid_idx);
    
    % --- 7) Compute Force Magnitudes ---
    % Clamp r (Your question: Yes, keep this for interpolation safety)
    r_clamped = r_val;
    r_clamped(r_clamped < pot_r_min) = pot_r_min;
    r_clamped(r_clamped > pot_r_max) = pot_r_max;
    
    Fij = H_interpolant(r_clamped);
    
    % Handle Out of bounds / NaN manually
    % (If your interpolant is robust, you can skip some of these)
    Fij(r_val >= pot_r_max) = 0;
    
    % Singularity protection (r < min) -> Cap force
    bad_low = isnan(Fij) & (r_val < pot_r_min);
    if any(bad_low)
        Fij(bad_low) = pot_F_min;
    end
    
    % Projection: F_vector = (F_scalar / r) * vec_r
    inv_r = 1 ./ r_val;
    force_factor = Fij .* inv_r;
    
    f_x_local = force_factor .* dx_val;
    f_y_local = force_factor .* dy_val;
    f_z_local = force_factor .* dz_val;
    
    % --- Accumulate Forces ---
    % We need to map linear indices (block) back to (i,j)
    % valid_idx are indices in the B x N matrix.
    % row (i) = rem(valid_idx-1, B) + 1
    % col (j) = floor((valid_idx-1)/B) + 1
    
    [row_sub, col_sub] = ind2sub(size(mask), valid_idx);
    
    % Map row_sub (1..B) back to real particle index (i_start...i_end)
    real_i = idx_i(row_sub)';
    
    % Accumulate forces on i
    % Since we have multiple j's acting on the same i, we sum them up.
    % For the block, we can just sum rows?
    % Actually, simpler:
    
    % Sum forces for each particle 'i' in this block
    % reshape to B x N (sparse-like) and sum over dim 2? No, mask is flat.
    
    % Fastest way in MATLAB without accumarray loop:
    % Initialize block forces
    bfx = zeros(length(idx_i), 1);
    bfy = zeros(length(idx_i), 1);
    bfz = zeros(length(idx_i), 1);
    
    % We iterate through the sparse list we just created?
    % Or, perform the operation on the masked matrices directly 
    % (preserving shape B x N).
    % Let's use the masked matrix approach (zero out invalid, sum rows).
    
    Fx_mat = zeros(size(mask));
    Fy_mat = zeros(size(mask));
    Fz_mat = zeros(size(mask));
    
    Fx_mat(valid_idx) = f_x_local;
    Fy_mat(valid_idx) = f_y_local;
    Fz_mat(valid_idx) = f_z_local;
    
    % Sum over j (columns) to get force on i
    fx_total(idx_i) = fx_total(idx_i) + sum(Fx_mat, 2);
    fy_total(idx_i) = fy_total(idx_i) + sum(Fy_mat, 2);
    fz_total(idx_i) = fz_total(idx_i) + sum(Fz_mat, 2);
    
    % --- Store Pairs (If requested) ---
    if return_pairs
        n_new = length(real_i);
        % Only store upper triangle pairs (i < j) to match original format
        real_j = col_sub; % Since we iterate 1:N
        
        % Filter strictly i < j for output list
        ut_mask = real_i < real_j;
        
        if any(ut_mask)
            u_i = real_i(ut_mask);
            u_j = real_j(ut_mask);
            u_d = [dx_val(ut_mask), dy_val(ut_mask), dz_val(ut_mask)];
            
            n_add = length(u_i);
            if pair_count + n_add > length(pairs_i)
                % Grow arrays
                pairs_i(end*2) = 0; 
                pairs_j(end*2) = 0; 
                d_mic_out(end*2,3) = 0;
            end
            
            pairs_i(pair_count+1 : pair_count+n_add) = u_i;
            pairs_j(pair_count+1 : pair_count+n_add) = u_j;
            d_mic_out(pair_count+1 : pair_count+n_add, :) = u_d;
            pair_count = pair_count + n_add;
        end
    end
end

% Trim pairs
if return_pairs
    pairs_i = pairs_i(1:pair_count);
    pairs_j = pairs_j(1:pair_count);
    d_mic_out = d_mic_out(1:pair_count, :);
else
    pairs_i = []; pairs_j = []; d_mic_out = [];
end

% Combine Forces
totalforces = [fx_total, fy_total, fz_total];

% Convert to displacements
potdisp = totalforces * force_scale;

% --- Displacement Clamping (Stability) ---
% Yes, this is needed if timestep is large or atoms overlap.
norms = vecnorm(potdisp, 2, 2);
overshoot = norms > maxstep;
if any(overshoot)
    potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
end

disppot = potdisp;

end