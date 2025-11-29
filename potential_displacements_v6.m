function [disppot, pairs_i, pairs_j, d_mic_out] = potential_displacements_v6(p, S, H, H_interpolant)
% V6: Block-processed + Early Exit Cascading Filter.
% Only calculates Y and Z components for pairs that pass the X check.

N = size(p,1);
L = 2*S.br;
invL = 1/L; 
rc = S.rc;
rc2 = rc^2;

% Potential constants
pot_r_min = H(1,1);
pot_r_max = H(end,1);
pot_F_min = H(1,2);

force_scale = (S.esdiff / S.kbT) * S.timestep;
maxstep = S.pot_clamp * sqrt(3) * S.stdx;

% Outputs
fx_total = zeros(N,1);
fy_total = zeros(N,1);
fz_total = zeros(N,1);

return_pairs = (nargout > 1);
if return_pairs
    % Dynamic allocation buffer
    est_pairs = N * 50; 
    pairs_i = zeros(est_pairs, 1);
    pairs_j = zeros(est_pairs, 1);
    d_mic_out = zeros(est_pairs, 3);
    pair_count = 0;
end

% -------------------------------------------------------------------------
% BLOCK PROCESSING
% -------------------------------------------------------------------------
% Target L3 Cache usage (e.g., 4 MB to be safe and conservative)
% 4 MB = 4 * 1024^2 bytes
TARGET_CACHE_BYTES = 4 * 1024^2; 

% Size of one double-precision number
BYTES_PER_DOUBLE = 8;

% Calculate optimal block size
% We want: BlockSize * N * 8 <= Target
ideal_block = floor(TARGET_CACHE_BYTES / (N * BYTES_PER_DOUBLE));

% CLAMPING
% Upper bound: No need to be larger than N
% Lower bound: Don't go below ~64 or MATLAB loop overhead dominates
BLOCK_SIZE = max(64, min(N, ideal_block));

% Extract columns for faster indexing
px = p(:,1);
py = p(:,2);
pz = p(:,3);

for i_start = 1:BLOCK_SIZE:N
    i_end = min(i_start + BLOCK_SIZE - 1, N);
    
    % Indices for this block (Global)
    idx_i_global = i_start:i_end; 
    
    % Coordinates for the block 'i'
    bx = px(idx_i_global);
    by = py(idx_i_global);
    bz = pz(idx_i_global);
    
    % =====================================================================
    % 1. X-COMPONENT CHECK (The Primary Filter)
    % =====================================================================
    % We compute X for the whole block because it's vectorized and contiguous.
    dx_block = bx - px.'; 
    
    % Fast MIC on X
    dx_block = dx_block - L .* round(dx_block .* invL);
    
    % Find valid pairs immediately.
    % This converts the dense BxN matrix into a sparse list of indices.
    % We stop processing the matrix structure here.
    [rows_sub, cols_sub] = find(abs(dx_block) <= rc);
    
    if isempty(rows_sub)
        continue;
    end
    
    % Extract the surviving dx values
    % (Indexing into dx_block uses linear indices)
    % Note: Using linear indexing is faster than sub2ind
    lin_idx = rows_sub + (cols_sub-1)*size(dx_block,1);
    dx_v = dx_block(lin_idx);
    
    % =====================================================================
    % 2. Y-COMPONENT CHECK (Sparse)
    % =====================================================================
    % Only calculate Y for the pairs that passed X.
    % We use the indices found above to pull specific values.
    
    % rows_sub maps to local block indices (1..BlockSize)
    % cols_sub maps to global indices (1..N)
    dy_v = by(rows_sub) - py(cols_sub);
    
    % Fast MIC on Y
    dy_v = dy_v - L .* round(dy_v .* invL);
    
    % Filter again
    mask_y = abs(dy_v) <= rc;
    
    % Prune our lists based on Y failure
    if ~any(mask_y), continue; end
    
    rows_sub = rows_sub(mask_y);
    cols_sub = cols_sub(mask_y);
    dx_v     = dx_v(mask_y);
    dy_v     = dy_v(mask_y);
    
    % =====================================================================
    % 3. Z-COMPONENT CHECK (Sparse)
    % =====================================================================
    dz_v = bz(rows_sub) - pz(cols_sub);
    dz_v = dz_v - L .* round(dz_v .* invL);
    
    mask_z = abs(dz_v) <= rc;
    
    if ~any(mask_z), continue; end
    
    % Prune lists based on Z failure
    rows_sub = rows_sub(mask_z);
    cols_sub = cols_sub(mask_z);
    dx_v     = dx_v(mask_z);
    dy_v     = dy_v(mask_z);
    dz_v     = dz_v(mask_z);

    % =====================================================================
    % 4. DISTANCE SQUARED
    % =====================================================================
    r2 = dx_v.^2 + dy_v.^2 + dz_v.^2;
    
    % Spherical cutoff check
    mask_r = (r2 <= rc2) & (r2 > 0); % r2>0 removes self-interaction (i==i)
    
    if ~any(mask_r), continue; end
    
    % Final pruning
    rows_sub = rows_sub(mask_r);
    cols_sub = cols_sub(mask_r);
    dx_v     = dx_v(mask_r);
    dy_v     = dy_v(mask_r);
    dz_v     = dz_v(mask_r);
    r2       = r2(mask_r);
    
    % =====================================================================
    % 5. FORCE CALCULATION
    % =====================================================================
    r = sqrt(r2);
    
    % Clamp r (Mathematical safety)
    r_clamped = r;
    r_clamped(r < pot_r_min) = pot_r_min;
    r_clamped(r > pot_r_max) = pot_r_max;
    
    Fmag = H_interpolant(r_clamped);
    
    % Handle Boundaries / NaN
    Fmag(r >= pot_r_max) = 0;
    
    bad_low = isnan(Fmag) & (r < pot_r_min);
    if any(bad_low)
        Fmag(bad_low) = pot_F_min;
    end
    
    % Vector Force: F_vec = (F_mag / r) * r_vec
    F_over_r = Fmag ./ r;
    
    fx_v = F_over_r .* dx_v;
    fy_v = F_over_r .* dy_v;
    fz_v = F_over_r .* dz_v;
    
    % =====================================================================
    % 6. ACCUMULATION
    % =====================================================================
    % rows_sub are local indices (1..BlockSize). 
    % We accumulate results into the 'i' block.
    
    % This sums up all forces acting on each particle in the block
    fx_total(idx_i_global) = fx_total(idx_i_global) + accumarray(rows_sub, fx_v, [length(idx_i_global), 1]);
    fy_total(idx_i_global) = fy_total(idx_i_global) + accumarray(rows_sub, fy_v, [length(idx_i_global), 1]);
    fz_total(idx_i_global) = fz_total(idx_i_global) + accumarray(rows_sub, fz_v, [length(idx_i_global), 1]);
    
    % =====================================================================
    % 7. PAIR STORAGE (Optional)
    % =====================================================================
    if return_pairs
        % Map local block index back to global i
        real_i = idx_i_global(rows_sub)';
        real_j = cols_sub;
        
        % Filter for Upper Triangle (i < j) for output consistency
        ut_mask = real_i < real_j;
        
        if any(ut_mask)
            n_add = sum(ut_mask);
            if pair_count + n_add > length(pairs_i)
                % Grow arrays
                new_size = length(pairs_i)*2 + n_add;
                pairs_i(new_size) = 0; pairs_j(new_size) = 0; d_mic_out(new_size,3) = 0;
            end
            
            rng = pair_count+1 : pair_count+n_add;
            pairs_i(rng) = real_i(ut_mask);
            pairs_j(rng) = real_j(ut_mask);
            d_mic_out(rng,:) = [dx_v(ut_mask), dy_v(ut_mask), dz_v(ut_mask)];
            pair_count = pair_count + n_add;
        end
    end
end

% Cleanup output arrays
if return_pairs
    pairs_i = pairs_i(1:pair_count);
    pairs_j = pairs_j(1:pair_count);
    d_mic_out = d_mic_out(1:pair_count, :);
else
    pairs_i=[]; pairs_j=[]; d_mic_out=[];
end

% Final Displacement Calculation
potdisp = [fx_total, fy_total, fz_total] * force_scale;

% Displacement Clamp (Simulation Stability)
norms = sqrt(sum(potdisp.^2, 2));
overshoot = norms > maxstep;
if any(overshoot)
    potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
end

disppot = potdisp;
end