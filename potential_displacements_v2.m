function disppot=potential_displacements_v2(p,N,cutoff, forcevector, interpolant,esdiff,clamp,kbT,stdx,timestep,ghostghost)

    % Number of real+ghost particles
    N_total = size(p,1);

    % ---------------------------------------------------------------------
    % 1. EFFICIENT NEIGHBOR SEARCH
    % ---------------------------------------------------------------------
    % This is the correct first step. It's extremely fast.
    [idx_list, dist_list] = rangesearch(p, p, cutoff);

    % Get the number of neighbors for each particle.
    num_neighbors = cellfun(@numel, idx_list);

    % Create the 'iy' vector by stacking all neighbor indices into one long column.
    iy = horzcat(idx_list{:})';

    
    % Create the 'ix' vector by repeating each particle's index by its number of neighbors.
    ix = repelem((1:N_total)', num_neighbors);
    
    % Create the 'r' (distance) vector by stacking all distances.
    % This is perfectly aligned with ix and iy. No 'ismember' needed!
    r = horzcat(dist_list{:})';

    % ---------------------------------------------------------------------
    % 3. VECTORIZED FILTERING OF PAIRS
    % ---------------------------------------------------------------------
    
    % Create a single logical mask to remove all unwanted pairs at once.
    % Condition 1: Keep only unique pairs where ix < iy. This automatically
    % removes self-pairs (ix == iy) and duplicate pairs (iy, ix).
    mask = ix < iy;
    
    % Condition 2 (Optional): Exclude ghost-ghost interactions.
    % A particle is a "ghost" if its index is > N.
    % We want to remove pairs where BOTH ix and iy are ghosts.
    if ghostghost == 0
        mask = mask & ~(ix > N & iy > N);
    end

    % Apply the mask to all vectors in a single operation.
    ix = ix(mask);
    iy = iy(mask);
    r  = r(mask);
    
    % ---------------------------------------------------------------------
    % 4. CALCULATE DISPLACEMENT VECTORS (dx, dy, dz)
    % ---------------------------------------------------------------------
    % Now that we have the final list of pairs, calculate the displacement
    % vectors for ONLY these pairs. This avoids creating giant NxN matrices.
    % Note on convention: rij = p(j) - p(i)
    rij = p(ix,:) - p(iy,:);
    dx = rij(:,1);
    dy = rij(:,2);
    dz = rij(:,3);
    
    % --- The rest of your code is already well-vectorized and correct ---
    
    % INTERPOLATION OF FORCE MAGNITUDES
    pot_r_min = forcevector(1,1);
    pot_r_max = forcevector(end,1);
    pot_F_min = forcevector(1,2);
    r_clamped = min(max(r, pot_r_min), pot_r_max);
    Fij = interpolant(r_clamped);
    Fij(r >= cutoff | r >= pot_r_max) = 0;
    Fij(isnan(Fij) & (r < pot_r_min)) = pot_F_min;
    Fij(isnan(Fij) & (r >= pot_r_max)) = 0;

    % FORCE COMPONENTS
    inv_r = 1 ./ r;
    fx = Fij .* dx .* inv_r;
    fy = Fij .* dy .* inv_r;
    fz = Fij .* dz .* inv_r;

    % ACCUMULATE FORCES
    Fx = accumarray(ix, fx, [N_total,1]) - accumarray(iy, fx, [N_total,1]);
    Fy = accumarray(ix, fy, [N_total,1]) - accumarray(iy, fy, [N_total,1]);
    Fz = accumarray(ix, fz, [N_total,1]) - accumarray(iy, fz, [N_total,1]);
    totalforces = [Fx, Fy, Fz];

    % CONVERT FORCES TO DISPLACEMENTS
    potdisp = (totalforces .* esdiff) / kbT * timestep;

    % Clamp overshoot
    overshoot = vecnorm(potdisp,2,2) > (clamp * sqrt(3) * stdx);
    if any(overshoot)
        potdisp(overshoot,:) = (potdisp(overshoot,:) ./ vecnorm(potdisp(overshoot,:),2,2)) .* (clamp * sqrt(3) * stdx);
    end

    % OUTPUT
    disppot = potdisp(1:N, :);
end
