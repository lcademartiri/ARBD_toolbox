function disppot=potential_displacements(p,N,cutoff,forcevector,esdiff,clamp,kbT,stdx,timestep,freeghosts)

    N_total = size(p,1);            

    % Pairwise distance components (NxN)
    pot.distvecx = p(:,1) - p(:,1)';
    pot.distvecy = p(:,2) - p(:,2)';
    pot.distvecz = p(:,3) - p(:,3)';

    % Squared distances
    distnorms = pot.distvecx.^2 + pot.distvecy.^2 + pot.distvecz.^2;

    % Mask of pairs within cutoff radius
    pot.idxpot = distnorms < cutoff^2;

    % Remove self-interactions
    pot.idxpot(1:N_total+1:end) = false;

    % 🔹 EXCLUDE ghost–ghost interactions
    % if freeghosts==0
        pot.idxpot(N+1:end, N+1:end) = false;
    % end

    % Get indices of interacting pairs (within cutoff)
    [ix, iy] = find(pot.idxpot);
    mask = ix < iy;
    ix = ix(mask);
    iy = iy(mask);

    % --- CALCULATE ACTUAL DISTANCES ONLY FOR INTERACTING PAIRS ---
    linear_indices_sq = sub2ind(size(distnorms), ix, iy);
    r = sqrt(distnorms(linear_indices_sq));

    % --- GET DISPLACEMENT COMPONENTS FOR INTERACTING PAIRS ---
    dx = pot.distvecx(linear_indices_sq);
    dy = pot.distvecy(linear_indices_sq);
    dz = pot.distvecz(linear_indices_sq);

    % --- SAFE INTERPOLATION of force magnitudes ---
    pot.r_min = forcevector(1,1);
    pot.r_max = forcevector(end,1);
    pot.F_min = forcevector(1,2);  % repulsive saturation value at smallest distance

    % Clamp query distances into [r_min, r_max]
    r_clamped = min(max(r, pot.r_min), pot.r_max);
    Fij = interp1(forcevector(:,1), forcevector(:,2), r_clamped, 'pchip');

    % Explicitly zero forces beyond r_max (or cutoff)
    Fij(r >= cutoff | r >= pot.r_max) = 0;
    Fij(isnan(Fij) & (r < pot.r_min)) = pot.F_min;
    Fij(isnan(Fij) & (r >= pot.r_max)) = 0;

    % Compute vector components of pairwise forces (Fij * unit vector)
    inv_r = 1 ./ r;
    fx = Fij .* dx .* inv_r;
    fy = Fij .* dy .* inv_r;
    fz = Fij .* dz .* inv_r;

    % --- ACCUMULATE FORCES ON EACH PARTICLE ---
    % (Newton’s 3rd law: i += Fij, j -= Fij)
    Fx = accumarray(ix, fx, [N_total,1], @sum, 0) - accumarray(iy, fx, [N_total,1], @sum, 0);
    Fy = accumarray(ix, fy, [N_total,1], @sum, 0) - accumarray(iy, fy, [N_total,1], @sum, 0);
    Fz = accumarray(ix, fz, [N_total,1], @sum, 0) - accumarray(iy, fz, [N_total,1], @sum, 0);

    pot.totalforces = [Fx, Fy, Fz];

    % --- DYNAMICS: convert forces to displacements (overdamped limit) ---
    pot.diffs = repmat(esdiff, N_total, 1);
    pot.potdisp = (pot.totalforces .* pot.diffs) / kbT * timestep;

    % Clamp deterministic displacements to avoid overshoot
    pot.idxovershoot = vecnorm(pot.potdisp,2,2) > (clamp * sqrt(3) * stdx);
    if any(pot.idxovershoot)
        pot.potdisp(pot.idxovershoot,:) = ...
            (pot.potdisp(pot.idxovershoot,:) ./ vecnorm(pot.potdisp(pot.idxovershoot,:),2,2)) .* ...
            (clamp * sqrt(3) * stdx);
    end

    % Output only the displacements of the real particles
    disppot = pot.potdisp(1:N, :);    
end