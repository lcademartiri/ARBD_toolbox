function clamp = mcdClamp(nodisp, DISP, esdiff, timestep, forcevector, kbT)
    % --- 1. Determine Sampling Range Dynamically ---
    r_contact = forcevector(1,1);       % (Start of table)
    r_cutoff  = forcevector(end,1);     % ~2.5*sigma or 3*sigma
    shell_width = r_cutoff - r_contact;
    
    % --- 2. MC array of random relative displacements ---
    % Sample uniformly in volume within the interaction shell
    % Range: [r_contact, r_contact + shell_width]
    mcd.az = rand(nodisp,1) * 2 * pi;
    mcd.el = acos(2 * rand(nodisp,1) - 1) - pi/2;
    
    % Generalized sampling shell
    u = rand(nodisp,1);
    r3_min = r_contact^3;
    r3_max = r_cutoff^3;
    mcd.rho = (r3_min + (r3_max - r3_min) .* u).^(1/3);
    
    [mcd.xi, mcd.yi, mcd.zi] = sph2cart(mcd.az, mcd.el, mcd.rho);
    
    % --- 3. Relative Brownian Displacements ---
    % DISP is usually N x 3. We need 'nodisp' samples.
    % If DISP is smaller than nodisp, we resample.
    n_lib = size(DISP,1);
    idx1 = randi(n_lib, nodisp, 1);
    idx2 = randi(n_lib, nodisp, 1);
    mcd.DISPrel = DISP(idx1,:) - DISP(idx2,:);
    
    % --- 4. Gauss-Legendre Quadrature Setup ---
    mcd.qnodes = 12;
    mcd.beta = 0.5 ./ sqrt(1-(2*(1:mcd.qnodes-1)).^(-2));
    mcd.T = diag(mcd.beta,1) + diag(mcd.beta,-1);
    [mcd.V, mcd.D] = eig(mcd.T);
    mcd.x = diag(mcd.D);          
    mcd.w = 2 * (mcd.V(1,:)').^2; 
    mcd.s_quad = 0.5 * (mcd.x + 1);
    mcd.w_quad = 0.5 * mcd.w;
    
    % --- 5. Path Integration ---
    % Reshape for vectorized ops: [nodisp x 1 x 3]
    mcd.r0_rep = reshape([mcd.xi, mcd.yi, mcd.zi], [nodisp,1,3]);
    mcd.rB_rep = reshape(mcd.DISPrel, [nodisp,1,3]);
    
    % Quadrature points: [1 x qnodes x 1]
    mcd.s_rep = reshape(mcd.s_quad, [1,mcd.qnodes,1]);
    
    % Paths: [nodisp x qnodes x 3]
    mcd.r_path = mcd.r0_rep + mcd.s_rep .* mcd.rB_rep; 
    mcd.rnorm_path = sqrt(sum(mcd.r_path.^2, 3)); 
    
    % --- 6. Force Lookup with Extrapolation Safety ---
    % We query the force table.
    mcd.r_query = mcd.rnorm_path;
    
    % Interpolate Force Magnitude
    % 'linear' is safer than 'pchip' for noisy tables, but pchip is fine if smooth
    mcd.Fmag_path = interp1(forcevector(:,1), forcevector(:,2), mcd.r_query, 'linear', NaN);
    
    % Handle Overlap (r < r_min) -> Use Max Force (Pot_F_min)
    % This assumes forcevector(:,2) is positive for repulsion? 
    % Usually LJ force is positive (repulsive) at short range.
    % forcevector(1,2) should be the highest repulsive force.
    max_repulsion = forcevector(1,2); 
    mcd.Fmag_path(mcd.r_query < r_contact) = max_repulsion;
    
    % Handle Cutoff (r > r_max) -> Zero
    mcd.Fmag_path(mcd.r_query > r_cutoff) = 0;
    mcd.Fmag_path(isnan(mcd.Fmag_path)) = 0; % Safety
    
    % Vector Force: [nodisp x qnodes x 3]
    mcd.nhat_path = mcd.r_path ./ mcd.rnorm_path;
    mcd.Fvec_path = mcd.Fmag_path .* mcd.nhat_path; 
    
    % --- 7. Integrate Displacement ---
    % Displacement = Mobility * Integral(Force dt)
    % Mobility = (D1 + D2) / kT
    mobility = (esdiff + esdiff) / kbT;
    
    % Weighted sum over quadrature points
    % weights: [1 x qnodes x 1]
    w_rep = reshape(mcd.w_quad, [1, mcd.qnodes, 1]);
    
    % Integral F(s) ds approx sum(w * F)
    % Note: The path parameter s goes 0->1. The Jacobian is |rB|.
    % But you are integrating over TIME (dt). 
    % The standard approximation is F_avg * dt.
    % F_avg = sum(w * F).
    F_avg = squeeze(sum(mcd.Fvec_path .* w_rep, 2)); % [nodisp x 3]
    
    mcd.det_disp_rel = F_avg * mobility * timestep;
    
    % --- 8. Calculate Ratios ---
    % Magnitude of deterministic vs stochastic displacement
    mag_det = vecnorm(mcd.det_disp_rel, 2, 2);
    mag_stoch = vecnorm(mcd.DISPrel, 2, 2);
    
    % Ratio
    mcd.forcedisps = mag_det ./ (mag_stoch + eps); % Avoid div/0
    
    % --- 9. Determine Clamp ---
    % Do NOT filter > 1. We want to know the tail statistics.
    % Filter only NaNs or Infs
    valid_ratios = mcd.forcedisps(isfinite(mcd.forcedisps) & mcd.forcedisps > 0);
    
    if isempty(valid_ratios)
        % Fallback if no interactions found (e.g. dilute)
        clamp = 1.0; 
        return;
    end
    
    % Fit Lognormal
    try
        pd = fitdist(valid_ratios, 'Lognormal');
        % 99th percentile (Safe limit)
        % If the ratio > clamp, we clip it.
        % A clamp of 5.0 means we allow potential moves 5x larger than thermal moves.
        clamp = icdf(pd, 0.99);
        
        % Hard safety limit: Don't let clamp be insane (e.g. > 100)
        % If physics says we need 1000x thermal step, timestep is too large.
        if clamp > 100
            fprintf('Warning: mcdClamp computed extremely high value (%.1f). Timestep may be too large.\n', clamp);
            clamp = 100;
        end
    catch
        % Fallback if fit fails
        clamp = quantile(valid_ratios, 0.99);
    end
end