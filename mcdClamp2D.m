function clamp = mcdClamp2D(nodisp, DISP, esdiff, timestep, forcevector, kbT)
    % --- 1. Determine Sampling Range Dynamically ---
    r_contact = forcevector(1,1);       % (Start of table)
    r_cutoff  = forcevector(end,1);     % ~2.5*sigma or 3*sigma
    stdx=std(DISP(:,1));
    
    % --- 2. MC array of random relative displacements (2D) ---
    % Sample uniformly in AREA within the interaction ring
    % Range: [r_contact, r_cutoff]
    
    % 2D Angle: 0 to 2pi
    mcd.theta = rand(nodisp,1) * 2 * pi;
    
    % Generalized sampling ring (Area scaling: r^2)
    u = rand(nodisp,1);
    r2_min = r_contact^2;
    r2_max = r_cutoff^2;
    
    % Inverse transform sampling for 2D area
    mcd.rho = sqrt(r2_min + (r2_max - r2_min) .* u);
    
    % Convert Polar to Cartesian (2D)
    [mcd.xi, mcd.yi] = pol2cart(mcd.theta, mcd.rho);
    
    % --- 3. Relative Brownian Displacements ---
    % DISP is expected to be N x 2. We need 'nodisp' samples.
    n_lib = size(DISP,1);
    
    % Safety check for dimensionality (Optional but recommended)
    if size(DISP, 2) ~= 2
        error('Input DISP must be N x 2 for 2D operations.');
    end

    idx1 = randi(n_lib, nodisp, 1);
    idx2 = randi(n_lib, nodisp, 1);
    mcd.DISPrel = DISP(idx1,:) - DISP(idx2,:);
    
    % --- 4. Gauss-Legendre Quadrature Setup ---
    % (This section is dimension-agnostic)
    mcd.qnodes = 12;
    mcd.beta = 0.5 ./ sqrt(1-(2*(1:mcd.qnodes-1)).^(-2));
    mcd.T = diag(mcd.beta,1) + diag(mcd.beta,-1);
    [mcd.V, mcd.D] = eig(mcd.T);
    mcd.x = diag(mcd.D);          
    mcd.w = 2 * (mcd.V(1,:)').^2; 
    mcd.s_quad = 0.5 * (mcd.x + 1);
    mcd.w_quad = 0.5 * mcd.w;
    
    % --- 5. Path Integration ---
    % Reshape for vectorized ops: [nodisp x 1 x 2] <--- CHANGED TO 2
    mcd.r0_rep = reshape([mcd.xi, mcd.yi], [nodisp,1,2]);
    mcd.rB_rep = reshape(mcd.DISPrel, [nodisp,1,2]);
    
    % Quadrature points: [1 x qnodes x 1]
    mcd.s_rep = reshape(mcd.s_quad, [1,mcd.qnodes,1]);
    
    % Paths: [nodisp x qnodes x 2] <--- CHANGED TO 2
    mcd.r_path = mcd.r0_rep + mcd.s_rep .* mcd.rB_rep; 
    mcd.rnorm_path = sqrt(sum(mcd.r_path.^2, 3)); 
    
    % --- 6. Force Lookup with Extrapolation Safety ---
    mcd.r_query = mcd.rnorm_path;
    
    % Interpolate Force Magnitude
    mcd.Fmag_path = interp1(forcevector(:,1), forcevector(:,2), mcd.r_query, 'linear', NaN);
    
    % Handle Overlap (r < r_min) -> Use Max Force
    max_repulsion = forcevector(1,2); 
    mcd.Fmag_path(mcd.r_query < r_contact) = max_repulsion;
    
    % Handle Cutoff (r > r_max) -> Zero
    mcd.Fmag_path(mcd.r_query > r_cutoff) = 0;
    mcd.Fmag_path(isnan(mcd.Fmag_path)) = 0; 
    
    % Vector Force: [nodisp x qnodes x 2] <--- CHANGED TO 2
    mcd.nhat_path = mcd.r_path ./ mcd.rnorm_path;
    mcd.Fvec_path = mcd.Fmag_path .* mcd.nhat_path; 
    
    % --- 7. Integrate Displacement ---
    % Mobility = (D1 + D2) / kT
    mobility = (esdiff + esdiff) / kbT;
    
    % Weighted sum over quadrature points
    w_rep = reshape(mcd.w_quad, [1, mcd.qnodes, 1]);
    
    % F_avg = sum(w * F) -> [nodisp x 2]
    F_avg = squeeze(sum(mcd.Fvec_path .* w_rep, 2)); 
    
    mcd.det_disp_rel = F_avg * mobility * timestep;
    
    % --- 8. Calculate Ratios ---
    % Magnitude of deterministic vs stochastic displacement
    % vecnorm handles the 2nd dimension automatically (columns)
    mag_det = vecnorm(mcd.det_disp_rel, 2, 2);
    mag_stoch = vecnorm(mcd.DISPrel, 2, 2);
    
    % Ratio
    mcd.forcedisps = mag_det ./ (mag_stoch + eps); 
    
    % --- 9. Determine Clamp ---
    % Filter only NaNs or Infs
    valid_ratios = mcd.forcedisps(isfinite(mcd.forcedisps) & mcd.forcedisps > 0);
    
    if isempty(valid_ratios)
        clamp = 1.0; 
        return;
    end
    
    % Fit Lognormal
    try
        pd = fitdist(valid_ratios, 'Lognormal');
        % 99th percentile (Safe limit)
        clamp = icdf(pd, 0.99);
        
        if clamp > stdx*10
            fprintf('Warning: mcdClamp computed extremely high value (%.1f). Timestep may be too large.\n', clamp);
            clamp = stdx*10;
        end
    catch
        % Fallback if fit fails
        clamp = quantile(valid_ratios, 0.99);
    end
end