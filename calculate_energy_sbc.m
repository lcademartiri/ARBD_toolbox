function [u_mean, U_total, U_particles_obs, Coords_spherical, real_indices, U_particles] = calculate_energy_sbc(pos, R_sphere, rc, U_interpolant)
    % INPUTS:
    % pos:           (N_total x 3) array of ALL particles (Real + Ghost)
    % R_sphere:      Radius defining the Real domain boundary
    % rc:            Cutoff radius (used for neighbor search)
    % U_interpolant: Function handle or griddedInterpolant object U(r)
    %                Must return potential energy given distance r.
    %
    % OUTPUTS:
    % u_mean:             Scalar mean energy per real particle
    % U_total:            Scalar total energy of the real system
    % U_particles:        (N_real x 1) Vector of individual particle energies
    % Coords_spherical:   (N_real x 3) Matrix of [az, el, r]
    % real_indices:       Indices of the particles considered "Real"

    % --- 1. Identify Real vs Ghost Particles ---
    r_from_origin = sqrt(sum(pos.^2, 2));
    
    is_real = r_from_origin < R_sphere;
    real_indices = find(is_real);
    N_real = length(real_indices);
	N_all = size(pos,1);
    
    % Initialize Output Arrays
    U_particles = zeros(N_all, 1);
    
    % --- 2. Calculate Spherical Coordinates for Real Particles ---
    % MATLAB's cart2sph returns: [azimuth, elevation, r]
    [az, el, r] = cart2sph(pos(real_indices,1), pos(real_indices,2), pos(real_indices,3));
    
    % Store in [az, el, r] format
    Coords_spherical = [az, el, r];

    % --- 3. Neighbor Search ---
    % Query = Real Particles, Reference = All Particles (Real + Ghost)
    [idx_list, dists] = rangesearch(pos, pos, rc);
    
    % --- 4. Loop over Particles ---
    for k = 1:N_all
        neighbors = idx_list{k};  % Indices in 'pos'
        d = dists{k};             % Distances
        
        % Filter self-interaction (d near 0)
        valid_mask = d > eps;
        neighbors = neighbors(valid_mask);
        d = d(valid_mask);
        
        if isempty(neighbors), continue; end
        
        % --- POTENTIAL CALCULATION (AGNOSTIC) ---
        % Replace hardcoded LJ with the interpolant call
        v_ij = U_interpolant(d);
        
        % --- ENERGY DISTRIBUTION LOGIC ---
        energy_contribution = sum(v_ij * 0.5);
                              
        U_particles(k) = energy_contribution;
    end
    U_particles_obs=U_particles(real_indices,:);
	
	
    % --- 6. Final Global Stats ---
    U_total = sum(U_particles_obs);
    u_mean = U_total / max(1, N_real); % Safety for empty sets
end