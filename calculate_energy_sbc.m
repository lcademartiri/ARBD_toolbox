function [u_mean, U_total, U_particles, Coords_spherical] = calculate_energy_sbc(pos, R_sphere, rc, epsilon, sigma)
    % INPUTS:
    % pos:       (N_total x 3) array of ALL particles (Real + Ghost)
    % R_sphere:  Radius defining the Real domain boundary
    % rc:        Cutoff radius (in same units as pos)
    % epsilon:   Energy scale (in your physical energy units)
    % sigma:     Particle diameter (in your physical length units)
    %
    % OUTPUTS:
    % u_mean:           Scalar mean energy per real particle
    % U_total:          Scalar total energy of the real system
    % U_particles:      (N_real x 1) Vector of individual particle energies
    % Coords_spherical: (N_real x 3) Matrix of [r, theta, phi]
    %                   r: Radial distance
    %                   theta: Polar angle (0 to pi, from +Z axis)
    %                   phi: Azimuthal angle (-pi to pi, in XY plane)

    % --- 1. Identify Real vs Ghost Particles ---
    % Calculate distance from origin for everyone
    r_from_origin = sqrt(sum(pos.^2, 2));
    
    is_real = r_from_origin < R_sphere;
    real_indices = find(is_real);
    N_real = length(real_indices);
    
    % Initialize Output Arrays
    U_particles = zeros(N_real, 1);
    
    % --- 2. Calculate Spherical Coordinates for Real Particles ---
    % MATLAB's cart2sph returns: [azimuth(phi), elevation, r]
    [az, el, r] = cart2sph(pos(real_indices,1), pos(real_indices,2), pos(real_indices,3));
    
    % Store in [r, theta, phi] format
    Coords_spherical = [az, el, r];

    % --- 3. Neighbor Search ---
    % Query = Real Particles, Reference = All Particles (Real + Ghost)
    [idx_list, dists] = rangesearch(pos, pos(real_indices, :), rc);
    
    % --- 4. Loop over Real Particles ---
    for k = 1:N_real
        neighbors = idx_list{k};  % Indices in 'pos'
        d = dists{k};             % Distances
        
        % Filter self-interaction (d=0)
        valid_mask = d > eps;
        neighbors = neighbors(valid_mask);
        d = d(valid_mask);
        
        if isempty(neighbors), continue; end
        
        % Calculate Pair Potentials (Lennard-Jones)
        inv_r2 = (sigma ./ d).^2;
        inv_r6 = inv_r2.^3;
        inv_r12 = inv_r6.^2;
        
        v_ij = 4 * epsilon * (inv_r12 - inv_r6);
        
        % --- ENERGY DISTRIBUTION LOGIC ---
        % Check if neighbors are Real or Ghost
        neighbor_is_real = is_real(neighbors);
        
        % Real-Real: 0.5 factor (Bond shared between two recorded particles)
        % Real-Ghost: 1.0 factor (Bond owned fully by the real particle)
        energy_contribution = sum(v_ij(neighbor_is_real) * 0.5) + ...
                              sum(v_ij(~neighbor_is_real) * 1.0);
                          
        U_particles(k) = energy_contribution;
    end
    
    % --- 5. Apply Tail Correction ---
    % Calculate the uniform background energy correction per particle
    V_sphere = (4/3) * pi * R_sphere^3;
    rho_real = N_real / V_sphere;
    rr = sigma / rc;
    
    % Tail correction per particle (Energy Density * Volume / N)
    u_tail_per_particle = (8/3) * pi * rho_real * epsilon * sigma^3 * ...
                          ( (1/3)*(rr^9) - (rr^3) );
    
    % Add this background field to every particle
    U_particles = U_particles + u_tail_per_particle;
    
    % --- 6. Final Global Stats ---
    U_total = sum(U_particles);
    u_mean = U_total / N_real;
end