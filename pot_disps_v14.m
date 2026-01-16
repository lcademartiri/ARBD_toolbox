function [potdisp, U_particles] = pot_disps_v14(p, S, H, H_interpolant, U_interpolant, ghostghost, cacheSizeMB)
% potential_displacements_v14
% Unified, optimized single-file implementation for SBC, Cubic PBC, FCC PBC, BB.
% Now includes OPTIONAL Potential Energy Output.
%
% Signature:
%   [potdisp, U_particles] = potential_displacements_v14(p, S, H_mat, H_int, U_int, ghostghost, cacheSizeMB)
%
% New Input:
%   U_interpolant   function handle: U = U_interpolant(r) (vectorized energy)
%                   If requested (nargout > 1), this is required.

    if nargin < 6, ghostghost = 0; end
    if nargin < 7 || isempty(cacheSizeMB)
        cacheSizeMB = 20;   % default L3 per core
    end
    
    % Check if energy output is requested to avoid overhead when not needed
    calc_energy = (nargout > 1);
    
    N = size(p,1);
    L = 2*S.br;
    invL=1/L;
    rc = S.rc;
    
    % Initialize Energy Array if requested
    if calc_energy
        U_particles = zeros(N, 1);
        if nargin < 5 || isempty(U_interpolant)
            error('U_interpolant required when requesting potential energy output.');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% === SBC PAIRS (original rangesearch + ghosts preserved) =
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if S.bc==1 || S.bc==4
    
        N_total = size(p,1);
        N = S.N;
        
        [idx_list, ~] = rangesearch(p, p, rc);
        
        num_neighbors = cellfun(@numel, idx_list);
        pairs_j = horzcat(idx_list{:})';
        pairs_i = repelem((1:N_total)', num_neighbors);
        
        mask = pairs_i < pairs_j;
        if ghostghost == 0
            mask = mask & ~(pairs_i > N & pairs_j > N);
        end
        
        pairs_i = pairs_i(mask);
        pairs_j = pairs_j(mask);
        
        % --- CASE 1: NO INTERACTIONS (Zero Force) ---
        if isempty(pairs_i)
            potdisp = zeros(N_total, 3);
            return; 
        end
        
        % raw displacement (no MIC in SBC)
        d_mic = p(pairs_i,:) - p(pairs_j,:);
        
        dx = d_mic(:,1);
        dy = d_mic(:,2);
        dz = d_mic(:,3);
        r  = sqrt(dx.^2 + dy.^2 + dz.^2);
        
        pot_r_min = H(1,1);
        pot_r_max = H(end,1);
        pot_F_min = H(1,2);
        
        % --- Force Calculation ---
        Fij = H_interpolant(r);
        
        Fij(r >= rc | r >= pot_r_max) = 0;
        Fij(isnan(Fij) & r < pot_r_min) = pot_F_min;
        Fij(isnan(Fij) & r >= pot_r_max) = 0;
        
        inv_r = 1./r;
        inv_r(isinf(inv_r))=0;
        fx = Fij .* dx .* inv_r;
        fy = Fij .* dy .* inv_r;
        fz = Fij .* dz .* inv_r;
        
        Fx = accumarray(pairs_i, fx, [N_total,1]) - accumarray(pairs_j, fx, [N_total,1]);
        Fy = accumarray(pairs_i, fy, [N_total,1]) - accumarray(pairs_j, fy, [N_total,1]);
        Fz = accumarray(pairs_i, fz, [N_total,1]) - accumarray(pairs_j, fz, [N_total,1]);
        
        totalforces = [Fx Fy Fz];
        
        % --- Optional Energy Calculation ---
        if calc_energy
            Uij = U_interpolant(r);
            % Apply same masking logic as Forces
            Uij(r >= rc | r >= pot_r_max) = 0;
            % For potentials, we often clamp very close overlaps to a max energy or linearize
            % Assuming U_interpolant handles r < pot_r_min via extrapolation or returning high val
            
            % Distribute 0.5 * Uij to each particle in the pair
            u_half = 0.5 * Uij;
            U_particles = accumarray(pairs_i, u_half, [N_total,1]) + ...
                          accumarray(pairs_j, u_half, [N_total,1]);
        end

        % convert to displacements
        potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;
        
        % clamp
        maxstep = S.pot_clamp * sqrt(3) * S.stdx;
        norms = vecnorm(potdisp,2,2);
        overshoot = norms > maxstep;
        if any(overshoot)
            potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
        end
    
    elseif S.bc==2 || S.bc==3
    
        % ... (Block size determination code is identical) ...
        bytesPerBlock = cacheSizeMB * 1048576 / 1.2; 
        overhead_N = 48 * N;
        if overhead_N < bytesPerBlock
            available = bytesPerBlock - overhead_N;
            B = floor(sqrt(available / 33));
        else
            B = floor(sqrt(bytesPerBlock / 33));
        end
        B = max(1, B);
        
        pot_r_min = H(1,1);
        pot_r_max = H(end,1);
        pot_F_min = H(1,2);
        
        px = p(:,1); 
        py = p(:,2); 
        pz = p(:,3);
        
        Fx = zeros(N,1);
        Fy = zeros(N,1);
        Fz = zeros(N,1);
        
        % BLOCK LOOP
        for ii = 1:B:N
            i_end = min(ii+B-1, N);
            pi_x = px(ii:i_end); pi_y = py(ii:i_end); pi_z = pz(ii:i_end);
            
            for jj = ii:B:N
                j_end = min(jj+B-1, N);
                pj_x = px(jj:j_end); pj_y = py(jj:j_end); pj_z = pz(jj:j_end);
                
                %% 1. Compute component deltas for block
                if S.bc == 2
                    cdiff = pi_x - pj_x.'; 
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = cdiff.^2;
                    
                    cdiff = pi_y - pj_y.'; 
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = r2+cdiff.^2;
            
                    cdiff = pi_z - pj_z.'; 
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = r2+cdiff.^2;
                    
                elseif S.bc == 3
                    dx = pi_x - pj_x.'; dy = pi_y - pj_y.'; dz = pi_z - pj_z.'; 
                    M = numel(dx);
                    D = [dx(:)'; dy(:)'; dz(:)'];
                    Sfrac = S.fcc.invA * D;
                    Sfrac = Sfrac - round(Sfrac);
                    Dwrap = S.fcc.A * Sfrac;
                    Dx = reshape(Dwrap(1,:), size(dx));
                    Dy = reshape(Dwrap(2,:), size(dx));
                    Dz = reshape(Dwrap(3,:), size(dx));
                    r2 = Dx.^2 + Dy.^2 + Dz.^2;
                end        
        
                % 2. Compute spherical mask and global coordinates
                mask=r2<rc^2;
                if ii == jj, mask = triu(mask, 1); end
                if ~any(mask(:)), continue; end
                
                [rows, cols] = find(mask);
                gI = (ii - 1) + rows;
                gJ = (jj - 1) + cols;
        
                % 4. Vector of masked distances
                r2 = sqrt(r2(mask));
                
                % 5. Interpolate forces
                Fij = H_interpolant(r2);
                Fij(r2 >= rc | r2 >= pot_r_max) = 0;
                if any(isnan(Fij))
                    bad = isnan(Fij);
                    Fij(bad & r2 < pot_r_min) = pot_F_min;
                    Fij(bad & r2 >= pot_r_max) = 0;
                end
                
                % 6. Projection factor: F / r
                r_inv_vec = 1./r2;
                r_inv_vec(isinf(r_inv_vec)) = 0;
                Fij = Fij .* r_inv_vec;
                
                % --- Optional Energy Interpolation ---
                if calc_energy
                    Uij = U_interpolant(r2);
                    Uij(r2 >= rc | r2 >= pot_r_max) = 0;
                    u_half = 0.5 * Uij;
                end
                
                % 7. Accumulate Forces & Energies
                if S.bc==2
                    % X
                    dx_val = pi_x(rows) - pj_x(cols);
                    dx_val = dx_val - L .* round(dx_val .* invL);
                    dx_val = Fij .* dx_val;
                    % Y
                    dy_val = pi_y(rows) - pj_y(cols);
                    dy_val = dy_val - L .* round(dy_val .* invL);
                    dy_val = Fij .* dy_val;
                    % Z
                    dz_val = pi_z(rows) - pj_z(cols);
                    dz_val = dz_val - L .* round(dz_val .* invL);
                    dz_val = Fij .* dz_val;
                    
                    n_pairs = length(rows);
                    for k = 1:n_pairs                        
                        iidx = gI(k); jidx = gJ(k);
                        
                        % Force Update
                        Fx(iidx) = Fx(iidx) + dx_val(k);
                        Fx(jidx) = Fx(jidx) - dx_val(k);
                        
                        Fy(iidx) = Fy(iidx) + dy_val(k); 
                        Fy(jidx) = Fy(jidx) - dy_val(k);
                        
                        Fz(iidx) = Fz(iidx) + dz_val(k); 
                        Fz(jidx) = Fz(jidx) - dz_val(k);

                        % Energy Update
                        if calc_energy
                            U_particles(iidx) = U_particles(iidx) + u_half(k);
                            U_particles(jidx) = U_particles(jidx) + u_half(k);
                        end
                    end
                    
                elseif S.bc == 3
                    % Reconstruct wrapped vectors
                    dx = pi_x(rows) - pj_x(cols);
                    dy = pi_y(rows) - pj_y(cols);
                    dz = pi_z(rows) - pj_z(cols);
                    D = [dx.'; dy.'; dz.'];
                    Sfrac = S.fcc.invA * D;
                    Sfrac = Sfrac - round(Sfrac);
                    Dwrap = S.fcc.A * Sfrac;
                    
                    dxw = Fij .* Dwrap(1,:).';
                    dyw = Fij .* Dwrap(2,:).';
                    dzw = Fij .* Dwrap(3,:).';
                    
                    n_pairs = length(rows);
                    for k = 1:n_pairs
                        iidx = gI(k); jidx = gJ(k);
                        
                        % Force Update
                        Fx(iidx) = Fx(iidx) + dxw(k);
                        Fx(jidx) = Fx(jidx) - dxw(k);
                        
                        Fy(iidx) = Fy(iidx) + dyw(k);
                        Fy(jidx) = Fy(jidx) - dyw(k);
                        
                        Fz(iidx) = Fz(iidx) + dzw(k);
                        Fz(jidx) = Fz(jidx) - dzw(k);

                        % Energy Update
                        if calc_energy
                            U_particles(iidx) = U_particles(iidx) + u_half(k);
                            U_particles(jidx) = U_particles(jidx) + u_half(k);
                        end
                    end
                end
            end
        end
        % Finalize Forces
        totalforces = [Fx Fy Fz];
        potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;
        
        % Clamp
        maxstep = S.pot_clamp * sqrt(3) * S.stdx;
        norms = vecnorm(potdisp, 2, 2);
        overshoot = norms > maxstep;
        if any(overshoot)
            potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
        end
    end
end