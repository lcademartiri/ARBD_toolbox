function potdisp = potential_displacements_2D_v13(p, S, H, H_interpolant, ghostghost, cacheSizeMB)
% potential_displacements_2D_v13
% Unified, optimized single-file implementation for 2D systems.
%
% Signature:
%   disppot = potential_displacements_2D_v13(p, S, H_mat, H_interpolant, ghostghost, cacheSizeMB)
%
% Inputs:
%   p              Nx2 positions (X, Y)
%   S              sim struct; required fields depend on S.bc:
%                    - S.bc: 1=SBC, 2=Square PBC, 3=Oblique/Hex PBC (General 2D), 4=BigBox
%                    - S.N (real particle count) for SBC/BB ghost filtering
%                    - S.br (half-box size) for PBCs
%                    - S.rc (cutoff)
%                    - For Lattice (bc=3): S.fcc.A and S.fcc.invA must be 2x2 matrices
%                    - S.esdiff, S.kbT, S.timestep, S.pot_clamp, S.stdx
%   H_mat          optional table [r, F] (used to infer pot bounds). May be []
%   H_interpolant  function handle: F = H_interpolant(r) (vectorized)
%   ghostghost     optional for SBC ghost filtering (0 or 1)
%   cacheSizeMB    optional, approx L3 cache in MB
%
% Outputs:
%   potdisp        N x 2 displacement vectors

    if nargin < 5, ghostghost = 0; end
    if nargin < 6 || isempty(cacheSizeMB)
        cacheSizeMB = 20;   % default L3 per core
    end
    
    N = size(p,1);
    L = 2*S.br;
    invL = 1/L;
    rc = S.rc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% === SBC PAIRS (original rangesearch + ghosts preserved) =
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if S.bc==1 || S.bc==4
    
        N_total = size(p,1);
        N = S.N;
        
        % rangesearch works in 2D automatically if p is Nx2
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
        
        % raw displacement (no MIC in SBC)
        d_mic = p(pairs_i,:) - p(pairs_j,:); 
        
        % --- CASE 1: NO INTERACTIONS (Zero Force) ---
        if isempty(pairs_i)
            potdisp = zeros(N, 2);
            return; 
        end
        
        rc = S.rc;
        
        dx = d_mic(:,1);
        dy = d_mic(:,2);
        % Z removed
        r  = sqrt(dx.^2 + dy.^2);
        
        pot_r_min = H(1,1);
        pot_r_max = H(end,1);
        pot_F_min = H(1,2);
        
        Fij = H_interpolant(r);
        
        Fij(r >= rc | r >= pot_r_max) = 0;
        Fij(isnan(Fij) & r < pot_r_min) = pot_F_min;
        Fij(isnan(Fij) & r >= pot_r_max) = 0;
        
        inv_r = 1./r;
        inv_r(isinf(inv_r))=0;
        
        % 2D components
        fx = Fij .* dx .* inv_r;
        fy = Fij .* dy .* inv_r;
        
        Fx = accumarray(pairs_i, fx, [N_total,1]) - accumarray(pairs_j, fx, [N_total,1]);
        Fy = accumarray(pairs_i, fy, [N_total,1]) - accumarray(pairs_j, fy, [N_total,1]);
        
        totalforces = [Fx Fy];
        
        % convert to displacements
        potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;
        
        % clamp (Note: S.pot_clamp usually scaled by sqrt(dim) -> sqrt(2) for 2D)
        maxstep = S.pot_clamp * sqrt(2) * S.stdx; 
        norms = vecnorm(potdisp,2,2);
        overshoot = norms > maxstep;
        if any(overshoot)
            potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
        end
    
    elseif S.bc==2 || S.bc==3
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BLOCK SIZE DETERMINATION %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjusted for 2D (4 doubles per particle instead of 6, etc)
        % Conservative estimate:
        % Bytes approx 33*B^2 + 48*B + 32*N (since 2D has fewer vectors)
        
        bytesPerBlock = cacheSizeMB * 1048576 / 1.2; 
        overhead_N = 32 * N; % Reduced from 48*N for 2D
        
        if overhead_N < bytesPerBlock
            available = bytesPerBlock - overhead_N;
            B = floor(sqrt(available / 33));
        else
            B = floor(sqrt(bytesPerBlock / 33));
        end
        B = max(1, B);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INITIALIZATION %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pot_r_min = H(1,1);
        pot_r_max = H(end,1);
        pot_F_min = H(1,2);
        
        px = p(:,1); 
        py = p(:,2); 
        
        Fx = zeros(N,1);
        Fy = zeros(N,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BLOCK LOOP %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii = 1:B:N
            i_end = min(ii+B-1, N);
        
            pi_x = px(ii:i_end);
            pi_y = py(ii:i_end);
        
            for jj = ii:B:N
                j_end = min(jj+B-1, N);
        
                pj_x = px(jj:j_end);
                pj_y = py(jj:j_end);
        
                %% 1. Compute component deltas for block
                if S.bc == 2
                    % -------- Square PBC MIC --------
                    cdiff = pi_x - pj_x.'; 
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = cdiff.^2;
        
                    cdiff = pi_y - pj_y.';   
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = r2+cdiff.^2;
                    
                    % No Z component
                
                elseif S.bc == 3
                    % -------- General 2D Lattice MIC --------
                    % Expects S.fcc.A and S.fcc.invA to be 2x2 matrices
                
                    dx = pi_x - pj_x.'; 
                    dy = pi_y - pj_y.'; 
                
                    % reshape into 2×M 
                    % (Stacking X and Y)
                    D = [dx(:)'; dy(:)'];
                
                    % convert to fractional -> wrap -> convert back
                    Sfrac = S.fcc.invA * D;
                    Sfrac = Sfrac - round(Sfrac);
                    Dwrap = S.fcc.A * Sfrac;
                
                    Dx = reshape(Dwrap(1,:), size(dx));
                    Dy = reshape(Dwrap(2,:), size(dx));
                    r2 = Dx.^2 + Dy.^2;
                end        
        
                % 2. Compute spherical mask and global coordinates
                mask = r2 < rc^2;
                if ii == jj
                     mask = triu(mask, 1);
                end
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
                r2 = 1./r2;
                r2(isinf(r2)) = 0;
                Fij = Fij .* r2;
                
                if S.bc==2
                    % 7. Calculate displacement components (X)
                    n_pairs = length(rows);
                    dx_val = pi_x(rows) - pj_x(cols);
                    dx_val = dx_val - L .* round(dx_val .* invL);
                    dx_val = Fij .* dx_val;
                    for k = 1:n_pairs                       
                        Fx(gI(k)) = Fx(gI(k)) + dx_val(k);                     
                        Fx(gJ(k)) = Fx(gJ(k)) - dx_val(k); 
                    end
                    
                    % 8. Calculate displacement components (Y)
                    dx_val = pi_y(rows) - pj_y(cols);
                    dx_val = dx_val - L .* round(dx_val .* invL);
                    dx_val = Fij .* dx_val;
                    for k = 1:n_pairs                      
                        Fy(gI(k)) = Fy(gI(k)) + dx_val(k);                       
                        Fy(gJ(k)) = Fy(gJ(k)) - dx_val(k); 
                    end
                    
                elseif S.bc == 3
                    n_pairs = length(rows);
                    
                    dx = pi_x(rows) - pj_x(cols);
                    dy = pi_y(rows) - pj_y(cols);
                    
                    % 2D lattice wrap for masked pairs
                    D = [dx.'; dy.'];       % 2 × n_pairs
                    Sfrac = S.fcc.invA * D;
                    Sfrac = Sfrac - round(Sfrac);
                    Dwrap = S.fcc.A * Sfrac;      % 2 × n_pairs
                    
                    dxw = Dwrap(1,:).';           
                    dyw = Dwrap(2,:).';
                    
                    dxw = Fij .* dxw;
                    dyw = Fij .* dyw;
                    
                    for k = 1:n_pairs
                        iidx = gI(k);
                        jidx = gJ(k);
                
                        Fx(iidx) = Fx(iidx) + dxw(k);
                        Fx(jidx) = Fx(jidx) - dxw(k);
                
                        Fy(iidx) = Fy(iidx) + dyw(k);
                        Fy(jidx) = Fy(jidx) - dyw(k);
                    end
                end
            end
        end
        % Finalize
        totalforces = [Fx Fy];
        potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;
        
        % 2D Clamp (sqrt(2))
        maxstep = S.pot_clamp * sqrt(2) * S.stdx;
        norms = vecnorm(potdisp, 2, 2);
        overshoot = norms > maxstep;
        if any(overshoot)
            potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
        end
    end
end