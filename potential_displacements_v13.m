function [potdisp,NNS] = potential_displacements_v13(p, S, H, H_interpolant, ghostghost, cacheSizeMB, locndens)
% potential_displacements_v13
% Unified, optimized single-file implementation for SBC, Cubic PBC, FCC PBC, BB.
% Design goals:
%  - SBC: use rangesearch + fully vectorized force accumulation (same as v2)
%  - PBC (cubic & FCC): blockified difference computation + global accumarray accumulation
%  - BB (big box): treated like SBC (real-space), vectorized
%  - Minimal per-pair interpreted loops; accumulate forces via accumarray once
%  - Robust masked r^2 computation (no full-matrix sqrt)
%
% Signature:
%   disppot = potential_displacements_v13(p, S, H_mat, H_interpolant, ghostghost, cacheSizeMB)
%
% Inputs:
%   p              Nx3 positions
%   S              sim struct; required fields depend on S.bc:
%                    - S.bc: 1=SBC, 2=Cubic PBC, 3=FCC PBC, 4=BigBox
%                    - S.N (real particle count) for SBC/BB ghost filtering
%                    - S.br (half-box size) for PBCs
%                    - S.rc (cutoff)
%                    - For FCC: S.fcc.A and S.fcc.invA
%                    - Optional block-sizing fields: S.cellsize, S.ncell (not required)
%                    - S.esdiff, S.kbT, S.timestep, S.pot_clamp, S.stdx for displacement conversion
%   H_mat          optional table [r, F] (used to infer pot bounds). May be []
%   H_interpolant  function handle: F = H_interpolant(r) (vectorized)
%   ghostghost     optional for SBC ghost filtering (0 or 1)
%   cacheSizeMB    optional, approx L3 cache in MB; used to pick block size (default 20)
%
% Outputs:
%   potdisp        N x 3 displacement vectors
%
% Notes:
%   - This file is a single self-contained implementation. It intentionally keeps
%     SBC/BB on a fully vectorized path and uses blockification only for PBCs.
%   - The function defers final force accumulation to three global accumarray calls
%     for optimal N^2 scaling while keeping memory bounded with blocks.

    if ~exist('ghostghost','var'), ghostghost = 0; end
    if ~exist('cacheSizeMB','var'), cacheSizeMB = 20; end
	if ~exist('locndens','var'), locndens = false; end
    
    N = size(p,1);
    L = 2*S.br;
    invL=1/L;
    rc = S.rc;
    
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
        
        % raw displacement (no MIC in SBC)
        d_mic = p(pairs_i,:) - p(pairs_j,:);% SBC uses real-space displacement 
        
        % --- CASE 1: NO INTERACTIONS (Zero Force) ---
        if isempty(pairs_i)
            % If pairs are empty, nobody is pushing anyone.
            % Return a zero displacement matrix of the correct size.           
            potdisp = zeros(N, 3);
            return; 
        end
        
        % N_total must cover all particles, even those with no collisions
        % N_total = max([max(pairs_i), max(pairs_j), N]);
        rc = S.rc;
        
        dx = d_mic(:,1);
        dy = d_mic(:,2);
        dz = d_mic(:,3);
        r  = sqrt(dx.^2 + dy.^2 + dz.^2);
		
		if locndens
			idxlocndens=r<S.pot_sigma*1.5;
			colliders=[pairs_i(idxlocndens);pairs_j(idxlocndens)];
			colliders(colliders>S.N)=0;
			nns1=histcounts(colliders,(-1:S.N)'+0.5)';
			idxlocndens=r<S.pot_sigma*2.5;
			colliders=[pairs_i(idxlocndens);pairs_j(idxlocndens)];
			colliders(colliders>S.N)=0;
			nns2=histcounts(colliders,(-1:S.N)'+0.5)';
			NNS=uint8([nns1(2:end,1),nns2(2:end,1)]);
		end
        
        pot_r_min = H(1,1);
        pot_r_max = H(end,1);
        pot_F_min = H(1,2);
        
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
        
        % convert to displacements
        potdisp = totalforces * (S.alpha*S.esdiff / S.kbT) * S.timestep;
        
        % clamp
        maxstep = S.pot_clamp * sqrt(3) * S.stdx;
        norms = vecnorm(potdisp,2,2);
        overshoot = norms > maxstep;
        if any(overshoot)
            potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
        end
    
    elseif S.bc==2 || S.bc==3
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BLOCK SIZE DETERMINATION BASED ON CPU CACHE %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Formula provided:
        % Doubles  = 6*N + 4*B^2 + 6*B
        % Logicals = 1*B^2
        %
        % Bytes = 8 * (6*N + 4*B^2 + 6*B) + 1 * (B^2)
        %       = 48*N + 32*B^2 + 48*B + 1*B^2
        %       = 33*B^2 + 48*B + 48*N
        
        bytesPerBlock = cacheSizeMB * 1048576 / 1.2; % Safety factor 1.2
        overhead_N = 48 * N;
        
        if overhead_N < bytesPerBlock
            % We have space for the N vectors in cache, solve for B
            % 33*B^2 approx (Bytes - 48*N)
            available = bytesPerBlock - overhead_N;
            B = floor(sqrt(available / 33));
        else
            % N vectors occupy more than cache (streaming mode)
            % We just ensure the BxB blocks fit
            B = floor(sqrt(bytesPerBlock / 33));
        end
        B = max(1, B);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INITIALIZATION %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pot_r_min = H(1,1);
        pot_r_max = H(end,1);
        pot_F_min = H(1,2);
        
        px = p(:,1); 
        py = p(:,2); 
        pz = p(:,3);
        
        Fx = zeros(N,1);
        Fy = zeros(N,1);
        Fz = zeros(N,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BLOCK LOOP %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii = 1:B:N
            i_end = min(ii+B-1, N);
        
            % positions for i block
            pi_x = px(ii:i_end);
            pi_y = py(ii:i_end);
            pi_z = pz(ii:i_end);
        
            for jj = ii:B:N
                j_end = min(jj+B-1, N);
        
                % positions for j block
                pj_x = px(jj:j_end);
                pj_y = py(jj:j_end);
                pj_z = pz(jj:j_end);
        
                %% 1. Compute component deltas for block
                if S.bc == 2
                    % -------- Cubic PBC MIC --------
                    cdiff = pi_x - pj_x.'; 
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = cdiff.^2;
        
                    cdiff = pi_y - pj_y.';   % Y component
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = r2+cdiff.^2;
            
                    cdiff = pi_z - pj_z.';  % Z component
                    cdiff = cdiff - L.*round(cdiff.*invL);
                    r2 = r2+cdiff.^2;
                
                elseif S.bc == 3
                    % -------- FCC MIC (minimal alloc) --------
                
                    % raw displacement for this axis
                    dx = pi_x - pj_x.'; 
                    dy = pi_y - pj_y.'; 
                    dz = pi_z - pj_z.'; 
                
                    % reshape into 3×M (smallest possible intermediate)
                    M = numel(dx);
                    D = [dx(:)'; dy(:)'; dz(:)'];
                
                    % convert to fractional → wrap → convert back
                    Sfrac = S.fcc.invA * D;
                    Sfrac = Sfrac - round(Sfrac);
                    Dwrap = S.fcc.A * Sfrac;
                
                    % extract current axis component only
                    Dx = reshape(Dwrap(1,:), size(dx));
                    Dy = reshape(Dwrap(2,:), size(dx));
                    Dz = reshape(Dwrap(3,:), size(dx));
                    r2 = Dx.^2 + Dy.^2 + Dz.^2;
                end        
        
                % 2. Compute spherical mask and global coordinates
                mask=r2<rc^2;
                if ii == jj
                     mask = triu(mask, 1);
                end
                if ~any(mask(:)), continue; end
                [rows, cols] = find(mask);
                gI = (ii - 1) + rows;
                gJ = (jj - 1) + cols;
        
                % 4. Vector of masked distances
                r2 = sqrt(r2(mask));
                
                % 5. Interpolate forces and clamp edge cases
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
                    % 7. Calculate displacement components due to potentials (X)
                    n_pairs = length(rows);
                    dx_val = pi_x(rows) - pj_x(cols);
                    dx_val = dx_val - L .* round(dx_val .* invL);
                    dx_val = Fij .* dx_val;
                    for k = 1:n_pairs                       
                        Fx(gI(k)) = Fx(gI(k)) + dx_val(k); % Update I                     
                        Fx(gJ(k)) = Fx(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
                    end
                    
                    % 8. Calculate displacement components due to potentials (Y)
                    dx_val = pi_y(rows) - pj_y(cols);
                    dx_val = dx_val - L .* round(dx_val .* invL);
                    dx_val = Fij .* dx_val;
                    for k = 1:n_pairs                      
                        Fy(gI(k)) = Fy(gI(k)) + dx_val(k);  % Update I                     
                        Fy(gJ(k)) = Fy(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
                    end
                    
                    % 9. Calculate displacement components due to potentials (Z)
                    dx_val = pi_z(rows) - pj_z(cols);
                    dx_val = dx_val - L .* round(dx_val .* invL);
                    dx_val = Fij .* dx_val;
                    for k = 1:n_pairs           
                        Fz(gI(k)) = Fz(gI(k)) + dx_val(k); % Update I                    
                        Fz(gJ(k)) = Fz(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
                    end
                elseif S.bc == 3
                    n_pairs = length(rows);
                    % 1. Raw displacement vectors (only for the masked pairs)
                    dx = pi_x(rows) - pj_x(cols);
                    dy = pi_y(rows) - pj_y(cols);
                    dz = pi_z(rows) - pj_z(cols);
                    % 2. Build minimal 3×n_pairs temporary
                    D = [dx.'; dy.'; dz.'];       % 3 × n_pairs
                    % 3. Convert to fractional coordinates
                    Sfrac = S.fcc.invA * D;
                    % 4. FCC minimum-image: wrap each fractional component to [-0.5, 0.5]
                    Sfrac = Sfrac - round(Sfrac);
                    % 5. Convert back to Cartesian (wrapped displacement vectors)
                    Dwrap = S.fcc.A * Sfrac;      % 3 × n_pairs
                    % 6. Extract displacement components
                    dxw = Dwrap(1,:).';           % n_pairs × 1
                    dyw = Dwrap(2,:).';
                    dzw = Dwrap(3,:).';
                    % 7. Multiply by Fij * (1/r)
                    dxw = Fij .* dxw;
                    dyw = Fij .* dyw;
                    dzw = Fij .* dzw;
                    % 8. Accumulate into global forces (Newton 3rd law)
                    for k = 1:n_pairs
                        iidx = gI(k);
                        jidx = gJ(k);
                
                        Fx(iidx) = Fx(iidx) + dxw(k);
                        Fx(jidx) = Fx(jidx) - dxw(k);
                
                        Fy(iidx) = Fy(iidx) + dyw(k);
                        Fy(jidx) = Fy(jidx) - dyw(k);
                
                        Fz(iidx) = Fz(iidx) + dzw(k);
                        Fz(jidx) = Fz(jidx) - dzw(k);
                    end
                end
            end
        end
        % Finalize
        totalforces = [Fx Fy Fz];
        potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;
        maxstep = S.pot_clamp * sqrt(3) * S.stdx;
        norms = vecnorm(potdisp, 2, 2);
        overshoot = norms > maxstep;
        if any(overshoot)
            potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
        end
    end

end
