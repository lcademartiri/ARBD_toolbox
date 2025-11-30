function potdisp = potential_displacements_v12(p, S, H_mat, H_interpolant, ghostghost, cacheSizeMB)
% Blockified version of potential_displacements_v3
% Automatically selects block size based on cache size
% Handles cubic PBC (S.bc == 2) only
%
% p: Nx3 positions
% S: simulation struct, needs fields:
%      - bc (must be 2)
%      - br
%      - rc
%      - esdiff, kbT, timestep, pot_clamp, stdx
%
% H_mat: potential table Nx2 [r F]
% H_interpolant: function handle returning force magnitude F(r)
% ghostghost: kept for compatibility
% cacheSizeMB: optional, override CPU L3 cache size (per core). Default = 20 MB.

if nargin < 5, ghostghost = 0; end
if nargin < 6 || isempty(cacheSizeMB)
    cacheSizeMB = 20;   % default L3 per core
end

N = size(p,1);
N_real=N;
L = 2*S.br;
invL=1/L;
rc = S.rc;

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
pot_r_min = H_mat(1,1);
pot_r_max = H_mat(end,1);
pot_F_min = H_mat(1,2);

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

        
        elseif S.bc == 1 || S.bc == 4
            % -------- SBC or Big Box: NO MIC --------
            cdiff = pi_x - pj_x.';   % X component
            r2 = cdiff.^2;
            
            cdiff = pi_y - pj_y.';   % Y component
            r2 = r2+cdiff.^2;
    
            cdiff = pi_z - pj_z.';  % Z component
            r2 = r2+cdiff.^2;
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

        % 3. ghost-ghost filter
        if S.bc==1 && ghostghost==0        
            keep = (gI <= N_real) | (gJ <= N_real);    
            if ~any(keep)
                continue;
            end
            rows = rows(keep);
            cols = cols(keep);
            gI   = gI(keep);
            gJ   = gJ(keep);
        end

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
        
        if S.bc==1 || S.bc==4
            n_pairs = length(rows);
            dx_val = pi_x(rows) - pj_x(cols);
            dx_val = Fij .* dx_val;
            for k = 1:n_pairs                       
                Fx(gI(k)) = Fx(gI(k)) + dx_val(k); % Update I                     
                Fx(gJ(k)) = Fx(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
            end
            dx_val = pi_y(rows) - pj_y(cols);
            dx_val = Fij .* dx_val;
            for k = 1:n_pairs                      
                Fy(gI(k)) = Fy(gI(k)) + dx_val(k);  % Update I                     
                Fy(gJ(k)) = Fy(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
            end
            dx_val = pi_z(rows) - pj_z(cols);
            dx_val = Fij .* dx_val;
            for k = 1:n_pairs           
                Fz(gI(k)) = Fz(gI(k)) + dx_val(k); % Update I                    
                Fz(gJ(k)) = Fz(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
            end
        elseif S.bc==2
            % 7. Calculate displacement components due to potentials (X)
            n_pairs = length(rows);
            dx_val = pi_x(rows) - pj_x(cols);
            dx_val = dx_val - L .* round(dx_val .* invL);
            dx_val = Fij .* dx_val;
            for k = 1:n_pairs                       
                Fx(gI) = Fx(gI(k)) + dx_val(k); % Update I                     
                Fx(gJ) = Fx(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
            end
            
            % 8. Calculate displacement components due to potentials (Y)
            dx_val = pi_y(rows) - pj_y(cols);
            dx_val = dx_val - L .* round(dx_val .* invL);
            dx_val = Fij .* dx_val;
            for k = 1:n_pairs                      
                Fy(gI) = Fy(gI(k)) + dx_val(k);  % Update I                     
                Fy(gJ) = Fy(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
            end
            
            % 9. Calculate displacement components due to potentials (Z)
            dx_val = pi_z(rows) - pj_z(cols);
            dx_val = dx_val - L .* round(dx_val .* invL);
            dx_val = Fij .* dx_val;
            for k = 1:n_pairs           
                Fz(gI) = Fz(gI(k)) + dx_val(k); % Update I                    
                Fz(gJ) = Fz(gJ(k)) - dx_val(k); % Update J (Newton's 3rd)
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
%% Finalize
totalforces = [Fx Fy Fz];
potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;
maxstep = S.pot_clamp * sqrt(3) * S.stdx;
norms = vecnorm(potdisp, 2, 2);
overshoot = norms > maxstep;
if any(overshoot)
    potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
end

end
