function disppot = potential_displacements_v11(p, S, H_mat, H_interpolant, ghostghost, cacheSizeMB)
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
if S.bc ~= 2
    error('Blocked version only supports cubic PBC (S.bc == 2).');
end
if nargin < 6 || isempty(cacheSizeMB)
    cacheSizeMB = 20;   % default L3 per core
end

N = size(p,1);
L = 2*S.br;
invL=1/L;
rc = S.rc;

%% 0. Choose block size based on cache
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

pot_r_min = H_mat(1,1);
pot_r_max = H_mat(end,1);
pot_F_min = H_mat(1,2);

px = p(:,1); 
py = p(:,2); 
pz = p(:,3);

Fx = zeros(N,1);
Fy = zeros(N,1);
Fz = zeros(N,1);

%% BLOCK LOOPS
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

        % BLOCK SIZE
        ni = i_end - ii + 1;
        nj = j_end - jj + 1;
        r2 = zeros(numel(pi_x),numel(pj_x));

        %% 1. Compute differences for block
        cdiff = pi_x - pj_x.';   % ni x nj        
        cdiff = cdiff - L.*round(cdiff.*invL);
        r2 = cdiff.^2;

        cdiff = pi_y - pj_y.';   % ni x nj
        cdiff = cdiff - L.*round(cdiff./L);
        r2 = r2+cdiff.^2;

        cdiff = pi_z - pj_z.';   % ni x nj
        cdiff = cdiff - L.*round(cdiff./L);
        r2 = r2+cdiff.^2;

        mask=r2<rc^2;
        if ii == jj
             mask = triu(mask, 1);
        end
        if ~any(mask(:)), continue; end
        [rows, cols] = find(mask);

        % =================================================================
        % PASS 2: Calculate Scalar Force Factor
        % =================================================================
        % We calculate F(r)/r. This depends only on r2.
        
        rr = sqrt(r2(mask));
        
        Fij = H_interpolant(rr);
        Fij(rr >= rc | rr >= pot_r_max) = 0;
        if any(isnan(Fij))
            bad = isnan(Fij);
            Fij(bad & rr < pot_r_min) = pot_F_min;
            Fij(bad & rr >= pot_r_max) = 0;
        end
        % Projection factor: F / r
        inv_r = 1./rr;
        inv_r(isinf(inv_r)) = 0;
        
        % Sparse Scalar Force stored in Dense Matrix shape
        % Using zeros() is faster than sparse() for temporary math
        f_scalar = Fij .* inv_r;

        % =================================================================
        % PASS 3: Vectorized Recalculation (Low Memory)
        % =================================================================
        % 1. Gather coordinates (pi_x(rows) is vector lookup)
        % 2. Compute component force vectors
        
        dx_val = pi_x(rows) - pj_x(cols);
        dx_val = dx_val - L .* round(dx_val .* invL);
        fx_val = f_scalar .* dx_val;
        
        dy_val = pi_y(rows) - pj_y(cols);
        dy_val = dy_val - L .* round(dy_val .* invL);
        fy_val = f_scalar .* dy_val;
        
        dz_val = pi_z(rows) - pj_z(cols);
        dz_val = dz_val - L .* round(dz_val .* invL);
        fz_val = f_scalar .* dz_val;

        % =================================================================
        % PASS 4: Direct Accumulation (The requested optimization)
        % =================================================================
        % This replaces 'accumarray'. It iterates only over the valid pairs.
        % Using global indices for direct update.
        
        Ii = (ii - 1) + rows;
        Jj = (jj - 1) + cols;
        n_pairs = length(Ii);
        
        for k = 1:n_pairs
            gI = Ii(k);
            gJ = Jj(k);
            
            % Update I
            Fx(gI) = Fx(gI) + fx_val(k);
            Fy(gI) = Fy(gI) + fy_val(k);
            Fz(gI) = Fz(gI) + fz_val(k);
            
            % Update J (Newton's 3rd)
            Fx(gJ) = Fx(gJ) - fx_val(k);
            Fy(gJ) = Fy(gJ) - fy_val(k);
            Fz(gJ) = Fz(gJ) - fz_val(k);
        end


        % =================================================================
        % Optional: Pairs
        % =================================================================
        if return_pairs
            if pair_count + n_pairs > length(I_list)
                growth = max(n_pairs, 100000);
                I_list(end+growth) = 0;
                J_list(end+growth) = 0;
                d_list(end+growth,3) = 0;
            end
            
            rng = pair_count+1 : pair_count+n_pairs;
            I_list(rng) = Ii;
            J_list(rng) = Jj;
            d_list(rng, :) = [dx_val, dy_val, dz_val];
            pair_count = pair_count + n_pairs;
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

disppot = potdisp;

if return_pairs
    pairs_i = uint32(I_list(1:pair_count));
    pairs_j = uint32(J_list(1:pair_count));
    d_mic_out = d_list(1:pair_count, :);
end
end
