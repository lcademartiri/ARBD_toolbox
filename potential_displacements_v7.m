function [disppot,pairs_i,pairs_j,d_mic_out] = potential_displacements_v7(p, S, H_mat, H_interpolant, ghostghost, cacheSizeMB)
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
rc = S.rc;

%% 0. Choose block size based on cache
% Memory per double = 8 bytes
% We must fit approx 3 matrices of size BxB into cache (xdiff, ydiff, zdiff)
% overhead factor = 1.2 for safety
% doubles=3*N... Fs
%     +3*N... ps
%     +3*N... p
%     +3*B... pis
%     +3*B... pjs
%     +3*B^2 ... diffs
%     +B^2 ... tmp_r2 (up to)
%     +3*B^2 ... r2+rr+invr
%     +2*B^2 ... ii_local & jj_local (up to)
%     +3*B^2; % d_loc (up to)
% doubles=9*N+12*B^2+6*B;
% 
% 
% 
% logicals=B^2 ... mask
%     + B^2; % mask_flat
% logicals=2*B^2;

bytesPerBlock = cacheSizeMB * 1048576 / 1.2;
if (72*N) < bytesPerBlock
    B = floor( sqrt( (bytesPerBlock - 72*N) / 98 ) );
else
    B = floor( sqrt( bytesPerBlock / 98 ) );
end
% B = max(200, min(B, 2000));  % enforce reasonable bounds

% For output
Fx = zeros(N,1);
Fy = zeros(N,1);
Fz = zeros(N,1);

% Before block loops, estimate max pairs
max_pairs = N * min(100, ceil(4/3*pi*rc^3 * N/(L^3)));
I_list = zeros(max_pairs, 1, 'uint32');
J_list = zeros(max_pairs, 1, 'uint32');
d_list = zeros(max_pairs, 3);
pair_count = 0;

pot_r_min = H_mat(1,1);
pot_r_max = H_mat(end,1);
pot_F_min = H_mat(1,2);

px = p(:,1); py = p(:,2); pz = p(:,3);

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

        %% 1. Compute differences for block
        xdiff = pi_x - pj_x.';   % ni x nj
        ydiff = pi_y - pj_y.';   % ni x nj
        zdiff = pi_z - pj_z.';   % ni x nj

        %% 2. MIC
        xdiff = xdiff - L.*round(xdiff./L);
        ydiff = ydiff - L.*round(ydiff./L);
        zdiff = zdiff - L.*round(zdiff./L);

        %% 3. Component mask
        mask = (abs(xdiff)<=rc) & (abs(ydiff)<=rc) & (abs(zdiff)<=rc);

        %% 4. r2 only for masked entries
        mask_flat = mask(:);
        if ~any(mask_flat), continue; end
        tmp_r2 = xdiff(mask_flat).^2 + ydiff(mask_flat).^2 + zdiff(mask_flat).^2;
        r2 = zeros(ni, nj);
        r2(mask) = tmp_r2;
        % radial cutoff
        mask = mask & (r2 <= rc^2);
        if ~any(mask(:)), continue; end

        %% 5. Extract valid (i_local, j_local)
        [ii_local, jj_local] = find(mask);

        if isempty(ii_local), continue; end

        % If diagonal block (jj==ii) filter to i_local < j_local
        if jj == ii
            ok = ( (ii_local - 1) < (jj_local - 1) );       % since both are offsets inside same block
            if ~any(ok), continue; end
            ii_local = ii_local(ok);
            jj_local = jj_local(ok);
            % FILTER ALL dependent vectors
            mask_idx = find(mask);       % linear indices of final mask
            mask_idx = mask_idx(ok);     % keep only ones consistent with diag-filter
        
            % rebuild dloc, rr, Fij inputs using filtered mask
            dloc = [xdiff(mask_idx), ydiff(mask_idx), zdiff(mask_idx)];
            rr   = sqrt(r2(mask_idx));
        else
            mask_idx = find(mask);
            dloc = [xdiff(mask_idx), ydiff(mask_idx), zdiff(mask_idx)];
            rr   = sqrt(r2(mask_idx));
        end

        % global indices
        Ii = ii - 1 + ii_local;
        Jj = jj - 1 + jj_local;

        %% 6. Force magnitude from interpolant
        Fij = H_interpolant(rr);

        Fij(rr>=rc | rr>=pot_r_max) = 0;
        bad = isnan(Fij);
        Fij(bad & rr<pot_r_min) = pot_F_min;
        Fij(bad & rr>=pot_r_max) = 0;

        inv_r = 1./rr;
        inv_r(isinf(inv_r)) = 0;

        fx = Fij .* dloc(:,1) .* inv_r;
        fy = Fij .* dloc(:,2) .* inv_r;
        fz = Fij .* dloc(:,3) .* inv_r;

        %% 7. Accumulate forces
        % Fx = Fx + accumarray(Ii, fx, [N 1]) - accumarray(Jj, fx, [N 1]);
        % Fy = Fy + accumarray(Ii, fy, [N 1]) - accumarray(Jj, fy, [N 1]);
        % Fz = Fz + accumarray(Ii, fz, [N 1]) - accumarray(Jj, fz, [N 1]);

        %% 7. Accumulate forces using sparse matrix trick
        np = length(Ii);  % number of pairs in this block

        % Build sparse accumulation matrices (only once per block)
        S_add = sparse([Ii; Jj], [1:np, 1:np], [ones(np,1); -ones(np,1)], N, np);

        % Accumulate all 3 dimensions at once
        F_block = S_add * [fx, fy, fz];  % N x 3

        Fx = Fx + F_block(:,1);
        Fy = Fy + F_block(:,2);
        Fz = Fz + F_block(:,3);

        % %% 7. Accumulate forces directly
        % for k = 1:length(Ii)
        %     Fx(Ii(k)) = Fx(Ii(k)) + fx(k);
        %     Fx(Jj(k)) = Fx(Jj(k)) - fx(k);
        % 
        %     Fy(Ii(k)) = Fy(Ii(k)) + fy(k);
        %     Fy(Jj(k)) = Fy(Jj(k)) - fy(k);
        % 
        %     Fz(Ii(k)) = Fz(Ii(k)) + fz(k);
        %     Fz(Jj(k)) = Fz(Jj(k)) - fz(k);
        % end

        % %% 7. Accumulate forces with combined indexing
        % np = length(Ii);
        % 
        % % Create combined index and value arrays
        % idx_combined = [Ii; Jj];           % length 2*np
        % vals_x = [fx; -fx];                 % length 2*np
        % vals_y = [fy; -fy];
        % vals_z = [fz; -fz];
        % 
        % % Single accumarray per dimension
        % Fx = Fx + accumarray(idx_combined, vals_x, [N 1]);
        % Fy = Fy + accumarray(idx_combined, vals_y, [N 1]);
        % Fz = Fz + accumarray(idx_combined, vals_z, [N 1]);

        %% record pairs (optional)
        % Inside block loops, replace concatenation:
        n_new = length(Ii);
        I_list(pair_count+1:pair_count+n_new) = Ii;
        J_list(pair_count+1:pair_count+n_new) = Jj;
        d_list(pair_count+1:pair_count+n_new, :) = dloc;
        pair_count = pair_count + n_new;
    end
end
% After loops, trim:
I_list = I_list(1:pair_count);
J_list = J_list(1:pair_count);
d_list = d_list(1:pair_count, :);

%% 8. Final forces → displacements
totalforces = [Fx Fy Fz];
potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;

maxstep = S.pot_clamp * sqrt(3) * S.stdx;
norms = vecnorm(potdisp,2,2);
overshoot = norms > maxstep;
if any(overshoot)
    potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
end

disppot = potdisp;

pairs_i = uint32(I_list);
pairs_j = uint32(J_list);
d_mic_out = d_list;
end
