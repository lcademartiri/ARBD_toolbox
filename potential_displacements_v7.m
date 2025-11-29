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
bytesPerBlock = cacheSizeMB * 1e6 / 1.2;
B = floor( sqrt(bytesPerBlock / (3*8)) );
B = max(200, min(B, 2000));  % enforce reasonable bounds

% For output
Fx = zeros(N,1);
Fy = zeros(N,1);
Fz = zeros(N,1);

% Keep lists of global pairs (optional but kept for compatibility)
I_list = [];
J_list = [];
d_list = [];

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

    for jj = ii+B:B:N
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
        r2 = zeros(ni,nj);
        mask_idx = mask(:);
        if any(mask_idx)
            tmp_r2 = xdiff(mask_idx).^2 + ydiff(mask_idx).^2 + zdiff(mask_idx).^2;
            mask(mask_idx) = mask(mask_idx) & (tmp_r2 <= rc^2);
            r2(mask) = tmp_r2;
        end

        %% 5. Extract valid (i_local, j_local)
        [ii_local, jj_local] = find(mask);

        if isempty(ii_local), continue; end

        % global indices
        Ii = ii - 1 + ii_local;
        Jj = jj - 1 + jj_local;

        % distances & displacement vectors
        rr = sqrt(r2(mask));
        dloc = [ xdiff(mask), ydiff(mask), zdiff(mask) ];

        %% 6. Force magnitude from interpolant
        r_clamped = min(max(rr,pot_r_min), pot_r_max);
        Fij = H_interpolant(r_clamped);

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
        Fx = Fx + accumarray(Ii, fx, [N 1]) - accumarray(Jj, fx, [N 1]);
        Fy = Fy + accumarray(Ii, fy, [N 1]) - accumarray(Jj, fy, [N 1]);
        Fz = Fz + accumarray(Ii, fz, [N 1]) - accumarray(Jj, fz, [N 1]);

        %% record pairs (optional)
        I_list = [I_list; Ii];
        J_list = [J_list; Jj];
        d_list = [d_list; dloc];
    end
end

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
