function [disppot,pairs_i,pairs_j,d_mic] = potential_displacements_v3(p, S, H, H_interpolant, ghostghost)
% Ultra-fast brute-force MIC + truncation version (cubic PBC only)
% Complexity: O(N^2) but heavily optimized vectorized.
% Best for N <= ~3000.

if nargin < 5, ghostghost = 0; end
if S.bc ~= 2
    error('This optimized version only supports cubic PBC (S.bc==2).');
end

N = size(p,1);

L = 2*S.br;
rc = S.rc;

%% ------------------------------------------------------------------------
% 1) Preallocate output force array
disppot = zeros(N,3);

%% ------------------------------------------------------------------------
% 2) Construct FULL difference matrices (vectorized, no loops)
%
% xdiff(i,j) = p(i,1) - p(j,1)
% Do this component-wise using broadcasting (bsxfun) or implicit expansion
px = p(:,1);  py = p(:,2);  pz = p(:,3);

xdiff = px - px.';     % NxN
ydiff = py - py.';     % NxN
zdiff = pz - pz.';     % NxN

%% ------------------------------------------------------------------------
% 3) Apply MIC once for all pairs
xdiff = xdiff - L .* round(xdiff ./ L);
ydiff = ydiff - L .* round(ydiff ./ L);
zdiff = zdiff - L .* round(zdiff ./ L);

%% ------------------------------------------------------------------------
% 4) Cheap prefilter BEFORE sqrt
% Keep only entries satisfying |dx|<=rc, |dy|<=rc, |dz|<=rc
mask = (abs(xdiff) <= rc) & (abs(ydiff) <= rc) & (abs(zdiff) <= rc);

% remove diagonal (i=j)
mask(1:N+1:end) = false;

%% ------------------------------------------------------------------------
% 5) Compute distances only where needed
r2 = zeros(size(xdiff));
r2(mask) = xdiff(mask).^2 + ydiff(mask).^2 + zdiff(mask).^2;
mask = mask & (r2 <= rc^2);

%% ------------------------------------------------------------------------
% 6) enforce i < j without using triu
mask = mask & ((1:N)' < 1:N);

[pairs_i, pairs_j] = find(mask);

% Extract displacements in same ordering as find(mask)
d_mic = [xdiff(mask), ydiff(mask), zdiff(mask)];
r     = sqrt(r2(mask));

K = numel(r);
if K == 0
    disppot = zeros(N,3);
    return;
end
%% ------------------------------------------------------------------------
% 7) Compute force magnitudes
pot_r_min = H(1,1);
pot_r_max = H(end,1);
pot_F_min = H(1,2);

Fij = H_interpolant(r);

Fij(r >= rc | r >= pot_r_max) = 0;
Fij(isnan(Fij) & r < pot_r_min) = pot_F_min;
Fij(isnan(Fij) & r >= pot_r_max) = 0;

inv_r = 1./r;
inv_r(isinf(inv_r))=0;
fx = Fij .* d_mic(:,1) .* inv_r;
fy = Fij .* d_mic(:,2) .* inv_r;
fz = Fij .* d_mic(:,3) .* inv_r;
N_total = max([max(pairs_i), max(pairs_j), N]);

Fx = accumarray(pairs_i, fx, [N_total,1]) - accumarray(pairs_j, fx, [N_total,1]);
Fy = accumarray(pairs_i, fy, [N_total,1]) - accumarray(pairs_j, fy, [N_total,1]);
Fz = accumarray(pairs_i, fz, [N_total,1]) - accumarray(pairs_j, fz, [N_total,1]);

totalforces = [Fx Fy Fz];

% convert to displacements
potdisp = totalforces * (S.esdiff / S.kbT) * S.timestep;

% clamp
maxstep = S.pot_clamp * sqrt(3) * S.stdx;
norms = vecnorm(potdisp,2,2);
overshoot = norms > maxstep;
if any(overshoot)
    potdisp(overshoot,:) = potdisp(overshoot,:) .* (maxstep ./ norms(overshoot));
end

disppot = potdisp;
end
