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

r2 = xdiff.^2 + ydiff.^2 + zdiff.^2;
r = sqrt(r2);

mask = mask & (r <= rc);

%% ------------------------------------------------------------------------
% 6) Upper-triangle only
maskUT = triu(mask,1);

[i_idx, j_idx] = find(maskUT);

pairs_i = uint32(i_idx);
pairs_j = uint32(j_idx);

% Extract displacements / distances exactly in this order
d_mic = [xdiff(maskUT), ydiff(maskUT), zdiff(maskUT)];
rr    = r(maskUT);

K = numel(rr);
if K == 0
    disppot = zeros(N,3);
    return;
end

%% ------------------------------------------------------------------------
% 7) Compute force magnitudes
if isa(H_interpolant,'function_handle')
    Fmag = H_interpolant(rr);
elseif isstruct(H) && isfield(H,'r') && isfield(H,'f')
    Fmag = interp1(H.r(:), H.f(:), rr, 'linear', 0);
else
    Fmag = (rr>0).* (1./(rr+eps));
end

% Vector forces
dirs = d_mic ./ rr;
Fvec = Fmag .* dirs;  % Kx3

%% ------------------------------------------------------------------------
% 8) Accumulate ±F
Fx = accumarray([i_idx; j_idx], [Fvec(:,1); -Fvec(:,1)], [N 1]);
Fy = accumarray([i_idx; j_idx], [Fvec(:,2); -Fvec(:,2)], [N 1]);
Fz = accumarray([i_idx; j_idx], [Fvec(:,3); -Fvec(:,3)], [N 1]);

disppot = [Fx Fy Fz];
end
