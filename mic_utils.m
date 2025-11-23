% === FILE: mic_utils.m ===
function dmat = mic_pairwise_distances(pos, S)
% Return condensed pairwise distance vector (same format as pdist output)
% pos: Nx3, S: structure with S.bc, S.fcc.A & S.fcc.invA or S.br for cubic
N = size(pos,1);
idx = 1;
dmat = zeros(N*(N-1)/2,1);
for i=1:N-1
    for j=i+1:N
        d = mic_min_image_disp(pos(i,:), pos(j,:), S);
        dmat(idx) = norm(d);
        idx = idx + 1;
    end
end
end

function D = mic_all_pair_displacements(pos, S)
% Return Mx3 displacements (both directions removed) for all unique pairs
N = size(pos,1);
M = N*(N-1)/2;
D = zeros(M,3);
idx = 1;
for i=1:N-1
    for j=i+1:N
        D(idx,:) = mic_min_image_disp(pos(i,:), pos(j,:), S);
        idx = idx + 1;
    end
end
end

function dmin = mic_min_image_disp(pi, pj, S)
% Compute minimum-image displacement vector from pj to pi (pi - pj)
d = pi(:)' - pj(:)';
if S.bc==2 % cubic PBC (orthorhombic)
    L = 2*S.br; % box side = diameter of sphere of same volume (existing convention)
    % component-wise MIC:
    dmin = d - L .* round(d ./ L);
elseif S.bc==3 % FCC/general Bravais
    A = S.fcc.A; invA = S.fcc.invA;
    frac = invA * d(:);            % fractional coords
    frac = frac - round(frac);     % wrap to [-0.5,0.5]
    dmin = (A * frac)';            % back to Cartesian row-vector
else
    % no PBC, plain displacement
    dmin = d;
end
end

function pos_wrapped = mic_wrap_positions(pos, S)
% Wrap Cartesian positions into canonical cell using MIC conventions
% pos: Nx3
pos_wrapped = pos;
if S.bc==2
    L = 2*S.br;
    pos_wrapped = mod(pos + 0.5*L, L) - 0.5*L; % center at origin convention
elseif S.bc==3
    A = S.fcc.A; invA = S.fcc.invA;
    % convert to fractional, wrap to [0,1), convert back
    f = (invA * pos')';
    f = f - floor(f); % in [0,1)
    pos_wrapped = (A * f')';
end
end
