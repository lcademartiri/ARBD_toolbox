function dmin = mic_min_image_disp(pi, pj, S)
% Minimum-image displacement vector: pi - pj (row vector)
% Supports:
%   S.bc==2  cubic PBC   (A = diag(L,L,L))
%   S.bc==3  FCC/Bravais (A = S.fcc.A)
%   Anything else -> no PBC

d = pi(:)' - pj(:)';

if S.bc == 2
    % cubic PBC
    L = 2 * S.br;
    dmin = d - L .* round(d ./ L);

elseif S.bc == 3
    % FCC / general Bravais
    A = S.fcc.A;
    invA = S.fcc.invA;
    frac = invA * d(:);   % fractional coords (col vector)
    frac = frac - round(frac);
    dmin = (A * frac).'; % back to row vector

else
    dmin = d;
end
end


function pos_wrapped = mic_wrap_positions(pos, S)
% Wraps Nx3 Cartesian positions into canonical cell
pos_wrapped = pos;

if S.bc == 2
    L = 2*S.br;
    pos_wrapped = mod(pos + 0.5*L, L) - 0.5*L;

elseif S.bc == 3
    A = S.fcc.A; invA = S.fcc.invA;
    f = (invA * pos.').';
    f = f - floor(f);
    pos_wrapped = (A * f.').';
end
end


function D = mic_all_pair_displacements(pos, S)
% Debug: compute all unique pair MIC displacements
N = size(pos,1);
M = N*(N-1)/2;
D = zeros(M,3);
idx = 1;
for i = 1:N-1
    for j = i+1:N
        D(idx,:) = mic_min_image_disp(pos(i,:), pos(j,:), S);
        idx = idx + 1;
    end
end
end


function shifts = mic_cartesian_shifts(S)
% Returns the standard 27 MIC image shifts for PBC neighbor searches
% For cubic: shifts are (-L,0,+L) in each direction
% For FCC: shifts are A * fractional_shifts

if S.bc == 2
    L = 2*S.br;
    S3 = [-1 0 1];
    [sx,sy,sz] = ndgrid(S3,S3,S3);
    shifts = [sx(:)*L, sy(:)*L, sz(:)*L];

elseif S.bc == 3
    A = S.fcc.A;
    S3 = [-1 0 1];
    [sx,sy,sz] = ndgrid(S3,S3,S3);
    frac = [sx(:), sy(:), sz(:)];
    shifts = (A * frac.').';

else
    shifts = zeros(1,3);
end
end
