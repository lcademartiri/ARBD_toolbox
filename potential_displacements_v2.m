function [disppot,pairs_i,pairs_j,d_mic] = potential_displacements_v2(p, S, H, H_interpolant, ghostghost)
% Unified SBC + PBC (Cubic/FCC) potential displacements with MIC
% p: N_total x 3
% S: structure containing
%       .bc           boundary condition (1=SBC, 2=Cubic PBC, 3=FCC PBC)
%       .N            number of real particles
%       .rc           cutoff
%       .br           for cubic box, L = 2*br
%       .fcc.A        lattice matrix for PBC case 3
%       .fcc.invA     inverse lattice matrix
% H, H_interpolant    force table interpolant + forcevector
% ghostghost          SBC ghost filtering flag

if S.bc == 1
    % ============= SBC branch ===========================
    [pairs_i, pairs_j, d_mic] = get_pairs_SBC(p, S, ghostghost);

elseif S.bc == 2
    % ============= Cubic PBC + MIC =======================
    [pairs_i, pairs_j, d_mic] = get_pairs_PBC_cubic_MIC(p, S);

elseif S.bc == 3
    % ============= FCC PBC + MIC =========================
    [pairs_i, pairs_j, d_mic] = get_pairs_PBC_FCC_MIC(p, S);

elseif S.bc == 4
    % ============= Big BOX ===============================
    [pairs_i, pairs_j, d_mic] = get_pairs_BB(p, S);
end
N=size(p,1);
% === Compute forces (shared for all BCs) ====================
disppot = compute_forces_from_pairs(pairs_i, pairs_j, d_mic, S, H, H_interpolant,N);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% === SBC PAIRS (original rangesearch + ghosts preserved) =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pairs_i, pairs_j, d_mic] = get_pairs_SBC(p, S, ghostghost)

N_total = size(p,1);
N = S.N;
rc = S.rc;

[idx_list, ~] = rangesearch(p, p, rc);

num_neighbors = cellfun(@numel, idx_list);
pairs_j = horzcat(idx_list{:})';
pairs_i = repelem((1:N_total)', num_neighbors);
% distances = horzcat(dist_list{:})';

mask = pairs_i < pairs_j;
if ghostghost == 0
    mask = mask & ~(pairs_i > N & pairs_j > N);
end

pairs_i = pairs_i(mask);
pairs_j = pairs_j(mask);

% raw displacement (no MIC in SBC)
d_mic = p(pairs_i,:) - p(pairs_j,:);% SBC uses real-space displacement 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% === Cubic PBC MIC pair generation (cell list) ===========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pairs_i, pairs_j, d_mic] = get_pairs_PBC_cubic_MIC(p, S)

N = size(p, 1);
L = 2*S.br;

% shift coordinates into [0, L)
p_shift = mod(p + 0.5*L, L);

cell_idx = floor(p_shift / S.cellsize) + 1;
cell_id = cell_idx(:,1) + S.ncell*(cell_idx(:,2)-1) + S.ncell*S.ncell*(cell_idx(:,3)-1);

% vectorized bucket: create a cell array where each cell contains indices
% of particles that fall in that spatial cell
tmp = accumarray(cell_id(:), (1:N)', [], @(x) {x});   % tmp is length max(cell_id)
cell_particles = cell(S.ncell^3,1);                     % preallocate full list
cell_particles(1:numel(tmp)) = tmp;                   % fill existing buckets
% any remaining cell_particles entries remain empty {}
pairs_i = [];
pairs_j = [];

for cid = 1:S.ncell^3
    plist = cell_particles{cid};
    if isempty(plist), continue; end
    % A. Self-Interaction (Particles within the same cell)
    % We need all unique pairs (i,j) where i < j within this list.
    n_p = length(plist);
    if n_p > 1
        % Generate pairs using nchoosek-like logic or simple meshgrid + mask
        % For speed in loops:
        % (If plist is 1:N, we want 1-2, 1-3... 2-3...)
        [I_self, J_self] = meshgrid(plist, plist);
        mask_self = I_self < J_self; % strictly less to avoid self and duplicates
        pairs_i = [pairs_i; I_self(mask_self)];
        pairs_j = [pairs_j; J_self(mask_self)];
    end

    % B. Neighbor Interactions (Half-Shell)
    % Check against the 13 forward neighbors
    for k = 1:13
        nid = S.neighbor_linear{k}(cid); % Ensure neighbor_linear is built for 13
        qlist = cell_particles{nid};
        if isempty(qlist), continue; end

        % Cross interaction: All p in Current vs All q in Neighbor
        % Since cells are distinct, we take ALL pairs (Cartesian product).
        % No need for i < j check because indices are in different cells.
        [I, J] = ndgrid(plist, qlist);
        pairs_i = [pairs_i; I(:)];
        pairs_j = [pairs_j; J(:)];
    end
end

% MIC displacement vectors
d_raw = p(pairs_i,:) - p(pairs_j,:);
d_mic = d_raw - L .* round(d_raw./L);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% === FCC PBC MIC pair generation (fractional cell list) ==
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pairs_i, pairs_j, d_mic] = get_pairs_PBC_FCC_MIC(p, S)

N = size(p, 1);
A = S.fcc.A;
invA = S.fcc.invA;

% 1. Convert to Fractional Coordinates [0,1)
f = (invA * p.').';
% Wrap to unit cell (periodic boundaries)
f = f - floor(f); 

% 2. Bin into Cell Lists
cell_idx = floor(f / S.cellsize) + 1;
% Clamp to be safe against 1.0
cell_idx = min(cell_idx, S.ncell); 

cell_id = cell_idx(:,1) + S.ncell*(cell_idx(:,2)-1) + S.ncell*S.ncell*(cell_idx(:,3)-1);

% Vectorized bucket creation
tmp = accumarray(cell_id(:), (1:N)', [], @(x) {x});
cell_particles = cell(S.ncell^3, 1);
cell_particles(1:numel(tmp)) = tmp;

pairs_i = [];
pairs_j = [];

% 3. Iterate over Cells
for cid = 1:S.ncell^3
    plist = cell_particles{cid};
    if isempty(plist), continue; end

    % A. Self-Interaction (Particles within the same cell)
    % Get unique pairs (i < j) inside this cell
    n_p = length(plist);
    if n_p > 1
        [I_self, J_self] = meshgrid(plist, plist);
        mask_self = I_self < J_self; 
        pairs_i = [pairs_i; I_self(mask_self)];
        pairs_j = [pairs_j; J_self(mask_self)];
    end

    % B. Neighbor Interactions (Half-Shell)
    % Check against the 13 forward neighbors defined in Setup
    for k = 1:13
        nid = S.neighbor_linear{k}(cid);
        qlist = cell_particles{nid};
        if isempty(qlist), continue; end

        % Cross interaction: All p vs All q
        [I, J] = ndgrid(plist, qlist);
        pairs_i = [pairs_i; I(:)];
        pairs_j = [pairs_j; J(:)];
    end
end

% 4. MIC Fractional Displacement (Final Calculation)
% (No 'mask' needed anymore because we eliminated duplicates logic-side)
d_raw = p(pairs_i,:) - p(pairs_j,:);

% Fractional displacement
frac = (invA * d_raw.').';
frac = frac - round(frac); % Apply MIC in fractional space

% Back to Cartesian
d_mic = (A * frac.').';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% === BB pair generation (Half-Shell Cell List) ===========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pairs_i, pairs_j, d_mic] = get_pairs_BB(p, S)

% --- BB setup ---
% L is the size of the Big Box (2*S.br)
L = 2 * S.br;
rc = S.rc;

% Determine grid size
ncell = max(1, floor(L / rc));
cellsize = L / ncell;

% 1. Shift positions to [0, L) frame
% Assuming 'p' is centered at 0 (range -br to +br)
p_shift = p + S.br; 

% 2. Calculate Cell Indices
cell_sub = floor(p_shift / cellsize) + 1;

% 3. Filter Valid Particles (Remove OOB)
% In BB, we simply ignore particles outside the grid.
valid_mask = all(cell_sub >= 1, 2) & all(cell_sub <= ncell, 2);

% Only proceed with valid particles
valid_ids = find(valid_mask);
if isempty(valid_ids)
    pairs_i = []; pairs_j = []; d_mic = []; return;
end

valid_sub = cell_sub(valid_mask, :);

% 4. Linear Cell ID
cell_id = valid_sub(:,1) + ...
          ncell*(valid_sub(:,2)-1) + ...
          ncell*ncell*(valid_sub(:,3)-1);

% 5. Vectorized Bucket Fill
% We map the *original* particle indices (valid_ids) into the cells
tmp = accumarray(cell_id, valid_ids, [ncell^3, 1], @(x) {x});
cell_particles = tmp;

% 6. Build Pairs (Half-Shell Logic)
pairs_i = [];
pairs_j = [];

% Iterate over non-empty cells
for cid = 1:ncell^3
    plist = cell_particles{cid};
    if isempty(plist), continue; end

    % A. Self-Interaction (Particles within the same cell)
    % Get unique pairs (i < j) inside this cell
    n_p = length(plist);
    if n_p > 1
        [I_self, J_self] = meshgrid(plist, plist);
        mask_self = I_self < J_self; 
        pairs_i = [pairs_i; I_self(mask_self)];
        pairs_j = [pairs_j; J_self(mask_self)];
    end

    % B. Neighbor Interactions (13 Forward Neighbors)
    % Check against the 13 forward neighbors defined in Setup
    for k = 1:13
        % Get the Neighbor ID from the precomputed table
        nid = S.neighbor_linear(k, cid); 
        
        % In BB, neighbor_linear has 0 for OOB neighbors (edges). Skip them.
        if nid == 0, continue; end
        
        qlist = cell_particles{nid};
        if isempty(qlist), continue; end

        % Cross interaction: All p vs All q
        % Since we only visit each pair of cells once (forward), we take all pairs.
        [I, J] = ndgrid(plist, qlist);
        pairs_i = [pairs_i; I(:)];
        pairs_j = [pairs_j; J(:)];
    end
end

% 7. Displacement (Euclidean, No MIC)
% No need to filter unique pairs anymore because Half-Shell logic guarantees uniqueness.
d_mic = p(pairs_i,:) - p(pairs_j,:);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% === Shared force module =================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disppot = compute_forces_from_pairs(pairs_i, pairs_j, d_mic, S, H, H_interpolant,N)

% --- CASE 1: NO INTERACTIONS (Zero Force) ---
if isempty(pairs_i)
    % If pairs are empty, nobody is pushing anyone.
    % Return a zero displacement matrix of the correct size.
    
    disppot = zeros(N, 3);
    return; 
end

% N_total must cover all particles, even those with no collisions
N_total = max([max(pairs_i), max(pairs_j), N]);
rc = S.rc;

dx = d_mic(:,1);
dy = d_mic(:,2);
dz = d_mic(:,3);
r  = sqrt(dx.^2 + dy.^2 + dz.^2);

pot_r_min = H(1,1);
pot_r_max = H(end,1);
pot_F_min = H(1,2);

r_clamped = min(max(r,pot_r_min), pot_r_max);
Fij = H_interpolant(r_clamped);

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

