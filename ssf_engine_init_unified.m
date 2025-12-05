function SSF = ssf_engine_init_unified(S, opts)
% SSF = ssf_engine_init_unified_safe(S, opts)
% Safe version that clamps the lattice enumeration to avoid 1e48 points.

if nargin<2, opts = struct(); end

% --- 1. DEFINE LIMITS ---

% A. The Thermal Limit (Highest resolution physically possible)
% This is used for the isotropic sampling to check for aliasing/noise
if isfield(S, 'stdx') && S.stdx > 1e-9 
    thermal_res = sqrt(3 * S.stdx^2);
    SSF.kmax_thermal = pi / thermal_res;
else
    % Fallback if stdx is zero or effectively zero
    fprintf('Warning: S.stdx is too small or missing. Defaulting to particle size.\n');
    SSF.kmax_thermal = 20 * (2*pi/S.rp); 
end

% B. The Lattice Limit (Computable limit)
% We restrict the EXPLICIT lattice summation to a multiple of the fundamental
% frequency to keep the grid size under ~1 million points.
% Standard SSF usually decays to 1 by k*sigma = 20.
if ~isfield(opts, 'n_lattice_shells'), opts.n_lattice_shells = 30; end

% --- 2. RECIPROCAL BASIS ---
% Formula: B = 2*pi * (A^-1)^T
B = 2 * pi * S.fcc.invA'; 
recip_vec_lengths = sqrt(sum(B.^2, 2));
SSF.k_fundamental = min(recip_vec_lengths);

% Define the Lattice Cutoff based on fundamental (Box) size
SSF.kmax_lattice = opts.n_lattice_shells * SSF.k_fundamental;


% --- 3. GENERATE LATTICE MODES (The "Allowed" Physics) ---
% We scan integers ONLY up to kmax_lattice (manageable size)

min_b = min(recip_vec_lengths);
maxn = ceil(SSF.kmax_lattice / min_b) + 1;

% Safety Clamp: Ensure we don't accidentally blow up memory
if maxn > 60
    fprintf('SSF Warning: Lattice index %d is too high. Clamping to 60.\n', maxn);
    maxn = 60;
end

% Check estimated size
est_points = (2*maxn)^3;
fprintf('SSF Init: Scanning lattice grid +/- %d (approx %d points)... ', maxn, est_points);

[nx, ny, nz] = ndgrid(-maxn:maxn, -maxn:maxn, -maxn:maxn);
n_triplets = [nx(:), ny(:), nz(:)];
n_triplets(all(n_triplets==0,2),:) = []; % Remove k=0

% Map to K-space
k_candidates = (B * n_triplets')'; 

% Filter
kmags = sqrt(sum(k_candidates.^2, 2));
mask = kmags <= SSF.kmax_lattice & kmags > 1e-6;

% Store and Sort
SSF.kvecs_lattice = k_candidates(mask, :);
SSF.kmag_lattice  = kmags(mask);
[SSF.kmag_lattice, sort_idx] = sort(SSF.kmag_lattice);
SSF.kvecs_lattice = SSF.kvecs_lattice(sort_idx, :);

fprintf('Kept %d modes.\n', numel(SSF.kmag_lattice));


% --- 4. GENERATE ISOTROPIC SAMPLING (The "Comparison" Probe) ---
% Here we go ALL THE WAY to the thermal limit to check high-k fluctuations.
% This does not explode because we choose a fixed number of shells.

if ~isfield(opts, 'num_shells'), opts.num_shells = 150; end
if ~isfield(opts, 'dirs_per_shell'), opts.dirs_per_shell = 120; end

% Create shells from k_fundamental out to kmax_thermal
k_rad = linspace(SSF.k_fundamental, SSF.kmax_thermal, opts.num_shells)';
num_dirs = opts.dirs_per_shell;
total_pts = numel(k_rad) * num_dirs;

k_iso = zeros(total_pts, 3);
km_iso = zeros(total_pts, 1);
idx = 1;

for i = 1:numel(k_rad)
    k_val = k_rad(i);
    dirs = fibonacci_sphere(num_dirs); 
    for j = 1:num_dirs
        k_iso(idx,:) = k_val * dirs(j,:);
        km_iso(idx)  = k_val;
        idx = idx + 1;
    end
end

SSF.kvecs_sampling = k_iso;
SSF.kmag_sampling  = km_iso;


% --- 5. ALLOCATE STORAGE ---
alloc = @(M) struct('sum_rho', zeros(1,M), ...     
                    'sum_abs2', zeros(1,M), ...   
                    'rho_ts', []);                 

SSF.lattice = alloc(size(SSF.kvecs_lattice, 1));
SSF.sampling = alloc(size(SSF.kvecs_sampling, 1));

if isfield(opts,'store_rho_ts') && opts.store_rho_ts
    SSF.lattice.rho_ts = complex(zeros(0, size(SSF.kvecs_lattice, 1)));
    SSF.sampling.rho_ts = complex(zeros(0, size(SSF.kvecs_sampling, 1)));
end

SSF.nsnap = 0;

end

% --- HELPER ---
function dirs = fibonacci_sphere(N)
    i = (0:N-1)';
    phi = acos(1 - 2*(i+0.5)/N);
    theta = pi*(1+sqrt(5))*(i+0.5);
    dirs = [sin(phi).*cos(theta), sin(phi).*sin(theta), cos(phi)];
end