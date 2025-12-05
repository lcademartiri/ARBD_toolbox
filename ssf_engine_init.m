function SSF = ssf_engine_init(S, opts)
% SSF = ssf_engine_init(S, opts)
% S: simulation struct (needs S.bc, S.fcc.A, S.fcc.invA, S.br, S.rp)
% opts: struct with fields (all optional)
%   .num_radial = 100;      % radial shells (for isotropic sampling)
%   .num_dirs_per_r = 64;   % directions per radial shell (for isotropic sampling)
%   .kmax_factor = 2*pi / S.rp;  % default k_max
%   .kmin_factor = 2*pi / (2*S.br); % default fundamental
%   .use_spherical_modes = false;  % compute j_l roots (slow)
%   .lmax = 8; .nperl = 10;        % spherical mode parameters
%   .store_rho_ts = false;        % store per-snapshot complex rho array (big)
%
% Returns SSF struct with fields:
%   .kvecs_pbc_cubic, .kvecs_pbc_fcc, .kvecs_isotropic, .kvecs_sphere (optional)
%   storage arrays: sum_rho_*, sum_abs2_*, nsnap, optionally rho_ts_* (if store)

if nargin<2, opts = struct(); end
% defaults
if ~isfield(opts,'num_radial'), opts.num_radial = 120; end
if ~isfield(opts,'num_dirs_per_r'), opts.num_dirs_per_r = 48; end
if ~isfield(opts,'kmax_factor'), opts.kmax_factor = 2*pi / sqrt(3*S.stdx^2); end % Nyquist limit as standard
if ~isfield(opts,'kmin_factor'), opts.kmin_factor = 2*pi / (2*S.br); end
if ~isfield(opts,'use_spherical_modes'), opts.use_spherical_modes = false; end
if ~isfield(opts,'lmax'), opts.lmax = 8; end
if ~isfield(opts,'nperl'), opts.nperl = 8; end
if ~isfield(opts,'store_rho_ts'), opts.store_rho_ts = false; end

SSF.opts = opts;
SSF.kmin = opts.kmin_factor;
SSF.kmax = opts.kmax_factor;

% ---------- PBC cubic lattice nodes ----------
% generate integer triplets (nx,ny,nz) up to num_radial steps until |k|<=kmax
maxn = ceil( (SSF.kmax/SSF.kmin) ) + 2;
[nx,ny,nz] = ndgrid(-maxn:maxn, -maxn:maxn, -maxn:maxn);
trip = [nx(:), ny(:), nz(:)];
trip(all(trip==0,2),:) = []; % remove 0 vector
k_candidates = SSF.kmin * trip;  % Cartesian cubic reciprocal vectors
mags = sqrt(sum(k_candidates.^2,2));
mask = mags <= SSF.kmax;
k_candidates = k_candidates(mask,:);
mags = mags(mask);
% sort by magnitude
[~,ord] = sort(mags);
SSF.kvecs_pbc_cubic = k_candidates(ord,:);
SSF.kmag_pbc_cubic = sqrt(sum(SSF.kvecs_pbc_cubic.^2,2));

% ---------- PBC FCC recip lattice (if S.bc==3 or always build) ----------
% if S.fcc.A defined use reciprocal basis B = 2*pi * inv(A)^T
B_matrix = 2*pi * S.fcc.invA';
% generate integer multipliers similar to above
trip_all = trip; % reuse
k_cand_fcc = (B_matrix * trip_all')';
mags_f = sqrt(sum(k_cand_fcc.^2,2));
maskf = mags_f <= SSF.kmax & mags_f>0;
k_cand_fcc = k_cand_fcc(maskf,:);
mags_f = mags_f(maskf);
[~,ordf] = sort(mags_f);
SSF.kvecs_pbc_fcc = k_cand_fcc(ordf,:);
SSF.kmag_pbc_fcc = sqrt(sum(SSF.kvecs_pbc_fcc.^2,2));

% ---------- Isotropic plane-wave sampling ----------
% radial grid
k_rad = linspace(SSF.kmin, SSF.kmax, opts.num_radial)';
num_dirs = opts.num_dirs_per_r;
k_isotropic = zeros(numel(k_rad)*num_dirs,3);
kmags_iso = zeros(numel(k_rad)*num_dirs,1);
idx = 1;
% use quasi-uniform directions: Fibonacci or random
for ii=1:numel(k_rad)
    k = k_rad(ii);
    % generate unit directions (random seeds deterministic could be used)
    dirs = fibonacci_sphere(num_dirs); % returns num_dirs x 3
    for jj = 1:num_dirs
        k_isotropic(idx,:) = k * dirs(jj,:);
        kmags_iso(idx) = k;
        idx = idx + 1;
    end
end
% trim unused
k_isotropic(idx:end,:) = [];
kmags_iso(idx:end) = [];
SSF.kvecs_isotropic = k_isotropic;
SSF.kmag_isotropic = kmags_iso;

% ---------- Optional: spherical eigenmodes (exact SBC basis) ----------
if opts.use_spherical_modes
    R = S.br; % sphere radius
    [k_sph, mode_info] = compute_spherical_bessel_modes(R, opts.lmax, opts.nperl);
    % k_sph: column vector of k values, mode_info parallel struct with l and n
    % build kvecs_sphere as a set of (k * Y_lm directions)? For spherical modes we store (k,l,m)
    SSF.k_sphere = k_sph;
    SSF.mode_info = mode_info;
else
    SSF.k_sphere = [];
    SSF.mode_info = [];
end

% ---------- Storage init ----------
SSF.nsnap = 0;

% choose which families to collect
SSF.collect = {'pbc_cubic','pbc_fcc','isotropic'}; 
if opts.use_spherical_modes, SSF.collect{end+1} = 'sphere'; end

% allocate running sums and optional time series
% helper function to allocate sums for M kvectors
alloc = @(M) struct('sum_rho', zeros(1,M), 'sum_abs2', zeros(1,M), 'count', 0, ...
                    'rho_ts', [] );

% cubic
M = size(SSF.kvecs_pbc_cubic,1);
SSF.pbc_cubic = alloc(M);
if opts.store_rho_ts, SSF.pbc_cubic.rho_ts = complex(zeros(0,M)); end

% fcc
M = size(SSF.kvecs_pbc_fcc,1);
SSF.pbc_fcc = alloc(M);
if opts.store_rho_ts, SSF.pbc_fcc.rho_ts = complex(zeros(0,M)); end

% isotropic
M = size(SSF.kvecs_isotropic,1);
SSF.isotropic = alloc(M);
if opts.store_rho_ts, SSF.isotropic.rho_ts = complex(zeros(0,M)); end

% sphere modes (not plane-wave rho but a_lm coefficients if computed)
if opts.use_spherical_modes
    M = numel(SSF.k_sphere);
    SSF.sphere = alloc(M);
    if opts.store_rho_ts, SSF.sphere.rho_ts = complex(zeros(0,M)); end
end

end


%%%% helper: fibonacci_sphere (returns num x 3 unit vectors)
function dirs = fibonacci_sphere(N)
    i = (0:N-1)';
    phi = acos(1 - 2*(i+0.5)/N);
    theta = pi*(1+sqrt(5))*(i+0.5);
    dirs = [sin(phi).*cos(theta), sin(phi).*sin(theta), cos(phi)];
end

%%%% helper: compute spherical bessel modes (slow) - optional
function [kvals, info] = compute_spherical_bessel_modes(R, lmax, nperl)
% Computes roots of derivative j_l'(kR)=0 for l=0..lmax, first nperl roots
% Returns flattened list kvals and info struct with l,n
    kvals = [];
    info = [];
    for l=0:lmax
        for n=1:nperl
            % initial guess: approx (n + l/2)*pi/R
            x0 = (n + l/2)*pi / R;
            f = @(x) spherical_besselj_derivative(l, x*R);
            % use fzero with bracket
            try
                root = fzero(f, x0);
            catch
                % try search bracket
                a = max(1e-6, x0-0.5*pi/R); b = x0+0.5*pi/R;
                root = fzero(f, [a,b]);
            end
            kvals(end+1,1) = root; %#ok<AGROW>
            info(end+1).l = l; info(end).n = n; %#ok<AGROW>
        end
    end
end

function dj = spherical_besselj_derivative(l, x)
% derivative wrt x of j_l(x). Using recurrence: j_l' = j_{l-1}/2 - j_{l+1}/2 ? use relation
    jl = sqrt(pi./(2*x)).*besselj(l+0.5, x);
    % numerical derivative (robust for small list)
    h = 1e-6 + 1e-6*abs(x);
    jlph = sqrt(pi./(2*(x+h))).*besselj(l+0.5, x+h);
    jlmh = sqrt(pi./(2*(x-h))).*besselj(l+0.5, x-h);
    dj = (jlph - jlmh)/(2*h);
end
