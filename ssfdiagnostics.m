% --- Post-Processing for Sanity Check ---
% 1. Normalize the raw sums
% S(k) = < |rho(k)|^2 > / N
N_part = size(p, 1); % Or S.N
if isfield(SSF, 'nsnap') && SSF.nsnap > 0
    norm_fac = 1 / (SSF.nsnap * N_part);
else
    error('No snapshots collected yet!');
end

% Extract Data
k_iso = SSF.kmag_sampling;
S_iso = SSF.sampling.sum_abs2 * norm_fac;

k_lat = SSF.kmag_lattice;
S_lat = SSF.lattice.sum_abs2 * norm_fac;

% --- PLOT 1: Radial S(k) check ---
figure(1); clf; hold on;
% Plot Isotropic Sampling (The continuous line probe)
plot(k_iso, S_iso, 'w.', 'MarkerSize', 1, 'DisplayName', 'Isotropic Sampling');

% Overplot Lattice Modes (The "Allowed" physics)
% We restrict x-axis because lattice modes might stop early due to your cutoff
plot(k_lat, S_lat, 'r', 'LineWidth', 1,'MarkerSize', 1, 'DisplayName', 'Lattice Modes');

% Add reference line at S=1
yline(1, 'b--', 'High-k Limit');

xlabel('k (wavenumber)'); ylabel('S(k)');
title(sprintf('Radial S(k) Check (N=%d, Snap=%d)', N_part, SSF.nsnap));
legend; grid on;

xscale log
yscale log

% --- CRITERIA ---
% PASS: The black dots form a noisy curve. The red circles sit roughly ON TOP of the black dots.
% PASS: At high k (right side), the data oscillates around 1.0.
% FAIL: Data is roughly equal to N (e.g. 1000). (Forgot /N)
% FAIL: Data is extremely small (e.g. 0.001). (Divided by N^2)

% --- PLOT 2 FIX: Auto-detect Peak and Scale ---

% 1. Find the k value where S(k) is maximum (The First Peak)
% We restrict search to the lattice modes or sampling modes that are non-zero
valid_idx = S_iso > 0;
[max_S, max_idx] = max(S_iso(valid_idx)); 
k_iso_valid = k_iso(valid_idx);
k_peak_detected = k_iso_valid(max_idx);

fprintf('Detected Peak S(k)=%.2f at k=%.2e\n', max_S, k_peak_detected);

% 2. Define a "Shell" around this peak
% We use a relative width (e.g., +/- 2% of the peak k)
width = 0.02 * k_peak_detected; 
mask = abs(SSF.kmag_sampling - k_peak_detected) < width;

% Check if we found points
if sum(mask) == 0
    warning('Shell is too thin. Widening...');
    width = 0.10 * k_peak_detected; % Try 10%
    mask = abs(SSF.kmag_sampling - k_peak_detected) < width;
end

k_vecs_shell = SSF.kvecs_sampling(mask, :);
S_vals_shell = S_iso(mask);

% 3. Plot
figure(2); clf;
% Scatter plot: x, y, z, size(40), Color(S_vals), filled circles
scatter3(k_vecs_shell(:,1), k_vecs_shell(:,2), k_vecs_shell(:,3), 60, S_vals_shell, 'filled');
axis equal; 
grid on;
colorbar;
title(sprintf('Anisotropy Check at k \\approx %.2e', k_peak_detected));
xlabel('kx'); ylabel('ky'); zlabel('kz');

% Make axis tight to see the shell
xlim([-1.2*k_peak_detected, 1.2*k_peak_detected]);
ylim([-1.2*k_peak_detected, 1.2*k_peak_detected]);
zlim([-1.2*k_peak_detected, 1.2*k_peak_detected]);