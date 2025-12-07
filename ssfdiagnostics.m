clear all
if exist('D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
elseif exist('G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
elseif exist('D:\GDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
end
if exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
elseif exist('C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
elseif exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
end
    
toolbox_folder = '..\ARBD_toolbox';
addpath(data_folder)
addpath(toolbox_folder)
addpath(output_folder)

load('SBCvsPBC_temp_29.mat');

%% --- Post-Processing for Sanity Check ---
% 1. Normalize the raw sums
% S(k) = < |rho(k)|^2 > / N
N_part = S.N; % Or S.N
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
figure; clf; hold on;
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
figure; clf;
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

% --- POST-PROCESSING: REMOVING GEOMETRY ---

% Normalization Factor
norm_fac = 1 / (S.N * SSF.nsnap); % Note: sum_rho needs /nsnap then | |^2

% === 1. PROCESS SBC DATA ===
% We use the Isotropic Sampling (White Dots) to show the continuum
k_sbc = SSF.kmag_sampling;
% Total (The Raw Data)
S_sbc_total = SSF.sampling.sum_abs2 * norm_fac;
% Static (The Sphere Form Factor)
S_sbc_static = (abs(SSF.sampling.sum_rho / SSF.nsnap).^2) / S.N;
% Fluctuation (The Pure Liquid)
S_sbc_fluct = S_sbc_total - S_sbc_static;


% === 2. PROCESS PBC DATA (Cubic or FCC) ===
% We use the Lattice Modes (Red Dots) because that's where the physics lives
k_pbc = SSF.kmag_lattice;
% Total
S_pbc_total = SSF.lattice.sum_abs2 * norm_fac;
% Static (The Artificial Crystal)
S_pbc_static = (abs(SSF.lattice.sum_rho / SSF.nsnap).^2) / S.N;
% Fluctuation
S_pbc_fluct = S_pbc_total - S_pbc_static;


% === 3. PLOTTING THE RESULT ===
figure; clf;

% -- Plot A: Decomposition --
subplot(1,2,1); hold on;
loglog(k_sbc, S_sbc_total, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Total (Raw)');
loglog(k_sbc, S_sbc_static, 'b-', 'LineWidth', 1, 'DisplayName', 'Static (Form Factor)');
loglog(k_sbc, S_sbc_fluct, 'r.', 'MarkerSize', 8, 'DisplayName', 'Pure Fluctuations');
xlabel('k'); ylabel('S(k)'); title('Removing the Sphere');
legend; axis tight;
% Interpretation: The Blue line will trace the huge oscillations. 
% The Red dots should be relatively flat (Isotropic Liquid).

% -- Plot B: Decomposition --
subplot(1,2,2); hold on;
stem(k_pbc, S_pbc_total, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Total (Raw)');
stem(k_pbc, S_pbc_static, 'b', 'LineWidth', 2, 'DisplayName', 'Static (Bragg Noise)');
stem(k_pbc, S_pbc_fluct, 'r', 'LineWidth', 2, 'DisplayName', 'Pure Fluctuations');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('k'); ylabel('S(k)'); title('Isolating Bragg Noise');
legend; axis tight;

hAxes = findobj(gcf,"Type","axes")
hAxes(1).YLim = [0.0100,96.7487]
hAxes(2).YLim = [0.0100,100]
hAxes(2).XScale = "log"
hAxes(2).YScale = "log"
StaticBraggNoise = findobj(gcf,"DisplayName","Static (Bragg Noise)")
StaticBraggNoise.LineStyle = "none"
StaticBraggNoise.MarkerSize = 2
StaticBraggNoise.MarkerEdgeColor = "none"
StaticBraggNoise.MarkerFaceColor = [0.2314,0.6667,0.1961]
TotalRaw = findobj(gcf,"DisplayName","Total (Raw)")
TotalRaw(1).LineStyle = "none"
PureFluctuations = findobj(gcf,"DisplayName","Pure Fluctuations")
PureFluctuations(1).LineStyle = "none"
PureFluctuations(1).MarkerSize = 2
PureFluctuations(1).MarkerFaceColor = [0.9608,0.4667,0.1608]
PureFluctuations(1).MarkerEdgeColor = "none"
hText = findall(gcf,"Type","text")
hText(1).String = "Isolating Bragg Noise"
PureFluctuations(2).Marker = "o"
PureFluctuations(2).MarkerSize = 2
PureFluctuations(2).MarkerFaceColor = [0.9608,0.4706,0.1608]
PureFluctuations(2).MarkerEdgeColor = "none"
TotalRaw(2).LineWidth = 1
StaticFormFactor = findobj(gcf,"DisplayName","Static (Form Factor)")
StaticFormFactor.Color = [0.2314,0.6667,0.1961]
% Interpretation: The Blue stems show how much the liquid is "Frozen".
% The Red stems show the actual liquid motion.