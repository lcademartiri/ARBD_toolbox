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

load('SBCvsPBC_25.mat','SSF');

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

% --- PLOT 4: Normalized Anisotropy at Structural Peak ---

% 1. Compute the Median Curve (Filter PBC spikes / SBC ringing)
% This collapses the "cloud" to a single representative curve
unique_k = unique(SSF.kmag_sampling);
median_S = zeros(size(unique_k));

for i = 1:numel(unique_k)
    k_val = unique_k(i);
    % Find all samples at this exact shell (within numerical tolerance)
    mask = abs(SSF.kmag_sampling - k_val) < 1e-9;
    median_S(i) = median(S_iso(mask));
end

% 2. Define Physical Search Window
% Theory: First peak is near k_theo = pi / Radius.
% We look for the maximum between 0.5 * k_theo and 2.0 * k_theo.
k_theo = pi / S.rp;
k_search_min = 0.5 * k_theo;
k_search_max = 2.0 * k_theo;

valid_search_mask = (unique_k >= k_search_min) & (unique_k <= k_search_max);

% 3. Find the Peak within the Window
if sum(valid_search_mask) == 0
    warning('No sampling shells found in the expected physical range!');
    [max_S_val, idx_in_unique] = max(median_S); % Fallback to global max
    k_peak_found = unique_k(idx_in_unique);
else
    S_subset = median_S(valid_search_mask);
    k_subset = unique_k(valid_search_mask);
    [max_S_val, idx_subset] = max(S_subset);
    k_peak_found = k_subset(idx_subset);
end

% 4. Extract Data for the Peak Shell
mask_peak = abs(SSF.kmag_sampling - k_peak_found) < 1e-9;
k_vecs_peak = SSF.kvecs_sampling(mask_peak, :);
S_raw_peak  = S_iso(mask_peak);

% 5. NORMALIZE (Divide by Median)
% A(k) = S(k) / <S(k)>
baseline_peak = median(S_raw_peak);
S_norm_peak = S_raw_peak / baseline_peak;

% Calculate Statistics for Title/Console
min_A = min(S_norm_peak);
max_A = max(S_norm_peak);
anisotropy_factor = max_A / min_A;

fprintf('\n--- Peak Anisotropy Analysis ---\n');
fprintf('Peak found at k = %.2e (S_raw ~ %.2f)\n', k_peak_found, baseline_peak);
fprintf('Anisotropy Range: %.2f to %.2f (Factor: %.2fx)\n', min_A, max_A, anisotropy_factor);

% 6. Plot
figure(4); clf;
scatter3(k_vecs_peak(:,1), k_vecs_peak(:,2), k_vecs_peak(:,3), 60, S_norm_peak, 'filled');
axis equal; grid on; colorbar;
xlabel('kx'); ylabel('ky'); zlabel('kz');

% Tight axis
xlim([-1.2*k_peak_found, 1.2*k_peak_found]);
ylim([-1.2*k_peak_found, 1.2*k_peak_found]);
zlim([-1.2*k_peak_found, 1.2*k_peak_found]);

% --- SMART COLOR SCALING ---

% Option 1: TIGHT FIT (Shows max contrast)
% This stretches the color map to fit exactly between min and max
% caxis([min(S_norm_peak), max(S_norm_peak)]);

% Option 2: SYMMETRIC FIT (Best for "Deviation from 1.0")
% This centers 1.0 as Green/Teal, and shows deviations equally
max_dev = max(abs(S_norm_peak - 1.0));
clim([1.0 - max_dev, 1.0 + max_dev]);

title({sprintf('Normalized Peak Anisotropy (k=%.2e)', k_peak_found), ...
       sprintf('Range: [%.2f, %.2f] (Ref=1.0)', min_A, max_A)});

% --- PLOT 3: Normalized Low-k Anisotropy ---

% 1. Target the Low-k Region (The "Force Field" Artifact)
% We pick a shell near the fundamental box mode or just a low value
k_target_low = SSF.k_fundamental; 

% 2. Find the Single Closest Shell (Smart Search)
unique_shells = unique(SSF.kmag_sampling);
[~, idx] = min(abs(unique_shells - k_target_low));
k_shell_low = unique_shells(idx);

% 3. Extract Data
mask = abs(SSF.kmag_sampling - k_shell_low) < 1e-9;
k_vecs = SSF.kvecs_sampling(mask, :);
S_raw  = S_iso(mask);

% 4. NORMALIZE (The "Subtraction" Equivalent)
% We divide by the median of THIS specific shell.
% This removes the global scaling of the Form Factor/Ringing.
baseline = median(S_raw);
S_norm = S_raw / baseline;

fprintf('Low-k Anisotropy Check at k=%.2e\n', k_shell_low);
fprintf('Baseline Signal: %.2e\n', baseline);
fprintf('Max Relative Deviation: %.2f x\n', max(S_norm));

% 5. Plot
figure(3); clf;
scatter3(k_vecs(:,1), k_vecs(:,2), k_vecs(:,3), 60, S_norm, 'filled');
axis equal; grid on; colorbar;
xlabel('kx'); ylabel('ky'); zlabel('kz');

% --- SMART COLOR SCALING ---

% Option 1: TIGHT FIT (Shows max contrast)
% This stretches the color map to fit exactly between min and max
% caxis([min(S_norm_peak), max(S_norm_peak)]);

% Option 2: SYMMETRIC FIT (Best for "Deviation from 1.0")
% This centers 1.0 as Green/Teal, and shows deviations equally
max_dev = max(abs(S_norm_peak - 1.0));
clim([1.0 - max_dev, 1.0 + max_dev]);
title(sprintf('Relative Anisotropy A(k) at k=%.2e', k_shell_low));

% Add annotation
text(min(xlim), min(ylim), min(zlim), ...
    sprintf('Mean: 1.0\nMax: %.2f', max(S_norm)), ...
    'BackgroundColor', 'white');

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

%%


trajectories=POS;
clear POS
% Results = batch_process_dynamics(S, trajectories, SSF)
%
% Memory-safe computation of F(k,t), Gamma, and Deff.
% 1. Holds Trajectories in RAM (approx 2.4 GB).
% 2. Processes K-vectors in chunks to prevent Rho matrix explosion.
% 3. Returns combined results for all K.

    % --- CONFIGURATION ---
    CHUNK_SIZE = 4000; % Processes 4000 k-vectors at a time (~6.4 GB RAM for Rho)
    dt_sim = 1.0;      % UPDATE THIS with your real time between snapshots
    
    % --- SETUP ---
    % Combine Lattice and Sampling vectors for one big pass (we split them later)
    all_kvecs = [SSF.kvecs_lattice; SSF.kvecs_sampling];
    num_lattice = size(SSF.kvecs_lattice, 1);
    num_total_k = size(all_kvecs, 1);
    
    % Check Trajectory Format
    if iscell(trajectories)
        n_snaps = numel(trajectories);
        % Convert cell to matrix if possible for speed (N_part x 3 x N_snaps)
        % (Only if you have RAM, otherwise stick to cell)
        fprintf('Converting cell array to matrix for speed...\n');
        P_all = cat(3, trajectories{:}); 
    else
        P_all = trajectories;
        n_snaps = size(P_all, 3);
    end
    [N_part, ~, ~] = size(P_all);
    
    fprintf('Starting Batch Analysis on %d K-vectors over %d snapshots.\n', num_total_k, n_snaps);
    fprintf('Trajectory Memory: %.2f GB\n', whos('P_all').bytes / 1024^3);
    
    % --- DEFINE LAGS ---
    % Logarithmic spacing for dynamics
    max_lag = floor(n_snaps / 10);
    lags = unique(round(logspace(0, log10(max_lag), 50)));
    t_vals = lags * dt_sim;
    
    % --- PRE-ALLOCATE RESULTS ---
    % We only store the summary stats, not the full time series
    Results.k_mags = sqrt(sum(all_kvecs.^2, 2));
    Results.S_k    = zeros(num_total_k, 1);
    Results.Gamma  = zeros(num_total_k, 1);
    Results.D_eff  = zeros(num_total_k, 1);
    Results.lags   = lags;
    
    % We CAN store F(k,t) because [NumK x 50 lags] is small (~90 MB)
    Results.F_kt   = zeros(numel(lags), num_total_k);

    % --- BATCH LOOP ---
    num_chunks = ceil(num_total_k / CHUNK_SIZE);
    
    total_timer = tic;
    
    for c = 1:num_chunks
        % Define range
        idx_start = (c-1)*CHUNK_SIZE + 1;
        idx_end   = min(c*CHUNK_SIZE, num_total_k);
        current_k_count = idx_end - idx_start + 1;
        
        k_chunk = all_kvecs(idx_start:idx_end, :);
        
        fprintf('Processing Chunk %d/%d (Vectors %d-%d)... ', c, num_chunks, idx_start, idx_end);
        
        % 1. EXTRACT RHO FOR THIS CHUNK
        % We need [Time x K] matrix
        Rho_Chunk = complex(zeros(n_snaps, current_k_count));
        
        % This inner loop is unavoidable, but P_all is in RAM so it's fast
        for t = 1:n_snaps
            p_t = P_all(:,:,t); % N_part x 3
            % Vectorized: sum(exp(-i * k * r))
            % (K_chunk x 3) * (3 x N_part) -> (K_chunk x N_part)
            phase = -1i * (k_chunk * p_t');
            Rho_Chunk(t, :) = sum(exp(phase), 2).'; 
        end
        
        % 2. CALCULATE DYNAMICS FOR THIS CHUNK
        % S(k)
        S_k_chunk = mean(abs(Rho_Chunk).^2, 1);
        Results.S_k(idx_start:idx_end) = S_k_chunk(:);
        
        % F(k,t) Autocorrelation
        for i = 1:numel(lags)
            tau = lags(i);
            % Vectorized correlation over time
            series_t0 = Rho_Chunk(1:end-tau, :);
            series_t1 = Rho_Chunk(1+tau:end, :);
            
            % Mean over time
            corr = mean(series_t0 .* conj(series_t1), 1);
            
            % Normalize
            F_val = real(corr) ./ S_k_chunk;
            Results.F_kt(i, idx_start:idx_end) = F_val;
        end
        
        % 3. FIT GAMMA AND D_EFF
        for k_local = 1:current_k_count
            k_global_idx = idx_start + k_local - 1;
            
            y = Results.F_kt(:, k_global_idx);
            k_mag = Results.k_mags(k_global_idx);
            
            % Fit initial decay (F > 0.3)
            valid_idx = y > 0.3;
            if sum(valid_idx) > 3
                p = polyfit(t_vals(valid_idx), log(y(valid_idx)), 1);
                gamma = -p(1);
            else
                gamma = NaN;
            end
            
            Results.Gamma(k_global_idx) = gamma;
            if ~isnan(gamma) && k_mag > 1e-6
                Results.D_eff(k_global_idx) = gamma / (k_mag^2);
            end
        end
        
        fprintf('Done. (%.1fs)\n', toc(total_timer));
        total_timer = tic; % Reset timer for next chunk
    end
    
    % --- SPLIT RESULTS BACK INTO LATTICE vs SAMPLING ---
    Results.Lattice.k     = Results.k_mags(1:num_lattice);
    Results.Lattice.S_k   = Results.S_k(1:num_lattice);
    Results.Lattice.D_eff = Results.D_eff(1:num_lattice);
    Results.Lattice.F_kt  = Results.F_kt(:, 1:num_lattice);
    
    Results.Sampling.k     = Results.k_mags(num_lattice+1:end);
    Results.Sampling.S_k   = Results.S_k(num_lattice+1:end);
    Results.Sampling.D_eff = Results.D_eff(num_lattice+1:end);
    Results.Sampling.F_kt  = Results.F_kt(:, num_lattice+1:end);
    
    fprintf('Batch Analysis Complete.\n');