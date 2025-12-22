

%% --- SSF ANALYSIS & VISUALIZATION ----------------------------------

% --- 0. SETUP & BASICS ---
fprintf('Starting SSF Analysis...\n');
if ~isfield(SSF, 'nsnap') || SSF.nsnap < 1
    error('SSF Error: No snapshots found in structure.');
end

N_part = S.N; 
norm_fac = 1 / (SSF.nsnap * N_part);

% Define High-Contrast Colors for Black Background
c_raw   = [0.9 0.9 0.9];       % White/Gray (Raw Data)
c_stat  = [0.2 1.0 0.2];       % Neon Green (Static/Bragg/FormFactor)
c_fluct = [1.0 0.2 1.0];       % Magenta (Pure Fluctuations)
c_ref   = [0.4 0.8 1.0 0.5];   % Cyan transparent (Reference)

% --- 1. PROCESS CARTESIAN PROBES (Common to all) ---
% A. "Random Directions" (Off-Lattice / Continuum)
% In SBC: Used for subtraction method. In PBC: Shows spectral leakage.
k_rand = SSF.kmag_sampling;
S_rand_tot  = SSF.sampling.sum_abs2 * norm_fac;
S_rand_stat = (abs(SSF.sampling.sum_rho / SSF.nsnap).^2) / N_part;
S_rand_fluc = S_rand_tot - S_rand_stat;

% B. "Reciprocal Lattice" (Cubic/FCC Grid)
% In SBC: Meaningless. In PBC: The correct Eigenmodes.
k_lat = SSF.kmag_lattice;
S_lat_tot  = SSF.lattice.sum_abs2 * norm_fac;
S_lat_stat = (abs(SSF.lattice.sum_rho / SSF.nsnap).^2) / N_part;
S_lat_fluc = S_lat_tot - S_lat_stat;

%% --- 2. BRANCH: SPHERICAL BOUNDARY CONDITIONS (SBC) ---
if S.bc == 1
    fprintf('Detected SBC. Computing Bessel Eigenmodes...\n');
    
    % --- BESSEL NORMALIZATION (Neumann + Volume + Degeneracy) ---
    R = S.br;
    V_sphere = (4/3)*pi*R^3;
    modes = SSF.sbc_modes.data; % [k, l, n]
    k_bes = modes(:,1);
    l_vals = modes(:,2);
    
    % Calculate Geometric Weights
    W_ln = zeros(size(k_bes));
    for i = 1:numel(k_bes)
        l = l_vals(i);
        x = k_bes(i) * R;
        if x < 1e-4
            W_ln(i) = 1.0;
        else
            val = sqrt(pi/(2*x)) * besselj(l + 0.5, x);
            % Neumann Factor: (1 - l(l+1)/x^2)
            geo_fac = (1 - (l*(l+1))/(x^2));
            W_ln(i) = ((R^3/2) * val^2 * geo_fac) / V_sphere;
        end
    end
    
    % Apply Normalizations
    raw_bes = SSF.sbc_modes.sum_power / SSF.nsnap;
    S_bes = (raw_bes ./ W_ln) / N_part;  % Geometric
    S_bes = S_bes ./ (2*l_vals + 1);     % Degeneracy Correction!

    % Separate Bulk (l>0) from Wall Layering (l=0)
    mask_bulk = l_vals > 0 & k_bes > 1e-4;
    mask_wall = l_vals == 0 & k_bes > 1e-4;

    % --- PLOT SBC ---
    figure('Color','k', 'Position', [50 50 1400 600]);
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % PANEL 1: CARTESIAN SUBTRACTION (Method A)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Raw
    plot(k_rand, S_rand_tot, 'o', 'Color', c_raw, 'MarkerSize', 2, ...
        'DisplayName', 'Raw (Includes Sphere Form Factor)');
    % Static Background
    plot(k_rand, S_rand_stat, '-', 'Color', c_stat, 'LineWidth', 1.5, ...
        'DisplayName', 'Static Profile (Form Factor)');
    % Result
    plot(k_rand, S_rand_fluc, '.', 'Color', c_fluct, 'MarkerSize', 10, ...
        'DisplayName', 'Cleaned S(k) (Subtraction)');
    
    title('Method A: Cartesian Basis (Subtraction)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w');
    ylim([0.01 10]); xlim([0 max(k_rand)]); set(gca, 'YScale', 'log');

    % PANEL 2: BESSEL EIGENMODES (Method B)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Bulk Modes (Rainbow)
    scatter(k_bes(mask_bulk), S_bes(mask_bulk), 30, l_vals(mask_bulk), 'filled', ...
        'DisplayName', 'Bulk Modes (l > 0)');
    colormap(gca, 'turbo'); 
    cb = colorbar; cb.Label.String = 'Angular Momentum (l)'; 
    cb.Color = 'w'; cb.Label.Color = 'w';
    
    % Wall Modes (White Outliers)
    plot(k_bes(mask_wall), S_bes(mask_wall), 'wo', 'MarkerSize', 6, 'LineWidth', 1.5, ...
        'DisplayName', 'Wall Layering (l = 0)');
    
    % Overlay Method A Result for comparison
    plot(k_rand, smooth(S_rand_fluc, 5), '-', 'Color', c_ref, 'LineWidth', 2, ...
        'DisplayName', 'Ref: Cartesian Cleaned');
    
    title('Method B: Bessel Basis (Natural Separation)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    yline(1, 'w--');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w', 'Location', 'Southeast');
    ylim([0.01 10]); xlim([0 max(k_bes)]); set(gca, 'YScale', 'log');
    
    sgtitle('SBC Analysis: Recovering Bulk S(k) in Confinement', 'Color', 'w', 'FontSize', 14);


%% --- 3. BRANCH: PERIODIC BOUNDARY CONDITIONS (PBC/FCC) ---
else
    fprintf('Detected PBC/FCC. Comparing Lattice vs Random Probes...\n');
    
    figure('Color','k', 'Position', [50 50 1400 600]);
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % PANEL 1: OFF-LATTICE (Showing Leakage)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Plot Random Vectors (Off-lattice)
    % These will look messy because they hit the sinc-function zeros of the Box
    plot(k_rand, S_rand_tot, 'o', 'Color', c_raw, 'MarkerSize', 3, ...
        'DisplayName', 'Off-Lattice (Random Directions)');
    
    % Overlay Lattice Points (To show where the peaks SHOULD be)
    plot(k_lat, S_lat_tot, 'r.', 'MarkerSize', 5, 'DisplayName', 'Lattice Points (Ref)');
    
    title('Probe A: Off-Lattice Vectors (Spectral Leakage)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w');
    ylim([0.01 100]); xlim([0 max(k_rand)]); set(gca, 'YScale', 'log', 'XScale', 'log');
    
    % PANEL 2: RECIPROCAL LATTICE (Correct Eigenmodes)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Use Stems for discrete lattice modes
    % Raw Total
    stem(k_lat, S_lat_tot, 'Color', [0.5 0.5 0.5], 'Marker', 'none', 'LineWidth', 1, ...
        'DisplayName', 'Total (Raw)');
    
    % Static (Bragg Peaks - Crystal Order)
    if max(S_lat_stat) > 0.1
        stem(k_lat, S_lat_stat, 'Color', c_stat, 'Marker', 's', 'MarkerSize', 4, ...
            'LineStyle', 'none', 'DisplayName', 'Static (Bragg Order)');
    end
    
    % Fluctuations (Liquid Physics)
    stem(k_lat, S_lat_fluc, 'Color', c_fluct, 'Marker', 'o', 'MarkerSize', 4, ...
        'LineWidth', 1.5, 'LineStyle', 'none', 'DisplayName', 'Pure Fluctuations');
        
    title('Probe B: Reciprocal Lattice (Eigenmodes)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    yline(1, 'w--');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w');
    ylim([0.01 100]); xlim([0 max(k_lat)]); set(gca, 'YScale', 'log', 'XScale', 'log');
    
    sgtitle('PBC Analysis: Importance of Lattice Commensurability', 'Color', 'w', 'FontSize', 14);
end