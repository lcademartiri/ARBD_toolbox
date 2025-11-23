function DCOMP = densityEngine(mode, DCOMP, S, P, p, data_folder)
%DENSITYENGINE Manages setup and accumulation for density analysis.
%
% This function operates in two modes: 'setup' and 'accumulate'.
% The 'setup' mode is computationally intensive and saves its results to a
% file to avoid re-running.
%
% --- SYNTAX ---
%
% 1. To initialize the data structures before the simulation starts:
%    DCOMP = densityEngine('setup', [], S, P, []);
%
% 2. To accumulate data at each timestep inside the simulation loop:
%    DCOMP = densityEngine('accumulate', DCOMP, S, P, p);
%
% --- INPUTS ---
%
%   mode  - String: 'setup' to initialize, or 'accumulate' to add data.
%   DCOMP - Struct: Pass an empty array [] for 'setup'. For 'accumulate',
%           pass the DCOMP struct from the previous timestep.
%   S     - Struct: Contains simulation settings (S.bc, S.br, S.rp, etc.).
%   P     - Struct: Contains simulation parameters (P.dens flag).
%   p     - Matrix: Particle data (positions). Required for 'accumulate'
%           mode only. Pass [] for 'setup'.
%
% --- OUTPUT ---
%
%   DCOMP - Struct: The initialized or updated density analysis data structure.

% --- Main Controller ---

% If density analysis is turned off, do nothing.
if ~isfield(P, 'dens') || ~P.dens
    return;
end

switch lower(mode)
    case 'setup'
        % The setup function handles its own file loading/saving
        DCOMP = setup_density_bins(S,data_folder);
    case 'accumulate'
        % The accumulation function updates the DCOMP struct in memory
        DCOMP = accumulate_density(DCOMP, S, p);
    otherwise
        error("Invalid mode specified. Use 'setup' or 'accumulate'.");
end

end


% --- Sub-function for Setup ---
function DCOMP = setup_density_bins(S,data_folder)
    
    fprintf('--- Setting up Density Analysis ---\n');

    filenamesaveseries = 'DMAP_%.0f_%.0e_%.0e_%.0f_%.0e.mat';
    filesave = sprintf(filenamesaveseries, S.bc, S.rp, S.phi, S.N, S.rc);

    % --- Attempt to Load Existing Setup File ---
    if exist(filesave, 'file')
        fprintf('Loading existing density setup file: %s\n', filesave);
        load(filesave, 'DCOMP');
        return; % Exit the function if loaded successfully
    end
    
    fprintf('No existing file found. Generating new density setup...\n');
    fprintf('(This may take a long time due to Monte Carlo simulations.)\n');

    % --- Common Parameters ---
    azbins = 90; 
    elbins = 45;

    % --- Common Monte Carlo Sphere Generation ---
    DCOMP.sphereN = 1e4;
    sphereaz = rand(DCOMP.sphereN, 1) * 2 * pi;
    sphereel = asin(2 * rand(DCOMP.sphereN, 1) - 1);
    sphererho = (rand(DCOMP.sphereN, 1).^(1/3)) .* S.rp;
    [x, y, z] = sph2cart(sphereaz, sphereel, sphererho);
    DCOMP.spherexyz = [x, y, z];
    clear x y z sphereaz sphereel sphererho;
    
    % --- BC-Specific Bin Definitions and Initialization ---
    DCOMP.edges.az = linspace(0, 2*pi, azbins)';
    DCOMP.edges.el = linspace(-pi/2, pi/2, elbins)';
    
    switch S.bc
        case 1 % SBC
            maxdist=S.br;
            BS=(maxdist:-0.02*S.rp:maxdist-S.rc)';
            IS=(maxdist-S.rc:-S.rp:0)';
            DCOMP.edges.rho=sort(unique([IS;BS]));
            DCOMP.counts_real.rho=zeros(numel(DCOMP.edges.rho)-1,1);
            DCOMP.counts_ghosts.rho=zeros(numel(DCOMP.edges.rho)-1,1);
            DCOMP.counts_real.az=zeros(azbins-1,1);
            DCOMP.counts_ghosts.az=zeros(azbins-1,1);
            DCOMP.counts_real.el=zeros(elbins-1,1);
            DCOMP.counts_ghosts.el=zeros(elbins-1,1);
            
        case {2, 4} % PBC CUBIC or BB
            paralleldist=S.br;
            BS=(paralleldist:-0.02*S.rp:(paralleldist-2*S.rc))';
            IS=linspace(0,(paralleldist-S.rc),ceil((paralleldist-S.rc)/S.rp))';
            DCOMP.edges.rho=sort(unique([IS;BS]));
            DCOMP.edges.slabs = sort(unique([IS;BS]));
            DCOMP.centers.slabs = DCOMP.edges.slabs(1:end-1,:) + diff(DCOMP.edges.slabs)./2;
            DCOMP.vols.slabs = DCOMP.edges.slabs(2:end,:).^3 - DCOMP.edges.slabs(1:end-1,:).^3;
            DCOMP.numshells.slabs = numel(DCOMP.centers.slabs);
            DCOMP.counts.slabs = zeros(DCOMP.numshells.slabs,1);
            
        case 3 % PBC FCC
            paralleldist=(2*S.br)/sqrt(6);
            BS=(paralleldist:-0.02*S.rp:(paralleldist-2*S.rc))';
            IS=linspace(0,(paralleldist-S.rc),ceil((paralleldist-S.rc)/S.rp))';
            DCOMP.edges.rho=sort(unique([IS;BS]));
            DCOMP.edges.slabs = sort(unique([IS;BS]));
            DCOMP.centers.slabs = DCOMP.edges.slabs(1:end-1,:) + diff(DCOMP.edges.slabs)./2;
            DCOMP.vols.slabs = (1/sqrt(2))*((DCOMP.edges.slabs(2:end,:)*2*sqrt(3/2)).^3-(DCOMP.edges.slabs(1:end-1,:)*2*sqrt(3/2)).^3);
            DCOMP.numshells.slabs = numel(DCOMP.centers.slabs);
            DCOMP.counts.slabs = zeros(DCOMP.numshells.slabs,1);
    end
    
    % --- Finalize Common Properties (Centers, Vols, Counts, etc.) ---
    DCOMP.centers.rho = DCOMP.edges.rho(1:end-1,:) + diff(DCOMP.edges.rho)./2;
    DCOMP.centers.az = DCOMP.edges.az(1:end-1,:) + (DCOMP.edges.az(2)-DCOMP.edges.az(1))/2;
    DCOMP.centers.el = DCOMP.edges.el(1:end-1,:) + (DCOMP.edges.el(2)-DCOMP.edges.el(1))/2;
    DCOMP.numshells.rho = numel(DCOMP.centers.rho);
    DCOMP.numshells.az = numel(DCOMP.centers.az);
    DCOMP.numshells.el = numel(DCOMP.centers.el);
    DCOMP.vols.rho = (4/3)*pi*(DCOMP.edges.rho(2:end).^3 - DCOMP.edges.rho(1:end-1).^3);
    DCOMP.vols.az = ones(DCOMP.numshells.az,1).*(S.bv/DCOMP.numshells.az);
    DCOMP.vols.el = ((2*pi*S.br^3)/3).*(sin(DCOMP.edges.el(2:end,1))-sin(DCOMP.edges.el(1:end-1,1))); % Corrected from cos to sin
    DCOMP.masks.core = DCOMP.centers.rho <= S.br;
    DCOMP.masks.outercore = DCOMP.masks.core & DCOMP.centers.rho > (S.br/2);
    DCOMP.masks.mittelcore = DCOMP.centers.rho > (S.br/4) & DCOMP.centers.rho < (3*S.br/4);
    DCOMP.masks.halo = DCOMP.centers.rho > S.br;
    DCOMP.counts.rho = zeros(DCOMP.numshells.rho,1);
    DCOMP.counts.az = zeros(DCOMP.numshells.az,1);
    DCOMP.counts.el = zeros(DCOMP.numshells.el,1);
    
    % --- Run the mass distribution calculation ---
    DCOMP = calculate_mass_distribution(DCOMP, S);

    % --- GENERATE AND SAVE DIAGNOSTIC FIGURES ---
    fprintf('Generating Diagnostic Figures...\n');
    
    % Strip extension from filesave to get base name
    [~, base_name, ~] = fileparts(filesave);
    full_prefix = fullfile(data_folder, base_name);
    
    % Call helper functions (figures are saved automatically)
    DCOMP=plot_density_diagnostics(DCOMP, S);
    
    fprintf('Diagnostics saved to: %s_*.fig\n', full_prefix);

    % --- Save the completed setup structure to the file ---
    fprintf('Saving newly generated density setup to: %s\n', filesave);
    save([data_folder,'\',filesave], 'DCOMP', 'S', '-v7.3'); % Use -v7.3 for potentially large structs
    DCOMP = rmfield(DCOMP, 'massdist');
end


% --- Helper function for the long mass distribution calculation ---
function DCOMP = calculate_mass_distribution(DCOMP, S)
    fprintf('Calculating mass distribution kernels...\n');
    
    % Pre-initialize massdist arrays here
    DCOMP.massdist.rho = zeros(DCOMP.numshells.rho, DCOMP.numshells.rho);
    DCOMP.massdist.az = zeros(DCOMP.numshells.az, DCOMP.numshells.rho, DCOMP.numshells.az);
    DCOMP.massdist.el = zeros(DCOMP.numshells.el, DCOMP.numshells.rho, DCOMP.numshells.el);
    if isfield(DCOMP, 'numshells') && isfield(DCOMP.numshells, 'slabs')
         DCOMP.massdist.slabs = (zeros(DCOMP.numshells.slabs, DCOMP.numshells.az, DCOMP.numshells.el, DCOMP.numshells.slabs));
    end

    switch S.bc
        case 1 % SBC
            fprintf('  SBC: rho...\n');
            for idcomp = 1:DCOMP.numshells.rho
                temp_pos = DCOMP.spherexyz;
                temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                temp_pos(vecnorm(temp_pos,2,2)>max(DCOMP.edges.rho),:)=[];
                counts = histcounts(vecnorm(temp_pos, 2, 2), DCOMP.edges.rho);
                DCOMP.massdist.rho(:,idcomp) = counts' / sum(counts); 
            end
            fprintf('  SBC: azimuth...\n');
            for idcomp = 1:DCOMP.numshells.rho % Can be parallelized
                for idcompaz = 1:DCOMP.numshells.az
                    temp_pos = DCOMP.spherexyz;
                    temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                    temp_pos(vecnorm(temp_pos,2,2)>(S.br+S.rc+S.rp),:)=[];
                    az = DCOMP.centers.az(idcompaz);
                    Raz = [cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                    temp_pos = (Raz * temp_pos')';
                    [az_out,~,~] = cart2sph(temp_pos(:,1), temp_pos(:,2), temp_pos(:,3));
                    counts = histcounts(az_out, DCOMP.edges.az);
                    DCOMP.massdist.az(:,idcomp,idcompaz) = counts' / sum(counts); % BUG FIX
                end
            end
            fprintf('  SBC: elevation...\n');
            for idcomp = 1:DCOMP.numshells.rho % Can be parallelized
                for idcompel = 1:DCOMP.numshells.el
                    temp_pos = DCOMP.spherexyz;
                    temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                    temp_pos(vecnorm(temp_pos,2,2)>(S.br+S.rc+S.rp),:)=[];
                    el = DCOMP.centers.el(idcompel);
                    Rel = [cos(el), 0, -sin(el); 0, 1, 0; sin(el), 0, cos(el)];
                    temp_pos = (Rel * temp_pos')';
                    [~,el_out,~] = cart2sph(temp_pos(:,1), temp_pos(:,2), temp_pos(:,3));
                    counts = histcounts(el_out, DCOMP.edges.el);
                    DCOMP.massdist.el(:,idcomp,idcompel) = counts' / sum(counts); % BUG FIX
                end
            end

        case {2, 4} % PBC Cubic or BB
             % --- Mass dist for spherical shells inside the cube ---
            fprintf('  CUBIC/BB: rho...\n');
            for idcomp = 1:DCOMP.numshells.rho
                temp_pos = DCOMP.spherexyz;
                temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                num_points_inside = size(temp_pos, 1);
                if num_points_inside > 0
                    counts = histcounts(vecnorm(temp_pos, 2, 2), DCOMP.edges.rho);
                    DCOMP.massdist.rho(:,idcomp) = counts' / num_points_inside; % BUG FIX
                end
            end
            fprintf('  CUBIC/BB: azimuth...\n');
            for idcomp = 1:DCOMP.numshells.rho % Can be parallelized
                for idcompaz = 1:DCOMP.numshells.az
                    temp_pos = DCOMP.spherexyz;
                    temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                    az = DCOMP.centers.az(idcompaz);
                    Raz = [cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                    temp_pos = (Raz * temp_pos')';
                    temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                    num_points_inside = size(temp_pos, 1);
                    if num_points_inside > 0 
                        [az_out,~,~] = cart2sph(temp_pos(:,1), temp_pos(:,2), temp_pos(:,3));
                        counts = histcounts(az_out, DCOMP.edges.az);
                        DCOMP.massdist.az(:,idcomp,idcompaz) = counts' / num_points_inside; % BUG FIX
                    end
                end
            end
            fprintf('  CUBIC/BB: elevation...\n');
            for idcomp = 1:DCOMP.numshells.rho % Can be parallelized
                for idcompel = 1:DCOMP.numshells.el
                    temp_pos = DCOMP.spherexyz;
                    temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                    el = DCOMP.centers.el(idcompel);
                    Rel = [cos(el), 0, -sin(el); 0, 1, 0; sin(el), 0, cos(el)];
                    temp_pos = (Rel * temp_pos')';
                    temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                    num_points_inside = size(temp_pos, 1);
                    if num_points_inside > 0 
                        [~,el_out,~] = cart2sph(temp_pos(:,1), temp_pos(:,2), temp_pos(:,3));
                        counts = histcounts(el_out, DCOMP.edges.el);
                        DCOMP.massdist.el(:,idcomp,idcompel) = counts' / num_points_inside; % BUG FIX
                    end
                end
            end
            % --- Mass dist for slabs ---
            fprintf('  CUBIC/BB: slabs...\n');
            for idcomp = 1:DCOMP.numshells.slabs
                for idcompaz = 1:DCOMP.numshells.az
                    for idcompel = 1:DCOMP.numshells.el
                        temp_pos = DCOMP.spherexyz;
                        temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.slabs(idcomp);
                        az = DCOMP.centers.az(idcompaz);
                        el = DCOMP.centers.el(idcompel);
                        Raz = [cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                        Rel = [cos(el), 0, -sin(el); 0, 1, 0; sin(el), 0, cos(el)];
                        Rtot = Rel * Raz;
                        temp_pos = (Rtot * temp_pos')';
                        temp_pos(max(abs(temp_pos),[],2)>S.br, :) = [];
                        num_points_inside = size(temp_pos, 1);
                        if num_points_inside > 0
                            slab_dists = max(abs(temp_pos), [], 2);
                            counts = histcounts(slab_dists, DCOMP.edges.slabs);
                            DCOMP.massdist.slabs(:,idcompaz,idcompel,idcomp) = counts' / num_points_inside; % BUG FIX
                        end
                    end
                end
            end
            
        case 3 % PBC FCC
            % --- Mass dist for spherical shells inside the cube ---
            fprintf('  FCC: rho...\n');
            for idcomp = 1:DCOMP.numshells.rho
                temp_pos = DCOMP.spherexyz;
                temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                num_points_inside = size(temp_pos, 1);
                if num_points_inside > 0
                    counts = histcounts(vecnorm(temp_pos, 2, 2), DCOMP.edges.rho);
                    DCOMP.massdist.rho(:,idcomp) = counts' / num_points_inside; % BUG FIX
                end
            end
            fprintf('  FCC: azimuth...\n');
            for idcomp = 1:DCOMP.numshells.rho % Can be parallelized
                for idcompaz = 1:DCOMP.numshells.az
                    temp_pos = DCOMP.spherexyz;
                    temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                    az = DCOMP.centers.az(idcompaz);
                    Raz = [cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                    temp_pos = (Raz * temp_pos')';
                    temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                    num_points_inside = size(temp_pos, 1);
                    if num_points_inside > 0
                        [az_out,~,~] = cart2sph(temp_pos(:,1), temp_pos(:,2), temp_pos(:,3));
                        counts = histcounts(az_out, DCOMP.edges.az);
                        DCOMP.massdist.az(:,idcomp,idcompaz) = counts' / num_points_inside; % BUG FIX
                    end
                end
            end
            fprintf('  FCC: elevation...\n');
            for idcomp = 1:DCOMP.numshells.rho % Can be parallelized
                for idcompel = 1:DCOMP.numshells.el
                    temp_pos = DCOMP.spherexyz;
                    temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.rho(idcomp);
                    el = DCOMP.centers.el(idcompel);
                    Rel = [cos(el), 0, -sin(el); 0, 1, 0; sin(el), 0, cos(el)];
                    temp_pos = (Rel * temp_pos')';
                    temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                    num_points_inside = size(temp_pos, 1);
                    if num_points_inside > 0
                        [~,el_out,~] = cart2sph(temp_pos(:,1), temp_pos(:,2), temp_pos(:,3));
                        counts = histcounts(el_out, DCOMP.edges.el);
                        DCOMP.massdist.el(:,idcomp,idcompel) = counts' / num_points_inside; % BUG FIX
                    end
                end
            end

            % --- Mass dist for slabs ---
            fprintf('  FCC: slabs...\n');
            for idcomp = 1:DCOMP.numshells.slabs
                for idcompaz = 1:DCOMP.numshells.az
                    for idcompel = 1:DCOMP.numshells.el
                        temp_pos = DCOMP.spherexyz;
                        % BUG FIX: Position particle at center of SLAB, not rho shell
                        temp_pos(:,1) = temp_pos(:,1) + DCOMP.centers.slabs(idcomp);
                        temp_pos(vecnorm(temp_pos,2,2) > S.br, :) = [];
                        az = DCOMP.centers.az(idcompaz);
                        el = DCOMP.centers.el(idcompel);
                        Raz = [cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                        Rel = [cos(el), 0, -sin(el); 0, 1, 0; sin(el), 0, cos(el)];
                        Rtot = Rel * Raz;
                        temp_pos = (Rtot * temp_pos')';
                        
                        % Find distance to wall for each point in the cloud
                        temp_frac = temp_pos * S.fcc.invA;
                        is_inside = all(temp_frac >= 0 & temp_frac < 1, 2);
                        temp_frac = temp_frac(is_inside, :);
                        
                        num_points_inside = size(temp_frac, 1);
                        if num_points_inside > 0
                            temp_frac_centered = temp_frac - round(temp_frac);
                            H_frac = max(abs(temp_frac_centered), [], 2);
                            a_prim = 2 * S.br; 
                            R_max = a_prim / sqrt(6);
                            H = H_frac * 2 * R_max;
                            counts = histcounts(H, DCOMP.edges.slabs);
                            DCOMP.massdist.slabs(:,idcompaz,idcompel,idcomp) = counts' / num_points_inside; % BUG FIX
                        end
                    end
                end
            end
    end
end


% --- Sub-function for Accumulation ---
function DCOMP = accumulate_density(DCOMP, S, p)

    if isempty(p)
        error('Particle data matrix "p" must be provided for accumulation mode.');
    end

    switch S.bc
        case 1 % SBC
            % Accumulate counts for all, real, and ghost particles separately
            [az, el, rho] = cart2sph(p(:,1), p(:,2), p(:,3));
            DCOMP.counts.rho = DCOMP.counts.rho + histcounts(rho, DCOMP.edges.rho)';
            DCOMP.counts.az  = DCOMP.counts.az  + histcounts(az, DCOMP.edges.az)';
            DCOMP.counts.el  = DCOMP.counts.el  + histcounts(el, DCOMP.edges.el)';
            
            [az_r, el_r, rho_r] = cart2sph(p(1:S.N,1), p(1:S.N,2), p(1:S.N,3));
            DCOMP.counts_real.rho = DCOMP.counts_real.rho + histcounts(rho_r, DCOMP.edges.rho)';
            DCOMP.counts_real.az  = DCOMP.counts_real.az  + histcounts(az_r, DCOMP.edges.az)';
            DCOMP.counts_real.el  = DCOMP.counts_real.el  + histcounts(el_r, DCOMP.edges.el)';
            
            if size(p, 1) > S.N
                [az_g, el_g, rho_g] = cart2sph(p(S.N+1:end,1), p(S.N+1:end,2), p(S.N+1:end,3));
                DCOMP.counts_ghosts.rho = DCOMP.counts_ghosts.rho + histcounts(rho_g, DCOMP.edges.rho)';
                DCOMP.counts_ghosts.az  = DCOMP.counts_ghosts.az  + histcounts(az_g, DCOMP.edges.az)';
                DCOMP.counts_ghosts.el  = DCOMP.counts_ghosts.el  + histcounts(el_g, DCOMP.edges.el)';
            end

        case {2, 4} % PBC cubic or BB
            % Spherical analysis within the cubic domain
            idx_cutoff = vecnorm(p(:,1:3), 2, 2) <= S.br;
            [az, el, rho] = cart2sph(p(idx_cutoff,1), p(idx_cutoff,2), p(idx_cutoff,3));
            DCOMP.counts.rho = DCOMP.counts.rho + histcounts(rho, DCOMP.edges.rho)';
            DCOMP.counts.az  = DCOMP.counts.az  + histcounts(az, DCOMP.edges.az)';
            DCOMP.counts.el  = DCOMP.counts.el  + histcounts(el, DCOMP.edges.el)';
            
            % Slab analysis
            densslabs = max(abs(p(1:S.N, 1:3)), [], 2);
            DCOMP.counts.slabs = DCOMP.counts.slabs + histcounts(densslabs, DCOMP.edges.slabs)';

        case 3 % PBC fcc
            % Spherical analysis within the domain
            idx_cutoff = vecnorm(p(:,1:3), 2, 2) <= S.br; % Note: This is an approximation
            [az, el, rho] = cart2sph(p(idx_cutoff,1), p(idx_cutoff,2), p(idx_cutoff,3));
            DCOMP.counts.rho = DCOMP.counts.rho + histcounts(rho, DCOMP.edges.rho)';
            DCOMP.counts.az  = DCOMP.counts.az  + histcounts(az, DCOMP.edges.az)';
            DCOMP.counts.el  = DCOMP.counts.el  + histcounts(el, DCOMP.edges.el)';
            
            % Slab analysis based on distance to Wigner-Seitz cell wall
            pos_frac_centered = p(1:S.N,1:3)*S.fcc.invA - round(p(1:S.N,1:3)*S.fcc.invA);
            H_frac = max(abs(pos_frac_centered), [], 2);
            densslabs = H_frac * 2 * ( (2 * S.br) / sqrt(6) );
            DCOMP.counts.slabs = DCOMP.counts.slabs + histcounts(densslabs, DCOMP.edges.slabs)';
    end
end