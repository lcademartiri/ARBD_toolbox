function EQUIP = equipartitionEngine(mode, EQUIP, S, P, p)
%EQUIPARTITIONENGINE Manages setup and accumulation for equipartition analysis.
%
% This function operates in two modes: 'setup' and 'accumulate'.
%
% --- SYNTAX ---
%
% 1. To initialize the data structures before the simulation starts:
%    EQUIP = equipartitionEngine('setup', [], S, P, []);
%
% 2. To accumulate data at each timestep inside the simulation loop:
%    EQUIP = equipartitionEngine('accumulate', EQUIP, S, P, p);
%
% --- INPUTS ---
%
%   mode  - String: 'setup' to initialize, or 'accumulate' to add data.
%   EQUIP - Struct: Pass an empty array [] for 'setup'. For 'accumulate',
%           pass the EQUIP struct from the previous timestep.
%   S     - Struct: Contains simulation settings (S.bc, S.br, S.rp, S.stdx, S.fcc.invA).
%   P     - Struct: Contains simulation parameters (P.equipartition flag).
%   p     - Matrix: Particle data [positions, ..., displacements]. Required
%           for 'accumulate' mode only. Pass [] for 'setup'.
%
% --- OUTPUT ---
%
%   EQUIP - Struct: The initialized or updated equipartition data structure.

% --- Main Controller ---

% If equipartition analysis is turned off, do nothing.
if ~isfield(P, 'equipartition') || ~P.equipartition
    return;
end

switch lower(mode)
    case 'setup'
        EQUIP = setup_bins(S);
    case 'accumulate'
        EQUIP = accumulate_data(EQUIP, S, p);
    otherwise
        error("Invalid mode specified. Use 'setup' or 'accumulate'.");
end

end


% --- Sub-function for Setup ---
function EQUIP = setup_bins(S)
    
    % --- Define Displacement Bins (common for all BCs) ---
    EQUIP.edges.disp = (0:0.02*S.stdx:10*S.stdx)';
    EQUIP.edges.disp = sort(unique([-EQUIP.edges.disp; EQUIP.edges.disp])); % Symmetric
    EQUIP.centers.disp = EQUIP.edges.disp(1:end-1) + 0.5*diff(EQUIP.edges.disp(1:2));
    EQUIP.numshells.disp = numel(EQUIP.centers.disp);

    % --- Define Spatial Bins (dependent on BC) ---
    switch S.bc
        case 1 % SBC (Spherical)
            EQUIP.bintype = 'radius';
            maxdist=S.br;
            BS=(maxdist:-0.02*S.rp:maxdist-S.rc)';
            IS=(maxdist-S.rc:-S.rp:0)';
            EQUIP.edges.rho=sort(unique([IS;BS]));
            
        case {2, 4} % PBC Cubic or BB (Big Box)
            EQUIP.bintype = 'max_component';
            maxdist=S.br;
            BS=(maxdist:-0.02*S.rp:maxdist-S.rc)';
            IS=(maxdist-S.rc:-S.rp:0)';
            EQUIP.edges.rho=sort(unique([IS;BS]));
            
        case 3 % PBC FCC
            EQUIP.bintype = 'dist_to_wall';
            % For clarity, the spatial bins for FCC are named 'slabs' internally
            L = 2 * S.br;
            maxdist = L / sqrt(6);
            BS=(maxdist:-0.02*S.rp:maxdist-S.rc)';
            IS=(maxdist-S.rc:-S.rp:0)';
            EQUIP.edges.slabs = sort(unique([IS;BS]));
            EQUIP.edges.rho = EQUIP.edges.slabs; % Maintain .rho for consistency if needed
            
        otherwise
            error('Unknown boundary condition S.bc = %d', S.bc);
    end

    % --- Finalize Spatial Bins and Initialize Histograms ---
    EQUIP.centers.rho  = EQUIP.edges.rho(1:end-1)  + 0.5*diff(EQUIP.edges.rho(1:2));
    EQUIP.numshells.rho  = numel(EQUIP.centers.rho);
    
    if S.bc == 1
        % x, y, z, and radial displacements for SBC
        EQUIP.histos = zeros(EQUIP.numshells.rho, EQUIP.numshells.disp, 4);
    else
        % x, y, z displacements for all other cases
        EQUIP.histos = zeros(EQUIP.numshells.rho, EQUIP.numshells.disp, 3);
    end

end


% --- Sub-function for Accumulation ---
function EQUIP = accumulate_data(EQUIP, S, p)

    if isempty(p)
        error('Particle data matrix "p" must be provided for accumulation mode.');
    end

    % --- Build intermediate quantities (common for all BCs) ---
    initial_positions = p(1:S.N, 1:3);
    displacements     = p(1:S.N, 8:10) - p(1:S.N, 1:3);
    
    % --- Bin displacements (common for all BCs) ---
    [~, ~, dx_bin] = histcounts(displacements(:, 1), EQUIP.edges.disp);
    [~, ~, dy_bin] = histcounts(displacements(:, 2), EQUIP.edges.disp);
    [~, ~, dz_bin] = histcounts(displacements(:, 3), EQUIP.edges.disp);

    % --- Perform binning and accumulation based on BC ---
    switch S.bc
        case 1 % SBC
            radius = vecnorm(initial_positions, 2, 2);
            [~, ~, r_bin] = histcounts(radius, EQUIP.edges.rho);
            
            unit_radial_vectors = initial_positions ./ radius;
            radial_displacement = sum(displacements .* unit_radial_vectors, 2);
            [~, ~, dr_bin] = histcounts(radial_displacement, EQUIP.edges.disp);
            
            valid_x = (r_bin > 0) & (dx_bin > 0);
            valid_y = (r_bin > 0) & (dy_bin > 0);
            valid_z = (r_bin > 0) & (dz_bin > 0);
            valid_r = (r_bin > 0) & (dr_bin > 0);
            
            if any(valid_x), EQUIP.histos(:,:,1) = EQUIP.histos(:,:,1) + accumarray([r_bin(valid_x), dx_bin(valid_x)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_y), EQUIP.histos(:,:,2) = EQUIP.histos(:,:,2) + accumarray([r_bin(valid_y), dy_bin(valid_y)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_z), EQUIP.histos(:,:,3) = EQUIP.histos(:,:,3) + accumarray([r_bin(valid_z), dz_bin(valid_z)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_r), EQUIP.histos(:,:,4) = EQUIP.histos(:,:,4) + accumarray([r_bin(valid_r), dr_bin(valid_r)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end

        case {2, 4} % PBC Cubic or BB
            dist_from_center = max(abs(initial_positions), [], 2);
            [~, ~, r_bin] = histcounts(dist_from_center, EQUIP.edges.rho);
            
            valid_x = (r_bin > 0) & (dx_bin > 0);
            valid_y = (r_bin > 0) & (dy_bin > 0);
            valid_z = (r_bin > 0) & (dz_bin > 0);

            if any(valid_x), EQUIP.histos(:,:,1) = EQUIP.histos(:,:,1) + accumarray([r_bin(valid_x), dx_bin(valid_x)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_y), EQUIP.histos(:,:,2) = EQUIP.histos(:,:,2) + accumarray([r_bin(valid_y), dy_bin(valid_y)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_z), EQUIP.histos(:,:,3) = EQUIP.histos(:,:,3) + accumarray([r_bin(valid_z), dz_bin(valid_z)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end

        case 3 % PBC FCC
            pos_frac_centered = initial_positions * S.fcc.invA - round(initial_positions * S.fcc.invA);
            H_frac = max(abs(pos_frac_centered), [], 2);
            dist_to_wall = H_frac * 2 * ( (2 * S.br) / sqrt(6) );
            [~, ~, slab_bin] = histcounts(dist_to_wall, EQUIP.edges.rho); % Use .rho for consistency

            valid_x = (slab_bin > 0) & (dx_bin > 0);
            valid_y = (slab_bin > 0) & (dy_bin > 0);
            valid_z = (slab_bin > 0) & (dz_bin > 0);

            if any(valid_x), EQUIP.histos(:,:,1) = EQUIP.histos(:,:,1) + accumarray([slab_bin(valid_x), dx_bin(valid_x)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_y), EQUIP.histos(:,:,2) = EQUIP.histos(:,:,2) + accumarray([slab_bin(valid_y), dy_bin(valid_y)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
            if any(valid_z), EQUIP.histos(:,:,3) = EQUIP.histos(:,:,3) + accumarray([slab_bin(valid_z), dz_bin(valid_z)], 1, [EQUIP.numshells.rho, EQUIP.numshells.disp]); end
    end
end