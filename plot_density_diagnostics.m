function DCOMP=plot_density_diagnostics(DCOMP, S)
    % PLOT_DENSITY_DIAGNOSTICS
    % Generates and saves 3 diagnostic figures for mass distribution:
    % 1. Radial (Distance vs Distance)
    % 2. Azimuthal (Radius vs Angle)
    % 3. Elevation (Radius vs Angle)
    %
    % Usage: plot_density_diagnostics(DCOMP, S, [filename_prefix])
    
    
    fprintf('Generating Density Diagnostics...\n');
    
    % --- 1. RADIAL DIAGNOSTIC ---
    try
        DCOMP=plot_radial(DCOMP, S);
    catch ME
        fprintf('Skipping Radial Plot: %s\n', ME.message);
    end
    
    % --- 2. AZIMUTHAL DIAGNOSTIC ---
    try
        DCOMP=plot_azimuth(DCOMP, S);
    catch ME
        fprintf('Skipping Azimuth Plot: %s\n', ME.message);
    end
    
    % --- 3. ELEVATION DIAGNOSTIC ---
    try
        DCOMP=plot_elevation(DCOMP, S);
    catch ME
        fprintf('Skipping Elevation Plot: %s\n', ME.message);
    end
end

% =========================================================================
% SUB-FUNCTIONS
% =========================================================================

function DCOMP=plot_radial(DCOMP, S)
    f_radial = setup_figure();
    
    % Data Extraction
    r = DCOMP.centers.rho;   
    mat = DCOMP.massdist.rho;     
    mat = check_orientation(mat, length(r), length(r));
    mat_norm = normalize_matrix(mat);
    
    % Plotting
    ax = axes('Parent', f_radial); hold(ax, 'on');
    
    step = 1; 
    for i = 1:step:length(r)
        x_pos = r(i);
        if x_pos > 1.05 * S.br, continue; end
        
        z_profile = mat_norm(:, i);
        x_vec = ones(size(r)) * x_pos;
        y_vec = r;
        plot_surface_hack(ax, x_vec, y_vec, z_profile);
    end
    
    % Annotations
    plot3(ax, r, r, ones(size(r))*1.05, 'w--', 'LineWidth', 2); 
    draw_wall_lines(ax, S.br, S.rc, 0, max(r));
    
    % Format
    format_axes(ax, 'Particle Center Position [m]', 'Mass Distribution r [m]', 'Radial Mass Distribution');
    xlim(ax, [0 S.br*1.05]); ylim(ax, [0 S.br*1.05]); zlim(ax, [0 1.1]);
    view(ax, -60, 30);

    DCOMP.massdist.f_radial=f_radial;
end

function DCOMP=plot_azimuth(DCOMP, S)
    f_az = setup_figure();
    
    theta_vec = DCOMP.centers.az(:)'; 
    r_vec = DCOMP.centers.rho;        
    val_3d = DCOMP.massdist.az;       
    
    % Slice Middle
    theta_idx = floor(length(theta_vec)/2);
    p_angle = theta_vec(theta_idx);
    
    % Extract Slice [Radius x Angle]
    [d1, d2, ~] = size(val_3d);
    if d1 == length(theta_vec) && d2 == length(r_vec)
        mat_slice = squeeze(val_3d(theta_idx, :, :)); 
    elseif d2 == length(theta_vec) && d1 == length(r_vec)
        mat_slice = squeeze(val_3d(:, theta_idx, :));
    else
        close(f_az); error('Dimension mismatch in Azimuth.');
    end
    
    mat_norm = normalize_rows(mat_slice);
    
    % Plotting
    ax = axes('Parent', f_az); hold(ax, 'on');
    
    for i = 1:length(r_vec)
        curr_r = r_vec(i);
        if curr_r > 1.05 * S.br, continue; end
        
        z_profile = mat_norm(i, :); 
        x_vec = theta_vec; 
        y_vec = ones(size(theta_vec)) * curr_r;
        
        plot_surface_hack(ax, x_vec, y_vec, z_profile);
    end
    
    % Annotations
    plot3(ax, [p_angle p_angle], [min(r_vec) max(r_vec)], [1.05 1.05], 'w--', 'LineWidth', 3);
    yline(ax, S.br, 'r-', 'LineWidth', 3);
    yline(ax, S.br-S.rc, 'g--', 'LineWidth', 3);
    
    format_axes(ax, 'Mass Angle \theta [rad]', 'Particle Radius r [m]', ...
        sprintf('Azimuthal Diagnostic (Slice @ %.2f rad)', p_angle));
    
    xlim(ax, [min(theta_vec) max(theta_vec)]); ylim(ax, [0 1.05*S.br]); zlim(ax, [0 1.1]);
    view(ax, -15, 60);

    DCOMP.massdist.f_az=f_az;
end

function DCOMP=plot_elevation(DCOMP, S)
    f_el = setup_figure();
    
    theta_vec = DCOMP.centers.el(:)'; 
    r_vec = DCOMP.centers.rho;        
    val_3d = DCOMP.massdist.el;       
    
    theta_idx = floor(length(theta_vec)/2);
    p_angle = theta_vec(theta_idx);
    
    [d1, d2, ~] = size(val_3d);
    if d1 == length(theta_vec) && d2 == length(r_vec)
        mat_slice = squeeze(val_3d(theta_idx, :, :)); 
    elseif d2 == length(theta_vec) && d1 == length(r_vec)
        mat_slice = squeeze(val_3d(:, theta_idx, :));
    else
        close(f_el); error('Dimension mismatch in Elevation.');
    end
    
    mat_norm = normalize_rows(mat_slice);
    
    ax = axes('Parent', f_el); hold(ax, 'on');
    
    for i = 1:length(r_vec)
        curr_r = r_vec(i);
        if curr_r > 1.05 * S.br, continue; end
        
        z_profile = mat_norm(i, :); 
        x_vec = theta_vec; 
        y_vec = ones(size(theta_vec)) * curr_r;
        
        plot_surface_hack(ax, x_vec, y_vec, z_profile);
    end
    
    plot3(ax, [p_angle p_angle], [min(r_vec) max(r_vec)], [1.05 1.05], 'w--', 'LineWidth', 3);
    yline(ax, S.br, 'r-', 'LineWidth', 3);
    yline(ax, S.br-S.rc, 'g--', 'LineWidth', 3);
    
    format_axes(ax, 'Mass Angle \phi [rad]', 'Particle Radius r [m]', ...
        sprintf('Elevation Diagnostic (Slice @ %.2f rad)', p_angle));
    
    xlim(ax, [min(theta_vec) max(theta_vec)]); ylim(ax, [0 1.05*S.br]); zlim(ax, [0 1.1]);
    view(ax, -15, 60);

    DCOMP.massdist.f_el=f_el;
end

% =========================================================================
% UTILITIES
% =========================================================================

function f = setup_figure()
    f = figure('Units','normalized','Position',[0.1 0.1 0.6 0.6], 'Color', 'k', 'Visible', 'off');
    colormap('jet'); 
end

function mat = check_orientation(mat, r1, r2)
    if size(mat, 1) ~= r1 || size(mat, 2) ~= r2
        if size(mat, 1) ~= r1, mat = mat'; end
    end
end

function mat_norm = normalize_matrix(mat)
    mat_norm = zeros(size(mat));
    for c = 1:size(mat, 2)
        col_max = max(mat(:, c));
        if col_max > 0, mat_norm(:, c) = mat(:, c) / col_max; end
    end
end

function mat_norm = normalize_rows(mat)
    mat_norm = zeros(size(mat));
    for r = 1:size(mat, 1)
        row_max = max(mat(r, :));
        if row_max > 0, mat_norm(r, :) = mat(r, :) / row_max; end
    end
end

function plot_surface_hack(ax, x_vec, y_vec, z_profile)
    % Ensure rows
    x_vec = x_vec(:)'; y_vec = y_vec(:)'; z_profile = z_profile(:)';
    
    xx = [x_vec; x_vec];
    yy = [y_vec; y_vec];
    zz = [z_profile; z_profile];
    cc = zz; 
    
    surface(xx, yy, zz, cc, ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.5, 'Parent', ax);
end

function draw_wall_lines(ax, br, rc, min_y, max_y)
    line([br br], [min_y max_y], [0 0], 'Color', 'r', 'LineWidth', 3, 'Parent', ax);
    line([min_y br], [br br], [0 0], 'Color', 'r', 'LineWidth', 3, 'Parent', ax);
    line([br-rc br-rc], [min_y max_y], [0 0], 'Color', 'g', 'LineWidth', 3, 'Parent', ax);
    line([min_y br-rc], [br-rc br-rc], [0 0], 'Color', 'g', 'LineWidth', 3, 'Parent', ax);
end

function format_axes(ax, xlab, ylab, titl)
    xlabel(ax, xlab, 'Color','w');
    ylabel(ax, ylab, 'Color','w');
    zlabel(ax, 'Normalized Intensity','Color','w');
    title(ax, titl, 'Color','w','FontSize',14,'FontWeight','bold');
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', ...
        'LineWidth', 3, 'FontSize', 12, 'FontWeight', 'bold');
    grid(ax, 'on'); ax.GridColor = [0.3 0.3 0.3];
    c = colorbar(ax); c.Color = 'w'; c.Label.String = 'Intensity'; c.Label.Color = 'w';
end