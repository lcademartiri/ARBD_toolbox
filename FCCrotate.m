function pr=FCCrotate(p,v)

    % Choose z-axis of new frame as v
    z_new = v;
    
    % Pick an arbitrary vector not parallel to z_new
    tmp = [0 0 1];
    if abs(dot(tmp,z_new)) > 0.9
        tmp = [0 1 0]; % fallback if nearly parallel
    end
    
    % Construct orthonormal basis (x_new, y_new, z_new)
    x_new = cross(tmp, z_new); x_new = x_new / norm(x_new);
    y_new = cross(z_new, x_new);
    
    % Rotation matrix
    R = [x_new; y_new; z_new];
    
    % Apply to displacement vectors (rows are dx,dy,dz)
    
    pr = (R * p')'; % rotated into FCC [111]-aligned frame


end