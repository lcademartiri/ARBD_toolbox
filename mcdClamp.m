function clamp=mcdClamp(nodisp,rp,DISP,esdiff,timestep,forcevector,kbT)

    % --- MC array of random relative displacements
    % generate disp random separations between two particles 
    mcd.az=rand(nodisp,1)*2*pi;
    mcd.el=acos(2*rand(nodisp,1)-1)-pi/2;
    mcd.rho=(rand(nodisp,1).^(1/3)).*(2e-9)+rp+rp;
    % convert them to cartesian
    [mcd.xi,mcd.yi,mcd.zi]=sph2cart(mcd.az,mcd.el,mcd.rho);
    % calculate nodisp relative displacements by difference of two shuffled versions of the DISP library
    mcd.DISPrel=DISP-DISP(randperm(nodisp),:);
    % calculate the step size for the relative motion of two brownian particles
    mcd.sigmaB_rel = sqrt(2*(esdiff+esdiff)*timestep);
    % ---
    
    % --- Gauss-Legendre
    % Constructs 12-point Gauss–Legendre nodes & weights on [0, 1]
    % This is used to integrate forces along a linear path between two positions.
    mcd.qnodes=12;
    mcd.beta = 0.5 ./ sqrt(1-(2*(1:mcd.qnodes-1)).^(-2));
    mcd.T = diag(mcd.beta,1) + diag(mcd.beta,-1);
    [mcd.V,mcd.D] = eig(mcd.T);
    mcd.x = diag(mcd.D);          % points in [-1,1]
    mcd.w = 2 * (mcd.V(1,:)').^2; % weights for [-1,1]
    mcd.s_quad = 0.5*(mcd.x+1);
    mcd.w_quad = 0.5*mcd.w;
    mcd.N = size(mcd.xi,1);
    % quadrature
    mcd.r0_rep = reshape([mcd.xi,mcd.yi,mcd.zi], [mcd.N,1,3]);          % Nx1x3
    mcd.rB_rep = reshape(mcd.DISPrel, [mcd.N,1,3]);      % Nx1x3
    mcd.s_rep = reshape(mcd.s_quad, [1,mcd.qnodes,1]);       % 1xMx1
    % for each random pair, it traces a straight path between the two configurations and samples the LJ force magnitude F(r) (from the lookup table H).
    mcd.r_path = mcd.r0_rep + mcd.s_rep .* mcd.rB_rep;      % NxMx3
    mcd.rnorm_path = sqrt(sum(mcd.r_path.^2,3));    % NxM
    mcd.r_query = mcd.rnorm_path;
    mcd.r_min = forcevector(1,1);
    mcd.r_max = forcevector(end,1);
    mcd.F_min = forcevector(1,2);
    % Clamp to interpolation domain
    mcd.r_query(mcd.r_query < mcd.r_min) = mcd.r_min;
    mcd.r_query(mcd.r_query > mcd.r_max) = mcd.r_max;
    mcd.Fmag_path = interp1(forcevector(:,1), forcevector(:,2), mcd.r_query, 'pchip', 0); % NxM
    mcd.Fmag_path(mcd.r_query < forcevector(1,1)) = mcd.F_min;  % manually saturate below contact
    mcd.nhat_path = mcd.r_path ./ mcd.rnorm_path(:,:,ones(1,1,3)); % NxMx3
    mcd.Fvec_path = mcd.Fmag_path(:,:,ones(1,1,3)) .* mcd.nhat_path; % NxMx3 % vector force along the path
    % integrate the force along the path (weighted by quadrature) to estimate the deterministic displacement due to LJ over one timestep for that random pair.
    mcd.det_disp_rel = squeeze(sum(mcd.Fvec_path .* mcd.w_quad(ones(mcd.N,1),:,ones(1,3)),2)) * ((esdiff + esdiff) / (kbT)) * timestep; % Nx3
    % ratio of "deterministic LJ displacement” to “random Brownian displacement” for each sample.
    % dimensionless measure of how strong LJ motion is compared to stochastic diffusion
    mcd.forcedisps=vecnorm(mcd.det_disp_rel,2,2)./vecnorm(mcd.DISPrel,2,2);
    % filters out absurd or negative ratios (probably from zero denominators or non-physical configurations)
    mcd.forcedisps(mcd.forcedisps>1,:)=[];
    mcd.forcedisps(mcd.forcedisps<=0,:)=[];
    % fits a log-normal PDF to the ratios
    mcd.pd = fitdist(mcd.forcedisps,'Lognormal');
    mcd.mu = mcd.pd.mu;
    mcd.sigma = mcd.pd.sigma;
    % 99th-percentile upper bound (since 2.576 ≈ z for 99%).
    clamp=exp(mcd.mu+mcd.sigma*2.576);
end