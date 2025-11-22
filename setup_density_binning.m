function dcomp=setup_density_binning(S,azbins,elbins)

%SETUP_DENSITY_BINNING  Initialize density comparison structure for all BCs.
%
%   dcomp = setup_density_binning(S)
%   Initializes binning edges, centers, volumes, masks, counts, and
%   Monte Carlo mass-distribution kernels for density diagnostics.
%
%   Required fields in S:
%       bc   - boundary condition type (1=SBC, 2=PBC cubic, 3=PBC FCC, 4=BB)
%       br   - domain radius or half box length
%       rc   - correction radius (SBC only)
%       rp   - particle radius
%       bv   - box volume
%       fcc.invA - (3x3) inverse lattice matrix, for FCC BC
%       S_br - (scalar) FCC lattice constant (for FCC slabs)

    if S.potential==1
        filenamesaveseries='DMAP_lj_%.0f_%.0e_%.0e_%.0f_%.0e.mat';
    elseif S.potential==2
        filenamesaveseries='DMAP_wca_%.0f_%.0e_%.0e_%.0f_%.0e.mat';
    elseif S.potential==0
        filenamesaveseries='DMAP_hs_%.0f_%.0e_%.0e_%.0f_%.0e.mat';
    end
    filesave=sprintf(filenamesaveseries,S.bc,S.rp,S.phi,S.N,S.rc);

    if exist(filesave,'file')
        load(filesave)
        return
    end

    if S.bc==1 % SBC
        % bins are calibrated on particle size because the behaviour at 
        % the boundary is dictated by the particle size. Binning must
        % reach where mass can reach so S.br+S.rc+S.rp since we include
        % reals and ghosts. In the case of PBCs we only need to go to
        % the boundary because we use MIC.
        dcomp.edges.rho=[0:0.02*S.rp:S.br+S.rc+S.rp]';
        dcomp.edges.az=linspace(0,2*pi,azbins)';
        dcomp.edges.el=linspace(-pi/2,pi/2,elbins)';
        dcomp.centers.rho=dcomp.edges.rho(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.centers.az=dcomp.edges.az(1:end-1,:)+(dcomp.edges.az(2)-dcomp.edges.az(1))/2; % centers of the azimuthal bins
        dcomp.centers.el=dcomp.edges.el(1:end-1,:)+(dcomp.edges.el(2)-dcomp.edges.el(1))/2; % centers of the elevation bins
        dcomp.vols.rho=(4/3)*pi*(dcomp.edges.rho(2:end).^3-dcomp.edges.rho(1:end-1).^3); % volumes of the onion shells
        dcomp.vols.az=ones(numel(dcomp.centers.az),1).*(S.bv/numel(dcomp.centers.az)); % volumes of the azimuthal onion shells
        dcomp.vols.el=((2*pi*S.br^3)/3).*(cos(dcomp.edges.el(2:end,1))-cos(dcomp.edges.el(1:end-1,1))); % volumes of the onion shells
        dcomp.numshells.rho=numel(dcomp.centers.rho); % number of radial shells
        dcomp.numshells.az=numel(dcomp.centers.az); % number of azimuthal shells
        dcomp.numshells.el=numel(dcomp.centers.el); % number of elevation shells
        dcomp.masks.core=dcomp.centers.rho<=S.br; % mask identifying core shells
        dcomp.masks.outercore=dcomp.masks.core & dcomp.centers.rho>(S.br/2); % mask identifying outer half of the domain
        dcomp.masks.mittelcore=dcomp.centers.rho>(S.br/4) & dcomp.centers.rho<(3*S.br/4); % mask identifying the central half of the domain
        dcomp.masks.halo=dcomp.centers.rho>S.br; % mask identifying halo shells
        % initializing the accumulation variables
        dcomp.counts.rho=zeros(dcomp.numshells.rho,1);
        dcomp.counts_real.rho=zeros(dcomp.numshells.rho,1);
        dcomp.counts_ghosts.rho=zeros(dcomp.numshells.rho,1);
        dcomp.counts.az=zeros(dcomp.numshells.az,1);
        dcomp.counts_real.az=zeros(dcomp.numshells.az,1);
        dcomp.counts_ghosts.az=zeros(dcomp.numshells.az,1);
        dcomp.counts.el=zeros(dcomp.numshells.el,1);
        dcomp.counts_real.el=zeros(dcomp.numshells.el,1);
        dcomp.counts_ghosts.el=zeros(dcomp.numshells.el,1);
        dcomp.massdist.rho=zeros(dcomp.numshells.rho,dcomp.numshells.rho);
        dcomp.massdist.az=zeros(dcomp.numshells.az,dcomp.numshells.rho,dcomp.numshells.az);
        dcomp.massdist.el=zeros(dcomp.numshells.el,dcomp.numshells.rho,dcomp.numshells.el);
        % generating spherical point cloud to MC simulate the mass
        % distribution. Sphere is radius S.rp
        dcomp.sphereN=1e4;
        dcomp.sphereaz=rand(dcomp.sphereN,1)*2*pi;
        dcomp.sphereel=asin(2*rand(dcomp.sphereN,1)-1);
        dcomp.sphererho=(rand(dcomp.sphereN,1).^(1/3)).*S.rp;
        [x,y,z]=sph2cart(dcomp.sphereaz,dcomp.sphereel,dcomp.sphererho);
        dcomp.spherexyz=[x,y,z];
        clear x y z
        % generating the massdist array which describes how a particle 
        % will distribute mass across the shells 
        fprintf('density - mass allocation matrix over rho \n')
        for idcomp=1:dcomp.numshells.rho
            dcomp.temp=dcomp.spherexyz;
            dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
            [dcomp.temp,dcomp.edges.rho]=histcounts(vecnorm(dcomp.temp,2,2),dcomp.edges.rho);
            dcomp.massdist.rho(:,idcomp)=dcomp.temp'/size(dcomp.temp,1);
        end
        fprintf('density - mass allocation matrix over rho/azimuth \n')
        for idcomp=1:dcomp.numshells.rho
            for idcompaz=1:dcomp.numshells.az
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                dcomp.temp=(Raz*dcomp.temp')';
                [az,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.az]=histcounts(az,dcomp.edges.az);
                dcomp.massdist.az(:,idcomp,idcompaz)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        fprintf('density - mass allocation matrix over rho/elevation \n')
        for idcomp=1:dcomp.numshells.rho
            for idcompel=1:dcomp.numshells.el
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                el=dcomp.centers.el(idcompel);
                Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                dcomp.temp=(Rel*dcomp.temp')';
                [el,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.el]=histcounts(el,dcomp.edges.el);
                dcomp.massdist.el(:,idcomp,idcompel)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
    elseif S.bc==2 % PBC CUBIC
        % bins are calibrated on particle size because the behaviour at 
        % the boundary is dictated by the particle size. Binning must
        % reach where mass can reach so S.br+S.rc+S.rp since we include
        % reals and ghosts. In the case of PBCs we only need to go to
        % the boundary because we use MIC.
        dcomp.edges.rho=[0:0.02*S.rp:S.br]';
        dcomp.edges.az=linspace(0,2*pi,azbins)';
        dcomp.edges.el=linspace(-pi/2,pi/2,elbins)';
        dcomp.centers.rho=dcomp.edges.rho(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.centers.az=dcomp.edges.az(1:end-1,:)+(dcomp.edges.az(2)-dcomp.edges.az(1))/2; % centers of the azimuthal bins
        dcomp.centers.el=dcomp.edges.el(1:end-1,:)+(dcomp.edges.el(2)-dcomp.edges.el(1))/2; % centers of the elevation bins
        dcomp.vols.rho=(4/3)*pi*(dcomp.edges.rho(2:end).^3-dcomp.edges.rho(1:end-1).^3); % volumes of the onion shells
        dcomp.vols.az=ones(numel(dcomp.centers.az),1).*(S.bv/numel(dcomp.centers.az)); % volumes of the azimuthal onion shells
        dcomp.vols.el=((2*pi*S.br^3)/3).*(cos(dcomp.edges.el(2:end,1))-cos(dcomp.edges.el(1:end-1,1))); % volumes of the onion shells
        dcomp.numshells.rho=numel(dcomp.centers.rho); % number of radial shells
        dcomp.numshells.az=numel(dcomp.centers.az); % number of azimuthal shells
        dcomp.numshells.el=numel(dcomp.centers.el); % number of elevation shells
        dcomp.masks.core=dcomp.centers.rho<=S.br; % mask identifying core shells
        dcomp.masks.outercore=dcomp.masks.core & dcomp.centers.rho>(S.br/2); % mask identifying outer half of the domain
        dcomp.masks.mittelcore=dcomp.centers.rho>(S.br/4) & dcomp.centers.rho<(3*S.br/4); % mask identifying the central half of the domain
        dcomp.masks.halo=dcomp.centers.rho>S.br; % mask identifying halo shells
        % initializing the accumulation variables
        dcomp.counts.rho=zeros(dcomp.numshells.rho,1);
        dcomp.counts.az=zeros(dcomp.numshells.az,1);
        dcomp.counts.el=zeros(dcomp.numshells.el,1);
        dcomp.massdist.rho=zeros(dcomp.numshells.rho,dcomp.numshells.rho);
        dcomp.massdist.az=zeros(dcomp.numshells.az,dcomp.numshells.rho,dcomp.numshells.az);
        dcomp.massdist.el=zeros(dcomp.numshells.el,dcomp.numshells.rho,dcomp.numshells.el);
        % generating spherical point cloud to MC simulate the mass
        % distribution. Sphere is radius S.rp
        dcomp.sphereN=1e4;
        dcomp.sphereaz=rand(dcomp.sphereN,1)*2*pi;
        dcomp.sphereel=asin(2*rand(dcomp.sphereN,1)-1);
        dcomp.sphererho=(rand(dcomp.sphereN,1).^(1/3)).*S.rp;
        [x,y,z]=sph2cart(dcomp.sphereaz,dcomp.sphereel,dcomp.sphererho);
        dcomp.spherexyz=[x,y,z];
        clear x y z
        % generating the massdist array which describes how a particle 
        % will distribute mass across the shells 
        for idcomp=1:dcomp.numshells.rho
            dcomp.temp=dcomp.spherexyz;
            dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
            dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
            [dcomp.temp,dcomp.edges.rho]=histcounts(vecnorm(dcomp.temp,2,2),dcomp.edges.rho);
            dcomp.massdist.rho(:,idcomp)=dcomp.temp'/size(dcomp.temp,1);
        end
        for idcomp=1:dcomp.numshells.rho
            for idcompaz=1:dcomp.numshells.az
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                dcomp.temp=(Raz*dcomp.temp')';
                dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                [az,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.az]=histcounts(az,dcomp.edges.az);
                dcomp.massdist.az(:,idcomp,idcompaz)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        for idcomp=1:dcomp.numshells.rho
            for idcompel=1:dcomp.numshells.el
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                el=dcomp.centers.el(idcompel);
                Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                dcomp.temp=(Rel*dcomp.temp')';
                dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                [el,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.el]=histcounts(el,dcomp.edges.el);
                dcomp.massdist.el(:,idcomp,idcompel)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        % slabs must be handled separately.
        dcomp.edges.slabs=[0:0.02*S.rp:S.br]';
        dcomp.centers.slabs=dcomp.edges.slabs(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.vols.slabs=dcomp.edges.slabs(2:end,:).^3-dcomp.edges.slabs(1:end-1,:).^3;
        dcomp.numshells.slabs=numel(dcomp.centers.slabs); % number of radial shells
        % initializing the accumulation variables
        dcomp.counts.slabs=zeros(dcomp.numshells.slabs,1);
        dcomp.massdist.slabs=zeros(dcomp.numshells.slabs,dcomp.numshells.az,dcomp.numshells.el,dcomp.numshells.slabs);
        for idcomp=1:dcomp.numshells.slabs
            for idcompaz=1:dcomp.numshells.az
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                for idcompel=1:dcomp.numshells.el
                    dcomp.temp=dcomp.spherexyz;
                    dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                    el=dcomp.centers.el(idcompel);
                    Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                    Rtot=Rel*Raz;
                    dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                    dcomp.temp=max(abs((Rtot*dcomp.temp')'),[],2);
                    [dcomp.temp,dcomp.edges.slabs]=histcounts(dcomp.temp,dcomp.edges.slabs);
                    dcomp.massdist.slabs(:,idcompaz,idcompel,idcomp)=dcomp.temp'/size(dcomp.temp,1);
                end
            end
        end            
    elseif S.bc==3 % PBC FCC
        % bins are calibrated on particle size because the behaviour at 
        % the boundary is dictated by the particle size. Binning must
        % reach where mass can reach so S.br+S.rc+S.rp since we include
        % reals and ghosts. In the case of PBCs we only need to go to
        % the boundary because we use MIC.
        dcomp.edges.rho=[0:0.02*S.rp:(2*S.br)/sqrt(6)]';
        dcomp.edges.az=linspace(0,2*pi,azbins)';
        dcomp.edges.el=linspace(-pi/2,pi/2,elbins)';
        dcomp.centers.rho=dcomp.edges.rho(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.centers.az=dcomp.edges.az(1:end-1,:)+(dcomp.edges.az(2)-dcomp.edges.az(1))/2; % centers of the azimuthal bins
        dcomp.centers.el=dcomp.edges.el(1:end-1,:)+(dcomp.edges.el(2)-dcomp.edges.el(1))/2; % centers of the elevation bins
        dcomp.vols.rho=(4/3)*pi*(dcomp.edges.rho(2:end).^3-dcomp.edges.rho(1:end-1).^3); % volumes of the onion shells
        dcomp.vols.az=ones(numel(dcomp.centers.az),1).*(S.bv/numel(dcomp.centers.az)); % volumes of the azimuthal onion shells
        dcomp.vols.el=((2*pi*S.br^3)/3).*(cos(dcomp.edges.el(2:end,1))-cos(dcomp.edges.el(1:end-1,1))); % volumes of the onion shells
        dcomp.numshells.rho=numel(dcomp.centers.rho); % number of radial shells
        dcomp.numshells.az=numel(dcomp.centers.az); % number of azimuthal shells
        dcomp.numshells.el=numel(dcomp.centers.el); % number of elevation shells
        dcomp.masks.core=dcomp.centers.rho<=S.br; % mask identifying core shells
        dcomp.masks.outercore=dcomp.masks.core & dcomp.centers.rho>(S.br/2); % mask identifying outer half of the domain
        dcomp.masks.mittelcore=dcomp.centers.rho>(S.br/4) & dcomp.centers.rho<(3*S.br/4); % mask identifying the central half of the domain
        dcomp.masks.halo=dcomp.centers.rho>S.br; % mask identifying halo shells
        % initializing the accumulation variables
        dcomp.counts.rho=zeros(dcomp.numshells.rho,1);
        dcomp.counts.az=zeros(dcomp.numshells.az,1);
        dcomp.counts.el=zeros(dcomp.numshells.el,1);
        dcomp.massdist.rho=zeros(dcomp.numshells.rho,dcomp.numshells.rho);
        dcomp.massdist.az=zeros(dcomp.numshells.az,dcomp.numshells.rho,dcomp.numshells.az);
        dcomp.massdist.el=zeros(dcomp.numshells.el,dcomp.numshells.rho,dcomp.numshells.el);
        % generating spherical point cloud to MC simulate the mass
        % distribution. Sphere is radius S.rp
        dcomp.sphereN=1e4;
        dcomp.sphereaz=rand(dcomp.sphereN,1)*2*pi;
        dcomp.sphereel=asin(2*rand(dcomp.sphereN,1)-1);
        dcomp.sphererho=(rand(dcomp.sphereN,1).^(1/3)).*S.rp;
        [x,y,z]=sph2cart(dcomp.sphereaz,dcomp.sphereel,dcomp.sphererho);
        dcomp.spherexyz=[x,y,z];
        clear x y z
        % generating the massdist array which describes how a particle 
        % will distribute mass across the shells 
        for idcomp=1:dcomp.numshells.rho
            dcomp.temp=dcomp.spherexyz;
            dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
            dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
            [dcomp.temp,dcomp.edges.rho]=histcounts(vecnorm(dcomp.temp,2,2),dcomp.edges.rho);
            dcomp.massdist.rho(:,idcomp)=dcomp.temp'/size(dcomp.temp,1);
        end
        for idcomp=1:dcomp.numshells.rho
            for idcompaz=1:dcomp.numshells.az
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                dcomp.temp=(Raz*dcomp.temp')';
                dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                [az,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.az]=histcounts(az,dcomp.edges.az);
                dcomp.massdist.az(:,idcomp,idcompaz)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        for idcomp=1:dcomp.numshells.rho
            for idcompel=1:dcomp.numshells.el
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                el=dcomp.centers.el(idcompel);
                Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                dcomp.temp=(Rel*dcomp.temp')';
                dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                [el,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.el]=histcounts(el,dcomp.edges.el);
                dcomp.massdist.el(:,idcomp,idcompel)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        % slabs must be handled separately.
        dcomp.edges.slabs=[0:0.02*S.rp:(2*S.br)/sqrt(6)]';
        dcomp.centers.slabs=dcomp.edges.slabs(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.vols.slabs=(1/sqrt(2))*((dcomp.edges.slabs(2:end,:)*2*sqrt(3/2)).^3-(dcomp.edges.slabs(1:end-1,:)*2*sqrt(3/2)).^3);
        dcomp.numshells.slabs=numel(dcomp.centers.slabs); % number of radial shells
        % initializing the accumulation variables
        dcomp.counts.slabs=zeros(dcomp.numshells.slabs,1);
        dcomp.massdist.slabs=zeros(dcomp.numshells.slabs,dcomp.numshells.az,dcomp.numshells.el,dcomp.numshells.slabs);
        for idcomp=1:dcomp.numshells.slabs
            for idcompaz=1:dcomp.numshells.az
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                for idcompel=1:dcomp.numshells.el
                    dcomp.temp=dcomp.spherexyz;
                    dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                    el=dcomp.centers.el(idcompel);
                    Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                    Rtot=Rel*Raz;
                    dcomp.temp=(Rtot*dcomp.temp')';
                    dcomp.tempfrac=dcomp.temp*S.fcc.invA;
                    is_inside = (dcomp.tempfrac >= 0) & (dcomp.tempfrac < 1);
                    rows_to_keep = all(is_inside, 2);
                    dcomp.tempfrac = dcomp.tempfrac(rows_to_keep, :);
                    dcomp.tempfrac_centered = dcomp.tempfrac - round(dcomp.tempfrac);
                    H_frac = max(abs(dcomp.tempfrac_centered), [], 2);
                    a_prim = 2 * S.br; 
                    R_max = a_prim / sqrt(6); % R_max corresponds to the H_frac = 0.5 boundary
                    H = H_frac * 2 * R_max; % H is the metric distance for each point (N x 1)
                    [dcomp.temp,dcomp.edges.slabs]=histcounts(H,dcomp.edges.slabs);
                    dcomp.massdist.slabs(:,idcompaz,idcompel,idcomp)=dcomp.temp'/size(dcomp.temp,1);
                end
            end
        end
    elseif S.bc==4 % BB
        % bins are calibrated on particle size because the behaviour at 
        % the boundary is dictated by the particle size. Binning must
        % reach where mass can reach so S.br+S.rc+S.rp since we include
        % reals and ghosts. In the case of PBCs we only need to go to
        % the boundary because we use MIC.
        dcomp.edges.rho=[0:0.02*S.rp:S.br]';
        dcomp.edges.az=linspace(0,2*pi,azbins)';
        dcomp.edges.el=linspace(-pi/2,pi/2,elbins)';
        dcomp.centers.rho=dcomp.edges.rho(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.centers.az=dcomp.edges.az(1:end-1,:)+(dcomp.edges.az(2)-dcomp.edges.az(1))/2; % centers of the azimuthal bins
        dcomp.centers.el=dcomp.edges.el(1:end-1,:)+(dcomp.edges.el(2)-dcomp.edges.el(1))/2; % centers of the elevation bins
        dcomp.vols.rho=(4/3)*pi*(dcomp.edges.rho(2:end).^3-dcomp.edges.rho(1:end-1).^3); % volumes of the onion shells
        dcomp.vols.az=ones(numel(dcomp.centers.az),1).*(S.bv/numel(dcomp.centers.az)); % volumes of the azimuthal onion shells
        dcomp.vols.el=((2*pi*S.br^3)/3).*(cos(dcomp.edges.el(2:end,1))-cos(dcomp.edges.el(1:end-1,1))); % volumes of the onion shells
        dcomp.numshells.rho=numel(dcomp.centers.rho); % number of radial shells
        dcomp.numshells.az=numel(dcomp.centers.az); % number of azimuthal shells
        dcomp.numshells.el=numel(dcomp.centers.el); % number of elevation shells
        dcomp.masks.core=dcomp.centers.rho<=S.br; % mask identifying core shells
        dcomp.masks.outercore=dcomp.masks.core & dcomp.centers.rho>(S.br/2); % mask identifying outer half of the domain
        dcomp.masks.mittelcore=dcomp.centers.rho>(S.br/4) & dcomp.centers.rho<(3*S.br/4); % mask identifying the central half of the domain
        dcomp.masks.halo=dcomp.centers.rho>S.br; % mask identifying halo shells
        % initializing the accumulation variables
        dcomp.counts.rho=zeros(dcomp.numshells.rho,1);
        dcomp.counts.az=zeros(dcomp.numshells.az,1);
        dcomp.counts.el=zeros(dcomp.numshells.el,1);
        dcomp.massdist.rho=zeros(dcomp.numshells.rho,dcomp.numshells.rho);
        dcomp.massdist.az=zeros(dcomp.numshells.az,dcomp.numshells.rho,dcomp.numshells.az);
        dcomp.massdist.el=zeros(dcomp.numshells.el,dcomp.numshells.rho,dcomp.numshells.el);
        % generating spherical point cloud to MC simulate the mass
        % distribution. Sphere is radius S.rp
        dcomp.sphereN=1e4;
        dcomp.sphereaz=rand(dcomp.sphereN,1)*2*pi;
        dcomp.sphereel=asin(2*rand(dcomp.sphereN,1)-1);
        dcomp.sphererho=(rand(dcomp.sphereN,1).^(1/3)).*S.rp;
        [x,y,z]=sph2cart(dcomp.sphereaz,dcomp.sphereel,dcomp.sphererho);
        dcomp.spherexyz=[x,y,z];
        clear x y z
        % generating the massdist array which describes how a particle 
        % will distribute mass across the shells 
        for idcomp=1:dcomp.numshells.rho
            dcomp.temp=dcomp.spherexyz;
            dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
            dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
            [dcomp.temp,dcomp.edges.rho]=histcounts(vecnorm(dcomp.temp,2,2),dcomp.edges.rho);
            dcomp.massdist.rho(:,idcomp)=dcomp.temp'/size(dcomp.temp,1);
        end
        for idcomp=1:dcomp.numshells.rho
            for idcompaz=1:dcomp.numshells.az
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                dcomp.temp=(Raz*dcomp.temp')';
                dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                [az,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.az]=histcounts(az,dcomp.edges.az);
                dcomp.massdist.az(:,idcomp,idcompaz)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        for idcomp=1:dcomp.numshells.rho
            for idcompel=1:dcomp.numshells.el
                dcomp.temp=dcomp.spherexyz;
                dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                el=dcomp.centers.el(idcompel);
                Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                dcomp.temp=(Rel*dcomp.temp')';
                dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                [el,~,~]=cart2sph(dcomp.temp(:,1),dcomp.temp(:,2),dcomp.temp(:,3));
                [dcomp.temp,dcomp.edges.el]=histcounts(el,dcomp.edges.el);
                dcomp.massdist.el(:,idcomp,idcompel)=dcomp.temp'/size(dcomp.temp,1);
            end
        end
        % slabs must be handled separately.
        dcomp.edges.slabs=[0:0.02*S.rp:S.br]';
        dcomp.centers.slabs=dcomp.edges.slabs(1:end-1,:)+0.01*S.rp; % centers of the bins
        dcomp.vols.slabs=dcomp.edges.slabs(2:end,:).^3-dcomp.edges.slabs(1:end-1,:).^3;
        dcomp.numshells.slabs=numel(dcomp.centers.slabs); % number of radial shells
        % initializing the accumulation variables
        dcomp.counts.slabs=zeros(dcomp.numshells.slabs,1);
        dcomp.massdist.slabs=zeros(dcomp.numshells.slabs,dcomp.numshells.az,dcomp.numshells.el,dcomp.numshells.slabs);
        for idcomp=1:dcomp.numshells.slabs
            for idcompaz=1:dcomp.numshells.az
                az=dcomp.centers.az(idcompaz);
                Raz=[cos(az), -sin(az), 0; sin(az), cos(az), 0; 0, 0, 1];
                for idcompel=1:dcomp.numshells.el
                    dcomp.temp=dcomp.spherexyz;
                    dcomp.temp(:,1)=dcomp.temp(:,1)+dcomp.centers.rho(idcomp);
                    el=dcomp.centers.el(idcompel);
                    Rel=[cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
                    Rtot=Raz*Rel;
                    dcomp.temp(vecnorm(dcomp.temp,2,2)>S.br,:)=[];
                    dcomp.temp=max(abs((Rtot*dcomp.temp')'),[],2);
                    [dcomp.temp,dcomp.edges.slabs]=histcounts(dcomp.temp,dcomp.edges.slabs);
                    dcomp.massdist.slabs(:,idcompaz,idcompel,idcomp)=dcomp.temp'/size(dcomp.temp,1);
                end
            end
        end
    end
    save(filesave,"dcomp","S")
end