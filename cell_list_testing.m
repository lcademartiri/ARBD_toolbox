clear all
%% cell list testing

S.rp=1e-8;
S.stdx=6.2793e-10;
S.pot_sigma=2*S.rp+0.25e-9;
S.rc=3*S.pot_sigma;

% parameters
phis=logspace(log10(1e-3),log10(1e-1),3)';
Ns=logspace(log10(1e3),log10(1e5),5)';
clms=(1:3)';

c=[];
q=1;
for i1=1:3
    for i2=1:5
        for i3=1:3
            phi=phis(i1);
            N=Ns(i2);
            clm=clms(i3);
            bv=(N*(4/3)*pi*S.rp^3)/phi;
            br=((3*bv)/(4*pi))^(1/3);
            cl=S.rc/clm+1e-16;
            c(q,:)=[S.rp,S.stdx,N,phi,br,bv,S.rc,cl,clm];
            q=q+1;
        end
    end
end
noconds=q-1;
clear phis Ns clms N phi clm br bv cl q i*

for ic=1:noconds
    S.N=c(ic,3);
    S.br=c(ic,5);
    S.cl=c(ic,8);
    S.icl=1/S.cl;
    S.clm=c(ic,9);
    DISP=build_noise_library(S.stdx,1e6);
    DISPcl=single(DISP./S.cl);
    % create cells padded to 2^n for morton indexing
    maxr=S.br+S.rc;
    cellcenters=(S.cl/2:S.cl:maxr)';
    if max(cellcenters-S.cl/2)<maxr
        cellcenters(end+1,:)=cellcenters(end,:)+S.cl;
    end
    padded_grid_size=2^nextpow2(numel(cellcenters)*2);
    padding_cells=(padded_grid_size/2)-numel(cellcenters);
    padding=(max(cellcenters):S.cl:padding_cells*S.cl+max(cellcenters))';
    cellcenters=[cellcenters;padding(2:end,1)];
    cellcenters=sort(unique([-cellcenters;cellcenters]));
    celledges=[cellcenters(1)-S.cl/2;cellcenters+S.cl/2];
    % morton indices
    morton_grid_dim=[padded_grid_size,padded_grid_size,padded_grid_size];
    num_morton_cells = prod(morton_grid_dim);
    [X, Y, Z] = ndgrid(0:morton_grid_dim(1)-1, 0:morton_grid_dim(2)-1, 0:morton_grid_dim(3)-1);
    morton_idx = zeros(size(X));
    for b = 0:9
        morton_idx = bitor(morton_idx, bitshift(bitget(X, b+1), 3*b + 0));
        morton_idx = bitor(morton_idx, bitshift(bitget(Y, b+1), 3*b + 1));
        morton_idx = bitor(morton_idx, bitshift(bitget(Z, b+1), 3*b + 2));
    end
    % neighbor list
    [dX, dY, dZ] = meshgrid(-S.clm:S.clm, -S.clm:S.clm, -S.clm:S.clm);
    offsets = [dX(:), dY(:), dZ(:)];
    no=size(offsets,1);
    NeighborTable = zeros(num_morton_cells, no, 'int32'); 
    for ix=padding_cells:padded_grid_size-padding_cells
        for iy=padding_cells:padded_grid_size-padding_cells
            for iz=padding_cells:padded_grid_size-padding_cells
                mi=morton_idx(ix,iy,iz);
                for io=1:no
                    nx = ix + offsets(io, 1);
                    ny = iy + offsets(io, 2);
                    nz = iz + offsets(io, 3);
                    neigh_morton_idx = morton_idx(nx, ny, nz);
                    NeighborTable(mi, io) = neigh_morton_idx;
                end                
            end
        end
    end
    vi=[find(abs(cellcenters-S.cl/2)<maxr,1)-1,find((cellcenters+S.cl/2)>maxr,1)]; % physical span of reachable cells
    morton_idx_phys=morton_idx(vi(1):vi(2),vi(1):vi(2),vi(1):vi(2));
    vmi=morton_idx_phys(:);
    % cleanup and organization
    cc_cl=single(cellcenters/S.cl); % cell center coordinates in single reduced dims
    ce_cl=single(celledges/S.cl); % cell center coordinates in single reduced dims
    vind=vi; % physical span of reachable cells
    vcc_cl=cc_cl(vind(1):vind(2),1); % cell center coordinates in single reduced dims
    vce_cl=ce_cl(vind(1):vind(2),1); % cell center coordinates in single reduced dims
    vmi=sort(vmi); % sorted list of physically reachable morton indices
    pmmi=morton_idx_phys; % map between coordinate indices and morton indices
    pmmig=size(pmmi,1);
    nl=NeighborTable; % neighbor list in morton indices;
    clear NeighborTable morton_idx_phys celledges cellcenters vi dX dY dZ X Y Z padding* io ix iy iz b neigh_morton_idx

    % generate starting particle configuration
    p=fillSphere(S.N,S.br,S.rp);
    p_cl=single(p.*S.icl);
    pcli=floor(p_cl);
    pcll=p_cl-pcli;
    offset_to_physical = 1+pmmig/2; % From your setup code
    ix = pcli(:,1) + offset_to_physical;
    iy = pcli(:,2) + offset_to_physical;
    iz = pcli(:,3) + offset_to_physical;
    lin_idx = ix + (iy-1)*pmmig + (iz-1)*pmmig^2;
    p_morton_id = pmmi(lin_idx);
end