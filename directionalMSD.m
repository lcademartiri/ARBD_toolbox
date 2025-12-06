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

mask = squeeze(POS(1,1,:) ~= 0);
POS = POS(:,:,mask);
snaps=size(POS,3);

% unwrap
if S.bc==2
    for in=1:size(POS,1)
        p=squeeze(POS(in,:,:))';
        for is=2:snaps
            displacement(is-1,:) = p(is,:) - p(is-1,:);
            displacement(is-1,:) = displacement(is-1,:) - 2*S.br * round(displacement(is-1,:) / (2*S.br)); % Minimum Image Convention for displacement
        end
        rho=vecnorm(displacement,2,2);
        ptemp=cumsum(displacement);
        p=[p(1,:);p(1,:)+ptemp];
        POS(in,:,:)=p';        
        disp(in)
    end
elseif S.bc==3
    [N_particles, n_dims, N_snaps] = size(POS);
    d = diff(POS, 1, 3);
    d_flat = reshape(permute(d, [1, 3, 2]), [], 3); 
    d_frac = d_flat * S.fcc.invA';
    d_frac = d_frac - round(d_frac);
    d_real = d_frac * S.fcc.A';
    d_corrected = permute(reshape(d_real, N_particles, N_snaps-1, 3), [1, 3, 2]);
    POS_unwrapped = zeros(size(POS));
    POS_unwrapped(:,:,1) = POS(:,:,1); % Initial positions
    POS_unwrapped(:,:,2:end) = POS(:,:,1) + cumsum(d_corrected, 3);
    POS = POS_unwrapped;
elseif S.bc==1
    for in=1:size(POS,1)        
        p=squeeze(POS(in,:,:))';
        for is=2:snaps
            displacement(is-1,:) = p(is,:) - p(is-1,:);
            if norm(displacement(is-1,:))>S.stdx*10
                displacement(is-1,:) = S.stdx*(displacement(is-1,:)./norm(displacement(is-1,:)));
            end            
        end
        rho=vecnorm(displacement,2,2);
        ptemp=cumsum(displacement);
        p=[p(1,:);p(1,:)+ptemp];
        POS(in,:,:)=p';    
        disp(in)
    end
end

% COM correction
for is=1:snaps
    ptemp=squeeze(POS(:,:,is));
    ptemp=ptemp-mean(ptemp);
    POS(:,:,is)=ptemp;
end
p0=squeeze(POS(:,:,1));
%% random directions

randu=500;
[x,y,z]=sph2cart(rand(randu,1)*2*pi-pi,asin(2*rand(randu,1)-1),1);
u=[x,y,z];
u(1,:)=[1,0,0];
u(2,:)=[0,1,0];
u(3,:)=[0,0,1];
u(4,:)=[1,1,0]./norm([1,1,0]);
u(5,:)=[0,1,1]./norm([0,1,1]);
u(6,:)=[1,0,1]./norm([1,0,1]);
u(7,:)=[1,1,1]./norm([1,1,1]);

%% MSD

for is=2:snaps
    dp=squeeze(POS(:,:,is))-p0;
    for iu=1:randu
        MSD(is-1,iu)=mean((dp*u(iu,:)').^2);
    end
    if mod(is,100)==0,disp(is); end
end

MSD100=mean(MSD(:,1:3),2);
MSD110=mean(MSD(:,4:6),2);
MSD111=mean(MSD(:,7),2);
MSDa=mean(MSD(:,8:end),2);
MSDn=[MSD100,MSD110,MSD111]./MSDa;
MSDrSTD=std(MSD(:,8:end),0,2)./MSDa;
figure
plot(MSDn(:,:),'DisplayName','MSDn(:,:)',YDataSource = 'MSDn(:,:)');
ylabel("MSDn(:,:)");
title("MSDn(:,:)");
legend("show");
figure
plot(MSDrSTD)