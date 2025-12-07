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

load('SBCvsPBC_temp_25.mat');

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

%% SVD 

%% ================================================================
%   INPUT: MSD is M×U (time × directions)
%   OUTPUT: tau_dom ± CI via block bootstrap
%
%   REQUIREMENTS:
%   - MSD: mean-square displacement projected on U directions
%   - times: M×1 vector of times (monotonic)
%
% =================================================================

%% --- USER INPUTS ------------------------------------------------
% MSD = ...          % your M×U array
% times = ...        % your M×1 vector

[M,U] = size(MSD);
fprintf('Loaded MSD matrix: %d snapshots × %d directions\n', M, U);

%% 1. REMOVE THE MEAN OVER DIRECTIONS (reference MSD)
%    This isolates anisotropy-induced fluctuations
ref_MSD = mean(MSD, 2);       % M × 1
D = MSD - ref_MSD;            % deviations: M × U

%% 2. CENTER OVER TIME (optional but helps SVD stability)
D = D - mean(D,1);            % subtract mean of each direction

%% 3. SVD DECOMPOSITION ------------------------------------------
%    D = Umat * S * V'
%    Left singular vectors (Umat) contain TEMPORAL MODES
[Umat,Smat,~] = svd(D, 'econ');

singvals = diag(Smat);
fprintf('Singular value spectrum (first five):\n');
disp(singvals(1:min(5,end)));

% Extract dominant temporal mode (with amplitude)
dominant_time_mode = Umat(:,1) * Smat(1,1);    % M×1

%% 4. FIT DECAY TIME τ ON DOMINANT TEMPORAL MODE ------------------
% Robust exponential fit of autocorrelation envelope

function tau = fit_tau(x)
    % Remove linear trend just in case
    x = detrend(x);

    Mloc = length(x);
    maxlag = floor(Mloc/4);

    % Normalized autocorrelation
    ac = xcov(x, maxlag, 'coeff');
    acpos = ac(maxlag+1:end);      % lag 0..maxlag
    lags = (0:maxlag)';

    % Use only positive part for log fit
    idx = find(acpos > 0.05);      
    if numel(idx) < 5
        tau = NaN; return;
    end
    idx = idx(2:end);   % skip lag=0

    % Fit log(acf) ~ -lag/tau
    p = polyfit(lags(idx), log(acpos(idx)), 1);
    tau = -1/p(1);
end

tau_dom = fit_tau(dominant_time_mode);
fprintf('tau (dominant mode) = %.4f\n', tau_dom);


%% 5. BLOCK BOOTSTRAP FOR CONFIDENCE INTERVAL --------------------
B = 1000;                     % bootstrap replicates
block = round(max(10, sqrt(M)));   % block length
nb   = ceil(M/block);         % total blocks
tau_boot = zeros(B,1);

rng(0);                       % reproducibility

for b = 1:B
    % --- Resample directions (with replacement)
    dir_idx = randi(U, U, 1);
    Dsample = D(:, dir_idx);   % M×U

    % --- Time block bootstrap (robust)
    idx_time = [];
    
    while numel(idx_time) < M
        bid = randi(nb);                  % pick a random block
        t1 = (bid-1)*block + 1;
        t2 = min(bid*block, M);
        idx_time = [idx_time, t1:t2];
    end
    
    idx_time = idx_time(1:M);             % now safe: idx_time is >= M


    Db = Dsample(idx_time, :);

    % --- SVD of bootstrap sample
    [Ub,Sb,~] = svd(Db - mean(Db,2), 'econ');
    dom_b = Ub(:,1) * Sb(1,1);

    % --- Fit tau
    tau_boot(b) = fit_tau(dom_b);
end

tau_CI = prctile(tau_boot, [2.5 50 97.5]);

fprintf('\nBootstrap Results (dominant temporal mode):\n');
fprintf('tau median = %.4f\n', tau_CI(2));
fprintf('tau 95%% CI = [%.4f , %.4f]\n', tau_CI(1), tau_CI(3));


%% 6. OPTIONAL: DIAGNOSTIC PLOTS --------------------------------
figure; 
subplot(2,1,1);
plot(times, dominant_time_mode, 'LineWidth',1.5);
xlabel('time'); ylabel('mode amplitude');
title('Dominant Temporal Mode (from SVD)');

subplot(2,1,2);
histogram(tau_boot, 40);
xlabel('\tau'); ylabel('count');
title('Bootstrap distribution of decay time \tau');


%%

% Assume:
% Umat, Smat, Vmat are from [Umat,Smat,Vmat] = svd(D, 'econ')
% times is Mx1
% directions is U x 3 (each row is a unit vector used to produce MSD columns)
% and MSD_original was MxU
[Umat,Smat,Vmat] = svd(D, 'econ');
singvals = diag(Smat);
% energy per singular value = sigma^2
energy = singvals.^2;
total_energy = sum(energy);
frac = energy / total_energy;
cumfrac = cumsum(frac);

fprintf('Explained variance fractions (first 10):\n');
disp(frac(1:min(10,end))');
fprintf('Cumulative (first 10):\n');
disp(cumfrac(1:min(10,end))');

% Plot spectrum
figure;
subplot(3,2,1);
semilogy(1:length(singvals), singvals, 'o-','LineWidth',1.2);
xlabel('mode index n'); ylabel('\sigma_n (log scale)');
title('Singular value spectrum');

subplot(3,2,2);
plot(1:length(frac), frac*100, 'o-','LineWidth',1.2); hold on;
plot(1:length(frac), cumfrac*100, 'r--','LineWidth',1);
xlabel('mode index n'); ylabel('percent energy');
legend('per-mode %','cumulative %','Location','best');
title('Explained variance per mode');

% Show first 3 temporal modes (left singular vectors * singular value)
nplot = min(3, size(Umat,2));
for n = 1:nplot
    mode_time = Umat(:,n) * Smat(n,n);  % amplitude-weighted time series
    subplot(3,2,2+n);
    plot(times, mode_time, 'LineWidth',1.2);
    xlabel('time'); ylabel(sprintf('mode %d', n));
    title(['Temporal mode ' num2str(n)]);
end

% Directional fingerprints (V columns) for first 3 modes
figure;
for n = 1:nplot
    subplot(nplot,1,n);
    bar(Vmat(:,n));
    xlabel('direction index'); ylabel(sprintf('V(:,%d)',n));
    title(['Directional fingerprint of mode ' num2str(n)]);
end

% Correlate mode fingerprints with canonical directions (100,110,111)
u100 = [1 0 0];
u110 = [1 1 0]/norm([1 1 0]);
u111 = [1 1 1]/norm([1 1 1]);
directions =u;
% compute projection of each mode's fingerprint on cos^2 with those axes
for n=1:nplot
    g = Vmat(:,n); % directional coefficients
    % form the projection of directions onto axis: (u dot dir)^2 - 1/3 (or raw)
    proj100 = (directions * u100').^2;
    proj110 = (directions * u110').^2;
    proj111 = (directions * u111').^2;
    corr100 = corr(g, proj100);
    corr110 = corr(g, proj110);
    corr111 = corr(g, proj111);
    fprintf('mode %d correlations: corr(g,proj100)=%.3f, proj110=%.3f, proj111=%.3f\n', n, corr100, corr110, corr111);
end
