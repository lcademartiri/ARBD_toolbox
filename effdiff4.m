function [DEFF,S,displacements]=effdiff4(CONDS,ic,P,C,ratio)


% in this function the alpha parameter (i.e., the correction to the diffusivity 
% that I need to use to model a brownian trajectory as a random walk) is being calculated on the basis of
% the physical parameters (phase of the particle, radius of the particle
% and solvent) and on the basis of the truncations of the VACF. More
% specifically the VACF gets calculated until a certain multiple of the
% tau_r (S.kmVACF, which will correspond to a S.ACFmaxlagVACF steps in the
% physical trajectory). The VACF is also of course calculated within the
% physical trajectory step in order to calculate the MSD associated with a
% single CT step.
% After that we calculate the autocorrelation coefficients which of course
% will be truncated at a maximum lag equal to S.ACFmaxlagVACF
% So, now I have a VACF that correlates up until S.ACFmaxlagVACF steps.
% This means that I have to precalculate S.ACFmaxlagVACF steps in order to
% establish the entire VACF, after which I can start taking the steps as
% representative of the autocorrelated trajectory. This autocorrelated
% physical trajectory will be then subsampled at a certain multiple of the
% tau_r (S.km, which will correspond to a S.ACFmaxlag steps in the
% physical trajectory). The idea is that I do know that I do need at least
% S.ACFmaxlag steps to completely vanish the correlation in the (even
% though this value was actually obtained using VACFs that were truncated
% at 5 tau_r!). The question is whether I need to calculate the VACF over
% the same number of steps over which the trajectory becomes
% decorrelated....

%% CONSOLIDATION OF PARAMETERS

S.rp=CONDS.rp(ic); % particle radius
S.rhoparticle=CONDS.rhop(ic); % particle density
S.mp=CONDS.mp(ic); % particle mass
S.ms=CONDS.ms(ic); % solvent mass
S.rs=CONDS.rs(ic); % solvent radius
S.rhosolvent=CONDS.rhos(ic); % solvent density
S.csol=CONDS.cs(ic); % speed of sound in the solvent calculated assuming 2000MPa modulus
S.nu=CONDS.eta(ic); % dynamic viscosity
S.nukin=S.nu/S.rhosolvent; % kinematic viscosity
mstar=((4/3)*pi*S.rp^3)*S.rhoparticle+0.5*((4/3)*pi*S.rp^3)*S.rhosolvent; % corrected particle mass
S.mpc=mstar; % corrected particle mass    
S.vtherm3D=sqrt(3*C.kbT/mstar); % 3D vRMS particle
S.vtherm=sqrt(S.vtherm3D^2/3); % mean thermal velocity in ONE dimension
S.diffusivity=C.kbT/(6*pi*S.nu*S.rp); % Einstein Stokes diffusivity      
% timescales
S.taur=(2/9)*(S.rhoparticle/S.nu)*S.rp^2; % relaxation time
S.tauc=(9.775e-11/sqrt(C.kbT))*((S.rs/S.rp)^2)*((sqrt(S.ms*S.mp))/(sqrt(S.ms)+sqrt(S.mp))); % collisional time
if S.taur/S.tauc>ratio
    S.tauc=S.taur/ratio;
end
S.tauv=S.rp^2/S.nukin; % vortex time    
S.taus=(S.rp)/S.csol; % sonic time 
S.gamma=1/S.taur; % correlation decay frequency
S.kmVACF=P.kuhnmultiplierVACF; % kuhn multiplier used for the calculation of the VACF and for the establishment of the autocorrelation
S.km=P.kuhnmultiplier; % kuhn multiplier used for the calculation of the VACF and for the establishment of the autocorrelation
S.dt=S.tauc; % CT timestep
S.kt=S.km*S.taur; % RF timestep
S.ksteps=ceil(S.kt/S.dt); % CT steps in each RF step.
S.rwsteps=P.rfstepsforeffdiffest; % number of random walk steps that I need to collect to extract the alpha parameter
S.nosteps=S.rwsteps*S.ksteps; % number of CT steps that need to be calculated in order to get enough RF steps to elaborate the Deff
S.ACFmaxlagVACF=S.kmVACF*ceil(S.taur/S.tauc); % this is the maximum lag (# of tauc steps) considered during the calculation of the CT VACF
S.ACFmaxlag=S.km*ceil(S.taur/S.tauc); % this is the maximum lag (# of tauc steps) considered during the calculation of the CT VACF
S.kbT=C.kbT;
S.CTsteps=1e6; % number of steps that get kept in memory during trajectory evaluation (this must be greater than S.ksteps)

filenameVACFgen='VACF_%d_%s_%s_%d.mat';
filenameVACFs=sprintf(filenameVACFgen,S.kmVACF,CONDS.solvent{ic},CONDS.phase{ic},CONDS.rp(ic));

%% VACF CALCULATION 

fprintf('%d:%d %s\n', ic, size(CONDS.rp,1), 'effdiff - VACF calculation');
% independent variable of the VACF, i.e., time interval
tlist=linspace(1,S.ACFmaxlagVACF,S.ACFmaxlagVACF)'.*S.dt; % time intervals at which we have to calculate the VACF because they correspond to the CT steps
tlistfull=linspace(0,1,1e4)'.*S.dt; % continuous time intervals between 0 and the CTtimestep at which to calculate the continuous VACF that is needed to extract the MSD of the continuous traj
% integration variable
x=logspace(log10(1e-50),log10(1e50),1e4)'; % integration variable - the 10^6 steps are EXTREMELY important to ensure the beta values are accurate

% shortcut variables
rhor=S.rhoparticle/S.rhosolvent; % density ratio - shortcut
sigma1=(1/9)*(7-4*rhor);
sigma2=(1/81)*(1+2*rhor)^2;
alpha1=1+S.rhosolvent/(2*S.rhoparticle);
alpha2=sqrt(1-0.25/(rhor^2));
c0     = 1/(1+2*rhor);
VACF=zeros(size(tlist,1),1);
% 3) Precompute denominator and sqrt(x)
den = 1 + sigma1*x + sigma2*(x.^2);
sqx = sqrt(x);
dx = diff(x);
% 4) Precompute short-time term st_all
tau_s = S.taus;
expDecay = exp(-alpha1 .* tlist ./ tau_s);
cosTerm   = cos(alpha2 .* tlist ./ tau_s);
sinTerm   = sin(alpha2 .* tlist ./ tau_s);
st_all    = c0 .* expDecay .* (cosTerm - (alpha1/alpha2).*sinTerm);
clear expDecay cosTerm sinTerm

N_t = numel(tlist);
VACF = zeros(N_t,1);

% 6) Loop over tlist and compute trapezoidal integral
for it = 1:N_t
    lambda = tlist(it) / S.tauv;
    % Evaluate integrand f(x)
    fx = (sqx .* exp(-lambda .* x)) ./ den;
    % Non-uniform trapezoidal rule:    
    fx_left  = fx(1:end-1);
    fx_right = fx(2:end);
    I_trap = sum( (fx_left + fx_right) .* dx ) * 0.5;
    % First term
    ft = (2*rhor)/(9*pi) * I_trap;
    % Total VACF
    VACF(it) = ft + st_all(it);
end

% Precompute short-time term st_all
tau_s = S.taus;
expDecay = exp(-alpha1 .* tlistfull ./ tau_s);
cosTerm   = cos(alpha2 .* tlistfull ./ tau_s);
sinTerm   = sin(alpha2 .* tlistfull ./ tau_s);
st_all    = c0 .* expDecay .* (cosTerm - (alpha1/alpha2).*sinTerm);
clear expDecay cosTerm sinTerm

VACFfull=zeros(size(tlistfull,1),1);
for it=1:size(tlistfull,1)
    lambda = tlistfull(it) / S.tauv;
    % Evaluate integrand f(x)
    fx = (sqx .* exp(-lambda .* x)) ./ den;
    % Non-uniform trapezoidal rule:    
    fx_left  = fx(1:end-1);
    fx_right = fx(2:end);
    I_trap = sum( (fx_left + fx_right) .* dx ) * 0.5;
    % First term
    ft = (2*rhor)/(9*pi) * I_trap;
    VACFfull(it,1)=ft+st_all(it);
end

S.VACF=VACF;
S.VACFF=VACFfull.*(S.kbT/S.mpc);
save(filenameVACFs,'S');
clear ft st num denom temp* alpha* sigma* x dx rhor VACFfull

%% MSD CALCULATION

fprintf('%d:%d %s\n', ic, size(CONDS.rp,1),'effdiff - MSD calculation')
for it=1:size(tlistfull,1)
    tempt=tlistfull(it);
    MSD(it,1)=2*sum((tempt-tlistfull(1:it)).*S.VACFF(1:it,1).*(tlistfull(2)-tlistfull(1)));
end
S.sigmax2=MSD(end);
clear MSD
sigmax2=S.sigmax2; % variance of the displacement x dist
sigmax=sqrt(sigmax2); % std of the displacement x dist

%% CALCULATE THE AUTOCORRELATION COEFFICIENTS

fprintf('%d:%d %s\n', ic, size(CONDS.rp,1),'effdiff - calculation of the autocorrelation coefficients')

tempVACF=[1;(S.VACF(1:S.ACFmaxlagVACF-1))];
ACV = S.VACF(1:S.ACFmaxlagVACF)';
beta=levinson(tempVACF,S.ACFmaxlagVACF);
beta=-beta';
beta(1,:)=[];
dbeta=[-1;diff(beta)];
idxbeta=dbeta>0;
idxbeta(1:ceil(length(dbeta)/2))=false;
beta(idxbeta,:)=0;
beta=[beta;0];
clear dbeta
sigma_noise=sqrt(1-ACV(1,1:numel(beta))*beta); % scaling term for the noise to preserve unitary variance 
clear ACM x VACFtemp VACFshort ACV
S.beta=beta;
beta3=pagemtimes(flip(beta),ones(1,3,10));

%% BUILD THE INITIAL SEQUENCE

fprintf('%d:%d %s\n', ic, size(CONDS.rp,1),'effdiff - initialize direction vectors and variables')
% initialize direction vectors and variables
o = zeros(S.ACFmaxlagVACF+1, 3,10);
o(1,:,:)=randn(1,3,10);
pk1=zeros(floor(S.rwsteps*1.05),3,10);

fprintf('%d:%d %s\n', ic, size(CONDS.rp,1),'effdiff - calculate precorrelated physical trajectories')
for istep=2:S.ACFmaxlagVACF  
    eta=sigma_noise*randn(1,3,10);
    o(istep,:,:)=sum(beta3(end-istep+2:end,:,:).*(o(1:istep-1,:,:)),1)+eta;
end
displacements=o(1:end-1,:,:).*sigmax;
p=[cumsum(displacements,1)];
p(end+1,:,:)=zeros(1,3,10);
pk1(1,:,:)=p(end-1,:,:);
kcount1=1;
clear displacements

fprintf('%d:%d %s\n', ic, size(CONDS.rp,1),'effdiff - calculate physical trajectory and extract random flight snapshots')
mcount=0;
istep=istep+1;
o=[o;zeros(S.CTsteps-size(o,1)+S.ACFmaxlag,3,10)];
p=[p;zeros(S.CTsteps-size(p,1)+S.ACFmaxlag,3,10)];
eta=sigma_noise.*randn(S.CTsteps+S.ACFmaxlag, 3,10);

while (kcount1-1)*30<=S.rwsteps
    fluxstep=istep-mcount*S.CTsteps;
    o(fluxstep,:,:)=sum(beta3.*(o(fluxstep-S.ACFmaxlagVACF:fluxstep-1,:,:)),1)+eta(fluxstep,:,:);
    p(fluxstep,:,:)=p(fluxstep-1,:,:)+o(fluxstep,:,:).*sigmax;        
    if mod(fluxstep+1,S.CTsteps+1+S.ACFmaxlag)==0
        mcount=mcount+1;
        o(1:S.ACFmaxlag,:,:)=o(end-S.ACFmaxlag+1:end,:,:);
        o(S.ACFmaxlag+1:end,:,:)=zeros(S.CTsteps,3,10);
        p(1:S.ACFmaxlag,:,:)=p(end-S.ACFmaxlag+1:end,:,:);
        p(S.ACFmaxlag+1:end,:,:)=zeros(S.CTsteps,3,10);
        eta=sigma_noise.*randn(S.CTsteps+S.ACFmaxlag, 3,10);
    end
    if mod(istep,floor(S.ksteps/1))==0
        kcount1=kcount1+1;
        tempp=p(fluxstep,:,:);
        pk1(kcount1,:,:)=tempp;
        
        pktemp=pk1;
        pktemp(pktemp(:,1)==0,:,:)=[];
        displacementstemp=diff(pktemp,1);
        vartemp=(std(displacementstemp,0,"all"))^2;
        DEFFtemp=vartemp/(2*S.kt);
        DEFFtemp=DEFFtemp/S.diffusivity;

        disp([sprintf('%.0f',S.rp*1e10),sprintf('%.0f',S.km),sprintf('%.0f',((kcount1-1)*30)),sprintf('%.5f',DEFFtemp)])
    end
    istep=istep+1;
end
pk1(pk1(:,1)==0,:,:)=[];
displacements=diff(pk1,1);
var=(std(displacements,0,"all"))^2;
DEFF=var/(2*S.kt);

