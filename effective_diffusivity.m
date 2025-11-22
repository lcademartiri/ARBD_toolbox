function CONDS=effective_diffusivity(data_folder,CONDS,P,C)

%% RECRUITING DATABASE

% finding all files with the right prefix in the data_folder
filestoconsider=dir(fullfile(data_folder,'effectivediffs_adaptive*'));

%% COLLATE CONDITIONS REQUIRING UNIQUE EFFECTIVE DIFFUSIVITIES

% unique diffusivities imply unique particle sizes, particle densities,
% solvent densities, temperatures, solvent viscosities, kuhn multipliers
templist=zeros(size(CONDS.rp,1),7);
for ic=1:size(CONDS.rp,1)
    templist(ic,:)=[CONDS.rp(ic),CONDS.rhop(ic),CONDS.rhos(ic),C.kbT,CONDS.eta(ic),P.kuhnmultiplier,P.kuhnmultiplierVACF]; 
end
[i1,i2,i3]=unique(templist,'rows');
CONDS.alpha=zeros(size(CONDS.rp,1),1);

%% IDENTIFY CONDITIONS THAT MUST BE CALCULATED

effdiffstocalculate=true(size(i2,1),1);
for i0=1:size(i2,1) % loop over all unique effective diffusivities
    ic=i2(i0);
    for ifile=1:size(filestoconsider,1) % loop over all relevant files in the database
        tempname=filestoconsider(ifile).name;
        % load the S struct
        Stemp=load(tempname,'S');
        % condition that checks whether the variables defining the
        % effective diffusivities are the ones required for the condition
        % under consideration
        if Stemp.S.rp==CONDS.rp(ic) && Stemp.S.rhoparticle==CONDS.rhop(ic) && Stemp.S.rhosolvent==CONDS.rhos(ic) && Stemp.S.kbT==C.kbT && Stemp.S.nu==CONDS.eta(ic) && Stemp.S.km==P.kuhnmultiplierVACF
            % if the file is correct, assign the alpha value to the CONDS
            % struct array and eliminate from list of "tocalculate" and
            % stop considering other files
            load(tempname,'DEFF');
            CONDS.alpha(i3==i3(ic),1)=DEFF(1,2);
            effdiffstocalculate(i0,1)=false;
            break
        else
            effdiffstocalculate(i0,1)=true;
        end            
    end
end

%% CALCULATE MISSING ALPHA VALUES

for i0=1:size(i2,1)
    ic=i2(i0);
    if effdiffstocalculate(i0,1)==true
        filenamedeffsgen='effectivediffs_adaptive_%d_%s_%s_%d_%d.mat';
        filenamedeffs=sprintf(filenamedeffsgen,P.kuhnmultiplierVACF,CONDS.solvent{ic},CONDS.phase{ic},CONDS.rp(ic),P.ratio);        
        [DEFF,S,rfsteps]=effdiff4(CONDS,ic,P,C,P.ratio);
        DEFF(1,2)=DEFF/S.diffusivity; % DEFF matrix whereby the first column is the corrected diffusivity while the second column is the ratio with the Einstein Stokes diffusivity
        CONDS.alpha(i3==i3(ic),1)=DEFF(1,2);
        save(filenamedeffs,'DEFF','S','rfsteps');
    end
end

end