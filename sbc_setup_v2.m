function [p,pgp,ASYMCORR,RADCOMPFIT,FLOWSbinned]=sbc_setup_v2(S,PDF,nsteps0)
    damping=0.01; % damping on the iterations corrections
    scalinggain=false; % switch for policy of scaling integral gain with SNR
    scalingsteps=false; % switch for policy of scaling nsteps with SNR
    maxsteps=1e6; % maximum number of steps per iteration
    integral_gain0=0.1;
    integral_gain=integral_gain0;
    gnumerator=[];
    nsteps=nsteps0;
    seriesname='damping1-ig10';
    if S.potential==1
        potname='lj';
    elseif S.potential==2
        potname='wca';
    elseif S.potential==0
        potname='hs';
    end
    filenamecorrectionseries=['ASYMCORR_',seriesname,'_',potname,'_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat'];
    filenamecorrectionseriestemp=['ASYMCORR_',seriesname,'_',potname,'_temp_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat'];
    filenamestartingconfigurationseries=['START_',potname,'_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat'];
    filenamepdfdenomseries='PDFdenom_%.0e_%.0e_%.0f.mat';
    filepdfdenom=sprintf(filenamepdfdenomseries,S.rp,S.phi,S.N);
    filecorrection=sprintf(filenamecorrectionseries,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma,nsteps0);
    filecorrectiontemp=sprintf(filenamecorrectionseriestemp,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma,nsteps0);
    filestartingconfiguration=sprintf(filenamestartingconfigurationseries,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma,nsteps0);
    debugging=false;
    graphing=true;
    if exist(filecorrection,'file') && exist(filestartingconfiguration,'file')
        load(filecorrection)
        load(filestartingconfiguration)
        return
    end

    %% --- CALCULATION PDF DENOMINATOR ------------------------------------
    if exist(filepdfdenom,'file')
        load(filepdfdenom,'gdenominator');
    else
        gdenominator=PDFdenom(S,PDF,nsteps0);
        save(filepdfdenom,'gdenominator','S')
    end
    % ---------------------------------------------------------------------

    % if exist(filecorrection,'file')
    %     load(filecorrection,"ASYMCORR","RADCOMPFIT","FLOWSbinned")
    % else
    %% --- CREATING THE INITIAL FCC LATTICE -----------------------------------
    disp('creating the initial fcc lattice')
    flag=1;
    if debugging
        rng(100);
    end
    while flag==1
        basis=[0,0.7071,0.7071;0.7071,0,0.7071;0.7071,0.7071,0].*(2.01*S.rp); % basis vectors of the fcc lattice
        maxsteps=2*ceil(((S.br*2)*sqrt(3))/(2*S.rp)); % calculate the maximum number of unit cells of the lattice we need to exceed the dimensions of the boundary
        templist=double(linspace(-maxsteps,maxsteps,2*maxsteps+1)'); % 1D coefficient set
        [x1,x2,x3] = meshgrid(templist,templist,templist); % create lattice of coefficients
        x1=x1(:); % create lattice of coefficients
        x2=x2(:); % create lattice of coefficients
        x3=x3(:); % create lattice of coefficients
        possiblepositions = x1.*basis(1,:)+x2.*basis(2,:)+x3.*basis(3,:); % multiplying the lattice of coefficients by the basis vectors to obtain the real lattice
        possiblepositions=bsxfun(@plus,possiblepositions,-S.rp*[1,1,1]./vecnorm([1,1,1],2)); % displacing the whole lattice to ensure maximum packing within the boundary
        % --- eliminate all lattice points outside 'br'-'r' ---
        tempnorms=vecnorm(possiblepositions,2,2);
        idxstart=tempnorms>(S.br-2*S.rp);   
        possiblepositions(idxstart,:)=[];
        % ---
        possiblepositions=possiblepositions(randperm(size(possiblepositions,1)),:); % randomize the order of the lattice nodes in the list
        p=possiblepositions; % pick the first N nodes in the list
        % --- double checking that particles in the nodes that were picked would not be touching ---
        D=pdist(p)';
        idxcoll=D<(2*S.rp);
        if sum(idxcoll)==0 % if all checks out tap out.
            flag=0;
        end
        % ---
    end 
    clear possiblepositions tempnorms x1 x2 x3 templist maxsteps theta u ct mct st basis rotmat D flag idxcoll idxstart
    p=p(1:S.N,:);    
    % -------------------------------------------------------------------------
    
    %% --- BROWNIAN EVOLUTION -------------------------------------------------
    disp('evolving initial configuration to determine starting configuration and radial displacement asymmetry')
    
    % --- INITIALIZATION ----------------------------------------------
    maxdist=S.br;
    OS=(maxdist:-0.02*S.rp:maxdist-S.rc)';
    IS=(maxdist-S.rc:-S.rp:0)';
    onionedges=sort(unique([IS;OS]));
    onionshells=size(onionedges,1)-1;
    edgeflows=linspace(0,10*S.stdx,1000)';
    edgeflows=sortrows(unique([-edgeflows;edgeflows]));
    FLOWS(:,1)=edgeflows(1:end-1)+mean(diff(edgeflows))/2;
    FLOWS(:,2:onionshells+1)=0;
    FLOWSID=FLOWS(:,1);
    FLOWS(:,1)=[];
    % ---------------------------------------------------------------------

    % DEBUGGING BLOCK ------------------------------------------------
    if debugging
        debug.edges=(0:0.02*S.rp:S.br+S.rc+4*S.rp)';
        debug.centers = debug.edges(1:end-1,:) + 0.01*S.rp;
        debug.counts=debug.centers.*0;
        debug.vols = (4/3)*pi*(debug.edges(2:end).^3 - debug.edges(1:end-1).^3);
        debug.ndens0=(S.N/S.bv);
    end
    % -----------------------------------------------------------------

    % GRAPHING SETUP --------------------------------------------------
    if graphing
        ndens.av_window=nsteps;
        ndens.edges=sort((S.br:-0.02*S.rp:0)');
        ndens.centers = ndens.edges(1:end-1,:) + 0.01*S.rp;
        ndens.counts=zeros(numel(ndens.centers),1);
        ndens.counts_storage=ndens.counts;
        ndens.vols = (4/3)*pi*(ndens.edges(2:end).^3 - ndens.edges(1:end-1).^3);
        ndens.ndens0=(S.N/S.bv);
        ndens.f1=figure;
        pdf.pre.counts=zeros(numel(PDF.pdfedges{3})-1,1);
        pdf.pre.counts_storage=pdf.pre.counts;
        pdf.pre.f1=figure;
    end
    % -----------------------------------------------------------------
    
    % ----- MAIN LOOP ---------------------------------------------------------
    qs=0;
    thermflag=0;
    r2_uniform = 3/5 * S.br^2;       % expected ⟨r²⟩ for uniform density in a sphere
    prho = vecnorm(p,2,2);
    % define initial ghosts as antipodal images
    pgp=p-(2*S.br).*(p./vecnorm(p,2,2)); 
    
        % DEBUGGING BLOCK ----------------------------------------------------
        if debugging
            POS{1,1}=p;
            POS{1,2}=pgp;
            DEBCOUNT=zeros(size(debug.edges,1)-1,1);        
        end
        % --------------------------------------------------------------------

    % ----- CALCULATE THE LJ FORCE --------------------------------------------
    if S.potential~=0
        H=pot_force(S.potential,S.rc,30000,S.pot_sigma,S.pot_epsilon);
        clamp=mcdClamp(1e6,S.rp,normrnd(0,S.stdx,1e6,3),S.esdiff,S.timestep,H,S.kbT);
        % Create the interpolant object once
        H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'pchip');
    end
    % -------------------------------------------------------------------------
    if debugging
        rng(100)
    end
    pdfmaxdev=10;
    qsm=0;
    iteration=1;
    reverseStr = '';
    while thermflag==0 | pdfmaxdev>0.01   
        qs=qs+1;
        iterationprogress=(qs-qsm)/nsteps;

        % --- assessing thermalization ------------------------------------
        if thermflag==0
            r2mean = mean(prho.^2);             % mean squared radius
            spread_ratio = r2mean / r2_uniform;
            SR(qs,1)=spread_ratio;
            if spread_ratio>0.99
                thermflag=1;
                qs=1;
                disp('configuration fully expanded - starting collection of radial displacement data')
                fprintf([seriesname,'_',potname,'\n']);
            end
            if qs>1000
                if mean(SR(qs-500:qs))<=mean(SR(qs-1000:qs-500))
                    thermflag=1;
                    qs=1;
                    disp('configuration fully expanded - starting collection of radial displacement data')
                    fprintf([seriesname,'_',potname,'\n']);
                end
            end
        end
        % -----------------------------------------------------------------

        % --- calculate radial displacement correction near the boundary ---
        if thermflag==1 & qs>nsteps0
            dr=zeros(size(p,1),3);
            % Mask particles within correctionwindow of the boundary
            mask = prho > (S.br - S.correctionwindow);
            if any(mask)
                % Interpolate correction value for each masked particle
                delta_r = interp1(ASYMCORR.correction(:,1), ASYMCORR.correction(:,2), prho(mask), 'pchip', 0);
                % Compute unit radial vectors
                rho_hat = p(mask,:) ./ prho(mask);
                % Apply correction only to the radial component of the displacement
                dr(mask,:) = dr(mask,:) + delta_r .* rho_hat;
            end
        end
        
        % --- attempt displacements -----
        if S.potential~=0
            idxgp=prho>(S.br-S.rc);
            ptemp=[p;pgp(idxgp,:)];
            % calculate the displacement components due to potentials for all particles (reals and active ghosts)
            potdisps=potential_displacements_v2(ptemp,size(ptemp,1),S.rc,H, H_interpolant,S.esdiff,clamp,S.kbT,S.stdx,S.timestep,1);
            % extract potential displacements for active ghosts
            potdispsgp=potdisps(S.N+1:end,:);
            % extract potential displacements for reals
            potdisps=potdisps(1:S.N,:);
            % calculate versors of initial ghost positions
            pgp_norm = vecnorm(pgp,2,2);
            pgp_norm(pgp_norm==0) = eps;
            pgp_dir = pgp ./ pgp_norm;
            % calculate brownian displacement for all ghosts
            v_rand = randn(S.N,3) * S.stdx;
            % add the potential displacements to the active ghosts
            v_rand(idxgp,:)=v_rand(idxgp,:)+potdispsgp;
            % extract radial component to the total ghost displacement
            v_rad_component = sum(v_rand .* pgp_dir, 2);          % scalar per particle
            % subtract radial component from the 3D displacement
            v_tan = v_rand - (v_rad_component .* pgp_dir);       % tangential component
            if debugging
                POS{qs,3}=v_tan;
            end
            % tangential move to the ghosts
            pgp2_temp = pgp + v_tan; 
            p2=p+randn(S.N,3)*S.stdx+potdisps;
            if debugging
                POS{qs,4}=p2;
            end
        else
            v_rand = randn(S.N,3)*S.stdx;
            % extract radial component to the total ghost displacement
            v_rad_component = sum(v_rand .* pgp_dir, 2);          % scalar per particle
            % subtract radial component from the 3D displacement
            v_tan = v_rand - (v_rad_component .* pgp_dir);       % tangential component
            % tangential move to the ghosts
            pgp2_temp = pgp + v_tan;
            % apply brownian displacement to reals
            p2=p+randn(S.N,3)*S.stdx;
        end
        if thermflag==1 & qs>nsteps0
            p2=p2+dr; % ADDING CORRECTION
        end
        % ---
        
        % project ghosts onto the proper distances from their reals
        p2rho = vecnorm(p2,2,2);
        % get norms of the tentative ghost positions
        pgp2_temp_norm = vecnorm(pgp2_temp,2,2);
        % (guard for tiny norms)
        pgp2_temp_norm(pgp2_temp_norm==0) = eps;
        % get versors of tentative ghost positions
        pgp2_dir = pgp2_temp ./ pgp2_temp_norm;
        % set distance of tentative ghosts to the tether distance from tentative real positions
        pgp2 = pgp2_dir .* (2*S.br-p2rho);
        if debugging
            POS{qs,5}=pgp2;
        end
        % after the antipodal setup the pgp list is updated  independently with the exception of the radial component
    
        % HS collision resolution through motion reset
        if S.potential==0
            distpp=pdist2(p2,[p2;pgp2],'squaredeuclidean');
            idxd=distpp>0 & distpp<(2*S.rp)^2;
            row=find(sum(idxd,2)>0);
            col=find(sum(idxd,1)>0);
            col=col(col>S.N)-S.N;
            resetters=[row;col'];
            p2(resetters,:)=p(resetters,:);
            pgp2(resetters,:)=pgp(resetters,:);
        end
        
        % --- characterization of radial displacements -----------------------
        if thermflag==1
            % Precompute norms once
            pnorms = vecnorm(p, 2, 2);
            % Assign each particle to a shell index (1..onionshells)
            shellIdx = discretize(pnorms, onionedges(:,1));
            % Keep only particles inside defined shells
            valid = ~isnan(shellIdx);
            shellIdx = shellIdx(valid);
            % Compute displacements
            bounddisps = p2(valid,:) - p(valid,:);
            % Compute unit vectors along radial direction
            pnorms2 = vecnorm(p(valid,:), 2, 2);
            pvers = p(valid,:) ./ pnorms2;
            % Radial components (positive inward)
            radcomp = -sum(bounddisps .* pvers, 2);
            % Bin radial components into flow bins
            binIdx = discretize(radcomp, edgeflows);
            % Keep only in-range bins
            valid2 = ~isnan(binIdx);
            binIdx = binIdx(valid2);
            shellIdx = shellIdx(valid2);
            % 2D accumulation: [radial-bin, shell]
            counts = accumarray([binIdx, shellIdx], 1, [numel(edgeflows)-1, onionshells]);
            % Update total flows
            FLOWS = FLOWS + counts;
        end
        % ---------------------------------------------------------------------

        % --- promotions and demotions ----------------------------------------
        % compare norms of tentative reals and tentative ghosts
        idxrgswap=p2rho>vecnorm(pgp2,2,2);
        if debugging
            POS{qs,6}=idxrgswap;
        end
        % initialize start position array for the next timestep 
        p=p2;
        pgp=pgp2;
        if debugging
            POS{qs,7}=[p,pgp];
        end
        % swap positions according to mask above
        pgp(idxrgswap,:)=p2(idxrgswap,:);
        p(idxrgswap,:)=pgp2(idxrgswap,:);
        % update norms and versors
        prho = vecnorm(p,2,2);
        % ---------------------------------------------------------------------

            % DEBUGGING BLOCK ------------------------------------------------
            if debugging
                pgpdens=pgp;
                pgpdens(vecnorm(p,2,2)<(S.br-S.rc),:)=[];
                if debugging
                    POS{qs,8}=pgpdens;
                end
                [tempcounts,debug.edges]=histcounts(vecnorm([p;pgpdens],2,2),debug.edges);
                DEBCOUNT(:,qs)=tempcounts';
            end
            % -----------------------------------------------------------------
        
        % --- characterization of number densities and PDF -----------------------
        if thermflag==1 && graphing && mod(qs,10)==0 
            [tempndc,ndens.edges]=histcounts(prho,ndens.edges);
            ndens.counts=ndens.counts+tempndc';
            pairdists=pdist(p);
            [temppdfcounts,PDF.pdfedges{3}]=histcounts(pairdists,PDF.pdfedges{3});
            pdf.pre.counts=pdf.pre.counts+temppdfcounts';
            pdf.pre.numerator=pdf.pre.counts/(mod(qs,ndens.av_window)/10+0.001);
            tempg=pdf.pre.numerator./gdenominator;
            tempg=[PDF.centers{3},tempg];
            pdf.pdf=tempg;
            pdfmaxdev=tempg(tempg(:,1)>S.br/2,2);
            pdfmaxdev=max(abs(pdfmaxdev(1:end-10,1)-1));

            if mod(qs,ndens.av_window)==0
                save(filecorrectiontemp,"ndens","pdf")
                ndens.counts_storage(:,end+1)=ndens.counts;
                ndens.counts=ndens.counts.*0;
                pdf.pre.counts_storage(:,end+1)=pdf.pre.counts;
                pdf.pre.counts=pdf.pre.counts.*0;
            end
        end
        % ---------------------------------------------------------------

        if qs>=nsteps0 && mod(qs-qsm,nsteps)==0
            fprintf('\n');
            qsm=qs;
            iteration=iteration+1;
            FLOWSforcorr=[FLOWSID,FLOWS];
            [ASYMCORR.update,SNR,RADCOMPFIT,FLOWSbinned]=correctAsymm(S,FLOWSforcorr,onionshells,onionedges);
            if qs~=nsteps0
                % Interpolate update correction onto the ASYMCORR radius grid
                ASYMCORR.interp = interp1(ASYMCORR.update(:,1), ASYMCORR.update(:,2), ASYMCORR.correction(:,1), 'pchip', 0);
                ASYMCORR.correction(:,2) = (1-damping).*ASYMCORR.correction(:,2) + integral_gain.*ASYMCORR.interp;
                ASYMCORR.storage(:,end+1)=ASYMCORR.correction(:,2);
                if scalingsteps & iteration>2
                    nsteps=nsteps*(SNR0/SNR)^2;
                    nsteps(nsteps>maxsteps)=maxsteps;
                end
                if scalinggain
                    integral_gain=integral_gain*SNR/SNR0;
                    integral_gain(integral_gain<2*damping)=2*damping;
                end
                save(filecorrectiontemp,"ASYMCORR","FLOWSID","FLOWS","onionshells","onionedges","-append")
            else
                if scalingsteps & iteration==2
                    SNR0=SNR;
                end
                ASYMCORR.storage=ASYMCORR.update;
                ASYMCORR.correction=ASYMCORR.update;
            end
            FLOWS(:,:)=0;
        end

        % --- counter ----
        if qs==1e3
            fprintf('iterative determination of radial displacement correction \n');
            fprintf('configuration expansion \n');
            
        end
        if mod(qs,1e3)==0
            if thermflag==1
                msg = sprintf('iteration %.0f - progress %.2f - iteration steps %.2e - max pdf deviation %.2f', iteration, iterationprogress*100, nsteps, pdfmaxdev*100);
                fprintf([reverseStr, msg]);
                reverseStr = repmat('\b', 1, length(msg));
                if graphing 
                    set(groot, 'CurrentFigure', ndens.f1);
                    ndens.av_counts=S.N.*(ndens.counts/sum(ndens.counts));
                    ndens.ndens=ndens.av_counts./ndens.vols;
                    ndens.ndensnorm=[ndens.centers,100.*(ndens.ndens./ndens.ndens0)];
                    plot(ndens.centers,100.*(ndens.ndens./ndens.ndens0));
                    xlim([S.br/10 S.br])
                    ylim([80 120])
                    title('Radial Number Density - uncorrected'); % <-- Creates the Title object
                    xlabel('distance');                         % <-- Creates the XLabel object
                    ylabel('number density');                         % <-- Creates the YLabel object
                    set(gcf, 'Color', 'k');
                    set(gca, 'Color', 'k');
                    set(gca, 'XColor', 'w');
                    set(gca, 'YColor', 'w');
                    set(get(gca, 'Title'), 'Color', 'w');
                    set(get(gca, 'XLabel'), 'Color', 'w');
                    set(get(gca, 'YLabel'), 'Color', 'w');
                    set(gca, 'LineWidth', 3);
                    set(findall(gca, 'Type', 'line'), 'LineWidth', 2);
                    xline(S.br-S.rp, '--w', 'LineWidth', 1);
                    xline(S.br-2*S.rp, '--w', 'LineWidth', 1);
                    xline(S.br-S.rc, '--w', 'LineWidth', 1);
                    yline(100, '--w', 'LineWidth', 1);
                    drawnow

                    
                    set(groot, 'CurrentFigure', pdf.pre.f1);
                    plot(tempg(:,1),tempg(:,2));
                    xlim([0 2*(S.br+S.rp)])
                    title('Pair Distribution Function - uncorrected'); % <-- Creates the Title object
                    xlabel('distance');                         % <-- Creates the XLabel object
                    ylabel('PDF');                         % <-- Creates the YLabel object
                    set(gcf, 'Color', 'k');
                    set(gca, 'Color', 'k');
                    set(gca, 'XColor', 'w');
                    set(gca, 'YColor', 'w');
                    set(get(gca, 'Title'), 'Color', 'w');
                    set(get(gca, 'XLabel'), 'Color', 'w');
                    set(get(gca, 'YLabel'), 'Color', 'w');
                    set(gca, 'LineWidth', 3);
                    set(findall(gca, 'Type', 'line'), 'LineWidth', 2);
                    yline(1, '--w', 'LineWidth', 1);
                    drawnow
                end
            end
            
        end
        % ----------------
    end
    if debugging
        save('debuggingdata-sbcsetup.mat',"S","DEBCOUNT","POS",'-v7.3')
        return
    end
    save(filecorrection,"ASYMCORR","RADCOMPFIT","FLOWSbinned","gdenominator","ndens","pdf","FLOWS","FLOWSID","onionshells","onionedges")
    save(filestartingconfiguration,"p","pgp","S")
end
% end