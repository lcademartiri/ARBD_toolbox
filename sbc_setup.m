function [p,pgp,ASYMCORR,ASYMFIT,RADCOMPFIT,FLOWSbinned]=sbc_setup(S,PDF,nsteps)
    
    if S.potential==1
        
        filenamecorrectionseries='ASYMCORR_lj_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat';
        filenamestartingconfigurationseries='START_lj_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat';
    elseif S.potential==2
        filenamepdfdenomseries='PDFdenom_wca_%.0e_%.0e_%.0f.mat';
        filenamecorrectionseries='ASYMCORR_wca_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat';
        filenamestartingconfigurationseries='START_wca_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat';
    elseif S.potential==0
        filenamepdfdenomseries='PDFdenom_hs_%.0e_%.0e_%.0f.mat';
        filenamecorrectionseries='ASYMCORR_hs_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat';
        filenamestartingconfigurationseries='START_hs_%.0e_%.0e_%.0f_%.1f_%.1e_%.0e.mat';
    end
    filenamepdfdenomseries='PDFdenom_%.0e_%.0e_%.0f.mat';
    filepdfdenom=sprintf(filenamepdfdenomseries,S.rp,S.phi,S.N);
    filecorrection=sprintf(filenamecorrectionseries,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma,nsteps);
    filestartingconfiguration=sprintf(filenamestartingconfigurationseries,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma,nsteps);
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
        gdenominator=PDFdenom(S,PDF,1e5);
        save(filepdfdenom,'gdenominator','S')
    end
    % ---------------------------------------------------------------------

    if exist(filecorrection,'file')
        load(filecorrection,"ASYMCORR","ASYMFIT","RADCOMPFIT","FLOWSbinned")
    else
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
        onionedges=[(S.br:-0.02*S.rp:0)';(S.br:0.02*S.rp:S.br+4*S.rp)'];
        onionedges=sort(unique(onionedges));
        onionshells=size(onionedges,1)-1;
        edgeflows=linspace(0,10*S.stdx,1000)';
        edgeflows=sortrows(unique([-edgeflows;edgeflows]));
        FLOWS(:,1)=edgeflows(1:end-1)+mean(diff(edgeflows))/2;
        FLOWS(:,2:onionshells+1)=0;
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
            ndens.av_window=1000;
            ndens.edges=sort((S.br:-0.02*S.rp:0)');
            ndens.centers = ndens.edges(1:end-1,:) + 0.01*S.rp;
            ndens.counts=zeros(numel(ndens.centers),ndens.av_window);
            ndens.vols = (4/3)*pi*(ndens.edges(2:end).^3 - ndens.edges(1:end-1).^3);
            ndens.ndens0=(S.N/S.bv);
            ndens.f1=figure;
            pdf.pre.counts=zeros(numel(PDF.pdfedges{3})-1,1);
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
        end
        % -------------------------------------------------------------------------
        if debugging
            rng(100)
        end
    
        while thermflag==0 | qs<nsteps    
            qs=qs+1;

            % --- assessing thermalization ------------------------------------
            if thermflag==0
                r2mean = mean(prho.^2);             % mean squared radius
                spread_ratio = r2mean / r2_uniform;
                SR(qs,1)=spread_ratio;
                if spread_ratio>0.99
                    thermflag=1;
                    qs=1;
                    disp('configuration fully expanded - starting collection of radial displacement data - restarting counter')
                end
                if qs>1000
                    if mean(SR(qs-500:qs))<=mean(SR(qs-1000:qs-500))
                        thermflag=1;
                        qs=1;
                        disp('configuration fully expanded - starting collection of radial displacement data - restarting counter')
                    end
                end
            end
            % -----------------------------------------------------------------
            
            % --- attempt displacements -----
            if S.potential~=0
                idxgp=prho>(S.br-S.rc);
                ptemp=[p;pgp(idxgp,:)];
                % calculate the displacement components due to potentials for all particles (reals and active ghosts)
                [potdisps,pairdists]=potential_displacements(ptemp,size(ptemp,1),S.rc,H,S.esdiff,clamp,S.kbT,S.stdx,S.timestep,1);
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
                FLOWS(:,2:end) = FLOWS(:,2:end) + counts;
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
            if thermflag==1 & graphing
                [tempndc,ndens.edges]=histcounts(prho,ndens.edges);
                ndens.counts(:,mod(qs,ndens.av_window)+1)=tempndc';
                [temppdfcounts,PDF.pdfedges{3}]=histcounts(nonzeros(pairdists(1:S.N,1:S.N)),PDF.pdfedges{3});
                pdf.pre.counts=pdf.pre.counts+temppdfcounts'./2;
            end
            % ---------------------------------------------------------------

            % --- counter ----
            if mod(qs,1e3)==0
                if thermflag==0
                    fprintf('configuration expansion: %.0e steps\n', qs);
                else
                    fprintf('configuration evolution and radial displacement asymmetry determination: %.0e steps of %.0e \n', qs, nsteps);
                    set(groot, 'CurrentFigure', ndens.f1);
                    ndens.av_counts=mean(ndens.counts,2);
                    ndens.ndens=ndens.av_counts./ndens.vols;
                    plot(ndens.centers,ndens.ndens);
                    xlim([S.br/10 S.br])
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
                    yline(ndens.ndens0, '--w', 'LineWidth', 1);
                    drawnow

                    pdf.pre.numerator=pdf.pre.counts/qs;
                    tempg=pdf.pre.numerator./gdenominator;
                    tempg=[PDF.centers{3},tempg];
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
                    drawnow
                end
                
            end
            % ----------------
        end
        if debugging
            save('debuggingdata-sbcsetup.mat',"S","DEBCOUNT","POS",'-v7.3')
            return
        end
        
        %% --- FITTING THE RADIAL DISPLACEMENT COMPONENTS
        disp('analyzing radial displacement statistics to elaborate correction')
        % --- Set up fittype and options ---
        ft = fittype( 'A*(mu*(2/pi)*(fwhm./(4*(x-xc).^2+fwhm.^2))+(1-mu)*(sqrt(4*log(2))/(sqrt(pi)*fwhm)).*exp(-(4*log(2)/fwhm^2)*(x-xc).^2))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [0 0 0 -Inf]; % order of parameters A fwhm mu xc
        opts.Robust = 'Bisquare';
        opts.StartPoint = [0.17 0.5 0.75 0]; % order of parameters A fwhm mu xc
        opts.Upper = [1 Inf 1 Inf]; % order of parameters A fwhm mu xc
        
        % --- main loop ---------
        binning=[1]; % binning of the distances to change statistics
        xData=FLOWS(:,1)./1e-9; % renormalization of the displacement to facilitate fits
        for i0=1:numel(binning)
            ib=binning(i0);
            % initialize binned FLOWS array
            FLOWStemp=[];
            % set first column as the centers of the radial displacement bins
            FLOWStemp(:,1)=FLOWS(:,1);
            % indices for binning columns together
            indices=[1:ib:size(FLOWS,2)];
            % miniloop to bin data together
            qi=1;
            for ibins=1:numel(indices)-1
                FLOWStemp(:,qi+1)=sum(FLOWS(:,indices(ibins)+1:indices(ibins+1)),2);
                qi=qi+1;
            end
            % define number of distance bins depending on the binning used
            % above
            onionshellstemp=onionshells/ib;
            % define edges of distance bins
            onionedgestemp=linspace(0,S.br+4*S.rp,onionshellstemp+1)';
            % define centers of distance bins
            onioncenterstemp=onionedges(1:end-1,:)+(onionedgestemp(2)-onionedgestemp(1))/2;
            % initialize FITRESULTS array with column 1 indicating shell
            % centers
            RADCOMPFIT{i0,1}(:,1)=onioncenterstemp;
            % main fitting loop
            for io=2:size(FLOWStemp,2)
                % define y data
                yData=FLOWStemp(:,io)./sum(FLOWStemp(:,io)); % normalize to a distribution        
                if sum(yData)>0 % if there are any counts
                    [fitresult, gof] = fit( xData, yData, ft, opts );
                    fitvalues=coeffvalues(fitresult);
                    fit99ci=confint(fitresult,0.99);
                    % FITRESULTS lists the 4 fit parameters sequentially by
                    % indicating first the -99CI value, the best fit value, and
                    % the +99CI value. Last column is the R square
                    RADCOMPFIT{i0,1}(io,2:14)=[fit99ci(1,1),fitvalues(1),fit99ci(2,1),fit99ci(1,2),fitvalues(2),fit99ci(2,2),fit99ci(1,3),fitvalues(3),fit99ci(2,3),fit99ci(1,4),fitvalues(4),fit99ci(2,4),gof.rsquare];            
                end
                % output counter
                if mod(io,10)==0
                    fprintf('fitting radial displacement distributions: binning %.0f, distance %.0f of %.0f \n', ib,io,size(FLOWStemp,2));
                end
            end
            % convert all data into distributions and compile into FLOWSbinned
            FLOWStemp(:,2:end)=FLOWStemp(:,2:end)./sum(FLOWStemp(:,2:end));
            FLOWSbinned{i0,1}=FLOWStemp;
        end
        
        %% --- FITTING THE RADIAL DISPLACEMENT ASYMMETRY -----
        
        disp('extracting and plotting fit of radial displacement asymmetry')
        % taking the unbinned data - evolution of the center of the
        % distribution from 2*S.rp from the boundary
        asym=[RADCOMPFIT{1,1}(:,[1,12])];
        S.correctionwindow=S.rc;
        
        idxcw2=asym(:,1)>(S.br-S.correctionwindow) & asym(:,1)<(S.br);
        asym(asym(:,1)>S.br,:)=[];
        idxcw=asym(:,1)>(S.br-S.correctionwindow) & asym(:,1)<(S.br);
        % asym(asym(:,1)<(S.br-S.correctionwindow) | asym(:,1)>(S.br),:)=[];
        s101=sgolayfilt(asym(:,2),3,101);
        s51=sgolayfilt(asym(:,2),3,51);
        s51=sgolayfilt(s51,3,51);
        s101=s101(idxcw,1);
        s51=s51(idxcw,1);
        asym=asym(idxcw,:);
        meancomp=mean(s101(asym(:,1)>(S.br-S.rc) & asym(:,1)<(S.br-2*S.rp),:));
        ramp=linspace(0,meancomp*2,sum(asym(:,1)>(S.br-S.rc) & asym(:,1)<(S.br-2*S.rp)))';
        f1=linspace(0,numel(s51)-1,numel(s51))'.*ramp(2);
        f2=s51;
        idxgap=asym(:,1)>(S.br-2*S.rp) & asym(:,1)<(S.br-1.5*S.rp);
        ra=find(idxgap==1,1);
        rb=find(idxgap==1,1,'last');
        t=(asym(:,1)-asym(ra,1))/(asym(rb,1)-asym(ra,1));
        t = max(0, min(1, t));   % clamp to [0,1]
        s = 3*t.^2 - 2*t.^3;
        F = (1 - s).*f1 + s.*f2;

        % asym=[RADCOMPFIT{1,1}(:,[1,11,12,13])];
        % asym(asym(:,1)<(S.br-S.correctionwindow) | asym(:,1)>(S.br),:)=[];
        % % renormalize the x data for facilitating fit
        % x=(asym(:,1)-asym(1,1))/1e-9;
        % 
        % % defining fit options - fit is a double exponential
        % ft = fittype( 'A*x*(1-exp((x-x0)/l1)).*exp((x-xc)/l2)', 'independent', 'x', 'dependent', 'y' );
        % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        % opts.Display = 'Off';
        % opts.Robust = 'Bisquare';
        % opts.Upper =        [max(asym(:,3))*2   75      75      max(x)          max(x)*10]; % parameter order: A l1 l2 x0 xc
        % opts.StartPoint =   [max(asym(:,3))/2   1       10      max(x)*(3/4)    max(x)]; % parameter order: A l1 l2 x0 xc
        % opts.Lower =        [1e-6               1e-6    1e-6    max(x)/2        1]; % parameter order: A l1 l2 x0 xc
        % 
        % % fit model to data.
        % [fitresult, gof] = fit( x, asym(:,3), ft, opts );
        % fittedtemp=fitresult(x);
        % asym(:,5)=fittedtemp;
        figure
        scatter(asym(:,1),asym(:,2),'filled')
        hold
        plot(asym(:,1),F)
        title('Radial Component of Displacements vs Distance from Origin'); 
        xlabel('distance from origin');                         
        ylabel('mean radial component of the displacement');                         
        set(gcf, 'Color', 'k');
        set(gca, 'Color', 'k');
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        set(get(gca, 'Title'), 'Color', 'w');
        set(get(gca, 'XLabel'), 'Color', 'w');
        set(get(gca, 'YLabel'), 'Color', 'w');
        set(gca, 'LineWidth', 3);
        set(findall(gca, 'Type', 'line'), 'LineWidth', 2);
        
        ASYMCORR=RADCOMPFIT{1,1}(:,1);
        ASYMCORR(:,2)=0;
        ASYMCORR(idxcw2,2)=F.*1e-9;
        % ASYMCORR(:,2)=ASYMCORR(:,2).*1e-9;
        % ASYMFIT.fit=fitresult;
        % ASYMFIT.gof=gof;
        
        save(filecorrection,"ASYMCORR","S","FLOWSbinned","RADCOMPFIT")
    end

    %% GENERATE STARTING THERMALIZED CONFIGURATION ------------------------
    
    flag=1;
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

    % ----- MAIN LOOP ---------------------------------------------------------
    qs=0;
    prho = vecnorm(p,2,2);
    % define initial ghosts as antipodal images
    pgp=p-(2*S.br).*(p./vecnorm(p,2,2)); 

    % ----- CALCULATE THE LJ FORCE --------------------------------------------
    if S.potential~=0
        H=pot_force(S.potential,S.rc,30000,S.pot_sigma,S.pot_epsilon);
        clamp=mcdClamp(1e6,S.rp,normrnd(0,S.stdx,1e6,3),S.esdiff,S.timestep,H,S.kbT);    
    end
    % -------------------------------------------------------------------------

    % GRAPHING SETUP --------------------------------------------------
    if graphing
        ndens.edges=sort((S.br:-0.02*S.rp:0)');
        ndens.centers = ndens.edges(1:end-1,:) + 0.01*S.rp;
        ndens.counts=zeros(numel(ndens.centers),1);
        ndens.vols = (4/3)*pi*(ndens.edges(2:end).^3 - ndens.edges(1:end-1).^3);
        ndens.ndens0=(S.N/S.bv);
        ndens.f1=figure;
        pdf.post.counts=zeros(numel(PDF.pdfedges{3})-1,1);
        pdf.post.f1=figure;
    end
    % -----------------------------------------------------------------

    while qs<nsteps
        qs=qs+1;

        % --- calculate radial displacement correction near the boundary ---
        dr=zeros(size(p,1),3);
        % Mask particles within correctionwindow of the boundary
        mask = prho > (S.br - S.correctionwindow);
        if any(mask)
            % Interpolate correction value for each masked particle
            delta_r = interp1(ASYMCORR(:,1), ASYMCORR(:,2), prho(mask), 'pchip', 0);
            % Compute unit radial vectors
            rho_hat = p(mask,:) ./ prho(mask);
            % Apply correction only to the radial component of the displacement
            dr(mask,:) = dr(mask,:) + delta_r .* rho_hat;
        end
                
        % --- attempt displacements -----
        if S.potential~=0
            idxgp=prho>(S.br-S.rc);
            ptemp=[p;pgp(idxgp,:)];
            % calculate the displacement components due to potentials for all particles (reals and active ghosts)
            [potdisps,pairdists]=potential_displacements(ptemp,size(ptemp,1),S.rc,H,S.esdiff,clamp,S.kbT,S.stdx,S.timestep,1);
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
            % tangential move to the ghosts
            pgp2_temp = pgp + v_tan; 
            p2=p+randn(S.N,3)*S.stdx+potdisps;
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
        p2=p2+dr; % ADDING CORRECTION
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
        % ---------------------------------------------------------------------
    
        % --- promotions and demotions ----------------------------------------
        % compare norms of tentative reals and tentative ghosts
        idxrgswap=p2rho>vecnorm(pgp2,2,2);
        % initialize start position array for the next timestep 
        p=p2;
        pgp=pgp2;
        % swap positions according to mask above
        pgp(idxrgswap,:)=p2(idxrgswap,:);
        p(idxrgswap,:)=pgp2(idxrgswap,:);
        % update norms and versors
        prho = vecnorm(p,2,2);
        % ---------------------------------------------------------------------

        % --- characterization of number densities -----------------------
        if graphing & qs>1e4
            [tempndc,ndens.edges]=histcounts(prho,ndens.edges);
            ndens.counts=ndens.counts+tempndc';
            [temppdfcounts,PDF.pdfedges{3}]=histcounts(nonzeros(pairdists(1:S.N,1:S.N)),PDF.pdfedges{3});
            pdf.post.counts=pdf.post.counts+temppdfcounts'./2;
        end
        % ---------------------------------------------------------------

        % --- counter ----
        if mod(qs,1e3)==0
            fprintf('starting configuration evolution and thermalization: %.0e steps of %.0e \n', qs, nsteps/10);
            set(groot, 'CurrentFigure', ndens.f1);
            ndens.av_counts=(ndens.counts)/(qs-1e4);
            ndens.ndens=ndens.av_counts./ndens.vols;
            plot(ndens.centers,ndens.ndens);
            xlim([S.br/10 S.br])
            title('Radial Number Density - corrected'); % <-- Creates the Title object
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
            yline(ndens.ndens0, '--w', 'LineWidth', 1);
            drawnow

            pdf.post.numerator=pdf.post.counts/(qs-1e4);
            tempg=pdf.post.numerator./gdenominator;
            tempg=[PDF.centers{3},tempg];
            set(groot, 'CurrentFigure', pdf.post.f1);
            plot(tempg(:,1),tempg(:,2));
            xlim([0 2*(S.br+S.rp)])
            title('Pair Distribution Function - corrected'); % <-- Creates the Title object
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
            drawnow
        end
        % ----------------
    end
    save(filestartingconfiguration,"p","pgp","S")

    % ---------------------------------------------------------------------
end