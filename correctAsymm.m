function [ASYMCORR,SNR,RADCOMPFIT,FLOWSbinned]=correctAsymm(S,FLOWS,onionshells,onionedges)

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
        
        xData=FLOWS(:,1)./1e-9; % renormalization of the displacement to facilitate fits
        onioncenters = onionedges(1:end-1,:) + 0.5.*diff(onionedges);
        RADCOMPFIT = zeros(onionshells, 14);
        RADCOMPFIT(:,1) = onioncenters;
        fprintf('fitting radial displacement distributions \n');
        for io=2:onionshells
            % define y data
            yData=FLOWS(:,io)./sum(FLOWS(:,io)); % normalize to a distribution        
            if sum(yData)>0 % if there are any counts
                [fitresult, gof] = fit( xData, yData, ft, opts );
                fitvalues=coeffvalues(fitresult);
                fit99ci=confint(fitresult,0.99);
                % FITRESULTS lists the 4 fit parameters sequentially by
                % indicating first the -99CI value, the best fit value, and
                % the +99CI value. Last column is the R square
                RADCOMPFIT(io,2:14)=[fit99ci(1,1),fitvalues(1),fit99ci(2,1),fit99ci(1,2),fitvalues(2),fit99ci(2,2),fit99ci(1,3),fitvalues(3),fit99ci(2,3),fit99ci(1,4),fitvalues(4),fit99ci(2,4),gof.rsquare];            
            end
            % output counter
            if mod(io,10)==0
                
            end
        end
        % convert all data into distributions and compile into FLOWSbinned
        FLOWSbinned=FLOWS;
        total_counts_per_shell = sum(FLOWS(:,2:end), 1);
        FLOWSbinned(:,2:end) = FLOWS(:,2:end) ./ total_counts_per_shell;
        
        %% --- FITTING THE RADIAL DISPLACEMENT ASYMMETRY -----
        
        disp('extracting and plotting fit of radial displacement asymmetry')
        % taking the unbinned data - evolution of the center of the
        % distribution from 2*S.rp from the boundary
        asym=[RADCOMPFIT(:,[1,12])];
        S.correctionwindow=S.rc;
        
        idxcw2=asym(:,1)>(S.br-S.correctionwindow) & asym(:,1)<(S.br);
        asym(asym(:,1)>S.br,:)=[];
        idxcw=asym(:,1)>(S.br-S.correctionwindow) & asym(:,1)<(S.br);
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
        yline(0, '--w', 'LineWidth', 1);
        set(findall(gca, 'Type', 'line'), 'LineWidth', 2);

        % SNR evaluation
        SNR=norm(F)/norm(asym(:,2)-F);
        
        ASYMCORR=RADCOMPFIT(:,1);
        ASYMCORR(:,2)=0;
        ASYMCORR(idxcw2,2)=F.*1e-9;
end