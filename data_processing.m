%% CLEANUP

clear all
close all
rng('shuffle')
warning('off', 'all');

%% SERIES NAME

filenamecollisionseries='SBCvsPBCtest_temp_%d.mat';
plottingenabled=1;

%% ANALYSIS

conditionstoprocess=[];

for ic=355:390
    filename=sprintf(filenamecollisionseries,ic);
    if exist(filename,'file')>0
        conditionstoprocess=[conditionstoprocess,ic];
    end
end

for ic=conditionstoprocess
    clearvars -except ISO conditionstoprocess ic filenamecollisionseries plottingenabled
    filename=sprintf(filenamecollisionseries,ic);
    load(filename,'C','P','V','S','EDGES','PDFT','SSF','DCOMP','DS','H','pdfedges','ASYMCORR','EQUIP')
    
    % --- ORGANIZE DATA ----------------------------------------------
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - organizing data\n', ic, S.bc, S.phi);
    % --- general purpose data
    N=S.N;
    globV=S.bv;
    dens=N/globV;
    kbT=1.38e-23*298.15;
    nosnap=size(PDFT,1); % number of snapshots in PDF
    % ---
    % --- sequences
    for i0=1:nosnap
        AZ(:,i0)=PDFT{i0,1};
        EL(:,i0)=PDFT{i0,2};
        RHO(:,i0)=PDFT{i0,3};
        AZS(:,i0)=PDFT{i0,4};
        ELS(:,i0)=PDFT{i0,5};
    end
    clear PDFT
    % ---
    % --- means
    MEAN_AZ=mean(AZ,2);
    MEAN_EL=mean(EL,2);
    MEAN_RHO=mean(RHO,2);
    MEAN_AZS=mean(AZS,2);
    MEAN_ELS=mean(ELS,2);
    MEAN_SSF100=mean(SSF.SSF100,1);
    MEAN_SSF110=mean(SSF.SSF110,1);
    MEAN_SSF111=mean(SSF.SSF111,1);
    % ---
    % --- moving averages
    MMEAN_AZ=movmean(AZ,100,2);
    MMEAN_AZ=MMEAN_AZ./sum(MMEAN_AZ);
    MMEAN_EL=movmean(EL,100,2);
    MMEAN_EL=MMEAN_EL./sum(MMEAN_EL);
    MMEAN_RHO=movmean(RHO,100,2);
    MMEAN_RHO=MMEAN_RHO./sum(MMEAN_RHO);
    MMEAN_AZS=movmean(AZS,100,2);
    MMEAN_AZS=MMEAN_AZS./sum(MMEAN_AZS);
    MMEAN_ELS=movmean(ELS,100,2);
    MMEAN_ELS=MMEAN_ELS./sum(MMEAN_ELS);
    MMEAN_SSF100=movmean(SSF.SSF100,100,1);
    MMEAN_SSF100=MMEAN_SSF100./sum(MMEAN_SSF100,2);
    MMEAN_SSF110=movmean(SSF.SSF110,100,1);
    MMEAN_SSF110=MMEAN_SSF110./sum(MMEAN_SSF110,2);
    MMEAN_SSF111=movmean(SSF.SSF111,100,1);
    MMEAN_SSF111=MMEAN_SSF111./sum(MMEAN_SSF111,2);
    % ---
    % --- variance of moving averages
    STD_MMEAN_AZ=std(MMEAN_AZ)';
    STD_MMEAN_EL=std(MMEAN_EL)';
    STD_MMEAN_RHO=std(MMEAN_RHO)';
    STD_MMEAN_AZS=std(MMEAN_AZS)';
    STD_MMEAN_ELS=std(MMEAN_ELS)';
    STD_MMEAN_SSF100=std(MMEAN_SSF100,0,2)';
    STD_MMEAN_SSF110=std(MMEAN_SSF110,0,2)';
    STD_MMEAN_SSF111=std(MMEAN_SSF111,0,2)';
    % --------------------------------------------------------------

    % --- VERIFY AND IDENTIFY THERMALIZATION -----------------------
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - verify and identify thermalization\n', ic, S.bc, S.phi);
    % --- find a good lognormal maximum using a coarse binning
    edges_azs=linspace(log10(min(STD_MMEAN_AZS)),log10(max(STD_MMEAN_AZS)),100)';
    [counts_azs,edges_azs]=histcounts(log10(STD_MMEAN_AZS),edges_azs);
    hist_azs=[edges_azs(1:end-1,:),counts_azs'./sum(counts_azs)];
    peakvalue=hist_azs(hist_azs(:,2)==max(hist_azs(:,2)),1);
    peakvalue=peakvalue(1);
    % ---
    % --- fine log10 histogram for windowed lognormal fit
    edges_azs=linspace(log10(min(STD_MMEAN_AZS)),log10(max(STD_MMEAN_AZS)),1000)';
    [counts_azs,edges_azs]=histcounts(log10(STD_MMEAN_AZS),edges_azs);
    hist_azs=[edges_azs(1:end-1,:),counts_azs'./sum(counts_azs)];
    [peakid,edges_azs]=histcounts((peakvalue),edges_azs);
    % ---
    % --- find peak bin
    peakid=[edges_azs(1:end-1,:),peakid'];
    peakid=find(peakid(:,2));
    peakA=hist_azs(peakid,2);
    peakA=peakA(1);
    % ---
    % --- define gaussian fit options
    ft = fittype('A*exp(-(x-mu)^2/(2*sigma^2))', ...
        'independent','x','coefficients',{'A','mu','sigma'});
    opts = fitoptions(ft);
    opts.Lower      = [0.5*peakA, 1.5*peakvalue,    0.01*abs(peakvalue)       ];
    opts.Upper      = [1.5*peakA, 0.8*peakvalue,    1*abs(peakvalue)     ];
    opts.StartPoint = [peakA,      peakvalue, 0.2*abs(peakvalue)];
    % ---
    % --- windowed fit
    xwins=[10:10:1000]';
    xall=hist_azs(:,1);
    yall=hist_azs(:,2);
    SStot = sum((yall - mean(yall)).^2);
    R2best=-1e6;
    bestfit=[];
    for iwin=1:numel(xwins)
        xwin=xwins(iwin);
        bottom=peakid-xwin;
        if bottom<1
            bottom=1;
        end
        top=peakid+xwin;
        if top>size(hist_azs,1)
            top=size(hist_azs,1);
        end
        x=hist_azs(bottom:top,1);
        y=hist_azs(bottom:top,2);
        [cf, gof] = fit(x, y, ft, opts);
        yfit_all  = cf.A * exp(-(xall - cf.mu).^2/(2*cf.sigma^2));
        SSres = sum((yall - yfit_all).^2);        
        R2(iwin) = 1 - SSres/SStot;
        if R2(iwin)>R2best
            bestfit.cf=cf;
            bestfit.gof=gof;
            R2best=R2(iwin);
        end
    end
    if plottingenabled==1
        figure
        plot(xall,yall)
        hold
        plot(xall,bestfit.cf.A * exp(-(xall - bestfit.cf.mu).^2/(2*bestfit.cf.sigma^2)))
        hold
    end
    confints=confint(bestfit.cf);
    thermalrange=[10^(confints(1,2)),10^(confints(2,2))];
    % ---
    % --- find thermalized data
    for idtherm=1:nosnap
        if STD_MMEAN_AZS(idtherm)>thermalrange(1) && STD_MMEAN_AZS(idtherm)<thermalrange(2)
            break
        end
    end
    % ---
    % --- cleanup
    clear peak* opts R2* bottom top x* y* edges* gof cf de iwin i0 ft count_azs SSres SStot
    % ---
    % -------------------------------------------------------------------

    % --- THERMALIZED DATA ----------------------------------------------
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculate aggregated, fully thermalized data\n', ic, S.bc, S.phi);
    % --- thermalized means
    MEAN_AZ=mean(AZ(:,idtherm:end),2);
    MEAN_AZ=MEAN_AZ./sum(MEAN_AZ);
    MEAN_EL=mean(EL(:,idtherm:end),2);
    MEAN_EL=MEAN_EL./sum(MEAN_EL);
    MEAN_RHO=mean(RHO(:,idtherm:end),2);
    MEAN_RHO=MEAN_RHO./sum(MEAN_RHO);
    MEAN_AZS=mean(AZS(:,idtherm:end),2);
    MEAN_AZS=MEAN_AZS./sum(MEAN_AZS);
    MEAN_ELS=mean(ELS(:,idtherm:end),2);
    MEAN_ELS=MEAN_ELS./sum(MEAN_ELS);
    MEAN_SSF100=mean(SSF.SSF100(idtherm:end,:),1);
    MEAN_SSF100=MEAN_SSF100./sum(MEAN_SSF100);
    MEAN_SSF110=mean(SSF.SSF110(idtherm:end,:),1);
    MEAN_SSF110=MEAN_SSF110./sum(MEAN_SSF110);
    MEAN_SSF111=mean(SSF.SSF111(idtherm:end,:),1);
    MEAN_SSF111=MEAN_SSF111./sum(MEAN_SSF111);
    % ---
    % --- plotting
    if plottingenabled==1
        figure
        plot(MEAN_AZS)
        hold
        plot(MEAN_AZ)
        figure
        plot(MEAN_ELS)
        hold
        plot(MEAN_EL)
    end
    % ---
    % -------------------------------------------------------------------

    
    % --- PDF AND SSF ANALYSIS --------------------------------------------------
       
    % --- determining the ideal gas distance distribution by Montecarlo ---
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculating pair distribution function\n', ic, S.bc, S.phi);
    gdenominator=zeros(size(pdfedges{3},1)-1,1);
    noreps=1e4;
    if S.bc==2
        ssf100_ig=zeros(size(krange,1),1);
        ssf111_ig=zeros(size(krange111,1),1);
        for iig=1:noreps
            temp=rand(S.N,3).*(2*S.br)-S.br;
            tempd=pdist(temp);
            [temphc,pdfedges{3}]=histcounts(tempd,pdfedges{3});
            gdenominator=gdenominator+temphc';
            if mod(100*iig/noreps,1)==0
                fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculating pair distribution denominator - percentage complete: %d\n', ic, S.bc, S.phi,100*iig/noreps);
            end
            ssf100=(1/S.N) * abs(sum(exp(-1i * SSF.kvec100 * temp'), 2)).^2;
            ssf010=(1/S.N) * abs(sum(exp(-1i * SSF.kvec010 * temp'), 2)).^2;
            ssf001=(1/S.N) * abs(sum(exp(-1i * SSF.kvec001 * temp'), 2)).^2;
            ssf111=(1/S.N) * abs(sum(exp(-1i * SSF.kvec111 * temp'), 2)).^2;
            ssf110=(1/S.N) * abs(sum(exp(-1i * SSF.kvec110 * temp'), 2)).^2;
            ssf011=(1/S.N) * abs(sum(exp(-1i * SSF.kvec011 * temp'), 2)).^2;
            ssf101=(1/S.N) * abs(sum(exp(-1i * SSF.kvec101 * temp'), 2)).^2;
            ssf100_ig = ssf100_ig+mean([ssf100,ssf010,ssf001],2);
            ssf110_ig = ssf110_ig+mean([ssf110,ssf011,ssf101],2);
            ssf111_ig = ssf111_ig+ssf111; 
        end
        ssf111_ig=ssf111_ig./noreps;
        ssf110_ig=ssf110_ig./noreps;
        ssf100_ig=ssf100_ig./noreps;
    elseif S.bc==1
        temprho=((rand(S.N,noreps)).^(1/3)).*S.br;
        tempaz=rand(S.N,noreps).*2*pi-pi;
        tempel=asin(2.*rand(S.N,noreps)-1);
        ssf100_ig=zeros(size(SSF.kvec100,1),1);
        ssf110_ig=zeros(size(SSF.kvec110,1),1);
        ssf111_ig=zeros(size(SSF.kvec111,1),1);

        for iig=1:noreps
            [tempx,tempy,tempz]=sph2cart(tempaz(:,iig),tempel(:,iig),temprho(:,iig));
            temp=[tempx,tempy,tempz];
            tempd=pdist(temp);
            [temphc,pdfedges{3}]=histcounts(tempd,pdfedges{3});
            gdenominator=gdenominator+temphc';
            if mod(100*iig/noreps,1)==0
                fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculating pair distribution denominator - percentage complete: %d\n', ic, S.bc, S.phi,100*iig/noreps);
            end
            ssf100=(1/S.N) * abs(sum(exp(-1i * SSF.kvec100 * temp'), 2)).^2;
            ssf010=(1/S.N) * abs(sum(exp(-1i * SSF.kvec010 * temp'), 2)).^2;
            ssf001=(1/S.N) * abs(sum(exp(-1i * SSF.kvec001 * temp'), 2)).^2;
            ssf111=(1/S.N) * abs(sum(exp(-1i * SSF.kvec111 * temp'), 2)).^2;
            ssf110=(1/S.N) * abs(sum(exp(-1i * SSF.kvec110 * temp'), 2)).^2;
            ssf011=(1/S.N) * abs(sum(exp(-1i * SSF.kvec011 * temp'), 2)).^2;
            ssf101=(1/S.N) * abs(sum(exp(-1i * SSF.kvec101 * temp'), 2)).^2;
            ssf100_ig = ssf100_ig+mean([ssf100,ssf010,ssf001],2);
            ssf110_ig = ssf110_ig+mean([ssf110,ssf011,ssf101],2);
            ssf111_ig = ssf111_ig+ssf111;
        end
        ssf111_ig=ssf111_ig./noreps;
        ssf110_ig=ssf110_ig./noreps;
        ssf100_ig=ssf100_ig./noreps;
    end
    gdenominator=gdenominator./noreps;
    clear temp tempd temprho tempaz tempel
    % ---
    % --- determine the g
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - average pairs numerator\n', ic, S.bc, S.phi);
    gnumerator=sum(RHO(:,idtherm:end),2)./(size(RHO,2)-idtherm); % average number of pairs in that bin
    g=gnumerator./gdenominator;
    g=[pdfedges{3}(1:end-1),g];
    % ---
    % --- plot g
    if plottingenabled==1
        figure
        plot(g(:,1),g(:,2))
        xlim([2*S.rp 2*S.br+2*S.rc]);
        yline(1)
        xline(2*S.br)
        xline(2*(S.br-S.rp))
    end
    % ---
    % -------------------------------------------------------------------

    % --- COLLATING DATA ------------------------------------------------
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - collating data\n', ic, S.bc, S.phi);
    ISO(ic).g=g;
    ISO(ic).pdfedges=pdfedges;
    ISO(ic).rr_dist_distrib=double([pdfedges{4}(1:end-1),double(DS),double(DS)./sum(double(DS))]);
    ISO(ic).thermalized_mAZ=[pdfedges{1}(1:end-1),MEAN_AZ];
    ISO(ic).thermalized_mAZS=[pdfedges{1}(1:end-1),MEAN_AZS];
    ISO(ic).thermalized_mEL=[pdfedges{2}(1:end-1),MEAN_EL];
    ISO(ic).thermalized_mELS=[pdfedges{2}(1:end-1),MEAN_ELS];
    ISO(ic).thermalized_mRHO=[pdfedges{3}(1:end-1),MEAN_RHO];
    ISO(ic).thermalized_mSSF=[krange,MEAN_SSF'];
    ISO(ic).thermalized_mSSF111=[krange111,MEAN_SSF111'];
    ISO(ic).stdmmAZS=[[1:numel(STD_MMEAN_AZS)]'.*S.kt,STD_MMEAN_AZS];
    ISO(ic).stdmmELS=[[1:numel(STD_MMEAN_ELS)]'.*S.kt,STD_MMEAN_ELS];
    ISO(ic).ssf100_ig=[krange(:,1),ssf100_ig];
    ISO(ic).ssf111_ig=[krange111(:,1),ssf111_ig];
    ISO(ic).thermalizationtime=idtherm*S.kt;
    ISO(ic).thermalizationfit=bestfit;
    ISO(ic).S=S;
    ISO(ic).P=P;
    ISO(ic).V=V;
    ISO(ic).C=C;
    % -------------------------------------------------------------------
    
end