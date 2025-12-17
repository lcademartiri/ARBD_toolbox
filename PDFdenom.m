function gdenominator=PDFdenom(S,PDF,noreps,data_folder)

    if S.bc==1
        filepdfdenom = sprintf('PDFdenom_SBC_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
    elseif S.bc==2
        filepdfdenom = sprintf('PDFdenom_PBCc_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
    elseif S.bc==3
        filepdfdenom = sprintf('PDFdenom_PBCFCC_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
    end

    gdenominator=zeros(size(PDF.pdfedges{3},1)-1,1);
    if S.bc==2
        for iig=1:noreps
            temp=rand(S.N,3).*(2*S.br)-S.br;
            tempd=mic_all_pair_displacements(temp, S);
            [temphc,PDF.pdfedges{3}]=histcounts(vecnorm(tempd,2,2),PDF.pdfedges{3});
            gdenominator=gdenominator+temphc'./2;
            if mod(100*iig/noreps,1)==0
                fprintf('boundary condition: %d -  phi: %.3f - calculating pair distribution denominator - percentage complete: %d\n', S.bc, S.phi,100*iig/noreps);
            end
        end
    elseif S.bc==1
        temprho=((rand(S.N,noreps)).^(1/3)).*S.br;
        tempaz=rand(S.N,noreps).*2*pi-pi;
        tempel=asin(2.*rand(S.N,noreps)-1);
        for iig=1:noreps
            [tempx,tempy,tempz]=sph2cart(tempaz(:,iig),tempel(:,iig),temprho(:,iig));
            temp=[tempx,tempy,tempz];
            tempd=pdist(temp);
            [temphc,PDF.pdfedges{3}]=histcounts(tempd,PDF.pdfedges{3});
            gdenominator=gdenominator+temphc';
            if mod(100*iig/noreps,1)==0
                fprintf('boundary condition: %d -  phi: %.3f - calculating pair distribution denominator - percentage complete: %d\n', S.bc, S.phi,100*iig/noreps);
            end
        end
    elseif S.bc==3
        for iig=1:noreps
            temp=rand(S.N,3);
            temp=temp*S.fcc.A;
            tempd=mic_all_pair_displacements(temp, S);
            [temphc,PDF.pdfedges{3}]=histcounts(vecnorm(tempd,2,2),PDF.pdfedges{3});
            gdenominator=gdenominator+temphc'./2;
            if mod(100*iig/noreps,1)==0
                fprintf('boundary condition: %d -  phi: %.3f - calculating pair distribution denominator - percentage complete: %d\n', S.bc, S.phi,100*iig/noreps);
            end
        end
    end
    
    gdenominator=gdenominator./noreps;

    save([data_folder,'/',filepdfdenom],"gdenominator",'S')

    clear temp tempd temprho tempaz tempel
end