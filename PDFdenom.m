function gdenominator=PDFdenom(S,PDF,noreps)

    gdenominator=zeros(size(PDF.pdfedges{3},1)-1,1);
    if S.bc==2
        for iig=1:noreps
            temp=rand(S.N,3).*(2*S.br)-S.br;
            tempd=pdist(temp);
            [temphc,PDF.pdfedges{3}]=histcounts(tempd,PDF.pdfedges{3});
            gdenominator=gdenominator+temphc';
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
    end
    gdenominator=gdenominator./noreps;
    clear temp tempd temprho tempaz tempel
end