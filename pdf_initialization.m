function PDF=pdf_initialization(S,P)
    % pdf edges for the azimut
    PDF.pdfedges{1}=linspace(-pi,pi,P.pdfbins(1)+1)';
    % pdf edges for the elevation
    PDF.pdfedges{2}=linspace(-pi/2,pi/2,P.pdfbins(2)+1)';            
    % pdf edges for the distances between reals
    if S.bc==1 % SBC
        maxdist=2*S.br;
        if 10*S.rp>(maxdist-2*S.rc)
            PDF.pdfedges{3}=sort((maxdist:-0.02*S.rp:0)'); % from 0 to box diameter
        else
            OS=sort((maxdist:-0.02*S.rp:(maxdist-2*S.rc))');
            IS=(0:0.02*S.rp:10*S.rp)';
            MS=linspace(max(IS),min(OS),ceil((min(OS)-max(IS))/S.rp))';
            PDF.pdfedges{3}=unique(sort([OS;MS;IS]));
        end
        
    elseif S.bc==2 % PBC cubic
        maxdist=sqrt(3)*S.br;
        paralleldist=S.br;
        if 10*S.rp>(paralleldist-2*S.rc)
            PDF.pdfedges{3}=sort((maxdist:-0.02*S.rp:0)'); % from 0 to 111 diagonal
        else
            OS=linspace(paralleldist,maxdist,ceil((maxdist-paralleldist)/S.rp))';
            IS=(0:0.02*S.rp:10*S.rp)';
            BS=(paralleldist:-0.02*S.rp:(paralleldist-2*S.rc))';
            MS=linspace(max(IS),min(BS),ceil((min(BS)-max(IS))/S.rp))';
            PDF.pdfedges{3}=unique(sort([OS;MS;BS;IS]));
        end
        
    elseif S.bc==3 % PBC FCC
        maxdist=sqrt(2)*S.br;
        paralleldist=S.br;
        if 10*S.rp>(paralleldist-2*S.rc)
            PDF.pdfedges{3}=sort((maxdist:-0.02*S.rp:0)'); % from 0 to long diagonal
        else
            OS=linspace(paralleldist,maxdist,ceil((maxdist-paralleldist)/S.rp))';
            IS=(0:0.02*S.rp:10*S.rp)';
            BS=(paralleldist:-0.02*S.rp:(paralleldist-2*S.rc))';
            MS=linspace(max(IS),min(BS),ceil((min(BS)-max(IS))/S.rp))';
            PDF.pdfedges{3}=unique(sort([OS;MS;BS;IS]));
        end
    end
    PDF.centers{1}=PDF.pdfedges{1}(1:end-1,1)+diff(PDF.pdfedges{1})./2;
    PDF.centers{2}=PDF.pdfedges{2}(1:end-1,1)+diff(PDF.pdfedges{2})./2;
    PDF.centers{3}=PDF.pdfedges{3}(1:end-1,1)+diff(PDF.pdfedges{3})./2;
end