function GPMAT = ghostparticlematrix()

% GHOST PARTICLE MATRIX 
br=10;
r=1;
GPM=[-1,-1,-1;-1,-1,0;-1,-1,1;-1,0,-1;-1,0,0;-1,0,1;-1,1,-1;-1,1,0;-1,1,1;0,-1,-1;0,-1,0;0,-1,1;0,0,-1;0,0,0;0,0,1;0,1,-1;0,1,0;0,1,1;1,-1,-1;1,-1,0;1,-1,1;1,0,-1;1,0,0;1,0,1;1,1,-1;1,1,0;1,1,1];
GPMsample=GPM.*(br-0.5*r);

for i0=1:size(GPMsample,1)
    clear temp temp2
    W=GPMsample(i0,:)+GPM.*2.*br;
    idxgpm=sum(abs(W)>(br+r),2)>0;    
    temp=GPM;
    temp(idxgpm,:)=[];
    idxgpm=sum(abs(temp),2)==0;
    temp(idxgpm,:)=[];
    temp2(1,:,:)=temp';
    GPMAT(i0,1:3,1:size(temp2,3))=temp2;
end

% the GPMAT matrix has 3dims: [numero della particella fantasma, coordinata
% cartesiana, condizione di uscita dal box della particella reale
% corrispondente]

for i1=1:size(GPMAT,1) % condition index
    for i3=1:size(GPMAT,3) % ghost particle index
        q=0;
        for i2=1:size(GPMAT,2) % cartesian index
            if GPMAT(i1,i2,i3)==0
                q=q+1;
            end
            if q==3
                GPMAT(i1,1,i3)=nan;
                GPMAT(i1,2,i3)=nan;
                GPMAT(i1,3,i3)=nan;
            end
        end
    end
end

end