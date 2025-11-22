function axv = excludedVolume(S,P,p)

% this function calculates the accessible volume in the simulation volume
% by creating a number of random points and looking at how many of them are
% in the excluded volume. The output is the accessible volume fraction.  

    if S.bc==1
        el=(pi/2)-acos(2.*rand(P.noevpoints,1)-1);
        az=(2*pi).*rand(P.noevpoints,1);
        rho=rand(P.noevpoints,1).*S.br;
        [x,y,z]=sph2cart(az,el,rho);
        D=((x-p(:,1)').^2+(y-p(:,2)').^2+(z-p(:,3)').^2);
        idxev=sum(D>(2*S.rp)^2,2)==2;
    elseif S.bc==2 || S.bc==3
        evp=rand(P.noevpoints,3).*(2*S.br)-S.br; 
        D=((evp(:,1)-p(:,1)').^2+(evp(:,2)-p(:,2)').^2+(evp(:,3)-p(:,3)').^2);
        idxev=sum(D>(2*S.rp)^2,2)==S.N;
    end
    axv=sum(idxev)/length(idxev);
    
end