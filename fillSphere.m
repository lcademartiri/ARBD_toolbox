function p=fillSphere(N,R,rp,overlapcheck)
    maxattempts=1e2;
    maxtries=1e1;
    rng('shuffle')
    p=zeros(N,3);
    tries=0;
    eps=1e-16;
    row=1;
    if overlapcheck==1
        while tries<maxtries && ~isempty(row)
            tries=tries+1;
            az=rand(N,1)*2*pi;
            el=asin(2*rand(N,1)-1);
            rho=R*(rand(N,1).^(1/3));
            [x,y,z]=sph2cart(az,el,rho);
            p=[x,y,z];
            dists=pdist2(p,p,'squaredeuclidean');
            [row,~]=find(dists>0 & dists<(2*rp+eps)^2);
            attempt=0;
            while ~isempty(row) & attempt<maxattempts
                attempt=attempt+1;
                p(row,:)=[];
                n=N-size(p,1);
                az=rand(n,1)*2*pi;
                el=asin(2*rand(n,1)-1);
                rho=R*(rand(n,1).^(1/3));
                [x,y,z]=sph2cart(az,el,rho);
                p=[p;x,y,z];
                dists=pdist2(p,p,'squaredeuclidean');
                [row,~]=find(dists>0 & dists<(2*rp+eps)^2);
            end
        end
        if tries==maxtries
            disp('failed to find non-overlapping configuration');
        end
    else
        az=rand(N,1)*2*pi;
        el=asin(2*rand(N,1)-1);
        rho=R*(rand(N,1).^(1/3));
        [x,y,z]=sph2cart(az,el,rho);
        p=[x,y,z];
    end
end