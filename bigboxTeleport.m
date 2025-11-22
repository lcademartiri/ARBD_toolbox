function p = bigboxTeleport(S,p)

% this function teleports any particle that hits the outer wall iun a BB
% condition to some other unoccupied spot outside of the center probe
% volume. Displacement vector in the generated particles are inherited from the
% particles that are eliminated. Furthermore it ensures that the particles that are generated
% are not in colliding conditions before or after the displacement.
    idxbb=sum(abs(p(:,8:10))>(S.br-S.rp),2)>0; % idx vector of all the particles whose destination position touches the boundary 
    if sum(idxbb,"all")>0 % if cndition triggered if any particle touches the boundary
        ptomove=p(idxbb,:); % list of the particles that are touching the boundary
        p(idxbb,:)=[]; % eliminate from the master list of particles the particles that are touching the boundary
        for ibb=1:size(ptomove,1) % for loop over all particles that are touching the boundary
            flag=1; % flag for the while loop that will be turned off when a safe landing position is found for the particle
            while flag==1
                ptemp1=rand(1,1)*2*(S.br-S.rp)-(S.br-S.rp);
                ptemp2=rand(1,1)*2*(S.br-S.rp)-(S.br-S.rp);
                ptemp3=rand(1,1)*(2*(S.br/S.bbm-S.rp))+(S.br-2*S.br/S.bbm+S.rp);
                randomsign=2*round(rand(1,1),0)-1;
                ptemp3=ptemp3*randomsign;
                ptemp=[ptemp1,ptemp2,ptemp3];
                ptemp=ptemp(randperm(3));
                % ptemp=rand(1,3)*2*(S.br-S.rp)-(S.br-S.rp); % random attempt at a position in the BB not touching the boundary
                ptemp=[ptemp-ptomove(ibb,5:7);ptemp]; % use existing displacement to backtrack the initial position
                idxcentercell=abs(ptemp)<(S.br-2*S.br/S.bbm+S.rp); % idx matrix identifying coordinates that are inside the probe volume of the BB
                if all(idxcentercell ~= 3) % in case the idx matrix find nothing (so the before and after coordinates are NOT in the probe volume)
                    idxDf=vecnorm(p(:,8:10)-ptemp(2,:),2,2)<(2*S.rp); % idx vector identifying whether there are existing particles in the master matrix p that are touching the randomly chosen position after displacement
                    idxDi=vecnorm(p(:,1:3)-ptemp(1,:),2,2)<(2*S.rp); % idx vector identifying whether there are existing particles in the master matrix p that are touching the randomly chosen position before displacement
                    if sum(idxDi)+sum(idxDf)==0 % if no existing particles are found to be touching the randomly chosen position then flag gets turned off, otherwise we try again.
                        flag=0;
                    end
                end
            end
            pmoved(ibb,:)=[ptemp(1,:),ptomove(ibb,4:7),ptemp(2,:)]; % compile the new position in the format of the master p matrix
            p=[p;pmoved(ibb,:)];% concatenate the new positions of the teleported particles
        end
         
        p=sortrows(p,4,"ascend"); % resort the p matrix on the basis of the ID.
    end

    % verify that, in case more than one particle
end