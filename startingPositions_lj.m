function p = startingPositions_lj(S,scaleflag)

        % starting positions are calculated by (1) taking a close-packed
        % lattice with minimum distance between nodes slightly bigger than
        % 2r; (2) displacing/rotating the lattice so that we maximize the number of
        % nodes in the boundary; (3) trimming the lattice so that all nodes
        % are within a radiu from the boundary; (4) picking randomly a number of nodes
        % equal to the number of particles desired a using those as
        % starting positions.

disp('assemblying starting positions')
if S.phi>=0.01
    flag=1;
    while flag==1
        if scaleflag
            lattice_scale = max(1, (0.50 / S.phi)^(1/3));
        else
            lattice_scale = 1;
        end
		basis=[0,0.7071,0.7071;0.7071,0,0.7071;0.7071,0.7071,0].*(2.01 * S.rp * lattice_scale);
        maxsteps=2*ceil(((S.br*2)*sqrt(3))/(2*S.rp)); % calculate the maximum number of unit cells of the lattice we need to exceed the dimensions of the boundary
        templist=double(linspace(-maxsteps,maxsteps,2*maxsteps+1)'); % 1D coefficient set
        [x1,x2,x3] = meshgrid(templist,templist,templist); % create lattice of coefficients
        x1=x1(:); % create lattice of coefficients
        x2=x2(:); % create lattice of coefficients
        x3=x3(:); % create lattice of coefficients
        possiblepositions = x1.*basis(1,:)+x2.*basis(2,:)+x3.*basis(3,:); % multiplying the lattice of coefficients by the basis vectors to obtain the real lattice
        possiblepositions=bsxfun(@plus,possiblepositions,-S.rp*[1,1,1]./vecnorm([1,1,1],2)); % displacing the whole lattice to ensure maximum packing within the boundary
        % --- eliminate all lattice points outside 'br'-'r' ---
        if S.bc==2 || S.bc==4
            idxstart=sum(abs(possiblepositions)>(S.br-S.rp),2)>0;
        elseif S.bc==1
            tempnorms=vecnorm(possiblepositions,2,2);
            idxstart=tempnorms>(S.br-S.rp);
        elseif S.bc==3
            A = [S.fcc.a1', S.fcc.a2', S.fcc.a3'];
            f = (A \ possiblepositions')';                  % fractional coords
            inside = all(f >= 0 & f <= 1, 2);
            d1 = min(f(:,1), 1-f(:,1)) * norm(S.fcc.a1);
            d2 = min(f(:,2), 1-f(:,2)) * norm(S.fcc.a2);
            d3 = min(f(:,3), 1-f(:,3)) * norm(S.fcc.a3);
            dist_to_boundary = min([d1,d2,d3],[],2);
            % Margin test
            margin_ok = dist_to_boundary >= S.rp;
            
            % Final mask
            keep_mask = inside & margin_ok;
            idxstart=~keep_mask;
        end
        possiblepositions(idxstart,:)=[];
        % ---
        possiblepositions=possiblepositions(randperm(size(possiblepositions,1)),:); % randomize the order of the lattice nodes in the list
        p=possiblepositions(1:S.N,:); % pick the first N nodes in the list
        % --- double checking that particles in the nodes that were picked would not be touching ---
        D=pdist(p)';
        idxcoll=D<(2*S.rp);
        if sum(idxcoll)==0 % if all checks out tap out.
            flag=0;
        end
        % ---
    end 
    clear possiblepositions tempnorms x1 x2 x3 templist maxsteps theta u ct mct st basis rotmat D
else
    flag=1;
    while flag==1
        if S.bc==1
            Nred=S.N*3;
            p=rand(Nred,3)*2*(S.br-S.rp)-(S.br-S.rp);
            idx=vecnorm(p,2,2)>(S.br-S.rp);
            p(idx,:)=[];
            idxcoll=sign(2*S.rp-pdist(p)')+1;
            Nred=size(p,1);
            if sum(idxcoll,'all')>0
                collpairs=find(idxcoll);
                bin=ceil(-0.5*sqrt(8*nchoosek(Nred,2)-8*collpairs+1)+Nred-0.5);
                coll1=bin;
                binedges=0.5*(Nred-bin)-0.5*(Nred-bin).^2+nchoosek(Nred,2);
                coll2=Nred-binedges+collpairs;
                p([coll1;coll2],:)=[];
                if size(p,1)>S.N
                    flag=0;
                    p=p(1:S.N,:);
                end
            else
                flag=0;
                p=p(1:S.N,:);
            end
        elseif S.bc==2 || S.bc==4
            Nred=S.N*2;
            p=rand(Nred,3)*2*(S.br-S.rp)-(S.br-S.rp);
            idxcoll=sign(2*S.rp-pdist(p)')+1;
            if sum(idxcoll,'all')>0
                collpairs=find(idxcoll);
                bin=ceil(-0.5*sqrt(8*nchoosek(Nred,2)-8*collpairs+1)+Nred-0.5);
                coll1=bin;
                binedges=0.5*(Nred-bin)-0.5*(Nred-bin).^2+nchoosek(Nred,2);
                coll2=Nred-binedges+collpairs;
                p([coll1;coll2],:)=[];
                if size(p,1)>S.N
                    flag=0;
                    p=p(1:S.N,:);
                end
            else
                flag=0;
                p=p(1:S.N,:);
            end
        elseif S.bc==3
            mindist=((3*(S.bv/S.N))/(4*pi))^(1/3);
            basis=[0,0.7071,0.7071;0.7071,0,0.7071;0.7071,0.7071,0].*(mindist); % basis vectors of the fcc lattice
            maxsteps=2*ceil(((S.br*2)*sqrt(3))/(mindist)); % calculate the maximum number of unit cells of the lattice we need to exceed the dimensions of the boundary
            templist=double(linspace(-maxsteps,maxsteps,2*maxsteps+1)'); % 1D coefficient set
            [x1,x2,x3] = meshgrid(templist,templist,templist); % create lattice of coefficients
            x1=x1(:); % create lattice of coefficients
            x2=x2(:); % create lattice of coefficients
            x3=x3(:); % create lattice of coefficients
            possiblepositions = x1.*basis(1,:)+x2.*basis(2,:)+x3.*basis(3,:); % multiplying the lattice of coefficients by the basis vectors to obtain the real lattice
            possiblepositions=bsxfun(@plus,possiblepositions,-S.rp*[1,1,1]./vecnorm([1,1,1],2)); % displacing the whole lattice to ensure maximum packing within the boundary
            a1 = 2*S.br.*S.fcc_unitvecs(1,:)';
            a2 = 2*S.br.*S.fcc_unitvecs(2,:)';
            a3 = 2*S.br.*S.fcc_unitvecs(3,:)';
            A = [a1, a2, a3];
            f = (A \ possiblepositions')';                  % fractional coords
            inside = all(f >= 0 & f <= 1, 2);
            d1 = min(f(:,1), 1-f(:,1)) * norm(a1');
            d2 = min(f(:,2), 1-f(:,2)) * norm(a2');
            d3 = min(f(:,3), 1-f(:,3)) * norm(a3');
            dist_to_boundary = min([d1,d2,d3],[],2);
            % Margin test
            margin_ok = dist_to_boundary >= S.rp;
            
            % Final mask
            keep_mask = inside & margin_ok;
            idxstart=~keep_mask;
            possiblepositions(idxstart,:)=[];
            % ---
            possiblepositions=possiblepositions(randperm(size(possiblepositions,1)),:); % randomize the order of the lattice nodes in the list
            p=possiblepositions(1:S.N,:); % pick the first N nodes in the list
            % --- double checking that particles in the nodes that were picked would not be touching ---
            D=pdist(p)';
            idxcoll=D<(2*S.rp);
            if sum(idxcoll)==0 % if all checks out tap out.
                flag=0;
            end
            % ---
        end
    end

end