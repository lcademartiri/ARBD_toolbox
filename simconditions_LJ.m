function [V,CONDS] = simconditions_LJ(V,C,P)

disp('compilation of all the conditions to analyze')

    % range of particle radii
    if V.r_range(4)==1
        V.r=logspace(log10(V.r_range(1)),log10(V.r_range(2)),V.r_range(3))';
    else
        V.r=linspace((V.r_range(1)),(V.r_range(2)),V.r_range(3))';
    end
    % range of particle numbers
    if V.N_range(4)==1
        V.N=ceil(logspace(log10(V.N_range(1)),log10(V.N_range(2)),V.N_range(3))');
    else
        V.N=ceil(linspace((V.N_range(1)),(V.N_range(2)),V.N_range(3))');
    end
    % range of volume fractions
    if V.phi_range(4)==1
        V.phi=(logspace(log10(V.phi_range(1)),log10(V.phi_range(2)),V.phi_range(3))');
    else
        V.phi=(linspace((V.phi_range(1)),(V.phi_range(2)),V.phi_range(3))');
    end
    % range of boundary conditions
    V.bc=linspace(V.bc_range(1),V.bc_range(2),V.bc_range(3))';
    % range of BBmultiplier conditions
    V.bbm=linspace(V.bbm_range(1),V.bbm_range(2),V.bbm_range(3))';
    % range of temperatures
    V.T=(logspace(log10(V.T_range(1)),log10(V.T_range(2)),V.T_range(3))');
    % range of LJ trigger
    V.LJpot=linspace(V.LJpot_range(1),V.LJpot_range(2),V.LJpot_range(3))';
        
    V.noconds=size(V.r,1)*size(V.phase_list,1)*size(V.N,1)*size(V.phi,1)*size(V.bc,1)*size(V.T,1)*size(V.solvent_list,1)*size(V.bbm,1);
    c=zeros(V.noconds,19);
    variableNames = {'condition #', 'boundary condition', '# of particles', 'volume fraction', 'phase', 'particle radius [m]', 'particle density [Kg/m^3]', 'particle mass [Kg]',...
    'collisional time [s]', 'relaxation time [s]', 'kuhn time [s]',  'box volume [m]', 'box radius [m]', 'solvent name', 'solvent radius [m]', ...
    'solvent mass [Kg]','solvent density [Kg/m^3]', 'solvent viscosity [Pa s]', 'solvent speed of sound [m/s]', 'BB multiplier','potential'};
    variableNamesShort = {'ID', 'bc', 'N', 'phi', 'phase', 'rp', 'rhop', 'mp',...
    'tauc', 'taur', 'kt', 'bv', 'br', 'solvent', 'rs', ...
    'ms','rhos', 'eta', 'cs', 'bbm','pot'};

    % --- fcc unit cell definition
    a1 = [0 1 1];
    a2 = [1 0 1];
    a3 = [1 1 0];
    u1=a1./norm(a1);
    u2=a2./norm(a2);
    u3=a3./norm(a3);
    vun=abs(dot(u1, cross(u2, u3)));
    V.fcc_unitvecs=[u1;u2;u3];
    V.fcc_unitvol=vun;
    % ---
    q=1;
    for ir=1:V.r_range(3)
        rp=V.r(ir);
        for ip=V.phase_list
            rhop=table2array(C.phasesproperties(ip,2));
            for is=V.solvent_list
                rs=table2array(C.solventproperties(is,4));
                ms=table2array(C.solventproperties(is,3));
                eta=table2array(C.solventproperties(is,6));
                rhos=table2array(C.solventproperties(is,5));
                cs=table2array(C.solventproperties(is,7));                
                for iN=1:V.N_range(3)
                    N=V.N(iN);
                    for iphi=1:V.phi_range(3)
                        phi=V.phi(iphi);
                        for iT=1:V.T_range(3)
                            T=V.T(iT);
                            for ibc=1:V.bc_range(3)
                                bc=V.bc(ibc);
                                for ibbm=1:V.bbm_range(3)
                                    bbm=V.bbm(ibbm);
                                    for ilj=1:V.LJpot_range(3)
                                        lj=V.LJpot(ilj);
                                        c(q,1)=q; % condition number
                                        c(q,2)=bc; % boundary condition
                                        if bc==4
                                            c(q,3)=N*bbm^3;
                                        else
                                            c(q,3)=N; % number of particles []
                                        end
                                        c(q,4)=phi; % volume fraction []
                                        c(q,5)=rp; % particle radius [m]
                                        c(q,6)=rhop; % particle density [m]
                                        mp=(4/3)*pi*rhop*c(q,5)^3;
                                        c(q,7)=mp; % mass of 1 particle [Kg]
                                        taur=((2*rhop)/(9*eta))*c(q,5)^2; % relaxation time [s]
                                        tauc=(9.775e-11/sqrt(C.kbT))*((rs/rp)^2)*((sqrt(ms*mp))/(sqrt(ms)+sqrt(mp))); % collisional time
                                        c(q,8)=tauc; % collisional time [s]
                                        c(q,9)=taur; % relaxation time [s]
                                        c(q,10)=c(q,9)*P.kuhnmultiplier; % time step extending step duration by a multiplier to ensure full relaxation.
                                        v=c(q,3)*(pi*(4/3))*c(q,5)^3; % particle volume
                                        if bc==4
                                            c(q,11)=((v/bbm^3)/c(q,4))*bbm^3; % BIG BOX volume
                                        else
                                            c(q,11)=v/c(q,4); % box volume
                                        end
                                        if c(q,2)==1 || c(q,2)==4 
                                            c(q,12)=((3*c(q,11))/(4*pi))^(1/3); % box radius
                                        elseif c(q,2)==2
                                            c(q,12)=(c(q,11)^(1/3))/2; % box radius
                                        elseif c(q,2)==3
                                            scalingfactor=(c(q,11)/V.fcc_unitvol)^(1/3);
                                            c(q,12)=scalingfactor/2; % box radius in FCC (this means half of the triclinic cell side)
                                        end
                                        temppn(q,1)=table2array(C.phasesproperties(ip,1));
                                        tempsn(q,1)=table2array(C.solventproperties(is,1));
                                        c(q,13)=rs; % solvent radius [m]
                                        c(q,14)=ms; % solvent molecular mass [Kg]
                                        c(q,15)=rhos; % solvent density [Kg/m^3]
                                        c(q,16)=eta; % solvent viscosity [Pa s]
                                        c(q,17)=cs; % solvent speed of sound [m/s]
                                        c(q,18)=bbm; % big box multiplier
                                        c(q,19)=lj; % LJ switch
                                        q=q+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    idx=c(:,2)~=4 & c(:,18)~=V.bbm_range(1);
    c(idx,:)=[];
    temppn(idx,:)=[];
    tempsn(idx,:)=[];
    c(c(:,2)~=3,18)=1;

    % idx=c(:,2)==4 & c(:,3)./(c(:,18).^3)>2;
    % c(idx,:)=[];
    % temppn(idx,:)=[];
    % tempsn(idx,:)=[];
   
    if V.assumptions(1)==1 % V.assumptions(1): exclude those where the relaxation time is smaller than 10 times the tauc
        idx=c(:,9)<10*c(:,8);
        tempsn(idx,:)=[];
        temppn(idx,:)=[];
        c(idx,:)=[];
    end
    if V.assumptions(2)==1 % V.assumptions(2): exclude those where the particle is smaller than the solvent
        idx=c(:,13)>c(:,5);
        tempsn(idx,:)=[];
        temppn(idx,:)=[];
        c(idx,:)=[];
    end
    if V.assumptions(3)==1 % V.assumptions(3): exclude those where the particle is lighter than the solvent
        idx=c(:,14)>c(:,7);
        tempsn(idx,:)=[];
        temppn(idx,:)=[];
        c(idx,:)=[];
    end
    V.noconds=size(c,1);
    c(:,1)=linspace(1,V.noconds,V.noconds)';
    c=array2table(c);
    insertposition1=5;
    insertposition2=12;
    leftPart = c(:, 1:insertposition1-1);
    centerPart = c(:, insertposition1:insertposition2);
    rightPart = c(:, insertposition2+1:end); 
    c=[leftPart, table(temppn),centerPart,table(tempsn), rightPart];
    c.Properties.VariableNames = variableNames;
    
    V.c=c;
    V.vn=variableNames;
    V.vns=variableNamesShort;

    % break up conditions table into separate variables into a struct array so
    % as to not having to edit the entire code everytime we add a variable into
    % the c array
    
    % Initialize a struct to store the columns
    CONDS = struct();
    
    % Loop through each column and store it in the struct
    for i0 = 1:length(V.vn)
        columnName = V.vn{i0};  % Get the column name
        variableName=V.vns{i0};
        columnData = V.c.(columnName);  % Extract the column data    
        % Store the column data in the struct
        CONDS.(variableName) = columnData;
    end
end