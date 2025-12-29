function H=pot_force(potentialtype,cutoff,points,eqdistance,welldepth)
	
	x=linspace(0,cutoff,points)';
	H(:,1)=x;
	sigma11=(eqdistance);
	xtilde=x/sigma11;
    if potentialtype==1 % Lennard Jones      
        Ftilde = 24 .* (2 ./ xtilde.^13 - 1 ./ xtilde.^7);
        H(:,2) = (welldepth ./ sigma11) .* Ftilde;
        H(1,:)=[];   
    elseif potentialtype==2 % Weeks Chandler Andersen
        Ftilde = 24 .* (2 ./ xtilde.^13 - 1 ./ xtilde.^7);  % correct derivative
        Ftilde(xtilde >= 2^(1/6)) = 0;  % truncate at r_min
        H(:,2) = (welldepth ./ sigma11) .* Ftilde;
        H(1,:) = []; 
    elseif potentialtype==3 % quadratic repulsion
		U = 0.5 * welldepth .* (1 - x./sigma11).^2;   
        dx = x(2) - x(1);
        F = -gradient(U, dx);   % proper derivative
        F(x >= sigma11) = 0;
        H(:,2)=F;
        H(1,:) = []; 		
	elseif potentialtype==4 % Hertzian repulsion
        U = welldepth .* (1 - x./sigma11).^(5/2);  
        dx = x(2) - x(1);
        F = -gradient(U, dx);   % proper derivative
        F(x >= sigma11) = 0;
        H(:,2)=F;
        H(1,:) = []; 
	end
end