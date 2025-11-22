function H=pot_force(potentialtype,cutoff,points,eqdistance,welldepth)

    if potentialtype==1 % Lennard Jones
        x=linspace(0,cutoff,points)';
        H(:,1)=x;
        sigma11=(eqdistance)/(2^(1/6));
        xtilde=x/sigma11;
        Ftilde = 24 .* (2 ./ xtilde.^13 - 1 ./ xtilde.^7);
        H(:,2) = (welldepth ./ sigma11) .* Ftilde;
        H(1,:)=[];   
    elseif potentialtype==2 % Weeks Chandler Andersen
        x = linspace(0, cutoff, points)';
        sigma11 = (eqdistance)/(2^(1/6));
        xtilde = x ./ sigma11;
        Ftilde = 24 .* (2 ./ xtilde.^13 - 1 ./ xtilde.^7) ./ xtilde;  % correct derivative
        Ftilde(xtilde >= 2^(1/6)) = 0;  % truncate at r_min
        H(:,1) = x;
        H(:,2) = (welldepth ./ sigma11) .* Ftilde;
        H(1,:) = []; 
    end
end