function p=MCshaker(p,S,shakes)	
	fprintf('Running Monte Carlo Shaker to erase memory...\n');
	for k = 1:shakes
		idx = randi(S.N);
		trial_pos = p(idx,:) + (rand(1,3)-0.5) * S.br; 
		% Simple Hard Sphere & Boundary Check
		if norm(trial_pos) < (S.br)
			% Only check neighbors if inside box
			dists = pdist2(trial_pos, p, 'squaredeuclidean');
			dists(idx) = inf; % Ignore self
			if min(dists) > (2*S.rp)^2
				p(idx,:) = trial_pos; % Accept Move
			end
		end
	end
end