function SSF = ssf_accumulation(p, SSF, index, N)
    % Ensure p is N x 3
    if size(p, 2) ~= 3, p = p'; end

    % --- NO ROTATION NEEDED ---
    % The k-vectors in SSF are already defined in the reciprocal basis 
    % corresponding to the unrotated p.
    
    % Note: Your previous code transposed p for calculation.
    % sum(exp(-1i * K * p'), 2) is standard if K is (NumK x 3) and p is (N x 3).
    
    % Helper for cleaner code
    function val = calc_sk(kvecs)
        % kvecs: M x 3
        % p: N x 3
        % Result: sum over N particles of exp(-i * k * p)
        % Argument for exp: (M x 3) * (3 x N) -> (M x N)
        phase = -1i * (kvecs * p'); 
        raw = sum(exp(phase), 2); % Sum over particles (dim 2) -> M x 1 vector
        val = raw; 
    end

    % Calculate Raw sums (Complex numbers)
    ssf100_raw = calc_sk(SSF.kvec100);
    ssf010_raw = calc_sk(SSF.kvec010);
    ssf001_raw = calc_sk(SSF.kvec001);
    
    ssf111_raw = calc_sk(SSF.kvec111);
    
    ssf110_raw = calc_sk(SSF.kvec110);
    ssf011_raw = calc_sk(SSF.kvec011);
    ssf101_raw = calc_sk(SSF.kvec101);
    
    % Calculate S(k) = (1/N) * |sum|^2
    % This is the Static Structure Factor
    
    % Store RAW averages (if you want the complex order parameter)
    % Note: Averaging [100], [010], [001] assumes cubic isotropy of the lattice
    SSF.SSF100_raw(index,:) = mean([abs(ssf100_raw).^2, abs(ssf010_raw).^2, abs(ssf001_raw).^2], 2) / N;
    SSF.SSF110_raw(index,:) = mean([abs(ssf110_raw).^2, abs(ssf011_raw).^2, abs(ssf101_raw).^2], 2) / N;
    SSF.SSF111_raw(index,:) = (abs(ssf111_raw).^2) / N;
    
    % Store S(k) (The final real number)
    % (You seem to be storing the same thing in _raw and non-raw here? 
    % Usually _raw might imply the complex sum before abs^2, but your code took abs^2.
    % I will maintain your logic of storing the magnitudes.)
    
    SSF.SSF100(index,:) = SSF.SSF100_raw(index,:);
    SSF.SSF110(index,:) = SSF.SSF110_raw(index,:);
    SSF.SSF111(index,:) = SSF.SSF111_raw(index,:);
end