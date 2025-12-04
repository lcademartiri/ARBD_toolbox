function SSF = ssf_accumulation(p, SSF, index, N)

    % Ensure p is N x 3
    if size(p,2) ~= 3, p = p'; end

    % Helper: compute rho_k = sum_j exp(-i k ⋅ r_j)
    function rho = calc_rho(kvecs)
        phase = -1i * (kvecs * p');    % (M x 3)*(3 x N) => (M x N)
        rho = sum(exp(phase), 2);      % M x 1 complex
    end

    % ---- Compute raw complex order parameters ----
    rho100 = calc_rho(SSF.kvec100);
    rho010 = calc_rho(SSF.kvec010);
    rho001 = calc_rho(SSF.kvec001);

    rho110 = calc_rho(SSF.kvec110);
    rho011 = calc_rho(SSF.kvec011);
    rho101 = calc_rho(SSF.kvec101);

    rho111 = calc_rho(SSF.kvec111);

    % ---- Store complex ρ_k snapshots (essential!) ----
    SSF.rho100(index,:) = rho100.';   % 1 x M complex
    SSF.rho010(index,:) = rho010.';
    SSF.rho001(index,:) = rho001.';

    SSF.rho110(index,:) = rho110.';
    SSF.rho011(index,:) = rho011.';
    SSF.rho101(index,:) = rho101.';

    SSF.rho111(index,:) = rho111.';

    % ---- Compute instantaneous S(k) = |rho_k|^2 / N ----
    SSF.S100(index,:) = abs(rho100.').^2 / N;
    SSF.S010(index,:) = abs(rho010.').^2 / N;
    SSF.S001(index,:) = abs(rho001.').^2 / N;

    SSF.S110(index,:) = abs(rho110.').^2 / N;
    SSF.S011(index,:) = abs(rho011.').^2 / N;
    SSF.S101(index,:) = abs(rho101.').^2 / N;

    SSF.S111(index,:) = abs(rho111.').^2 / N;

end
