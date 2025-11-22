function SSF=ssf_initialization(S)
    SSF.numk = 1000; 
    SSF.k_min = 2 * pi / (2*S.br);
    SSF.k_max = 2 * pi / S.rp;
    SSF.n = (1:SSF.numk)';
    % --- [100] Direction ---
    kvec100_candidates = SSF.k_min * [SSF.n, zeros(SSF.numk,1), zeros(SSF.numk,1)]; % Vectors are of the form (n*k_min, n*k_min, 0)
    k_magnitudes_100 = SSF.k_min * SSF.n;
    idx = k_magnitudes_100 <= SSF.k_max;
    SSF.kvec100 = kvec100_candidates(idx, :);
    % --- [010] Direction ---
    kvec010_candidates = SSF.k_min * [zeros(SSF.numk,1), SSF.n, zeros(SSF.numk,1)]; % Vectors are of the form (0, n*k_min, n*k_min)
    k_magnitudes_010 = SSF.k_min * SSF.n;
    idx = k_magnitudes_010 <= SSF.k_max;
    SSF.kvec010 = kvec010_candidates(idx, :);
    % --- [001] Direction ---
    kvec001_candidates = SSF.k_min * [zeros(SSF.numk,1), zeros(SSF.numk,1), SSF.n]; % Vectors are of the form (SSF.n*k_min, 0, SSF.n*k_min)
    k_magnitudes_001 = SSF.k_min * SSF.n;
    idx = k_magnitudes_001 <= SSF.k_max;
    SSF.kvec001 = kvec001_candidates(idx, :);
    % --- [111] Direction ---
    kvec111_candidates = SSF.k_min * [SSF.n, SSF.n, SSF.n]; % Each component must be an integer multiple of SSF.k_min. This creates (SSF.k_min,SSF.k_min,SSF.k_min), (2*SSF.k_min,2*SSF.k_min,2*SSF.k_min), etc.
    k_magnitudes_111 = SSF.k_min * SSF.n * sqrt(3);
    idx = k_magnitudes_111 <= SSF.k_max;
    SSF.kvec111 = kvec111_candidates(idx, :);         
    % --- [110] Direction ---
    kvec110_candidates = SSF.k_min * [SSF.n, SSF.n, zeros(SSF.numk,1)]; % Vectors are of the form (SSF.n*SSF.k_min, SSF.n*SSF.k_min, 0)
    k_magnitudes_110 = SSF.k_min * SSF.n * sqrt(2);
    idx = k_magnitudes_110 <= SSF.k_max;
    SSF.kvec110 = kvec110_candidates(idx, :);
    % --- [011] Direction ---
    kvec011_candidates = SSF.k_min * [zeros(SSF.numk,1), SSF.n, SSF.n]; % Vectors are of the form (0, SSF.n*SSF.k_min, SSF.n*SSF.k_min)
    k_magnitudes_011 = SSF.k_min * SSF.n * sqrt(2);
    idx = k_magnitudes_011 <= SSF.k_max;
    SSF.kvec011 = kvec011_candidates(idx, :);
    % --- [101] Direction ---
    kvec101_candidates = SSF.k_min * [SSF.n, zeros(SSF.numk,1), SSF.n]; % Vectors are of the form (SSF.n*SSF.k_min, 0, SSF.n*SSF.k_min)
    k_magnitudes_101 = SSF.k_min * SSF.n * sqrt(2);
    idx = k_magnitudes_101 <= SSF.k_max;
    SSF.kvec101 = kvec101_candidates(idx, :);
end