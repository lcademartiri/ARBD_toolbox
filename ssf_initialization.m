function SSF = ssf_initialization(S)
    SSF.numk = 10000;
    
    % --- 1. Determine Fundamental k-scale based on BC ---
    if S.bc == 3 % FCC PBC
        % For FCC, periodicity is defined by the reciprocal lattice vectors.
        % S.fcc.A contains row vectors of the real-space lattice.
        % Reciprocal lattice B = 2*pi * inv(A).
        % The valid k-vectors are integer combinations of B's rows.
        % Note: Your code defines S.fcc.A with basis vectors as rows?
        % Let's check MICimplementation_v4.txt:
        % S.fcc.A = [a1; a2; a3]'; So columns are a1, a2, a3.
        % S.fcc.invA = inv(S.fcc.A);
        
        % Reciprocal basis vectors b1, b2, b3 are columns of:
        % B_matrix = 2 * pi * S.fcc.invA'; 
        % (Using standard definition: b_i . a_j = 2pi delta_ij)
        
        B_matrix = 2 * pi * S.fcc.invA';
        
        % We can still define a "scalar" k_min for magnitude limiting,
        % but vector generation must use B_matrix.
        SSF.k_min = norm(B_matrix(:,1)); % Magnitude of first reciprocal vector
        
    else % SBC (1), Cubic PBC (2), Big Box (4)
        % Box size L = 2 * S.br
        L = 2 * S.br;
        SSF.k_min = 2 * pi / L;
        
        % For Cubic, the basis is orthogonal
        B_matrix = diag([SSF.k_min, SSF.k_min, SSF.k_min]);
    end

    SSF.k_max = 2 * pi / S.rp;
    n_vec = (1:SSF.numk)';
    
    % --- Helper to generate and filter vectors ---
    function k_out = get_k_vectors(n_multipliers)
        % n_multipliers is a [NumK x 3] matrix of integers (n1, n2, n3)
        
        if S.bc == 3
            % General Lattice: k = n1*b1 + n2*b2 + n3*b3
            % B_matrix columns are b1, b2, b3
            k_candidates = (B_matrix * n_multipliers')';
        else
            % Cubic/SBC: Simple scalar multiplication
            k_candidates = SSF.k_min * n_multipliers;
        end
        
        % Filter by magnitude
        mags = sqrt(sum(k_candidates.^2, 2));
        valid = mags <= SSF.k_max & mags > 0;
        k_out = k_candidates(valid, :);
    end

    % --- [100] Direction (Along first lattice vector) ---
    % n * (1, 0, 0)
    SSF.kvec100 = get_k_vectors([n_vec, zeros(SSF.numk,1), zeros(SSF.numk,1)]);

    % --- [010] Direction ---
    % n * (0, 1, 0)
    SSF.kvec010 = get_k_vectors([zeros(SSF.numk,1), n_vec, zeros(SSF.numk,1)]);

    % --- [001] Direction ---
    % n * (0, 0, 1)
    SSF.kvec001 = get_k_vectors([zeros(SSF.numk,1), zeros(SSF.numk,1), n_vec]);

    % --- [110] Direction ---
    % n * (1, 1, 0)
    SSF.kvec110 = get_k_vectors([n_vec, n_vec, zeros(SSF.numk,1)]);

    % --- [101] Direction ---
    % n * (1, 0, 1)
    SSF.kvec101 = get_k_vectors([n_vec, zeros(SSF.numk,1), n_vec]);

    % --- [011] Direction ---
    % n * (0, 1, 1)
    SSF.kvec011 = get_k_vectors([zeros(SSF.numk,1), n_vec, n_vec]);

    % --- [111] Direction ---
    % n * (1, 1, 1)
    SSF.kvec111 = get_k_vectors([n_vec, n_vec, n_vec]);
end