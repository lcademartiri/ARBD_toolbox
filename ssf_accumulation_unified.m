function SSF = ssf_accumulation_unified(p, SSF,cacheSizeMB)
% SSF = ssf_accumulation_unified(p, SSF)
%
% Accumulates data for both Lattice and Sampling modes.
% 1. Updates running sums (Mean SSF, Mean Density) -> Zero RAM growth.
% 2. Appends to time-series (rho_ts) ONLY if initialized in 'init'.
% 3. Uses chunking to prevent calculation memory spikes.

    % Ensure p is N x 3
    if size(p, 2) ~= 3, p = p'; end
    N = size(p, 1);
    
    % Update global snapshot counter
    SSF.nsnap = SSF.nsnap + 1;
    index = SSF.nsnap;

    % --- PROCESS LATTICE MODES ---
    SSF.lattice = accumulate_set(p, SSF.kvecs_lattice, SSF.lattice, index, N,cacheSizeMB);

    % --- PROCESS SAMPLING MODES ---
    SSF.sampling = accumulate_set(p, SSF.kvecs_sampling, SSF.sampling, index, N,cacheSizeMB);

end


function DataStruct = accumulate_set(p, kvecs, DataStruct, index, N,cacheSizeMB)
% Helper to process a specific set of k-vectors (Lattice or Sampling)

    num_k = size(kvecs, 1);
    
    % --- CHUNKING PARAMETERS ---
    % To calculate exp(-i*k*r), we need a matrix of size [NumK x NumPart].
    % If NumK is 100,000 and NumPart is 1,000, that's 1.6 GB complex double.
    % We process in chunks of 'chunk_size' k-vectors to stay in CPU cache.
    bytes_per_element = 16;
    
    % We use a safety factor (0.5) to leave room for the input coordinates, 
    % result vectors, and OS overhead/instructions.
    safety_factor = 0.5; 
    
    % Target memory in bytes
    target_mem_bytes = (cacheSizeMB * 1024 * 1024) * safety_factor;
    
    % Calculate how many K-vectors we can handle per chunk given N particles
    mem_per_k_row = N * bytes_per_element;
    
    chunk_size = floor(target_mem_bytes / mem_per_k_row);
    
    % Clamping:
    % 1. Must be at least 1.
    % 2. Cap at 10000 to maintain temporal locality and prevent OS paging
    %    issues even if cache is reported as huge.
    chunk_size = max(1, chunk_size);
    chunk_size = min(10000, chunk_size);
    
    % Prepare temporary arrays for this snapshot
    % (We need the full vector for this snapshot to store it, 
    %  but we calculate it pieces)
    current_rho = complex(zeros(1, num_k));
    
    for start_idx = 1:chunk_size:num_k
        end_idx = min(start_idx + chunk_size - 1, num_k);
        k_chunk = kvecs(start_idx:end_idx, :);
        
        % Vectorized Phase Calculation: (N_k_chunk x 3) * (3 x N_part)
        phase = -1i * (k_chunk * p'); 
        
        % Sum over particles to get rho(k) for this chunk
        % Result is (N_k_chunk x 1)
        rho_chunk = sum(exp(phase), 2).'; 
        
        % Store in the full snapshot vector
        current_rho(start_idx:end_idx) = rho_chunk;
    end
    
    % --- ACCUMULATE STATISTICS (Zero RAM Cost) ---
    
    % 1. Accumulate Complex Sum (for Mean Density / Order Parameter)
    DataStruct.sum_rho = DataStruct.sum_rho + current_rho;
    
    % 2. Accumulate |rho|^2 (for Structure Factor S(k))
    % S(k) = < |rho|^2 > / N
    % We sum |rho|^2 here, and divide by (nsnap * N) later during analysis.
    DataStruct.sum_abs2 = DataStruct.sum_abs2 + (abs(current_rho).^2);
    
    % --- OPTIONAL TIME SERIES STORAGE ---
    % Only store if the array was pre-allocated in the init function.
    % (Checking isempty is fast)
    if ~isempty(DataStruct.rho_ts)
        DataStruct.rho_ts(index, :) = current_rho;
    end
    
end