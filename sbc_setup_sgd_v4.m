function [p,pgp,sgd_correction,sgd_edges]=sbc_setup_sgd_v4(S,PDF,nsteps)

    % --- SETUP PARAMETERS ------------------------------------------------
    
    gCS=(1-S.phi/2)/(1-S.phi)^3;
    diffE=S.esdiff*S.alpha/gCS;
    tau_alpha=(S.rp^2)/(6*diffE);
    relaxsteps=tau_alpha/S.timestep;
    
    sgd_base_gain = 0.1;  
    sgd_gain = sgd_base_gain; 
    sgd_smooth_win = 5; 
    sgd_cap = 0.003 * S.rp;
    potdepth=2*S.rc;
    if 2*S.br-2*potdepth<(10*S.rp)
        potdepth=2*S.br-10*S.rp;
    end
    if 2*S.br-2*S.rc<(10*S.rp)
        potdepth=S.rc;
    end
    
    % Extended to capture structural correlations
    sgd_edges = sort((S.br:-0.02*S.rp:S.br - potdepth)');
    metric_smoothing_param=0.8;
    bold=true;
    batch=true;
    aggro=false;
    grace=true;
    graceperiod=20e3;
    patience0=2;
    
    % Automatic Batch Sizing
    base_batch_size = ceil(10*relaxsteps);
    sgd_batch_size = base_batch_size;
    
    seriesname='sgd_msp03_ig01_batch1e3_grace';

    if S.potential==1, potname='lj';
    elseif S.potential==2, potname='wca';
    elseif S.potential==0, potname='hs';
    end
    
    filenamecorrection=sprintf(['ASYMCORR_',seriesname,'_%s_%.0e_%.0e_%.0f_%.1f_%.1e.mat'],...
        potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
    filestartingconfiguration=sprintf(['START_',potname,'_%.0e_%.0e_%.0f_%.1f_%.1e.mat'],...
        potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
    filepdfdenom=sprintf('PDFdenom_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);

    debugging=false;
    graphing=true;

    if exist(filenamecorrection,'file') && exist(filestartingconfiguration,'file')
        load(filenamecorrection); load(filestartingconfiguration);
        return
    end

    % --- PRE-CALCULATIONS ------------------------------------------------
    if exist(filepdfdenom,'file')
        load(filepdfdenom,'gdenominator');
    else
        gdenominator=PDFdenom(S,PDF,nsteps);
        save(filepdfdenom,'gdenominator','S')
    end

    if S.potential~=0
        H=pot_force(S.potential,S.rc,30000,S.pot_sigma,S.pot_epsilon);
        clamp=mcdClamp(1e6,S.rp,normrnd(0,S.stdx,1e6,3),S.esdiff,S.timestep,H,S.kbT);
        H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'linear', 'nearest');
    end

    % --- INITIAL PLACEMENT (FCC) -----------------------------------------
    disp('Creating initial FCC lattice...')
    flag=1;
    if debugging, rng(100); end
    
    while flag==1
        basis=[0,0.7071,0.7071;0.7071,0,0.7071;0.7071,0.7071,0].*(2.01*S.rp);
        maxsteps=2*ceil(((S.br*2)*sqrt(3))/(2*S.rp));
        templist=double(linspace(-maxsteps,maxsteps,2*maxsteps+1)');
        [x1,x2,x3] = meshgrid(templist,templist,templist);
        possiblepositions = x1(:).*basis(1,:)+x2(:).*basis(2,:)+x3(:).*basis(3,:);
        possiblepositions=bsxfun(@plus,possiblepositions,-S.rp*[1,1,1]./vecnorm([1,1,1],2));
        
        tempnorms=vecnorm(possiblepositions,2,2);
        possiblepositions(tempnorms>(S.br-2*S.rp),:)=[];
        possiblepositions=possiblepositions(randperm(size(possiblepositions,1)),:);
        p=possiblepositions(1:S.N,:);
        
        D=pdist(p)';
        if sum(D<(2*S.rp))==0, flag=0; end
    end 
    clear possiblepositions tempnorms x1 x2 x3 templist maxsteps basis D

    % --- SGD INIT --------------------------------------------------------
    sgd_bins = numel(sgd_edges)-1;
    sgd_centers = sgd_edges(1:end-1) + diff(sgd_edges)/2;
    sgd_correction = zeros(sgd_bins, 1);
    F_corr_interp = griddedInterpolant(sgd_centers, sgd_correction, 'linear', 'nearest');
    
    % Accumulators
    batch_sum_drift = zeros(sgd_bins, 1);
    batch_counts = zeros(sgd_bins, 1);
    steps_in_batch = 0;

    % --- CONTROLLER STATE ---
    best_pdf_metric = inf;
    best_correction = sgd_correction;
    best_p_snapshot = p; 
    best_pgp_snapshot = [];
    
    patience_counter = 0;
    pdf_metric = 0;
    
    PHASE_SEARCH = 1;
    PHASE_REFINE = 2;
    current_phase = PHASE_SEARCH;
    
    fprintf('SGD Init. Phase: SEARCH. Batch: %d.\n', sgd_batch_size);

    % --- GRAPHING INIT ---------------------------------------------------
    if graphing
        ndens.av_window=1000; 
        ndens.edges=sort((S.br:-0.02*S.rp:0)');
        ndens.centers = ndens.edges(1:end-1,:) + 0.01*S.rp;
        ndens.counts=zeros(numel(ndens.centers),1);
        ndens.counts_storage=ndens.counts;
        ndens.vols = (4/3)*pi*(ndens.edges(2:end).^3 - ndens.edges(1:end-1).^3);
        ndens.ndens0=(S.N/S.bv);
        
        pdf.pre.counts=zeros(numel(PDF.pdfedges{3})-1,1);
        
        history.steps = [];
        history.pdf_dev = [];
        history.pdf_smooth = [];
        history.max_corr = [];
        history.gain = [];
        history.batch_size = [];
        
        f_fig = figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
        ax_dens = subplot(2,2,1);
        ax_pdf  = subplot(2,2,2);
        ax_conv = subplot(2,2,3);
        ax_ctrl = subplot(2,2,4);
    end

    % --- MAIN LOOP -------------------------------------------------------
    qs=0;
    thermflag=0;
    r2_uniform = 3/5 * S.br^2;
    
    pgp=p-(2*S.br).*(p./vecnorm(p,2,2));
    reverseStr = '';
    pdf_metric_val = 10;
    
    disp('Starting SGD Evolution...')

    while thermflag==0 || pdf_metric_val>0.002 
        qs=qs+1;
        
        prho = vecnorm(p, 2, 2);
        pvers = p ./ (prho + eps);
        idxgp = prho > (S.br - S.rc);
        
        if thermflag==0
            spread_ratio = mean(prho.^2) / r2_uniform;
            if spread_ratio > 0.99 || (qs>2000 && abs(spread_ratio-1)<0.01)
                thermflag=1;
                qs=1; 
                disp('--- Thermalization Complete ---');
            end
        end

        % 3. Base Displacement
        if S.potential ~= 0
            ptemp = [p; pgp(idxgp,:)];
            all_potdisps = potential_displacements_v2(ptemp, size(ptemp,1), S.rc, ...
                H, H_interpolant, S.esdiff, clamp, S.kbT, S.stdx, S.timestep, 1);
            potdisps = all_potdisps(1:S.N, :);
            potdispsgp = all_potdisps(S.N+1:end, :);
            v_rand = randn(S.N, 3) * S.stdx;
            base_disp = v_rand + potdisps;
        else
            v_rand = randn(S.N, 3) * S.stdx;
            base_disp = v_rand;
            potdispsgp = zeros(sum(idxgp), 3);
        end

        % --- SGD ACCUMULATION ----------------------------------------
        if thermflag == 1
            % 1. Accumulate Drift
            dr_raw = sum(base_disp .* pvers, 2);
            [~, bin_idx] = histc(prho, sgd_edges);
            valid_mask = bin_idx > 0 & bin_idx <= sgd_bins;
            
            if any(valid_mask)
                new_sums = accumarray(bin_idx(valid_mask), dr_raw(valid_mask), [sgd_bins, 1]);
                new_counts = accumarray(bin_idx(valid_mask), 1, [sgd_bins, 1]);
                batch_sum_drift = batch_sum_drift + new_sums;
                batch_counts = batch_counts + new_counts;
            end
            
            % 2. Accumulate Diagnostics (Every step)
            if graphing
                [hc, ~] = histcounts(vecnorm(p,2,2), ndens.edges);
                ndens.counts = ndens.counts + hc';
                pairdists = pdist(p);
                [hc_pdf, ~] = histcounts(pairdists, PDF.pdfedges{3});
                pdf.pre.counts = pdf.pre.counts + hc_pdf';
            end
            
            steps_in_batch = steps_in_batch + 1;
            
            % --- BATCH UPDATE & CONTROL TRIGGER ----------------------
            if steps_in_batch >= sgd_batch_size
                
                % A. SGD Update
                has_data = batch_counts > 10;
                if any(has_data)
                    batch_mean_drift = batch_sum_drift(has_data) ./ batch_counts(has_data);
                    sgd_correction(has_data) = sgd_correction(has_data) - sgd_gain * batch_mean_drift;
                    sgd_correction = max(min(sgd_correction, sgd_cap), -sgd_cap);
                    sgd_correction = smoothdata(sgd_correction, 'movmean', sgd_smooth_win);
                    F_corr_interp.Values = sgd_correction;
                end
                
                % B. Diagnostics & Control
                if graphing
                    % Normalization uses EXACT accumulated steps
                    w_count = steps_in_batch;
                    curr_g = (pdf.pre.counts / w_count) ./ gdenominator;
                    
                    valid_mask = PDF.centers{3} > 2*(S.br-potdepth) & PDF.centers{3} < 2*(S.br-0.5*S.rp);
                    
                    % Auto-Target
                    expected_counts = gdenominator(valid_mask) * w_count;
                    expected_counts(expected_counts == 0) = inf; 
                    bin_noise_sigma = 1 ./ sqrt(expected_counts);
                    convergence_target = sqrt(mean(bin_noise_sigma.^2));
                    
                    if any(valid_mask)
                        residuals = curr_g(valid_mask) - 1;
                        weights = gdenominator(valid_mask);
                        raw_pdf_metric = sqrt(sum(weights .* (residuals.^2)) / sum(weights));
                    else
                        raw_pdf_metric = 1;
                    end
                    
                    % Smoothing
                    if pdf_metric == 0, pdf_metric = raw_pdf_metric; end 
                    pdf_metric = (1-metric_smoothing_param) * pdf_metric + metric_smoothing_param * raw_pdf_metric;
                    pdf_metric_val = pdf_metric;
                    
                    % --- CONTROLLER LOGIC ---
                    if pdf_metric < best_pdf_metric
                        best_pdf_metric = pdf_metric;
                        best_correction = sgd_correction;
                        best_p_snapshot = p; 
                        best_pgp_snapshot = pgp;
                        patience_counter = 0; 
                        
                        if bold && current_phase == PHASE_SEARCH && pdf_metric > 2*convergence_target
                             sgd_base_gain = sgd_base_gain * 1.05;
                        end
                        
                        if current_phase == PHASE_SEARCH && pdf_metric < (convergence_target * 2.0) && qs>graceperiod
                             fprintf('\n=== ENTERING REFINEMENT PHASE ===\n');
                             sgd_correction = best_correction;
                             F_corr_interp.Values = sgd_correction;
                             p = best_p_snapshot;
                             pgp = best_pgp_snapshot;
                             
                             sgd_batch_size = 10 * base_batch_size; 
                             sgd_base_gain = sgd_base_gain * sqrt(0.1);
                             sgd_gain = sgd_base_gain;
                             
                             % Clear and skip plot to start fresh batch
                             batch_sum_drift(:) = 0; batch_counts(:) = 0; 
                             ndens.counts(:) = 0; pdf.pre.counts(:) = 0;
                             steps_in_batch = 0;
                             current_phase = PHASE_REFINE;
                             continue; 
                        end
                        
                    elseif pdf_metric > (best_pdf_metric * 1.3) && grace && qs > graceperiod
                        patience_counter = patience_counter + 1;
                        if patience_counter > patience0
                            if sgd_base_gain > 1e-4
                                sgd_base_gain = sgd_base_gain * 0.5;
                                fprintf('--- Gain slashed to %.1e ---\n', sgd_base_gain);
                            end
                            patience_counter = 0;
                            best_pdf_metric = pdf_metric; 
                        end
                    end
                    
                    % Gain Scaling
                    if current_phase == PHASE_SEARCH
                        target_ratio = pdf_metric / 0.05;
                        sgd_gain = min(sgd_base_gain, sgd_base_gain * target_ratio);
                        sgd_gain = max(sgd_gain, 1e-5);
                    else
                        sgd_gain = sgd_base_gain;
                    end
                    
                    % Plotting
                    history.steps(end+1) = qs;
                    history.pdf_dev(end+1) = raw_pdf_metric;
                    history.pdf_smooth(end+1) = pdf_metric;
                    history.max_corr(end+1) = max(abs(sgd_correction));
                    history.gain(end+1) = sgd_gain;
                    history.batch_size(end+1) = sgd_batch_size;
                    
                    
                    set(0, 'CurrentFigure', f_fig);
                    curr_ndens_norm = 100 * (((ndens.counts/w_count)./ndens.vols) ./ ndens.ndens0);
    
                    subplot(ax_dens);
                    plot(ndens.centers, curr_ndens_norm, 'LineWidth', 2, 'Color', 'w');
                    xline(S.br-S.rp, '--w'); xline(S.br, '-r');
                    ylim([80 120]); xlim([0 S.br]);
                    title(sprintf('Density Dev: %.2f%%', max(abs(curr_ndens_norm-100))));
                    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
                    set(gca, 'LineWidth', 3);
    
                    subplot(ax_pdf);
                    plot(PDF.centers{3}, curr_g, 'LineWidth', 2, 'Color', 'y');
                    yline(1, '--w');
                    xlim([0 2.1*S.br]); ylim([0.5 1.5]);
                    title(sprintf('PDF RMS: %.4f', raw_pdf_metric));
                    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
                    set(gca, 'LineWidth', 3);
    
                    % --- MODIFIED BOTTOM PLOTS ---
                    subplot(ax_conv);
                    yyaxis left
                    plot(history.steps, history.pdf_dev, '-', 'Color', [1 1 0 0.4], 'LineWidth', 0.5);
                    hold on;
                    plot(history.steps, history.pdf_smooth, '-y', 'LineWidth', 2); hold off;
                    ylabel('PDF RMS'); set(gca, 'YColor', 'y', 'YScale', 'log');
                    yyaxis right
                    plot(history.steps, history.max_corr, '-c', 'LineWidth', 1.5);
                    ylabel('Max Corr [m]'); set(gca, 'YColor', 'c', 'YScale', 'linear');
                    title('Convergence Metrics');
                    set(gca, 'Color', 'k', 'XColor', 'w');
                    set(gca, 'LineWidth', 3);
                    set(gcf, 'Color', 'k');
    
                    subplot(ax_ctrl);
                    yyaxis left
                    plot(history.steps, history.gain, '-m', 'LineWidth', 2);
                    ylabel('Gain'); set(gca, 'YColor', 'm', 'YScale', 'log');
                    yyaxis right
                    plot(history.steps, history.batch_size, '-g', 'LineWidth', 2);
                    ylabel('Batch Size'); set(gca, 'YColor', 'g', 'YScale', 'linear');
                    title('Control State');
                    set(gcf, 'Color', 'k');
                    set(gca, 'Color', 'k', 'XColor', 'w');
                    set(gca, 'LineWidth', 3);
                    % -----------------------------
    
                    
                    set(gca, 'XColor', 'w');
                    set(gca, 'YColor', 'w');
                    set(get(gca, 'Title'), 'Color', 'w');
                    set(get(gca, 'XLabel'), 'Color', 'w');
                    set(get(gca, 'YLabel'), 'Color', 'w');
                    set(gca, 'LineWidth', 3);
                    
                    drawnow;
                    
                    msg = sprintf('Step %d | RMS: %.4f (Target: %.4f) | Gain: %.1e', ...
                        qs, pdf_metric, convergence_target, sgd_gain);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat('\b', 1, length(msg));
                    
                    if pdf_metric <= convergence_target && qs > 20000 && current_phase == PHASE_REFINE
                        disp('--- CONVERGED (Hit Noise Floor) ---');
                        break;
                    end
                end
                
                % Reset Accumulators
                batch_sum_drift(:) = 0;
                batch_counts(:) = 0;
                ndens.counts(:) = 0;
                pdf.pre.counts(:) = 0;
                steps_in_batch = 0;
            end
        end
        % -----------------------------------------------------------------
        
        % 4. Apply Correction
        dr_corr_mag = F_corr_interp(prho);
        core_mask = prho < (S.br - potdepth);
        dr_corr_mag(core_mask) = 0;
        total_disp = base_disp + (dr_corr_mag .* pvers);
        
        % 5. Position Update
        p2 = p + total_disp;
        
        % 6. Ghost Update
        v_rand_gp = randn(S.N, 3) * S.stdx;
        if S.potential ~= 0
            v_rand_gp(idxgp,:) = v_rand_gp(idxgp,:) + potdispsgp;
        end
        pgp_norm = vecnorm(pgp, 2, 2) + eps;
        pgp_dir  = pgp ./ pgp_norm;
        v_rad_comp = sum(v_rand_gp .* pgp_dir, 2);
        v_tan_gp   = v_rand_gp - (v_rad_comp .* pgp_dir);
        pgp_temp = pgp + v_tan_gp;
        pgp_temp_norm = vecnorm(pgp_temp, 2, 2) + eps;
        pgp_next_dir  = pgp_temp ./ pgp_temp_norm;
        p2rho = vecnorm(p2, 2, 2);
        pgp2  = pgp_next_dir .* (2*S.br - p2rho);
        
        % 7. Reset
        if S.potential == 0
            distpp = pdist2(p2, [p2; pgp2], 'squaredeuclidean');
            idxd = distpp > 0 & distpp < (2*S.rp)^2; 
            [r_idx, c_idx] = find(idxd);
            c_idx(c_idx > S.N) = c_idx(c_idx > S.N) - S.N;
            resetters = unique([r_idx; c_idx]);
            if ~isempty(resetters)
                p2(resetters, :) = p(resetters, :);
                pgp2(resetters, :) = pgp(resetters, :);
            end
        end
        
        % 8. Handover
        p2rho = vecnorm(p2,2,2);
        pgp2rho = vecnorm(pgp2,2,2);
        idx_swap = p2rho > pgp2rho;
        p_next = p2; pgp_next = pgp2;
        p_next(idx_swap, :) = pgp2(idx_swap, :);
        pgp_next(idx_swap, :) = p2(idx_swap, :);
        p = p_next; pgp = pgp_next;
        
        if qs > nsteps
            break
        end
    end

    ASYMCORR.correction = [sgd_centers, sgd_correction];
    save(filenamecorrection, 'ASYMCORR', 'sgd_edges');
    save(filestartingconfiguration, 'p', 'pgp', 'S');
    disp('SGD Setup Complete.');
end