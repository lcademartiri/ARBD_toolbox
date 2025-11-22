function [p,pgp,sgd_correction,sgd_edges,history] = sbc_setup_sgd_v5(S,PDF,nsteps,opts)
% SBC_SETUP_SGD_V5
% Upgraded SGD control for learning radial boundary correction under SBC.
% - Progressive multi-step refinement (3x multipliers)
% - Gain scaling ~ 1/sqrt(batch_multiplier)
% - Warm-up after each refinement stage
% - Tightened thermalization: wait 10 * relaxsteps after expansion
% - Plotting: black background, thick axes, bold fonts
%
% Usage:
%   [p,pgp,sgd_correction,sgd_edges,history] = sbc_setup_sgd_v5(S,PDF,nsteps)
% or with options:
%   opts.base_gain, opts.batch_multiplier_sequence, opts.ramp_factor, ...
%
% Based on your previous v4 implementation (reference used during edit). :contentReference[oaicite:1]{index=1}

if nargin < 4, opts = struct(); end

% -------------------------- Defaults & derived params --------------------
% Physical / derived
gCS = (1 - S.phi/2) / (1 - S.phi)^3;
diffE = S.esdiff * S.alpha / gCS;
tau_alpha = (S.rp^2) / (6 * diffE);
relaxsteps = ceil(tau_alpha / S.timestep);

% Tightened thermalization: require at least 10 * relaxsteps after expansion
min_therm_steps = max( ceil(10 * relaxsteps), 2000 ); % safety lower bound

% SGD defaults (tunable via opts)
if isfield(opts,'base_batch_mult')
    base_batch_mult = opts.base_batch_mult;
else
    base_batch_mult = 10;
end

if isfield(opts,'batch_mult_sequence')
    batch_mult_sequence = opts.batch_mult_sequence;
else
    batch_mult_sequence = [3,3,3];
end
base_batch_size = ceil(10 * relaxsteps); % original design: 10 * relax
sgd_batch_size = base_batch_size;             % will be updated during phases

if isfield(opts,'base_gain')
    sgd_base_gain = opts.base_gain;
else
    sgd_base_gain = 0.1;
end

if isfield(opts,'gain_ramp')
    gain_ramp = opts.gain_ramp;
else
    gain_ramp = true;
end

if isfield(opts,'gain_ramp_steps')
    gain_ramp_steps = opts.gain_ramp_steps;
else
    gain_ramp_steps = 3;
end

if isfield(opts,'sgd_smooth_win')
    sgd_smooth_win = opts.sgd_smooth_win;
else
    sgd_smooth_win = 5;
end

if isfield(opts,'sgd_cap')
    sgd_cap = opts.sgd_cap;
else
    sgd_cap = 0.003 * S.rp;
end

potdepth = 2 * S.rc;
if 2*S.br - 2*potdepth < (10 * S.rp)
    potdepth = 2*S.br - 10*S.rp;
end
if 2*S.br - 2*S.rc < (10 * S.rp)
    potdepth = S.rc;
end

% controller flags
if isfield(opts,'metric_smoothing_param')
    metric_smoothing_param = opts.metric_smoothing_param;
else
    metric_smoothing_param = 0.8;
end

metric_smoothing_reset_on_stage = true; % we reset smoothed metric when stage changes
bold = true;
grace = true;
if isfield(opts,'graceperiod')
    graceperiod = opts.graceperiod;
else
    graceperiod = 20000;
end

if isfield(opts,'patience0')
    patience0 = opts.patience0;
else
    patience0 = 2;
end

% plotting / debugging
if isfield(opts,'debugging')
    debugging = opts.debugging;
else
    debugging = false;
end

if isfield(opts,'graphing')
    graphing = opts.graphing;
else
    graphing = true;
end

% series/file names (kept from v4 for compatibility)
seriesname = 'sgd_msp03_v5';
if S.potential==1, potname='lj'; elseif S.potential==2, potname='wca'; else potname='hs'; end
filenamecorrection = sprintf(['ASYMCORR_',seriesname,'_%s_%.0e_%.0e_%.0f_%.1f_%.1e.mat'],...
    potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
filestartingconfiguration = sprintf(['START_',potname,'_%.0e_%.0e_%.0f_%.1f_%.1e.mat'],...
    potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
filepdfdenom = sprintf('PDFdenom_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);

% quick load check (same behavior as before)
if exist(filenamecorrection,'file') && exist(filestartingconfiguration,'file')
    load(filenamecorrection); load(filestartingconfiguration);
    % return loaded variables (keep signature)
    history = [];
    return
end

% -------------------------- Pre-calculations ------------------------------
if exist(filepdfdenom,'file')
    load(filepdfdenom,'gdenominator');
else
    gdenominator = PDFdenom(S,PDF,nsteps); %#ok<NASGU>
    save(filepdfdenom,'gdenominator','S');
end

if S.potential~=0
    H = pot_force(S.potential,S.rc,30000,S.pot_sigma,S.pot_epsilon);
    clamp = mcdClamp(1e6,S.rp,normrnd(0,S.stdx,1e6,3),S.esdiff,S.timestep,H,S.kbT);
    H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'linear', 'nearest');
end

% -------------------------- Initial placement (FCC) -----------------------
disp('Creating initial FCC lattice.')
flag = 1;
if debugging, rng(100); end

while flag==1
    basis=[0,0.7071,0.7071;0.7071,0,0.7071;0.7071,0.7071,0].*(2.01*S.rp);
    maxsteps=2*ceil(((S.br*2)*sqrt(3))/(2*S.rp));
    templist=double(linspace(-maxsteps,maxsteps,2*maxsteps+1)');
    [x1,x2,x3] = meshgrid(templist,templist,templist);
    possiblepositions = x1(:).*basis(1,:)+x2(:).*basis(2,:)+x3(:).*basis(3,:);
    possiblepositions = bsxfun(@plus,possiblepositions,-S.rp*[1,1,1]./vecnorm([1,1,1],2));
    tempnorms = vecnorm(possiblepositions,2,2);
    possiblepositions(tempnorms > (S.br - 2*S.rp), :) = [];
    possiblepositions = possiblepositions(randperm(size(possiblepositions,1)), :);
    p = possiblepositions(1:S.N, :);
    D = pdist(p)';
    if sum(D < (2*S.rp)) == 0, flag = 0; end
end
clear possiblepositions tempnorms x1 x2 x3 templist maxsteps basis D

% -------------------------- SGD INIT -------------------------------------
% edges and bins (unchanged logic, but make centers robust)
sgd_edges = sort((S.br:-0.02*S.rp:S.br - potdepth)');
sgd_bins = numel(sgd_edges) - 1;
sgd_centers = sgd_edges(1:end-1) + diff(sgd_edges)/2;

% correction vector and interpolant
sgd_correction = zeros(sgd_bins, 1);
F_corr_interp = griddedInterpolant(sgd_centers, sgd_correction, 'linear', 'nearest');

% accumulators
batch_sum_drift = zeros(sgd_bins,1);
batch_counts = zeros(sgd_bins,1);
steps_in_batch = 0;

% controller state
best_pdf_metric = inf;
best_correction = sgd_correction;
best_p_snapshot = p;
best_pgp_snapshot = [];

patience_counter = 0;
pdf_metric = 0;

PHASE_SEARCH = 1;
PHASE_REFINE = 2;
current_phase = PHASE_SEARCH;

% Stage management (progressive refinement)
stage_index = 0;
max_stages = numel(batch_mult_sequence);
stage_warmup = false;         % true while filling the new big batch
stage_steps_filled = 0;
stage_batch_multiplier = 1;   % current multiplier
stage_target_for_1e5 = [];    % will compute expected noise floor for 1e5 steps

% compute the expected noise-floor for 1e5 timesteps (so we can stop when reached)
w_count_1e5 = 1e5;
valid_mask_global = PDF.centers{3} > 2*(S.br - potdepth) & PDF.centers{3} < 2*(S.br - 0.5*S.rp);
if any(valid_mask_global)
    expected_counts_1e5 = gdenominator(valid_mask_global) * w_count_1e5;
    expected_counts_1e5(expected_counts_1e5==0) = inf;
    bin_noise_sigma_1e5 = 1 ./ sqrt(expected_counts_1e5);
    target_1e5 = sqrt(mean(bin_noise_sigma_1e5.^2));
else
    target_1e5 = 1e-3; % fallback small number
end

% initial batch size & gain
sgd_batch_size = base_batch_size;
sgd_base_gain = sgd_base_gain;
sgd_gain = sgd_base_gain;

fprintf('SGD Init. Phase: SEARCH. Batch: %d. MinThermSteps: %d\n', sgd_batch_size, min_therm_steps);

% -------------------------- GRAPHING INIT --------------------------------
history = struct('steps',[],'pdf_dev',[],'pdf_smooth',[],'max_corr',[],'gain',[],'batch_size',[],'stage',[]);
if graphing
    ndens.av_window = 1000;
    ndens.edges = sort((S.br:-0.02*S.rp:0)');
    ndens.centers = ndens.edges(1:end-1) + diff(ndens.edges)/2;
    ndens.counts = zeros(numel(ndens.centers),1);
    ndens.counts_storage = ndens.counts;
    ndens.vols = (4/3)*pi*(ndens.edges(2:end).^3 - ndens.edges(1:end-1).^3);
    ndens.ndens0 = (S.N / S.bv);
    pdf.pre.counts = zeros(numel(PDF.pdfedges{3})-1,1);

    f_fig = figure('Units','normalized','Position',[0.05 0.05 0.8 0.85]);
    % style: black bg, thick axes later
    set(gcf,'Color','k');

    ax_dens = subplot(2,2,1);
    ax_pdf  = subplot(2,2,2);
    ax_conv = subplot(2,2,3);
    ax_ctrl = subplot(2,2,4);

    % set default axes style
    axs = [ax_dens, ax_pdf, ax_conv, ax_ctrl];
    for a = axs
        set(a,'Color','k','XColor','w','YColor','w','LineWidth',3);
        set(get(a,'Title'),'FontWeight','bold');
        set(get(a,'XLabel'),'FontWeight','bold');
        set(get(a,'YLabel'),'FontWeight','bold');
        set(a,'FontWeight','bold');
    end
end

% -------------------------- MAIN LOOP ------------------------------------
qs = 0;
thermflag = 0;
r2_uniform = 3/5 * S.br^2;
pgp = p - (2*S.br).*(p ./ vecnorm(p,2,2));
reverseStr = '';
pdf_metric_val = 10;

disp('Starting SGD Evolution...')

% make sure variables that might be used later are defined
idxgp = false(S.N,1);

while true
    qs = qs + 1;

    prho = vecnorm(p, 2, 2);
    pvers = p ./ (prho + eps);
    idxgp = prho > (S.br - S.rc);

    % --- THERMALIZATION CHECK (tightened) ------------------------------
    if thermflag == 0
        spread_ratio = mean(prho.^2) / r2_uniform;
        is_expanded = spread_ratio > 0.99 || (qs > 2000 && abs(spread_ratio - 1) < 0.01);
        is_relaxed = qs > min_therm_steps;
        if is_expanded && is_relaxed
            thermflag = 1;
            qs = 1; % reset internal step counter used for batches
            disp('--- Thermalization Complete (tight criterion) ---');
        end
    end

    % 3. Base Displacement (potentials handled same as before)
    if S.potential ~= 0
        ptemp = [p; pgp(idxgp,:)];
        all_potdisps = potential_displacements_v2(ptemp, size(ptemp,1), S.rc, H, H_interpolant, S.esdiff, clamp, S.kbT, S.stdx, S.timestep, 1);
        potdisps = all_potdisps(1:S.N, :);
        potdispsgp = all_potdisps(S.N+1:end, :);
        v_rand = randn(S.N, 3) * S.stdx;
        base_disp = v_rand + potdisps;
    else
        v_rand = randn(S.N, 3) * S.stdx;
        base_disp = v_rand;
        potdispsgp = zeros(sum(idxgp), 3);
    end

    % --- SGD ACCUMULATION ----------------------------------------------
    if thermflag == 1
        % 1. Accumulate drift into bins (use discretize, robust)
        dr_raw = sum(base_disp .* pvers, 2);

        % bin indices via discretize (more robust than histc)
        bin_idx = discretize(prho, sgd_edges);
        valid_mask = bin_idx > 0 & bin_idx <= sgd_bins;

        if any(valid_mask)
            % accumarray needs indices >=1; patch accordingly
            idxs = bin_idx(valid_mask);
            new_sums = accumarray(idxs, dr_raw(valid_mask), [sgd_bins, 1]);
            new_counts = accumarray(idxs, 1, [sgd_bins, 1]);
            batch_sum_drift = batch_sum_drift + new_sums;
            batch_counts = batch_counts + new_counts;
        end

        % 2. Accumulate Diagnostics (Every step)
        if graphing
            [hc, ~] = histcounts(vecnorm(p,2,2), ndens.edges);
            ndens.counts = ndens.counts + hc';
            pairdists = pdist(p);
            [hc_pdf, ~] = histcounts(pairdists, PDF.pdfedges{3});
            % ensure same length
            if numel(hc_pdf) == numel(pdf.pre.counts)
                pdf.pre.counts = pdf.pre.counts + hc_pdf';
            else
                % fallback: reinitialize
                pdf.pre.counts = zeros(numel(PDF.pdfedges{3})-1,1);
            end
        end

        steps_in_batch = steps_in_batch + 1;

        % --- BATCH UPDATE & CONTROL TRIGGER --------------------------
        if steps_in_batch >= sgd_batch_size
            % A. SGD Update (only update bins with enough counts)
            has_data = batch_counts > 10;
            if any(has_data)
                batch_mean_drift = batch_sum_drift(has_data) ./ batch_counts(has_data);
                sgd_correction(has_data) = sgd_correction(has_data) - sgd_gain * batch_mean_drift;
                sgd_correction = max(min(sgd_correction, sgd_cap), -sgd_cap);
                sgd_correction = smoothdata(sgd_correction, 'movmean', sgd_smooth_win);
                F_corr_interp.Values = sgd_correction;
            end

            % B. Diagnostics & Control (compute auto-target and pdf metric)
            if graphing
                w_count = steps_in_batch;
                curr_g = (pdf.pre.counts / w_count) ./ gdenominator;

                valid_mask = PDF.centers{3} > 2*(S.br - potdepth) & PDF.centers{3} < 2*(S.br - 0.5*S.rp);

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

                % Smoothing (optionally reset on stage changes)
                if pdf_metric == 0 || metric_smoothing_reset_on_stage && stage_warmup
                    pdf_metric = raw_pdf_metric;
                else
                    pdf_metric = (1 - metric_smoothing_param) * pdf_metric + metric_smoothing_param * raw_pdf_metric;
                end
                pdf_metric_val = pdf_metric;
            else
                raw_pdf_metric = 1;
                convergence_target = 1;
            end

            % ----------------- Controller logic (Search -> Stage progression) -----------------
            % If metric improved, keep best snapshot
            if pdf_metric < best_pdf_metric
                best_pdf_metric = pdf_metric;
                best_correction = sgd_correction;
                best_p_snapshot = p;
                best_pgp_snapshot = pgp;
                patience_counter = 0;

                % Bold driver: only in search phases and if far from noise
                if bold && current_phase == PHASE_SEARCH && pdf_metric > 2 * convergence_target
                    sgd_base_gain = sgd_base_gain * 1.05;
                end

                % Trigger stage entry: enter refinement/stage progression only if
                %   - currently in SEARCH, 
                %   - pdf_metric < 2 * convergence_target (coarse trigger)
                %   - and we passed graceperiod
                if current_phase == PHASE_SEARCH && pdf_metric < (convergence_target * 2.0) && qs > graceperiod
                    % start multi-stage refinement sequence:
                    stage_index = stage_index + 1;
                    if stage_index <= max_stages
                        fprintf('\n=== ENTERING REFINEMENT STAGE %d ===\n', stage_index);
                        % restore best correction/state
                        sgd_correction = best_correction;
                        F_corr_interp.Values = sgd_correction;
                        p = best_p_snapshot;
                        pgp = best_pgp_snapshot;

                        % update batch multiplier progressively
                        stage_batch_multiplier = batch_mult_sequence(stage_index);
                        sgd_batch_size = stage_batch_multiplier * base_batch_size;

                        % scale gain by 1/sqrt(batch_multiplier) (less aggressive)
                        stage_scale = 1 / sqrt(stage_batch_multiplier);
                        if gain_ramp
                            % we will ramp in multiple steps rather than instantaneous
                            ramp_factors = linspace(1, stage_scale, gain_ramp_steps+1);
                            % keep first entry = 1 (no change), ramp later inside stage
                            ramp_factors = ramp_factors(2:end);
                            sgd_ramp_queue = ramp_factors;
                        else
                            sgd_base_gain = sgd_base_gain * stage_scale;
                            sgd_ramp_queue = [];
                        end
                        sgd_gain = sgd_base_gain;

                        % Clear accumulators and set warmup flag
                        batch_sum_drift(:) = 0; batch_counts(:) = 0;
                        ndens.counts(:) = 0; pdf.pre.counts(:) = 0;
                        steps_in_batch = 0;
                        stage_warmup = true;
                        stage_steps_filled = 0;
                        current_phase = PHASE_REFINE;
                        continue; % skip plotting until warm
                    else
                        % if beyond configured stages, just go to a single refine step
                        fprintf('\n=== ENTERING FINAL REFINEMENT ===\n');
                        sgd_correction = best_correction;
                        F_corr_interp.Values = sgd_correction;
                        p = best_p_snapshot;
                        pgp = best_pgp_snapshot;

                        % big refine to target_1e5 (cap multiplier)
                        target_mult = ceil((w_count_1e5 / base_batch_size));
                        % but to be safe limit multiplier
                        target_mult = min(target_mult, 1000);
                        sgd_batch_size = target_mult * base_batch_size;
                        sgd_base_gain = sgd_base_gain * (1 / sqrt(target_mult));
                        sgd_gain = sgd_base_gain;

                        batch_sum_drift(:) = 0; batch_counts(:) = 0;
                        ndens.counts(:) = 0; pdf.pre.counts(:) = 0;
                        steps_in_batch = 0;
                        stage_warmup = true;
                        stage_steps_filled = 0;
                        current_phase = PHASE_REFINE;
                        continue;
                    end
                end
            else
                % no improvement: patience logic (gain cut)
                if pdf_metric > (best_pdf_metric * 1.3) && grace && qs > graceperiod
                    patience_counter = patience_counter + 1;
                    if patience_counter > patience0
                        if sgd_base_gain > 1e-6
                            sgd_base_gain = sgd_base_gain * 0.5;
                            fprintf('--- Gain slashed to %.1e ---\n', sgd_base_gain);
                        end
                        patience_counter = 0;
                        best_pdf_metric = pdf_metric;
                    end
                end
            end

            % ----------------- Gain scaling inside phase ------------------
            if current_phase == PHASE_SEARCH
                % keep adaptive scaling, but no less than floor
                target_ratio = pdf_metric / 0.05;
                sgd_gain = min(sgd_base_gain, sgd_base_gain * target_ratio);
                sgd_gain = max(sgd_gain, 1e-6);
            else
                % in REFINE: use sgd_base_gain as set, but if we have ramp entries, consume them
                if exist('sgd_ramp_queue','var') && ~isempty(sgd_ramp_queue)
                    % apply one ramp step per filled stage batch
                    if stage_warmup == false && stage_steps_filled >= sgd_batch_size
                        next_scale = sgd_ramp_queue(1);
                        sgd_base_gain = sgd_base_gain * next_scale;
                        sgd_ramp_queue(1) = [];
                        fprintf('Ramping gain: new base_gain = %.2e\n', sgd_base_gain);
                    end
                end
                sgd_gain = sgd_base_gain;
            end

            % ----------------- Plotting & History -------------------------
            history.steps(end+1) = qs;
            history.pdf_dev(end+1) = raw_pdf_metric;
            history.pdf_smooth(end+1) = pdf_metric;
            history.max_corr(end+1) = max(abs(sgd_correction));
            history.gain(end+1) = sgd_gain;
            history.batch_size(end+1) = sgd_batch_size;
            history.stage(end+1) = stage_index;

            if graphing
                set(0, 'CurrentFigure', f_fig);
                curr_ndens_norm = 100 * (((ndens.counts / w_count)./ndens.vols) ./ ndens.ndens0);

                subplot(ax_dens);
                plot(ndens.centers, curr_ndens_norm, 'LineWidth', 2, 'Color', 'w');
                xline(S.br - S.rp, '--w'); xline(S.br, '-r');
                ylim([80 120]); xlim([0 S.br]);
                title(sprintf('Density Dev: %.2f%%', max(abs(curr_ndens_norm-100))));
                set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 3);
                set(gca,'FontWeight','bold');

                subplot(ax_pdf);
                plot(PDF.centers{3}, curr_g, 'LineWidth', 2, 'Color', 'y');
                yline(1, '--w');
                xlim([0 2.1*S.br]); ylim([0.5 1.5]);
                title(sprintf('PDF RMS: %.4f', raw_pdf_metric));
                set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 3);
                set(gca,'FontWeight','bold');

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
                set(gca, 'Color', 'k', 'XColor', 'w', 'LineWidth', 3);
                set(gca,'FontWeight','bold');
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
                set(gca, 'Color', 'k', 'XColor', 'w', 'LineWidth', 3);
                set(gca,'FontWeight','bold');

                drawnow;
                msg = sprintf('Step %d | RMS: %.4f (Target: %.4f) | Gain: %.1e | Batch: %d', qs, pdf_metric, convergence_target, sgd_gain, sgd_batch_size);
                fprintf([reverseStr, msg]);
                reverseStr = repmat('\b', 1, length(msg));
            end

            % ----------------- CONVERGENCE CHECK (stop conditions) -----------
            % If we have reached the target defined by 1e5 expected noise floor AND we
            % are in a refine stage, stop
            if pdf_metric <= target_1e5 && current_phase == PHASE_REFINE
                disp('--- CONVERGED: reached 1e5 noise-floor target ---');
                break;
            end
            % also a fallback: if pdf_metric <= convergence_target and large enough qs
            if pdf_metric <= convergence_target && qs > 20000 && current_phase == PHASE_REFINE
                disp('--- CONVERGED (Hit local noise floor) ---');
                break;
            end
            % Reset accumulators (end of batch)
            batch_sum_drift(:) = 0;
            batch_counts(:) = 0;
            ndens.counts(:) = 0;
            pdf.pre.counts(:) = 0;
            steps_in_batch = 0;

            % If we were in warmup, clear warmup now (we have a filled stage batch)
            if stage_warmup
                stage_warmup = false;
                stage_steps_filled = 0;
                % reset smoothed metric history to avoid long-memory conflict
                if metric_smoothing_reset_on_stage
                    pdf_metric = 0;
                end
            end
        end
    end

    % ----------------- 4. Apply Correction & Integrate ------------------------
    dr_corr_mag = F_corr_interp(prho);
    core_mask = prho < (S.br - potdepth);
    dr_corr_mag(core_mask) = 0;
    total_disp = base_disp + (dr_corr_mag .* pvers);

    % 5. Position Update
    p2 = p + total_disp;

    % 6. Ghost update (same as before)
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

    % 7. Reset on overlaps if HS
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

    % 8. Handover (same as before)
    p2rho = vecnorm(p2,2,2);
    pgp2rho = vecnorm(pgp2,2,2);
    idx_swap = p2rho > pgp2rho;
    p_next = p2; pgp_next = pgp2;
    p_next(idx_swap, :) = pgp2(idx_swap, :);
    pgp_next(idx_swap, :) = p2(idx_swap, :);
    p = p_next; pgp = pgp_next;

    % termination guard (safety)
    if qs > nsteps
        disp('Reached nsteps limit; exiting.');
        break;
    end
end

% Save correction and starting configuration
ASYMCORR.correction = [sgd_centers, sgd_correction];
save(filenamecorrection, 'ASYMCORR', 'sgd_edges');
save(filestartingconfiguration, 'p', 'pgp', 'S');

disp('SGD Setup Complete (v5).');
end
