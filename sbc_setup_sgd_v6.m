function [p,pgp,sgd_correction,sgd_edges,history] = sbc_setup_sgd_v6(S,PDF,nsteps,opts)
% SBC_SETUP_SGD_V5
% Final v5 with:
% - SNR gating (soft/hard)
% - Plateau detection & plateau protection (no plateau from inactive updates)
% - SNR-starvation safety Mode C (adaptive relaxation only during refine)
% - Per-update clipping (relative & absolute)
% - Reset SNR thresholds at start of each refinement stage (Option A)
% - Diagnostics & black-figure plotting with bright colored curves; axes frames/titles/labels white
%
% Usage:
%   [p,pgp,sgd_correction,sgd_edges,history] = sbc_setup_sgd_v5(S,PDF,nsteps,opts)
%
% opts fields and defaults (selected):
%   batch_mult_sequence = [3,3,3]
%   base_gain = 0.1
%   gain_ramp = true
%   gain_ramp_steps = 3
%   sgd_smooth_win = 5
%   sgd_cap = 0.003 * S.rp
%   metric_smoothing_param = 0.8
%   graceperiod = 20000
%   patience0 = 2
%   debugging = false
%   graphing = true
%   snr_th_search = 3.0
%   snr_th_refine = 1.5
%   snr_target = 5.0
%   n_min = 20
%   clip_frac = 0.2
%   abs_cap_frac = 0.005
%   plateau_window_batches = 5
%   plateau_tol = 5e-3
%   plateau_patience = 3
%   use_soft_snr = true
%   starvation_relax_trigger_batches = 3
%   starvation_fallback_batches = 6
%
% Notes: relies on helper functions present in your environment.

if nargin < 4, opts = struct(); end

% -------------------- Option parsing (MATLAB-safe) -----------------------
if isfield(opts,'batch_mult_sequence'),  batch_mult_sequence = opts.batch_mult_sequence; else batch_mult_sequence = [3,3,3]; end
if isfield(opts,'base_gain'),            sgd_base_gain_default = opts.base_gain; else sgd_base_gain_default = 0.2; end
if isfield(opts,'gain_ramp'),            gain_ramp = opts.gain_ramp; else gain_ramp = true; end
if isfield(opts,'gain_ramp_steps'),      gain_ramp_steps = opts.gain_ramp_steps; else gain_ramp_steps = 3; end
if isfield(opts,'sgd_smooth_win'),       sgd_smooth_win = opts.sgd_smooth_win; else sgd_smooth_win = 5; end
if isfield(opts,'sgd_cap'),              sgd_cap = opts.sgd_cap; else sgd_cap = 0.0015 * S.rp; end
if isfield(opts,'metric_smoothing_param'), metric_smoothing_param = opts.metric_smoothing_param; else metric_smoothing_param = 0.8; end
if isfield(opts,'graceperiod'),          graceperiod = opts.graceperiod; else graceperiod = 20000; end
if isfield(opts,'patience0'),            patience0 = opts.patience0; else patience0 = 2; end
if isfield(opts,'debugging'),            debugging = opts.debugging; else debugging = false; end
if isfield(opts,'graphing'),             graphing = opts.graphing; else graphing = true; end
if isfield(opts,'enable_io'),            enable_io = opts.enable_io; else enable_io = true; end

% SNR gating & starvation params
if isfield(opts,'snr_th_search'),        snr_th_search_default = opts.snr_th_search; else snr_th_search_default = 3.0; end
if isfield(opts,'snr_th_refine'),        snr_th_refine_default = opts.snr_th_refine; else snr_th_refine_default = 1.5; end
if isfield(opts,'snr_target'),           snr_target_default = opts.snr_target; else snr_target_default = 5.0; end
if isfield(opts,'n_min'),                n_min = opts.n_min; else n_min = 20; end
if isfield(opts,'clip_frac'),            clip_frac = opts.clip_frac; else clip_frac = 0.3; end
if isfield(opts,'abs_cap_frac'),         abs_cap_frac = opts.abs_cap_frac; else abs_cap_frac = 0.005; end
if isfield(opts,'plateau_window_batches'), plateau_window_batches = opts.plateau_window_batches; else plateau_window_batches = 5; end
if isfield(opts,'plateau_tol'),          plateau_tol = opts.plateau_tol; else plateau_tol = 5e-3; end
if isfield(opts,'plateau_patience'),     plateau_patience = opts.plateau_patience; else plateau_patience = 3; end
if isfield(opts,'use_soft_snr'),         use_soft_snr = opts.use_soft_snr; else use_soft_snr = true; end
if isfield(opts,'starvation_relax_trigger_batches'), starvation_relax_trigger_batches = opts.starvation_relax_trigger_batches; else starvation_relax_trigger_batches = 3; end
if isfield(opts,'starvation_fallback_batches'), starvation_fallback_batches = opts.starvation_fallback_batches; else starvation_fallback_batches = 6; end

% -------------------- Derived physical params ----------------------------
gCS = (1 - S.phi/2) / (1 - S.phi)^3;
diffE = S.esdiff * S.alpha / gCS;
tau_alpha = (S.rp^2) / (6 * diffE);
relaxsteps = ceil(tau_alpha / S.timestep);
noise_prefactor = 2.0;
% tightened thermalization: require 10 * relaxsteps or minimum 2000
min_therm_steps = max( ceil(10 * relaxsteps), 2000 );

% -------------------- filenames / quick loads ----------------------------
if S.potential==1, potname='lj'; elseif S.potential==2, potname='wca'; else potname='hs'; end
seriesname = 'sgd_msp03_v5_final';
filenamecorrection = sprintf(['ASYMCORR_',seriesname,'_%s_%.0e_%.0e_%.0f_%.1f_%.1e.mat'],...
    potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
filestartingconfiguration = sprintf(['START_',potname,'_%.0e_%.0e_%.0f_%.1f_%.1e.mat'],...
    potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
filepdfdenom = sprintf('PDFdenom_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);

if enable_io && exist(filenamecorrection,'file') && exist(filestartingconfiguration,'file')
    load(filenamecorrection); load(filestartingconfiguration);
    history = [];
    return
end

% -------------------- Pre-calc helpers -----------------------------------
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

% -------------------- initial placement (FCC-like) -----------------------
disp('Creating initial FCC-like lattice...');
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

% -------------------- SGD init ------------------------------------------
potdepth = 2 * S.rc;
if 2*S.br - 2*potdepth < (10 * S.rp)
    potdepth = 2*S.br - 10*S.rp;
end

sgd_edges = sort((S.br:-0.02*S.rp:S.br - potdepth)');
sgd_bins = numel(sgd_edges) - 1;
sgd_centers = sgd_edges(1:end-1) + diff(sgd_edges)/2;

sgd_correction = zeros(sgd_bins, 1);
F_corr_interp = griddedInterpolant(sgd_centers, sgd_correction, 'linear', 'nearest');

batch_sum_drift = zeros(sgd_bins,1);
batch_sum_drift_sq = zeros(sgd_bins,1);
batch_counts = zeros(sgd_bins,1);
steps_in_batch = 0;

% controller state
best_pdf_metric = inf;
best_correction = sgd_correction;
best_p_snapshot = p;
best_pgp_snapshot = [];
patience_counter = 0;
pdf_metric = 0;

% phases & stages
PHASE_SEARCH = 1; PHASE_REFINE = 2;
current_phase = PHASE_SEARCH;
stage_index = 0;
max_stages = numel(batch_mult_sequence);
stage_warmup = false;
stage_steps_filled = 0;
stage_batch_multiplier = 1;

% SNR thresholds: will be reset at stage start (Option A)
snr_th_search = snr_th_search_default;
snr_th_refine = snr_th_refine_default;
snr_target = snr_target_default;
sgd_base_gain = sgd_base_gain_default;
sgd_gain = sgd_base_gain;

% starvation counters
no_update_counter = 0;
starvation_relax_count = 0;
no_update_history = [];

% plateau tracking
max_corr_queue = [];
plateau_counter = 0;

% expected target for 1e5 timesteps
w_count_1e5 = 1e5;
valid_mask_global = PDF.centers{3} > 2*(S.br - potdepth) & PDF.centers{3} < 2*(S.br - 0.5*S.rp);
if any(valid_mask_global)
    expected_counts_1e5 = gdenominator(valid_mask_global) * w_count_1e5;
    expected_counts_1e5(expected_counts_1e5==0) = inf;
    bin_noise_sigma_1e5 = 1 ./ sqrt(expected_counts_1e5);
    target_1e5 = sqrt(mean(bin_noise_sigma_1e5.^2));
else
    target_1e5 = 1e-3;
end

base_batch_size = ceil(10 * relaxsteps);
sgd_batch_size = base_batch_size;

fprintf('SGD init. Phase: SEARCH. Batch size: %d. MinThermSteps: %d\n', sgd_batch_size, min_therm_steps);

% -------------------- plotting & diagnostics init -----------------------
history = struct('steps',[],'pdf_dev',[],'pdf_smooth',[],'max_corr',[],'gain',[],...
    'batch_size',[],'stage',[],'fraction_updated',[],'avg_snr',[],'median_snr',[],...
    'bins_with_data',[],'bins_updated',[],'no_update_counter',[],'starvation_relax_count',[]);

if graphing
    ndens.edges = sort((S.br:-0.02*S.rp:0)');
    ndens.centers = ndens.edges(1:end-1) + diff(ndens.edges)/2;
    ndens.counts = zeros(numel(ndens.centers),1);
    ndens.vols = (4/3)*pi*(ndens.edges(2:end).^3 - ndens.edges(1:end-1).^3);
    ndens.ndens0 = (S.N / S.bv);
    pdf.pre.counts = zeros(numel(PDF.pdfedges{3})-1,1);

    f_fig = figure('Units','normalized','Position',[0.05 0.05 0.85 0.85]);
    set(gcf,'Color','k');

    ax_dens = subplot(3,2,1);
    ax_pdf  = subplot(3,2,2);
    ax_conv = subplot(3,2,3);
    ax_ctrl = subplot(3,2,4);
    ax_diag = subplot(3,2,5);
    ax_snr  = subplot(3,2,6);

    axs = [ax_dens, ax_pdf, ax_conv, ax_ctrl, ax_diag, ax_snr];
    for a = axs
        set(a,'Color','k','XColor','w','YColor','w','LineWidth',3);
        set(get(a,'Title'),'FontWeight','bold','Color','w');
        set(get(a,'XLabel'),'FontWeight','bold','Color','w');
        set(get(a,'YLabel'),'FontWeight','bold','Color','w');
        set(a,'FontWeight','bold');
        set(a,'FontSize',10);
    end
end

% -------------------- main loop -----------------------------------------
qs = 0;
thermflag = 0;
r2_uniform = 3/5 * S.br^2;
pgp = p - (2*S.br).*(p ./ (vecnorm(p,2,2) + eps));
reverseStr = '';
pdf_metric_val = 10;

disp('Starting SGD evolution...');

while true
    qs = qs + 1;

    prho = vecnorm(p, 2, 2);
    pvers = p ./ (prho + eps);
    idxgp = prho > (S.br - S.rc);

    % --- thermalisation check (tight) -----------------------------------
    if thermflag == 0
        spread_ratio = mean(prho.^2) / r2_uniform;
        is_expanded = spread_ratio > 0.99 || (qs > 2000 && abs(spread_ratio - 1) < 0.01);
        is_relaxed = qs > min_therm_steps;
        if is_expanded && is_relaxed
            thermflag = 1;
            qs = 1;
            disp('--- Thermalization complete (tight) ---');
        end
    end

    % --- base displacement & potential -----------------------------------
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

    % --- SGD accumulation -----------------------------------------------
    if thermflag == 1
        dr_raw = sum(base_disp .* pvers, 2);

        bin_idx = discretize(prho, sgd_edges);
        valid_mask = bin_idx > 0 & bin_idx <= sgd_bins;

        if any(valid_mask)
            idxs = bin_idx(valid_mask);
            values = dr_raw(valid_mask);
            new_sums = accumarray(idxs, values, [sgd_bins, 1]);
            new_sums_sq = accumarray(idxs, values.^2, [sgd_bins, 1]);
            new_counts = accumarray(idxs, 1, [sgd_bins, 1]);

            batch_sum_drift = batch_sum_drift + new_sums;
            batch_sum_drift_sq = batch_sum_drift_sq + new_sums_sq;
            batch_counts = batch_counts + new_counts;
        end

        if graphing
            [hc, ~] = histcounts(vecnorm(p,2,2), ndens.edges);
            ndens.counts = ndens.counts + hc';
            pairdists = pdist(p);
            [hc_pdf, ~] = histcounts(pairdists, PDF.pdfedges{3});
            if numel(hc_pdf) == numel(pdf.pre.counts)
                pdf.pre.counts = pdf.pre.counts + hc_pdf';
            else
                pdf.pre.counts = zeros(numel(PDF.pdfedges{3})-1,1);
            end
        end

        steps_in_batch = steps_in_batch + 1;

        % --- when batch full: process ----------------------------------
        if steps_in_batch >= sgd_batch_size
            % per-bin stats
            has_data = batch_counts > 0;
            bin_mean = zeros(sgd_bins,1);
            bin_var = zeros(sgd_bins,1);
            stderr = zeros(sgd_bins,1);
            snr = zeros(sgd_bins,1);
            n_i = batch_counts;

            idxs = find(has_data);
            for ii = idxs'
                ni = n_i(ii);
                mu = batch_sum_drift(ii) / ni;
                ss = batch_sum_drift_sq(ii);
                if ni > 1
                    v = max(0, (ss - ni * mu^2) / (ni - 1));
                else
                    v = 0;
                end
                se = sqrt(v) / sqrt(max(1, ni));
                s = abs(mu) / (se + eps);

                bin_mean(ii) = mu;
                bin_var(ii) = v;
                stderr(ii) = se;
                snr(ii) = s;
            end

            % choose SNR thresholds depending on phase
            if current_phase == PHASE_SEARCH
                snr_th = snr_th_search;
            else
                snr_th = snr_th_refine;
            end

            % per-bin adaptive gain & gating
            bins_to_update = false(sgd_bins,1);
            delta = zeros(sgd_bins,1);
            for ii = 1:sgd_bins
                if n_i(ii) < n_min
                    lambda = 100;
                    if n_i(ii) == 0
                        mu_shrunk = 0;
                    else
                        mu_shrunk = (n_i(ii)/(n_i(ii)+lambda)) * bin_mean(ii);
                    end
                    mu_eff = mu_shrunk;
                    s_eff = snr(ii) * (n_i(ii)/(n_i(ii)+lambda));
                else
                    mu_eff = bin_mean(ii);
                    s_eff = snr(ii);
                end

                if use_soft_snr
                    per_bin_factor = min(1, s_eff / snr_target);
                else
                    per_bin_factor = double(s_eff >= snr_th);
                end

                if per_bin_factor > 0
                    eta_i = sgd_gain * per_bin_factor;
                    delta_i = - eta_i * mu_eff;
                    delta(ii) = delta_i;
                    bins_to_update(ii) = true;
                end
            end

            % clipping caps
            % max_current_corr = max(abs(sgd_correction)) + eps;
            % max_delta_rel = clip_frac * max_current_corr;
            
            % Clipping caps
            % FIX: Add a 'floor' to the relative clipping to prevent stall at machine epsilon.
            % Using 0.1% of the hard cap (sgd_cap) allows the correction to grow linearly 
            % from zero until it establishes a physical magnitude.
            min_clip_floor = 5e-2 * sgd_cap; 
            
            % The valid base for relative clipping is the larger of:
            % 1. The actual current max correction
            % 2. The minimum floor (approx 3e-14 in your case)
            max_current_corr = max(max(abs(sgd_correction)), min_clip_floor);
            
            max_delta_rel = clip_frac * max_current_corr;
            abs_cap = abs_cap_frac * S.br;

            bins_updated = find(bins_to_update);

            % apply updates if any, else increment no_update counter
            if ~isempty(bins_updated)
                no_update_counter = 0;
                starvation_relax_count = 0;
                no_update_history = [no_update_history, 0];
            else
                no_update_counter = no_update_counter + 1;
                no_update_history = [no_update_history, 1];
            end

            for ii = bins_updated'
                d = delta(ii);
                d = sign(d) * min(abs(d), max_delta_rel);
                d = sign(d) * min(abs(d), abs_cap);
                sgd_correction(ii) = sgd_correction(ii) + d;
            end

            sgd_correction = max(min(sgd_correction, sgd_cap), -sgd_cap);
            sgd_correction = smoothdata(sgd_correction, 'movmean', sgd_smooth_win);
            F_corr_interp.Values = sgd_correction;

            % diagnostics & PDF metric
            if graphing
                w_count = steps_in_batch;
                curr_g = (pdf.pre.counts / max(1,w_count)) ./ gdenominator;

                valid_pdf_mask = PDF.centers{3} > 2*(S.br - potdepth) & PDF.centers{3} < 2*(S.br - 0.5*S.rp);
                expected_counts = gdenominator(valid_pdf_mask) * w_count;
                expected_counts(expected_counts==0) = inf;
                bin_noise_sigma = 1 ./ sqrt(expected_counts);
                convergence_target = noise_prefactor*sqrt(mean(bin_noise_sigma.^2));

                if any(valid_pdf_mask)
                    residuals = curr_g(valid_pdf_mask) - 1;
                    weights = gdenominator(valid_pdf_mask);
                    raw_pdf_metric = sqrt(sum(weights .* (residuals.^2)) / sum(weights));
                else
                    raw_pdf_metric = 1;
                end

                if pdf_metric == 0 || (stage_warmup && true)
                    pdf_metric = raw_pdf_metric;
                else
                    pdf_metric = (1 - metric_smoothing_param) * pdf_metric + metric_smoothing_param * raw_pdf_metric;
                end
            else
                raw_pdf_metric = 1;
                convergence_target = 1;
            end

            % update best snapshot & bold driver
            if pdf_metric < best_pdf_metric
                best_pdf_metric = pdf_metric;
                best_correction = sgd_correction;
                best_p_snapshot = p;
                best_pgp_snapshot = pgp;
                patience_counter = 0;

                if current_phase == PHASE_SEARCH
                    if pdf_metric > 2 * convergence_target && gain_ramp
                        sgd_base_gain = sgd_base_gain * 1.05;
                    end
                end

                if current_phase == PHASE_SEARCH && pdf_metric < (convergence_target * 2.0) && qs > graceperiod
                    stage_candidate_ready = true;
                else
                    stage_candidate_ready = false;
                end
            else
                if pdf_metric > (best_pdf_metric * 1.3) && qs > graceperiod
                    patience_counter = patience_counter + 1;
                    if patience_counter > patience0
                        if sgd_base_gain > 1e-7
                            sgd_base_gain = sgd_base_gain * 0.5;
                            fprintf('--- Gain slashed to %.1e ---\n', sgd_base_gain);
                        end
                        patience_counter = 0;
                        best_pdf_metric = pdf_metric;
                    end
                end
            end

            % Plateau detection (only consider plateau if bins were updated recently)
            max_corr_now = max(abs(sgd_correction));
            if ~isempty(bins_updated)
                max_corr_queue = [max_corr_queue, max_corr_now];
                if numel(max_corr_queue) > plateau_window_batches
                    max_corr_queue(1) = [];
                end
                if numel(max_corr_queue) == plateau_window_batches
                    delta_rel = (max(max_corr_queue) - min(max_corr_queue)) / (mean(max_corr_queue) + eps);
                    if delta_rel < plateau_tol
                        plateau_counter = plateau_counter + 1;
                    else
                        plateau_counter = 0;
                    end
                end
            else
                % If no bins updated this batch, do NOT advance plateau counter
                % but maintain the queue (we want plateau only from active updates)
                % shift queue with same value to avoid false-positive plateau from inactivity
                if numel(max_corr_queue) >= plateau_window_batches
                    % leave queue as is
                else
                    max_corr_queue = [max_corr_queue, max_corr_now];
                end
            end
            is_capped = max_corr_now >= (0.99 * sgd_cap); % Check if we are hugging the ceiling
            plateau_flag = (plateau_counter >= plateau_patience) && ~is_capped;
            
            if is_capped
                fprintf('WARNING: SGD hit the cap (%.2e). Plateau detection suspended.\n', sgd_cap);
                % Optional: Auto-expand the cap?
                % sgd_cap = sgd_cap * 1.5; 
            end
            % Stage transition: need both pdf_metric condition and plateau_flag (and actual updates)
            if exist('stage_candidate_ready','var') && stage_candidate_ready
                if plateau_flag && ~isempty(bins_updated)
                    % advance stage
                    stage_index = stage_index + 1;
                    if stage_index <= max_stages
                        fprintf('\n=== ENTERING REFINEMENT STAGE %d ===\n', stage_index);
                        sgd_correction = best_correction;
                        F_corr_interp.Values = sgd_correction;
                        p = best_p_snapshot; pgp = best_pgp_snapshot;

                        stage_batch_multiplier = batch_mult_sequence(stage_index);
                        sgd_batch_size = stage_batch_multiplier * base_batch_size;
                        stage_scale = 1 / sqrt(stage_batch_multiplier);

                        % reset SNR thresholds at each stage start (Option A)
                        snr_th_search = snr_th_search_default;
                        snr_th_refine = snr_th_refine_default;
                        snr_target = snr_target_default;

                        if gain_ramp
                            ramp_factors = linspace(1, stage_scale, gain_ramp_steps+1);
                            ramp_factors = ramp_factors(2:end);
                            sgd_ramp_queue = ramp_factors;
                        else
                            sgd_base_gain = sgd_base_gain * stage_scale;
                            sgd_ramp_queue = [];
                        end
                        sgd_gain = sgd_base_gain;

                        batch_sum_drift(:) = 0; batch_sum_drift_sq(:) = 0; batch_counts(:) = 0;
                        ndens.counts(:) = 0; pdf.pre.counts(:) = 0;
                        steps_in_batch = 0;
                        stage_warmup = true;
                        stage_steps_filled = 0;
                        current_phase = PHASE_REFINE;

                        max_corr_queue = [];
                        plateau_counter = 0;
                        no_update_counter = 0;
                        starvation_relax_count = 0;
                        continue;
                    else
                        % final big refine toward 1e5
                        fprintf('\n=== ENTERING FINAL REFINEMENT to 1e5 target ===\n');
                        sgd_correction = best_correction;
                        F_corr_interp.Values = sgd_correction;
                        p = best_p_snapshot; pgp = best_pgp_snapshot;

                        target_mult = ceil((w_count_1e5 / base_batch_size));
                        target_mult = min(target_mult, 1000);
                        sgd_batch_size = target_mult * base_batch_size;
                        sgd_base_gain = sgd_base_gain * (1 / sqrt(target_mult));
                        sgd_gain = sgd_base_gain;

                        batch_sum_drift(:) = 0; batch_sum_drift_sq(:) = 0; batch_counts(:) = 0;
                        ndens.counts(:) = 0; pdf.pre.counts(:) = 0;
                        steps_in_batch = 0;
                        stage_warmup = true;
                        stage_steps_filled = 0;
                        current_phase = PHASE_REFINE;

                        max_corr_queue = [];
                        plateau_counter = 0;
                        no_update_counter = 0;
                        starvation_relax_count = 0;
                        continue;
                    end
                else
                    % not plateaued or no active updates: keep collecting
                end
            end

            % Starvation handling (Mode C: adaptive relax only in refine / large batches)
            if isempty(bins_updated)
                % if current batch size is large enough (phase-aware), consider relaxing
                large_batch_condition = (sgd_batch_size >= 3 * base_batch_size) || (current_phase == PHASE_REFINE);
                if large_batch_condition
                    if no_update_counter >= starvation_relax_trigger_batches
                        % relax multiplicatively
                        snr_th_refine = snr_th_refine * 0.8;
                        snr_target = snr_target * 0.8;
                        starvation_relax_count = starvation_relax_count + 1;
                        fprintf('SNR starvation: relaxing thresholds (count=%d). New snr_th_refine=%.3f, snr_target=%.3f\n', starvation_relax_count, snr_th_refine, snr_target);
                        no_update_counter = 0;
                    end
                    if no_update_counter >= starvation_fallback_batches
                        % fallback tiny unbiased shrink update to nudge system
                        fallback_scale = 0.01; % very small
                        fprintf('SNR starvation persistent: applying fallback tiny update (scale=%.3g)\n', fallback_scale);
                        % apply small shrink update toward zero of the mean
                        for ii = 1:sgd_bins
                            if batch_counts(ii) > 0
                                d_fb = - fallback_scale * (batch_sum_drift(ii) / max(1, batch_counts(ii)));
                                d_fb = sign(d_fb) * min(abs(d_fb), max_delta_rel);
                                d_fb = sign(d_fb) * min(abs(d_fb), abs_cap);
                                sgd_correction(ii) = sgd_correction(ii) + d_fb;
                            end
                        end
                        sgd_correction = max(min(sgd_correction, sgd_cap), -sgd_cap);
                        sgd_correction = smoothdata(sgd_correction, 'movmean', sgd_smooth_win);
                        F_corr_interp.Values = sgd_correction;
                        no_update_counter = 0;
                    end
                end
            end

            % consume ramp queue if present
            if exist('sgd_ramp_queue','var') && ~isempty(sgd_ramp_queue)
                if ~stage_warmup
                    next_scale = sgd_ramp_queue(1);
                    sgd_base_gain = sgd_base_gain * next_scale;
                    sgd_ramp_queue(1) = [];
                    fprintf('Ramping gain: new base_gain = %.2e\n', sgd_base_gain);
                end
            end
            sgd_gain = sgd_base_gain;

            % record diagnostics
            history.steps(end+1) = qs;
            history.pdf_dev(end+1) = raw_pdf_metric;
            history.pdf_smooth(end+1) = pdf_metric;
            history.max_corr(end+1) = max_corr_now;
            history.gain(end+1) = sgd_gain;
            history.batch_size(end+1) = sgd_batch_size;
            history.stage(end+1) = stage_index;
            fraction_updated = numel(bins_updated) / max(1, sgd_bins);
            history.fraction_updated(end+1) = fraction_updated;
            history.avg_snr(end+1) = mean(snr(has_data));
            history.median_snr(end+1) = median(snr(has_data));
            history.bins_with_data(end+1) = sum(has_data);
            history.bins_updated(end+1) = numel(bins_updated);
            history.no_update_counter(end+1) = no_update_counter;
            history.starvation_relax_count(end+1) = starvation_relax_count;

            % plotting (black bg, axes white; curves bright and match corresponding axis colors)
            if graphing
                set(0, 'CurrentFigure', f_fig);
                curr_ndens_norm = 100 * (((ndens.counts / max(1,w_count))./ndens.vols) ./ ndens.ndens0);

                subplot(ax_dens);
                plot(ndens.centers, curr_ndens_norm, 'LineWidth', 2, 'Color', 'w');
                xline(S.br - S.rp, '--w'); xline(S.br, '-r');
                ylim([80 120]); xlim([0 S.br]);
                title('Density Deviation [%]','Color','w');
                set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 3);

                subplot(ax_pdf);
                plot(PDF.centers{3}, curr_g, 'LineWidth', 2, 'Color', [1 1 0]); % yellow
                yline(1, '--w');
                xlim([0 2.1*S.br]); ylim([0.5 1.5]);
                title('Pair Distribution Function (g(r))','Color','w');
                set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 3);

                subplot(ax_conv);
                yyaxis left
                ax = gca;
                ax.YColor = [1 1 0]; % left axis color = yellow
                plot(history.steps, history.pdf_smooth, '-','Color',[1 1 0],'LineWidth',2); hold on;
                plot(history.steps, history.pdf_dev, '-','Color',[1 1 0 0.4],'LineWidth',0.8); hold off;
                ylabel('PDF RMS (smoothed)','Color',[1 1 0]);
                yyaxis right
                ax = gca;
                ax.YColor = [0 1 1]; % right axis color = cyan
                plot(history.steps, history.max_corr, '-','Color',[0 1 1],'LineWidth',1.5);
                ylabel('Max Correction','Color',[0 1 1]);
                title('Convergence Metrics','Color','w');
                set(gca,'Color','k','XColor','w','LineWidth',3);

                subplot(ax_ctrl);
                yyaxis left
                ax = gca;
                ax.YColor = [1 0 1]; % magenta
                plot(history.steps, history.gain, '-','Color',[1 0 1],'LineWidth',2);
                ylabel('Gain','Color',[1 0 1]);
                yyaxis right
                ax = gca;
                ax.YColor = [0 1 0]; % bright green
                plot(history.steps, history.batch_size, '-','Color',[0 1 0],'LineWidth',2);
                ylabel('Batch Size','Color',[0 1 0]);
                title('Control State','Color','w');
                set(gca,'Color','k','XColor','w','LineWidth',3);

                subplot(ax_diag);
                plot(history.steps, history.fraction_updated, '-','Color','w','LineWidth',1.5); hold on;
                plot(history.steps, history.bins_with_data ./ sgd_bins, '--','Color','w','LineWidth',1.5); hold off;
                ylim([0 1]);
                title('Fractionity: updated bins / bins with data','Color','w');
                set(gca,'Color','k','XColor','w','LineWidth',3);

                subplot(ax_snr);
                plot(history.steps, history.avg_snr, '-','Color',[0 1 1],'LineWidth',1.5); hold on;
                plot(history.steps, history.median_snr, '--','Color',[0 1 0],'LineWidth',1); hold off;
                title('SNR (avg, median)','Color','w');
                set(gca,'Color','k','XColor','w','LineWidth',3);

                drawnow;
                msg = sprintf('Step %d | RMS: %.4f (Target: %.4f) | Gain: %.1e | Batch: %d | Updated: %.2f', qs, pdf_metric, convergence_target, sgd_gain, sgd_batch_size, fraction_updated);
                fprintf([reverseStr, msg]);
                reverseStr = repmat('\b', 1, length(msg));
            end

            % convergence checks & reset batch accumulators
            if pdf_metric <= target_1e5 && current_phase == PHASE_REFINE
                disp('--- CONVERGED: reached 1e5 noise-floor target ---');
                break;
            end
            if pdf_metric <= convergence_target && qs > 20000 && current_phase == PHASE_REFINE
                disp('--- CONVERGED (local noise floor) ---');
                break;
            end

            batch_sum_drift(:) = 0;
            batch_sum_drift_sq(:) = 0;
            batch_counts(:) = 0;
            ndens.counts(:) = 0;
            pdf.pre.counts(:) = 0;
            steps_in_batch = 0;

            if stage_warmup
                stage_warmup = false;
                stage_steps_filled = 0;
                pdf_metric = 0;
            end
        end
    end

    % --- apply correction & integrate -----------------------------------
    dr_corr_mag = F_corr_interp(prho);
    core_mask = prho < (S.br - potdepth);
    dr_corr_mag(core_mask) = 0;
    total_disp = base_disp + (dr_corr_mag .* pvers);

    p2 = p + total_disp;

    % ghost update
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

    p2rho = vecnorm(p2,2,2);
    pgp2rho = vecnorm(pgp2,2,2);
    idx_swap = p2rho > pgp2rho;
    p_next = p2; pgp_next = pgp2;
    p_next(idx_swap, :) = pgp2(idx_swap, :);
    pgp_next(idx_swap, :) = p2(idx_swap, :);
    p = p_next; pgp = pgp_next;

    if qs > nsteps
        disp('Reached nsteps limit; exiting.');
        break;
    end
end

% save correction & starting configuration
ASYMCORR.correction = [sgd_centers, sgd_correction];
if enable_io
    save(filenamecorrection, 'ASYMCORR', 'sgd_edges');
    save(filestartingconfiguration, 'p', 'pgp', 'S');
end

disp('SGD Setup Complete (v5 final).');
end
