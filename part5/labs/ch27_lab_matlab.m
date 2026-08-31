%% Chapter 27 Lab: Experimental Design and Task fMRI (MATLAB)
% This lab accompanies Chapter 27, "Experimental Design and Task fMRI".
% You will build blocked, dense event-related, and sparse jittered
% event-related designs, score their efficiency for a target contrast,
% diagnose collinearity with variance inflation factors (VIFs), and map
% the tradeoff between contrast detection and HRF shape estimation.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path, plus the
% Statistics and Machine Learning Toolbox.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab, canlab.github.io:
% "Efficiency, colinearity, and variance inflation" and the
% FMRI_simulations design-optimization tools).
%
% Runtime: about a minute. No scanner data needed -- efficiency depends
% only on the design matrix, so we score designs before any data exist.

which create_random_er_design   % if this errors, CanlabCore is not on your path

%% 1. Build three competing designs
% All designs share: two conditions (A, B), 2-s events, 480-s run,
% minimum ISI of 4 s, TR = 1 s, 128-s high-pass filter. Only the
% arrangement of events in time differs.

TR = 1;                 % repetition time (s)
scanLength = 480;       % run length (s)
ISI = 4;                % minimum inter-stimulus interval (s)
eventduration = 2;      % event duration (s)
HPlength = 128;         % high-pass filter length (s)
dononlin = 0;           % nonlinear saturation model off (see Section 6)

% Blocked design: 2 conditions in alternating 16-s blocks
create_figure('block design');
[Xb, eb] = create_block_design(scanLength, TR, 2, 16, HPlength, dononlin);

% Dense fixed-ISI event-related: every 4-s slot holds an event (50% A, 50% B)
create_figure('dense ER design');
[Xd, ed] = create_random_er_design(TR, ISI, eventduration, [.5 .5], ...
    HPlength, dononlin, 'scanLength', scanLength);

% Sparse jittered event-related: 25% A, 25% B, 50% rest slots ("jitter")
create_figure('sparse ER design');
[Xs, es, ons_s] = create_random_er_design(TR, ISI, eventduration, [.25 .25], ...
    HPlength, dononlin, 'scanLength', scanLength);

fprintf('Overall efficiency (A-optimality): block %3.1f, dense %3.1f, sparse %3.1f\n', ...
    eb, ed, es);

%% 2. Efficiency of the contrasts you actually plan to test
% Overall efficiency averages across regressors, but designs should be
% scored on the specific contrasts of interest. The design-related
% variance of contrast c is c * inv(X'X) * c', and efficiency is its
% inverse:  e_c = 1 / (c * inv(X'X) * c').
% We use pinv(X) * pinv(X)' as a numerically stable form of inv(X'X).

c_diff = [1 -1 0];      % A - B: do the conditions differ?
c_base = [1  0 0];      % A vs. implicit baseline: does A activate at all?

eff = @(X, c) 1 ./ (c * pinv(X) * pinv(X)' * c');

designs = {Xb Xd Xs};
names = {'Blocked' 'Dense ER' 'Sparse jittered'};

e_diff = cellfun(@(X) eff(X, c_diff), designs);
e_base = cellfun(@(X) eff(X, c_base), designs);

disp(table(names', e_diff', e_base', ...
    'VariableNames', {'Design' 'Eff_AminusB' 'Eff_AvsBaseline'}))

create_figure('contrast efficiency', 1, 2);
subplot(1, 2, 1); bar(e_diff); set(gca, 'XTickLabel', names);
title('Efficiency of [A - B]'); ylabel('Efficiency (higher = better)');
subplot(1, 2, 2); bar(e_base); set(gca, 'XTickLabel', names);
title('Efficiency of [A vs. baseline]');

% The blocked design should dominate for A - B (compare Figure 27.1B in
% the book). For A vs. baseline, the dense design suffers: with no rest,
% its regressors hover near a plateau and are confounded with the
% intercept. Jitter creates the peaks AND valleys needed to compare
% events with an implicit resting baseline.

%% 3. Collinearity: correlations among regressors
% Efficiency drops when regressors are correlated -- the model cannot
% assign credit among predictors that make similar predictions. Look at
% the correlations among the task regressors in our three designs.

create_figure('correlations', 1, 3);
for i = 1:3
    subplot(1, 3, i);
    imagesc(corr(designs{i}(:, 1:2)), [-1 1]); colorbar;
    title(names{i}); set(gca, 'XTick', 1:2, 'YTick', 1:2);
end

for i = 1:3
    r = corr(designs{i}(:, 1), designs{i}(:, 2));
    fprintf('%-16s corr(A, B) = %5.2f\n', names{i}, r);
end

% Block and dense designs have strongly NEGATIVE A-B correlations: with no
% rest, whenever A is on, B is off. That is fine for A - B but terrible
% for estimating either condition against baseline.

%% 4. Variance inflation factors
% Pairwise correlations miss the worst case: a regressor predictable from
% a COMBINATION of the others (including the intercept). The variance
% inflation factor for regressor j is VIF_j = 1 / (1 - R^2_j), where
% R^2_j comes from regressing predictor j on all remaining predictors.
% VIF = 1 is perfect (orthogonal); ~2+ deserves attention; ~4-8+ is
% serious trouble.

for i = 1:3
    create_figure(['vifs ' names{i}]);
    v = getvif(designs{i}, false, 'plot');   % CanlabCore
    title(sprintf('VIFs: %s', names{i}));
    fprintf('%-16s VIF(A) = %6.2f   VIF(B) = %6.2f\n', names{i}, v(1), v(2));
end

% The no-rest trap: block and dense designs show VIFs in the tens. Each
% regressor is almost perfectly predictable from the other plus the
% intercept, so neither beta is well estimated on its own -- yet the
% A - B contrast is estimated precisely (Section 2). VIF is a
% PER-PARAMETER diagnostic; specific contrasts can escape collinearity.
% The sparse jittered design keeps VIFs near 1.

%% 4b. The fixed-lag trap: the most common collinearity mistake in tasks
% A very natural way to program a task is: cue, then stimulus exactly 2 s
% later, every trial. The two convolved regressors then look nearly
% identical, and the model cannot tell cue-related from stimulus-related
% signal.

cue  = ons_s{1}(:, 1);              % cue onsets from the sparse design
dur  = eventduration * ones(size(cue));
lags = [2 6 10];                    % jittered cue -> stimulus intervals
lag  = lags(randi(3, numel(cue), 1))';

Xfix = onsets2fmridesign({[cue dur] [cue + 2   dur]}, TR, scanLength);
Xjit = onsets2fmridesign({[cue dur] [cue + lag dur]}, TR, scanLength);

vf = getvif(Xfix);  vj = getvif(Xjit);
fprintf('cue-stim, fixed 2-s lag: corr %5.2f  VIF %5.2f  e(A-B) %6.1f\n', ...
    corr(Xfix(:, 1), Xfix(:, 2)), vf(1), eff(Xfix, c_diff));
fprintf('cue-stim, jittered lag : corr %5.2f  VIF %5.2f  e(A-B) %6.1f\n', ...
    corr(Xjit(:, 1), Xjit(:, 2)), vj(1), eff(Xjit, c_diff));

create_figure('cue-stimulus lag', 1, 2);
subplot(1, 2, 1); plot(Xfix(:, 1:2)); title('Fixed 2-s cue-stimulus lag');
xlabel('Time (TRs)'); legend({'Cue' 'Stimulus'});
subplot(1, 2, 2); plot(Xjit(:, 1:2)); title('Jittered cue-stimulus lag');
xlabel('Time (TRs)');

% Jittering the cue -> stimulus interval decorrelates the regressors and
% restores efficiency. Jitter is the cure.

%% 5. Detection vs. estimation efficiency: the fundamental tradeoff
% Detection efficiency assumes a canonical HRF and asks how precisely we
% can estimate the A - B contrast. Estimation efficiency uses a finite
% impulse response (FIR) model -- one free parameter per post-stimulus
% time point -- and asks how precisely we can estimate the HRF's SHAPE.
% Both use the same formula with different design matrices.
%
% We score random designs with varying amounts of rest (jitter), plus
% the blocked design, on both criteria -- a home-made Figure 27.4.

n_lags = 16;            % FIR: estimate response at lags 0..15 s
nT = scanLength / TR;

% Build an FIR design matrix from onset times (in sec) by shifting a
% stick function; intercept is the last column.
fir_design = @(onsets_cell) local_fir(onsets_cell, nT, TR, n_lags);

rest_props = [0 .125 .25 .375 .5 .625];
n_reps = 8;
[det_eff, est_eff, jit] = deal([]);

for i = 1:numel(rest_props)
    p_ev = (1 - rest_props(i)) / 2;
    for j = 1:n_reps
        create_figure('scratch'); % create_random_er_design plots into current fig
        [Xr, ~, ons_r] = create_random_er_design(TR, ISI, eventduration, ...
            [p_ev p_ev], HPlength, dononlin, 'scanLength', scanLength);
        det_eff(end+1) = eff(Xr, c_diff);                     %#ok<*SAGROW>
        Xf = fir_design(ons_r);
        v = diag(pinv(Xf' * Xf));
        est_eff(end+1) = (2 * n_lags) / sum(v(1:end-1));      % FIR cols only
        jit(end+1) = rest_props(i);
    end
end
close all

% Score the blocked design the same way. Block onsets: events every 4 s
% within alternating 16-s blocks.
onsA = []; onsB = [];
for bstart = 0:32:scanLength-32
    onsA = [onsA, bstart + (0:4:12)];
    onsB = [onsB, bstart + 16 + (0:4:12)];
end
Xf_block = fir_design({onsA' onsB'});
v = diag(pinv(Xf_block' * Xf_block));
est_block = (2 * n_lags) / sum(v(1:end-1));
det_block = eff(Xb, c_diff);

create_figure('tradeoff');
scatter(est_eff, det_eff, 40, jit, 'filled'); hold on
plot(est_block, det_block, 'p', 'MarkerSize', 18, ...
    'MarkerFaceColor', [1 .6 0], 'MarkerEdgeColor', 'k');
xlabel('HRF shape estimation efficiency (FIR model)');
ylabel('Contrast detection efficiency (canonical HRF)');
title('Detection vs. estimation: no design wins both');
hcb = colorbar; ylabel(hcb, 'Proportion rest (jitter)');
legend({'Random ER designs' 'Blocked (16 s)'});

% The blocked design sits high on detection but far left on estimation:
% its events are packed into predictable runs, so shifted FIR regressors
% are nearly collinear and the response shape is unrecoverable. Adding
% jitter improves shape estimation at some cost in detection.

%% 6. Going further: genetic algorithms and design optimization
% Random search cannot cover the astronomical space of possible designs.
% The CANlab OptimizeDesign toolbox (in CanlabCore/OptimizeDesign11)
% implements a genetic algorithm, optimizeGA, that recombines pieces of
% the best designs across "generations". You specify a GA structure with
% your conditions, contrasts and their weights, high-pass filter, and
% the balance between contrast detection and HRF estimation fitness --
% see ga_example_script.m and Genetic_Algorithm_readme.rtf in that
% folder. A typical call is:
%
%   M = optimizeGA(GA);        % GA = structure of specifications
%
% GA-optimized designs beat random event-related designs on BOTH
% detection and estimation (Figure 27.4 in the book), approaching the
% theoretical limit. m-sequences (see the M-sequence subfolder) are the
% estimation-efficiency extreme.
%
% Also try: create_random_er_design(..., dononlin = 1) applies a
% saturation model. Rapid designs look far less efficient once BOLD
% nonlinearity (Section 27.2.4) is taken into account -- linear-model
% efficiency flatters designs with events packed closer than ~5 s.

%% Questions to answer
% 1. Re-run Section 1 with 8-s and 32-s blocks. Where does A - B
%    efficiency peak? Why are very long blocks risky given the 128-s
%    high-pass filter?
% 2. In Section 2, why does the dense design do so poorly for
%    [A vs. baseline] despite having the most trials?
% 3. How high do the VIFs get in Section 4? What correlation between two
%    regressors corresponds to a VIF of 2?
% 4. In Section 4b, how much does the fixed 2-s cue-stimulus lag cost you
%    relative to the jittered version? What lag distribution would you
%    use in a real task, and what psychological cost does jitter carry?
% 5. Set dononlin = 1 in Section 1 and re-score the three designs. Which
%    design's efficiency drops most, and why?

%% Local function: FIR design matrix from onsets
function Xf = local_fir(onsets_cell, nT, TR, n_lags)
% One shifted stick-function column per lag per condition; intercept last.
% onsets_cell{i} may be [n x 1] onsets or [n x 2] (onset, duration) in sec.
Xf = ones(nT, 2 * n_lags + 1);
for i = 1:2
    ons_i = onsets_cell{i}(:, 1);
    stick = zeros(nT, 1);
    stick(round(ons_i ./ TR) + 1) = 1;
    for lag = 0:n_lags-1
        col = [zeros(lag, 1); stick(1:nT-lag)];
        Xf(:, (i - 1) * n_lags + lag + 1) = col;
    end
end
end
