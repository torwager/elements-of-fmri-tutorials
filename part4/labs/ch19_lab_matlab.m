%% Chapter 19 Lab: GLM Design Specification (MATLAB)
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part4/ch19-glm-design-specification
%
% This lab accompanies Chapter 19, "GLM Design Specification". You will
% build design matrices with variable-duration events and parametric
% modulators, measure collinearity with variance inflation factors (VIFs),
% see exactly what orthogonalization does to parameter estimates, and add
% nuisance covariates (motion, spikes) and a cosine drift set that
% implements high-pass filtering inside the GLM.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples, incl. variable_dur_and_param_mod_example_code.mlx).
%
% Runtime: under a minute. All data are simulated.

%% 1. Event vs. epoch regressors, and variable durations
% The first modeling choice is how to represent neural activity: brief
% events, sustained epochs, or events whose duration varies trial to trial
% (like reaction time). onsets2fmridesign accepts a cell array with one
% cell per condition; adding a second column of durations (in sec) to each
% cell produces duration-modulated regressors.

TR = 2;                                      % repetition time (s)
runlen = 360;                                % run length in seconds

onsets = {(10:24:340)'};                     % one event type, 14 events
rt = [1.2 3.4 0.8 4.5 2.0 3.1 0.6 5.0 2.6 1.5 4.0 1.0 3.6 2.2]';  % RTs (s)

% Fixed (mean) duration model vs. variable-duration model
ons_fixed = {[onsets{1} repmat(mean(rt), 14, 1)]};   % [onset dur] columns
ons_vardur = {[onsets{1} rt]};

Xfixed = onsets2fmridesign(ons_fixed, TR, runlen, 'hrf');
Xvar   = onsets2fmridesign(ons_vardur, TR, runlen, 'hrf');

create_figure('durations', 2, 1);
subplot(2, 1, 1); plotDesign(ons_fixed, [], TR, 'samefig');
title('Fixed (mean) duration'); set(gca, 'XLim', [0 runlen]);
subplot(2, 1, 2); plotDesign(ons_vardur, [], TR, 'samefig');
title('Variable duration (RT)'); set(gca, 'XLim', [0 runlen]);

% Question: In which regressor do long-RT trials produce larger predicted
% responses? Why does duration modulation diverge from amplitude
% modulation as durations exceed ~2-3 sec?

%% 2. Parametric modulation
% A parametric modulator asks: does the response amplitude scale with a
% trial-by-trial variable (here, pain ratings)? 'parametric_standard'
% builds, for each condition, an average-response regressor plus a
% modulator regressor built from MEAN-CENTERED modulator values (it
% centers, but does not orthogonalize -- see Section 4). Columns:
% [Avg, Modulator, intercept].

ratings = [3 5 2 7 4 6 1 8 5 3 7 2 6 4]';    % one rating per event

X = onsets2fmridesign(onsets, TR, runlen, 'hrf', 'parametric_standard', ...
    {ratings});

create_figure('parametric modulation');
h = plot_matrix_cols(X(:, 1:2), 'horiz');
set(h, 'LineWidth', 2); set(h(2), 'LineStyle', '--');
set(gca, 'YTick', [1 2], 'YTickLabel', {'Avg response' 'Rating modulator'});
xlabel('Time (TRs)'); title('Average + parametric modulator regressors');

%% 3. Collinearity: correlations and variance inflation factors
% To see the collinearity problem clearly, build a NON-orthogonalized
% modulator by convolving a delta function scaled by raw (uncentered)
% ratings, and compare correlations and VIFs for the raw vs. the
% mean-centered version.

% Build a raw modulated regressor by hand: onsets weighted by ratings
Xavg = onsets2fmridesign(onsets, TR, runlen, 'hrf');    % [Avg, intercept]
x1 = Xavg(:, 1);                                        % average regressor

% 'singletrial' gives one column per event; weight columns by ratings
Xtrial = onsets2fmridesign(onsets, TR, runlen, 'hrf', 'singletrial');
xmod_raw      = Xtrial(:, 1:14) * ratings;              % raw modulator
xmod_centered = Xtrial(:, 1:14) * (ratings - mean(ratings));

fprintf('corr(avg, raw modulator)      = %3.2f\n', corr(x1, xmod_raw));
fprintf('corr(avg, centered modulator) = %3.2f\n', corr(x1, xmod_centered));

% VIFs (getvif adds an intercept internally; pass task columns only)
vif_raw      = getvif([x1 xmod_raw]);
vif_centered = getvif([x1 xmod_centered]);
disp(table(vif_raw', vif_centered', 'VariableNames', ...
    {'VIF_raw' 'VIF_centered'}, 'RowNames', {'Avg' 'Modulator'}))

% With evenly spaced onsets like these, mean-centering removes essentially
% all of the collinearity (VIF -> 1); with jittered or clustered onsets some
% usually remains. Rule of thumb: VIF < 2.5 is comfortable; > 4 warrants
% concern; > 8 is bad.

%% 4. The orthogonalization experiment
% Simulate a "voxel" whose true response includes both an average effect
% and a modulator effect, then fit the pair as-is vs. with the modulator
% orthogonalized -- using both the RAW and the mean-centered modulator.
% Watch which beta changes and which does not.

rng(42);                                     % seed for reproducibility
beta_true = [1.0 0.5];                        % avg amplitude, modulator slope
y = 100 + beta_true(1) * x1 + beta_true(2) * xmod_centered ...
    + 2 * randn(size(x1));

int = ones(size(x1));                        % intercept column
b = zeros(4, 2);                              % rows: avg/mod x as-is/orth

for i = 1:2
    if i == 1, xmod = xmod_raw; else, xmod = xmod_centered; end

    % resid(a, b, true) returns b with everything a explains removed
    xmod_orth = resid(x1, xmod, true);        % orthogonalize wrt average

    b_asis = [x1 xmod      int] \ y;          % as-is (correlated) model
    b_orth = [x1 xmod_orth int] \ y;          % orthogonalized model

    b(:, i) = [b_asis(1); b_orth(1); b_asis(2); b_orth(2)];
end

disp(table(b(:, 1), b(:, 2), 'VariableNames', {'Raw_mod' 'Centered_mod'}, ...
    'RowNames', {'Avg_asis' 'Avg_orth' 'Mod_asis' 'Mod_orth'}))
fprintf('True values: average = %3.1f, modulator = %3.1f\n', beta_true);

% Key result: the MODULATOR beta is identical with and without
% orthogonalization -- it already reflects the unique (partial) effect.
% What changes is the AVERAGE regressor's beta, which absorbs all the
% shared variance. With the raw modulator (r = .89 with the average) that
% shift is dramatic; with the centered modulator there is almost nothing
% left to reassign. Orthogonalization reassigns credit; it does not create
% information -- and mean-centering, not orthogonalization, is the real fix.

%% 5. High-pass filtering inside the GLM: cosine drift set
% Slow scanner drift dominates low frequencies. Instead of pre-filtering,
% we can add a discrete cosine basis set spanning frequencies below the
% cutoff (here 1/128 s, the SPM default) as nuisance columns in X.

n = runlen / TR;                              % number of volumes
HPlength = 128;                               % cutoff period (s)
k = fix(2 * runlen / HPlength + 1);           % number of DCT functions
drift = spm_dctmtx(n, k);                     % col 1 is constant

create_figure('drift basis');
plot(drift(:, 2:end), 'LineWidth', 1.5);
xlabel('Time (TRs)'); title('Cosine drift regressors (128-s cutoff)');

% Simulate two kinds of drift -- one unrelated to the task, and one that
% shares the task's low-frequency structure -- and fit each with and
% without the drift set, tracking the task beta, its SE, and its t value.

low = resid(ones(n, 1), drift * pinv(drift) * x1);   % task's low-freq part
driftsigs = {4 * sin(2 * pi * (1:n)' / n * 1.5) + linspace(0, -5, n)', ...
             5 * low ./ std(low)};
driftnames = {'unrelated to task', 'correlated w/ task'};
models = {[x1 ones(n, 1)], [x1 drift]};
modelnames = {'no drift model', 'with drift set'};

for d = 1:2
    y2 = 100 + 1.0 * x1 + driftsigs{d} + 2 * randn(n, 1);

    for m = 1:2
        Xi = models{m};
        bi = Xi \ y2;
        r = y2 - Xi * bi;
        dfe = n - size(Xi, 2);
        s2 = (r' * r) / dfe;
        se = sqrt(s2 * diag(inv(Xi' * Xi)));
        fprintf('drift %-19s %-15s beta = %5.3f (true 1.00), SE = %5.3f, t = %6.2f, dfe = %d\n', ...
            driftnames{d}, modelnames{m}, bi(1), se(1), bi(1) / se(1), dfe);
    end
end

% Unmodeled drift always inflates residual variance -- the SE roughly
% doubles and t is roughly halved. When the drift also shares structure
% with the task, the task beta is biased as well, because the drift has
% nowhere to go but the task regressor.

% Also try the CANlab pre-filtering route, which removes the same variance:
% y2f = hpfilter(y2, TR, HPlength, n);   % filter data (and X!) before OLS

%% 5b. What does the filter cost you? Check the CONTRAST
% Filtering is not free: task variance below the cutoff is removed too. In
% event-related designs the variance is spread across frequencies, and what
% matters is the variance of the CONTRAST you plan to test, not just the
% individual regressors. Here we build a randomized two-condition design
% and measure how much variance an 80-s filter removes from each.

rng(3);                                      % seed for this random design
allons = sort(8 + (runlen - 28) * rand(28, 1));   % 28 jittered events
isA = false(28, 1); isA(randperm(28, 14)) = true; % random A/B labels
Xab = onsets2fmridesign({allons(isA) allons(~isA)}, TR, runlen, 'hrf');

D = spm_dctmtx(n, fix(2 * runlen / 80 + 1));      % 80-s high-pass filter
P = D * pinv(D);                                  % projection onto drift space
pctrem = @(x) 100 * (1 - sum((x - P * x).^2) ./ sum((x - mean(x)).^2));

fprintf('variance removed from A:              %4.1f%%\n', pctrem(Xab(:, 1)));
fprintf('variance removed from B:              %4.1f%%\n', pctrem(Xab(:, 2)));
fprintf('variance removed from A - B contrast: %4.1f%%\n', ...
    pctrem(Xab(:, 1) - Xab(:, 2)));

% Re-run this cell with different rng seeds. The cost varies a lot across
% random orderings, and for roughly one design in five the CONTRAST loses
% more variance than either regressor -- which is why you check the filter
% cost on your contrasts before you scan (see Chapter 27 on design
% efficiency).

%% 6. Motion covariates, spike regressors, and task VIFs
% Add six simulated realignment parameters (random-walk motion) and spike
% regressors for two "bad" volumes, then compare task VIFs with and
% without the nuisance set. Correlations AMONG nuisance regressors are
% harmless; their correlation with TASK regressors is what inflates task
% VIFs and (rightly) reduces significance.

mot = cumsum(randn(n, 6) * 0.02);             % simulated motion params
spikes = intercept_model(n, [40 41]);         % spike regs for volumes 40-41
spikes = spikes(:, 2:end);                    % drop the intercept column

Xtask = [x1 xmod_centered];
Xfull = [Xtask mot spikes drift];             % full first-level design

vif_task  = getvif(Xtask);
vif_full  = getvif(Xfull);
fprintf('Task VIFs alone:        Avg %3.2f, Modulator %3.2f\n', vif_task(1:2));
fprintf('Task VIFs full model:   Avg %3.2f, Modulator %3.2f\n', vif_full(1:2));

create_figure('full design');
imagesc(zscore([Xtask mot spikes drift(:, 2:end)])); colormap gray; colorbar
xlabel('Regressor'); ylabel('Time (TRs)');
title('Full design: task, motion, spikes, cosine drift');

% Questions:
% 1. Re-simulate motion so it correlates with the task (add 0.1*x1 to one
%    motion column). What happens to the task VIFs? To significance?
% 2. Double the number of drift regressors (halve HPlength). How many
%    error degrees of freedom do you lose? When would that matter most?
% 3. Why must nuisance regressors added AFTER pre-filtering be filtered
%    with the same filter before entering the GLM?
