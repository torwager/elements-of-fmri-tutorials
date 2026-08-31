%% Chapter 21 Lab: Group Analysis (MATLAB)
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part4/ch21-group-analysis
%
% This lab accompanies Chapter 21, "Group Analysis". You will simulate
% hierarchical (multi-subject) data with known ground truth, compare fixed
% effects (FFX) and random effects (RFX) inference, fit mixed effects
% models, and see how robust regression protects group results from
% outlier subjects.
%
% Requirements: Statistics and Machine Learning Toolbox (ttest, fitlme,
% robustfit). The mixed-model section also demonstrates CANlab's
% glmfit_multilevel from CanlabCore (https://github.com/canlab/CanlabCore);
% that section is optional if CanlabCore is not on your path.
% Code adapted from CANlab tutorials (github.com/canlab,
% CANlab_help_examples/canlab_mixed_effects_matlab_demo1) and the
% Computational Foundations mixed-effects tutorials.
%
% Runtime: about a minute. All data are simulated.

%% 1. Simulate hierarchical (multi-subject) data
% Group fMRI data have two levels. For one voxel's contrast:
%
%   Level 1 (within-subject):   y_ij = beta_i + eps_ij,  eps ~ N(0, sigma_w^2)
%   Level 2 (between-subject):  beta_i = beta_G + eta_i, eta ~ N(0, sigma_b^2)
%
% Trials scatter around each subject's true effect beta_i, and the subject
% effects scatter around the population effect beta_G. The two variance
% components -- within-subject noise and true between-subject variability
% -- drive everything in this lab.

rng(2021);          % seed the random number generator for reproducibility

n_subj   = 20;      % subjects
n_trials = 40;      % trials (1st-level observations) per subject
beta_G   = 0.5;     % true population effect
sigma_b  = 0.5;     % SD of true subject effects (eta_i)
sigma_w  = 1.0;     % SD of trial-level noise (eps_ij)

subj_fx = beta_G + sigma_b .* randn(1, n_subj);              % true beta_i
Y = repmat(subj_fx, n_trials, 1) + sigma_w .* randn(n_trials, n_subj);

subj_means = mean(Y)';       % first-level summary statistic per subject

% Plot trials, subject means, and true effects, sorted by true effect
[~, order] = sort(subj_fx);

figure; hold on;
for k = 1:n_subj
    i = order(k);
    plot(k * ones(n_trials, 1), Y(:, i), '.', 'Color', [.75 .75 .75], 'MarkerSize', 4);
end
plot(1:n_subj, subj_means(order), 'o', 'Color', [0 .35 .75], 'MarkerFaceColor', [0 .35 .75]);
plot(1:n_subj, subj_fx(order), '_r', 'LineWidth', 2, 'MarkerSize', 10);
yline(beta_G, '--k'); yline(0, 'k');
xlabel('Subject (sorted by true effect)'); ylabel('Response (a.u.)');
title('Trials (gray), subject means (blue), true subject effects (red)');

%% 2. Decompose the variance of subject means
% The subject means vary because of BOTH real individual differences and
% leftover trial noise:  var(means) ~ sigma_b^2 + sigma_w^2 / n_trials.
% With 40 trials, between-subject variance dominates.

var_obs      = var(subj_means);
var_expected = sigma_b^2 + sigma_w^2 / n_trials;

fprintf('Variance of subject means, observed:  %.3f\n', var_obs);
fprintf('Expected sigma_b^2 + sigma_w^2/n:     %.3f + %.3f = %.3f\n', ...
    sigma_b^2, sigma_w^2 / n_trials, var_expected);

%% 3. Fixed effects vs. random effects inference on one dataset
% FFX ("super subject"): concatenate all trials from all subjects into one
% t-test. Its error term contains only within-subject noise (df ~ 799).
% RFX (summary statistics): one-sample t-test on the 20 subject means;
% variation ACROSS SUBJECTS is the error term (df = 19).

[~, p_ffx, ~, stat_ffx] = ttest(Y(:));          % FFX: pool 800 trials
[~, p_rfx, ~, stat_rfx] = ttest(subj_means);    % RFX: subject means

fprintf('\n%-28s  est = %.3f  t(%4d) = %6.2f  p = %.2g\n', ...
    'FFX (pool all trials):', mean(Y(:)), stat_ffx.df, stat_ffx.tstat, p_ffx);
fprintf('%-28s  est = %.3f  t(%4d) = %6.2f  p = %.2g\n', ...
    'RFX (subject means):', mean(subj_means), stat_rfx.df, stat_rfx.tstat, p_rfx);

% Same group estimate -- but FFX standard errors are far too small, because
% between-subject variance is missing from its error term.

%% 4. The false-positive experiment: FFX inflates FPR under the null
% Set the TRUE group effect to zero but keep real between-subject
% variability: some subjects activate, others deactivate, population mean
% is nil. A valid test should reject ~5% of the time. Repeat with
% sigma_b = 0 (identical subjects), the only world where FFX is valid.

n_iter = 1000;                 % simulated experiments; more -> stabler FPR estimates
sigma_b_levels = [0 0.5];      % between-subject SDs: identical subjects vs. real differences
fpr = zeros(numel(sigma_b_levels), 2);            % rows: sigma_b, cols: FFX RFX

for s = 1:numel(sigma_b_levels)
    sb = sigma_b_levels(s);
    p1 = zeros(n_iter, 1); p2 = zeros(n_iter, 1);
    for it = 1:n_iter
        sfx = sb .* randn(1, n_subj);                                 % beta_G = 0
        Yi  = repmat(sfx, n_trials, 1) + sigma_w .* randn(n_trials, n_subj);
        [~, p1(it)] = ttest(Yi(:));                                   % FFX
        [~, p2(it)] = ttest(mean(Yi)');                               % RFX
    end
    fpr(s, :) = [mean(p1 < .05) mean(p2 < .05)];
end

t = array2table(fpr, 'VariableNames', {'FFX', 'RFX'}, ...
    'RowNames', {'sigma_b = 0', 'sigma_b = 0.5'});
disp('False positive rates at alpha = .05 (nominal: 0.05):'); disp(t)

figure;
bar(fpr); hold on; yline(0.05, '--k', 'nominal \alpha = .05');
set(gca, 'XTickLabel', {'\sigma_b = 0 (identical)', '\sigma_b = 0.5 (subjects differ)'});
legend({'FFX (super subject)', 'RFX (summary statistics)'}, 'Location', 'northwest');
ylabel('False positive rate');
title('FFX inflates false positives when between-subject variance exists');

% With realistic individual differences, FFX rejects the true null most of
% the time. Valid population inference requires between-subject variation
% in the error term -- which is exactly what RFX provides.

%% 5. A genuine mixed effects model with fitlme
% Now a within-subject design: each subject completes trials in two
% conditions (effects coded -0.5 / +0.5), with random intercepts AND random
% slopes across subjects. We fit the full mixed model to trial-level data
% and compare it with the summary statistics approach.

rng(33);                               % seed for reproducibility
n_subj = 24; n_per_cond = 30;          % subjects; trials per condition
betaG_0 = 0.2; betaG_1 = 0.6;          % population intercept, condition effect
sd_int = 0.4; sd_slope = 0.5;          % between-subject SDs of intercept, slope
sigma_w = 1.0;                         % SD of trial-level (within-subject) noise

subject = []; cond = []; y = [];
for i = 1:n_subj
    b0 = betaG_0 + sd_int * randn;
    b1 = betaG_1 + sd_slope * randn;
    ci = [-0.5 * ones(n_per_cond, 1); 0.5 * ones(n_per_cond, 1)];
    yi = b0 + b1 .* ci + sigma_w .* randn(2 * n_per_cond, 1);
    subject = [subject; i * ones(2 * n_per_cond, 1)]; %#ok<AGROW>
    cond = [cond; ci]; y = [y; yi];                   %#ok<AGROW>
end
tbl = table(y, cond, categorical(subject), 'VariableNames', {'y', 'cond', 'subject'});

% Wilkinson notation: fixed effect of cond, random intercept + slope by subject
lme = fitlme(tbl, 'y ~ cond + (cond | subject)', 'FitMethod', 'REML');
disp(lme)

% Compare: summary statistics approach (t-test on per-subject differences)
diffs = zeros(n_subj, 1);
for i = 1:n_subj
    wh = subject == i;
    diffs(i) = mean(y(wh & cond > 0)) - mean(y(wh & cond < 0));
end
[~, p_ss, ~, stat_ss] = ttest(diffs);
fprintf('\nSummary statistics: est = %.3f, SE = %.3f, t(%d) = %.2f, p = %.2g\n', ...
    mean(diffs), std(diffs) / sqrt(n_subj), stat_ss.df, stat_ss.tstat, p_ss);

% With a balanced design and homogeneous error, the two agree closely: the
% one-sample t-test on contrasts IS a simplified mixed effects analysis.
% The mixed model additionally returns the variance components (compare the
% random-effects SDs in the fitlme output with sd_int/sd_slope above).

%% 6. CANlab glmfit_multilevel: precision-weighted summary statistics
% CANlab's glmfit_multilevel fits a first-level model per subject, then a
% precision-weighted second-level model -- a fast mixed-effects compromise
% used in the CANlab GLM and mediation toolboxes (igls.m, mediation.m).
% Skip this cell if CanlabCore is not on your path.
% Adapted from CANlab tutorials (github.com/canlab).

if ~isempty(which('glmfit_multilevel'))
    YY = cell(1, n_subj); XX = cell(1, n_subj);
    for i = 1:n_subj
        wh = subject == i;
        YY{i} = y(wh);
        XX{i} = cond(wh);       % intercept is added automatically (first column)
    end

    stats_ml = glmfit_multilevel(YY, XX, [], 'weighted', 'noverbose', ...
        'names', {'Intercept' 'Condition'});

    fprintf('glmfit_multilevel: Condition effect = %.3f, t(%.1f) = %.2f, p = %.2g\n', ...
        stats_ml.beta(1, 2), stats_ml.dfe(2), stats_ml.t(1, 2), stats_ml.p(1, 2));
else
    disp('CanlabCore not found -- skipping glmfit_multilevel demo.');
end

%% 7. Outlier subjects: OLS vs. robust (IRLS) regression
% Second-level scenario: one contrast (COPE) value per subject, regressed
% on a behavioral covariate. There is NO true relationship -- until one
% outlier subject creates one. Robust iteratively reweighted least squares
% (IRLS) down-weights points far from the central mass of subjects.
% CANlab's robfit.m applies this voxelwise to contrast images.

rng(11);                               % seed for reproducibility
n2 = 30;                               % subjects contributing one COPE each
perf = randn(n2, 1);                   % behavioral covariate (mean-centered)
cope = 0.3 .* randn(n2, 1);            % contrast values: true slope = 0

perf_o = perf; cope_o = cope;
perf_o(end) = 4; cope_o(end) = 3;      % one extreme outlier subject

[b_clean, ~, st_clean] = glmfit(perf, cope);            % OLS, clean
[b_ols,   ~, st_ols]   = glmfit(perf_o, cope_o);        % OLS, with outlier
[b_rob, st_rob]        = robustfit(perf_o, cope_o);     % IRLS (bisquare)

fprintf('Slope, clean OLS:          %6.3f  (p = %.3f)\n', b_clean(2), st_clean.p(2));
fprintf('Slope, OLS with outlier:   %6.3f  (p = %.3f)  <- spurious!\n', b_ols(2), st_ols.p(2));
fprintf('Slope, robust IRLS:        %6.3f  (p = %.3f)\n', b_rob(2), st_rob.p(2));
fprintf('IRLS weight for the outlier subject: %.3f\n', st_rob.w(end));

xg = linspace(-2.5, 4.5, 50)';
figure; hold on;
scatter(perf_o, cope_o, 40, 1 - st_rob.w, 'filled', 'MarkerEdgeColor', 'k');
colormap(flipud(gray)); cb = colorbar; ylabel(cb, '1 - IRLS weight (darker = kept)');
plot(xg, b_ols(2) .* xg + b_ols(1), 'r-', 'LineWidth', 1.5);
plot(xg, b_rob(2) .* xg + b_rob(1), 'k-', 'LineWidth', 1.5);
legend({'subjects (shaded by IRLS weight)', 'OLS fit', 'robust IRLS fit'}, ...
    'Location', 'northwest');
xlabel('Performance (covariate)'); ylabel('Contrast value (COPE)');
title('One outlier subject flips the OLS slope; IRLS down-weights it');

%% 8. Wrap-up
% - Multi-subject data are hierarchical; valid population inference puts
%   between-subject variability in the error term.
% - FFX omits that variance component: valid only if subjects are
%   identical, badly anticonservative otherwise.
% - The summary statistics approach implicitly treats subject as a random
%   effect and is fully efficient when first-level precision is homogeneous.
% - Mixed models (fitlme, glmfit_multilevel, igls, FLAME, 3dMEMA) estimate
%   variance components explicitly and precision-weight subjects; they help
%   most with unbalanced, missing, or variable-quality data.
% - Robust IRLS regression (robustfit; robfit.m for images) is cheap
%   insurance against outlier subjects at the second level.
%
% Next: Chapter 22 -- running this group test at 100,000 voxels at once
% raises the multiple comparisons problem.
