%% Chapter 29 Lab: Statistical Power and Sample Size (MATLAB)
% This lab accompanies Chapter 29, "Statistical Power and Sample Size".
% You will build a power calculator from the noncentral t distribution,
% quantify the cost of multiple comparisons correction, find the minimum
% detectable effect size for a planned sample, simulate the "winner's
% curse" (effect size inflation at a statistical threshold), and compare
% univariate vs. multivariate power in the spirit of the BWAS debate
% (Marek et al., 2022).
%
% Requirements: Statistics and Machine Learning Toolbox (tinv, nctcdf,
% normcdf, ttest). CanlabCore is not required for this lab.
% Adapted from the CANlab power utilities (power_calc.m) and the
% FMRI_simulations scripts behind the book's power figures
% (github.com/canlab).
%
% Runtime: under a minute. All data are simulated.

%% 1. Power for a one-sample t-test, analytically
% With true effect size d = mu/sigma and sample size N, the one-sample
% t statistic follows a noncentral t distribution with N - 1 degrees of
% freedom and noncentrality delta = d*sqrt(N). Power is the probability
% of exceeding the critical value (two-tailed for a planned test):
%
%   power = 1 - nctcdf(tcrit, N - 1, d*sqrt(N))
%
% (Local functions implementing this are defined at the end of the file.)
%
% One convention to fix up front: a PLANNED test is usually two-tailed,
% whereas thresholds applied to fMRI statistic maps ("p < .001") are
% conventionally DIRECTIONAL (one tail), because we threshold positive
% activations. The helper functions take a 'tails' argument so we can
% follow that convention -- as CANlab's power_calc.m does, where the
% alpha input is one-tailed.

dvals = [0.2 0.5 0.8];   % Cohen's benchmarks: small, medium, large
alpha = 0.05;

for i = 1:length(dvals)
    ncrit = n_for_power(@power_1samp, dvals(i), alpha);
    fprintf('d = %.1f: N = %4d for 80%% power (one-sample, p < .05 two-tailed)\n', ...
        dvals(i), ncrit);
end

% A "medium" effect (d = 0.5) needs about N = 34, matching standard power
% software such as G*Power.

%% 2. Power curves as a function of sample size
% The classic planning plot: power vs. N for several effect sizes, with
% the 80% power criterion marked. Drop-lines show the N needed for 80%.

dvals = [0.2 0.3 0.4 0.5 0.8];
N = 3:250;

figure('Color', 'w'); hold on
colors = copper(length(dvals) + 2);

for i = 1:length(dvals)
    pow = power_1samp(dvals(i), N, alpha);
    plot(N, pow, 'LineWidth', 2, 'Color', colors(i + 1, :), ...
        'DisplayName', sprintf('d = %.1f', dvals(i)));

    ncrit = n_for_power(@power_1samp, dvals(i), alpha);
    if ~isempty(ncrit) && ncrit <= max(N)
        plot([ncrit ncrit], [0 .8], ':', 'Color', [.5 .5 .5], ...
            'HandleVisibility', 'off');
    end
end

plot([min(N) max(N)], [.8 .8], 'k--', 'HandleVisibility', 'off')
xlabel('Sample size (N)'); ylabel('Power')
title('One-sample t-test power, p < .05 two-tailed')
legend('Location', 'southeast'); ylim([0 1.02])

%% 3. Sanity check: verify analytic power by simulation
% Power is just the long-run fraction of significant experiments, so we
% can verify the formula by brute force: simulate many experiments with
% true d = 0.5 and count how often p < .05.

rng(29)
d_true = 0.5;  n_sims = 2000;

fprintf('%4s %10s %10s\n', 'N', 'analytic', 'simulated')
for n = [10 20 34 50 80]
    dat = d_true + randn(n, n_sims);           % subjects x experiments
    [~, p] = ttest(dat);                       % one-sample t-test each column
    fprintf('%4d %10.3f %10.3f\n', n, power_1samp(d_true, n, alpha), mean(p < .05));
end

%% 4. The cost of multiple comparisons
% Mass-univariate mapping cannot use alpha = .05 per voxel. Typical
% effective per-test thresholds:
%   p < .05        one pre-registered ROI test
%   p < .001       uncorrected mapping; often approximates FDR q < .05
%   p < .05/1000   Bonferroni over ~1,000 parcels/tests
%   p < 4.26e-6    whole-brain FWER (permutation-based, empirical average)
% How much do required sample sizes grow?

alphas     = [0.05  0.001  0.05/1000  4.26e-6];
tails      = [2     1      1          1];
alphanames = {'p < .05, two-tailed (ROI)', 'p < .001 (~FDR q < .05)', ...
              'Bonferroni, 1000 tests', 'FWER whole brain'};
d_grid = [0.2 0.3 0.5 0.8];

fprintf('\nN for 80%% power, ONE-SAMPLE test\n');
fprintf('%-26s', 'threshold');
fprintf('  d=%-5.1f', d_grid); fprintf('\n');
for i = 1:length(alphas)
    fprintf('%-26s', alphanames{i});
    for j = 1:length(d_grid)
        fprintf('  %-7d', n_for_power(@power_1samp, d_grid(j), alphas(i), tails(i)));
    end
    fprintf('\n');
end

fprintf('\nTOTAL N for 80%% power, TWO-GROUP comparison (patients vs. controls)\n');
fprintf('%-26s', 'threshold');
fprintf('  d=%-5.1f', d_grid); fprintf('\n');
for i = 1:length(alphas)
    fprintf('%-26s', alphanames{i});
    for j = 1:length(d_grid)
        fprintf('  %-7d', 2 * n_for_power(@power_2samp, d_grid(j), alphas(i), tails(i)));
    end
    fprintf('\n');
end

% Correction multiplies the required sample by ~3-4x: a d = 0.5 effect
% needs 34 participants for one ROI test but ~122 under FWER correction.
% A two-group comparison needs ~4x the total sample: 128 at p < .05 and
% ~460 with FWER correction. (Book planning values: 34, 121, 130, 466.)

%% 5. Correlations: reproducing the book's Figure 29.2 reference values
% For brain-behavior correlations we use the Fisher z approximation:
% atanh(r) is ~normal with SE = 1/sqrt(N-3). Two-group comparisons need
% roughly 4x the total sample of one-sample tests (power_2samp below).

r_grid = [0.1 0.2 0.3 0.4 0.5];

fprintf('\nN for 80%% power to detect a correlation:\n');
fprintf('%5s %8s %8s %8s\n', 'r', 'p<.05', 'p<.001', 'FWER');
for i = 1:length(r_grid)
    n05   = n_for_power(@power_corr, r_grid(i), 0.05,    2);
    n001  = n_for_power(@power_corr, r_grid(i), 0.001,   1);
    nfwer = n_for_power(@power_corr, r_grid(i), 4.26e-6, 1);
    fprintf('%5.1f %8d %8d %8d\n', r_grid(i), n05, n001, nfwer);
end

% These reproduce the book's Figure 29.2 values to within a participant
% or two: r = 0.5 needs ~30 (single ROI), ~55 (p < .001), ~96 (FWER);
% r = 0.1 needs ~780, ~1,540, and ~2,790. Small effects plus whole-brain
% search is a brutal combination.

%% 6. Minimum detectable effect size for a planned sample
% Invert the power analysis: given the N you can afford, what is the
% smallest effect you have an 80% chance of detecting? This "minimum
% detectable effect size" is an honest summary of a study's sensitivity.
% (Mirrors power_min_detectable_effect_size.m from FMRI_simulations.)

n_grid = [30 50 100 200 500 1000];
d_search = 0.05:0.005:3;
r_search = 0.02:0.002:0.99;

fprintf('\nMinimum detectable effect with 80%% power:\n');
fprintf('%6s %10s %10s %10s %10s\n', 'N', 'd, p<.05', 'd, FWER', 'r, p<.05', 'r, FWER');
for i = 1:length(n_grid)
    d05   = min(d_search(power_1samp(d_search, n_grid(i), 0.05,    2) >= .8));
    dfwer = min(d_search(power_1samp(d_search, n_grid(i), 4.26e-6, 1) >= .8));
    r05   = min(r_search(power_corr(r_search, n_grid(i), 0.05,     2) >= .8));
    rfwer = min(r_search(power_corr(r_search, n_grid(i), 4.26e-6,  1) >= .8));
    fprintf('%6d %10.2f %10.2f %10.2f %10.2f\n', n_grid(i), d05, dfwer, r05, rfwer);
end

% This reconstructs Figure 29.2C/D from scratch: N = 30 -> d ~ 1.17,
% N = 100 -> d ~ 0.56, N = 1000 -> d ~ 0.17. With FWER correction, an
% N = 30 study is powered only for very large effects -- far larger than
% the d ~ 0.5 typical of task effects in individual voxels (based on HCP
% reference activations).

%% 7. The winner's curse: effect size inflation at a threshold
% Simulate 20,000 voxels that ALL share the same true effect (d = 0.5,
% as in the book's Figure 29.1), run a group t-test, threshold, and
% compare post hoc effect sizes in significant voxels with the truth.
% Concept from effect_size_inflation_example_sim.m (github.com/canlab).

rng(29)
n_sub = 30;  n_vox = 20000;  d_true = 0.5;

dat = d_true + randn(n_sub, n_vox);      % subjects x voxels
[~, ~, ~, st] = ttest(dat);
p_dir = 1 - tcdf(st.tstat, n_sub - 1);   % directional p, as in mapping
d_hat = st.tstat ./ sqrt(n_sub);         % observed effect size per voxel

sig = p_dir < .001;
fprintf('\nTrue effect size:                 d = %.2f\n', d_true);
fprintf('Mean estimate, ALL voxels:        d = %.2f  (unbiased)\n', mean(d_hat));
fprintf('Mean estimate, significant only:  d = %.2f  (%.0f%% inflated)\n', ...
    mean(d_hat(sig)), 100 * (mean(d_hat(sig)) / d_true - 1));
fprintf('Voxels significant at p < .001:   %.1f%% (analytic power: %.1f%%)\n', ...
    100 * mean(sig), 100 * power_1samp(d_true, n_sub, .001, 1));

figure('Color', 'w'); hold on
edges = -0.2:0.03:1.4;
histogram(d_hat, edges, 'FaceColor', [.7 .7 .7], 'DisplayName', 'all voxels');
histogram(d_hat(sig), edges, 'FaceColor', [.8 .2 .3], 'DisplayName', 'significant (p < .001)');
plot([d_true d_true], ylim, 'k:', 'LineWidth', 2, 'DisplayName', 'true d = 0.5');
xlabel('Estimated effect size'); ylabel('Number of voxels')
title('Winner''s curse: selection inflates post hoc effect sizes')
legend

%% 8. Inflation depends on threshold stringency and sample size
% Stricter thresholds select luckier noise and inflate more; larger
% samples shrink the bias. Paradoxically, correcting for multiple
% comparisons makes false positives rarer but post hoc effect sizes MORE
% inflated -- so power analyses should never be based on the significant
% voxels of the same map.

thresholds   = [0.05 0.005 0.001 4.26e-6];
sample_sizes = [15 30 60 120];

fprintf('\nMean d-hat in significant voxels (true d = %.1f everywhere):\n', d_true);
fprintf('%5s', 'N');
fprintf('  p<%-8.2g', thresholds); fprintf('\n');
for i = 1:length(sample_sizes)
    n = sample_sizes(i);
    dat = d_true + randn(n, n_vox);
    [~, ~, ~, st] = ttest(dat);
    p_dir = 1 - tcdf(st.tstat, n - 1);
    d_hat = st.tstat ./ sqrt(n);
    fprintf('%5d', n);
    for j = 1:length(thresholds)
        s = p_dir < thresholds(j);
        if sum(s) >= 10
            fprintf('  %-10.2f', mean(d_hat(s)));
        else
            fprintf('  %-10s', '--');
        end
    end
    fprintf('\n');
end

%% 9. The BWAS debate: univariate vs. multivariate power
% Marek et al. (2022): the largest univariate brain-behavior correlations
% are ~r = 0.1, but multivariate models reach ~r = 0.4 in the same data.
% Compare power curves and required samples for both regimes.
% (Based on univ_vs_multivar_power_based_on_marek2022.mlx.)

r_effects = [0.095 0.11 0.18 0.39];
labels = {'best univariate edge (r = 0.095)', 'median multivariate (r = 0.11)', ...
          '75th pct multivariate (r = 0.18)', 'best multivariate (r = 0.39)'};

N = 4:2500;
figure('Color', 'w'); hold on
colors = lines(length(r_effects));
for i = 1:length(r_effects)
    plot(N, power_corr(r_effects(i), N, 0.05), 'LineWidth', 2, ...
        'Color', colors(i, :), 'DisplayName', labels{i});
end
plot([min(N) max(N)], [.8 .8], 'k--', 'HandleVisibility', 'off')
xlabel('Sample size (N)'); ylabel('Power')
title('Power for BWAS-scale effects, p < .05 (cf. Marek et al. 2022)')
legend('Location', 'southeast')

n_uni   = n_for_power(@power_corr, 0.095, 0.05, 2);
n_multi = n_for_power(@power_corr, 0.39,  0.05, 2);
fprintf('\np < .05:  univariate r = 0.095 -> N = %d; multivariate r = 0.39 -> N = %d (%.0f-fold)\n', ...
    n_uni, n_multi, n_uni / n_multi);

% Same comparison at a p < .001 mapping threshold, using the round values
% quoted in the literature (univariate r = 0.1 vs. multivariate r = 0.4)
n_uni_001   = n_for_power(@power_corr, 0.10, 0.001, 1);
n_multi_001 = n_for_power(@power_corr, 0.40, 0.001, 1);
fprintf('p < .001: univariate r = 0.10  -> N = %d; multivariate r = 0.40 -> N = %d (%.0f-fold)\n', ...
    n_uni_001, n_multi_001, n_uni_001 / n_multi_001);

% The best univariate effects need thousands of participants; multivariate
% effects of r ~ 0.4 are detectable with N in the tens to low hundreds --
% and a multivariate model yields ONE test, so no multiple comparisons
% correction is needed. Task-evoked patterns can be stronger still (d > 3
% for some validated signatures), detectable in very small samples.

%% Local functions
% Power formulas used throughout. power_1samp and power_2samp use the
% noncentral t distribution (cf. power_calc.m); power_corr uses the
% Fisher z approximation for correlations.

function pow = power_1samp(d, n, alpha, tails)
% Power of a one-sample t-test with true effect size d.
% tails = 2 for a two-sided planned test (default); tails = 1 for the
% directional thresholds conventionally applied to fMRI statistic maps.
if nargin < 4, tails = 2; end
tcrit = tinv(1 - alpha/tails, n - 1);
pow = 1 - nctcdf(tcrit, n - 1, d .* sqrt(n));
end

function pow = power_2samp(d, n_per_group, alpha, tails)
% Power of a two-sample t-test, balanced groups. Total N is
% 2*n_per_group -- roughly 4x the one-sample requirement.
if nargin < 4, tails = 2; end
df = 2 .* n_per_group - 2;
tcrit = tinv(1 - alpha/tails, df);
pow = 1 - nctcdf(tcrit, df, d .* sqrt(n_per_group ./ 2));
end

function pow = power_corr(r, n, alpha, tails)
% Power to detect a correlation r, Fisher z approximation.
if nargin < 4, tails = 2; end
zcrit = norminv(1 - alpha/tails);
pow = normcdf(sqrt(n - 3) .* atanh(r) - zcrit);
end

function ncrit = n_for_power(power_fn, effect, alpha, tails)
% Smallest N (or n per group) achieving 80% power; empty if > 6000.
if nargin < 4, tails = 2; end
nvals = 3:6000;
pow = power_fn(effect, nvals, alpha, tails);
ncrit = nvals(find(pow >= .80, 1));
end
