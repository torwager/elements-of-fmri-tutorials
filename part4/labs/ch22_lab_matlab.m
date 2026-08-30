%% Chapter 22 Lab — Multiple Comparisons
% In this lab you will *experience* the multiple comparisons problem by
% simulating thousands of statistical tests with known ground truth. We
% compare uncorrected thresholds, Bonferroni (FWER) correction,
% Benjamini-Hochberg FDR correction, and a permutation-based max-statistic
% threshold.
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox.
% (With CANlab Core tools on your path you can also use FDR.m; see below.)
% Adapted in part from fdr_sims_playground.m in the FMRI_simulations
% repository (Tor Wager, github.com/canlab/FMRI_simulations).

%% Part 1 — The multiple comparisons problem: 10,000 null tests
% A "brain" of 10,000 voxels in which nothing is truly active: 30 subjects
% of pure noise, one-sample t-test at every voxel. How many reach p < .05?

rng(42);
n = 30;          % subjects
k = 10000;       % tests (voxels)

dat = randn(n, k);            % pure noise: the null is true everywhere
[~, p] = ttest(dat);          % one-sample t-test at each voxel

nsig = sum(p < .05);
fprintf('%d of %d null tests are "significant" at p < .05 (expected ~%d)\n', ...
    nsig, k, .05 * k);

figure('Color', 'w');
subplot(1, 2, 1)
histogram(p, 40); hold on
xline(.05, 'r--', 'p = .05');
xlabel('p-value'); ylabel('count'); title('Null p-values are uniform');

subplot(1, 2, 2)
imagesc(reshape(p < .05, 100, 100)); axis image
colormap(flipud(gray))
title(sprintf('%d false positives at p < .05', nsig));
set(gca, 'XTick', [], 'YTick', []);

% Every "significant" voxel here is a false positive, scattered
% salt-and-pepper across the image. A real 2-mm brain mask has ~240,000
% voxels, predicting ~12,000 false positives per uncorrected map.

%% Part 2 — Bonferroni and Benjamini-Hochberg under the null
% Bonferroni controls the family-wise error rate (any false positive) by
% testing each voxel at alpha/m. Benjamini-Hochberg (BH) controls the
% false discovery rate: rank p-values p(1) <= ... <= p(m) and find the
% largest rank r with p(r) <= (r/m)*q, rejecting all tests with p <= p(r).
% Under the global null, both should find (almost) nothing.

alpha = 0.05;

n_unc  = sum(p < alpha);
n_bonf = sum(p < alpha / k);
p_fdr  = bh_threshold(p, alpha);       % local function, bottom of script
n_fdr  = sum(p <= p_fdr);
% With CANlab Core on your path, equivalently: pt = FDR(p, .05);

fprintf('Pure noise, %d tests:\n', k);
fprintf('  Uncorrected p < .05 : %5d significant (all false positives!)\n', n_unc);
fprintf('  Bonferroni          : %5d significant\n', n_bonf);
fprintf('  FDR (BH)            : %5d significant\n', n_fdr);

%% Part 3 — With signal present: the sensitivity/specificity tradeoff
% Plant a true effect of d = 0.5 in 10% of the 10,000 tests, n = 50
% subjects, and score each method against ground truth:
%   TPR (sensitivity), FPR, and the observed false discovery rate.

rng(7);
n = 50;  k = 10000;  d = 0.5;
numtrue = k / 10;                       % 1,000 truly active tests
istrue = (1:k) <= numtrue;
mu = [d * ones(1, numtrue), zeros(1, k - numtrue)];

dat = mu + randn(n, k);
[~, p] = ttest(dat);

p_fdr = bh_threshold(p, alpha);
thresholds = [alpha, alpha / k, p_fdr];
names = {'Uncorrected', 'Bonferroni', 'FDR (BH)'};

fprintf('%-14s%7s%8s%9s%10s\n', 'method', 'n sig', 'TPR', 'FPR', 'obs FDR');
for i = 1:3
    sig = p <= thresholds(i);
    tpr = sum(sig & istrue) / numtrue;
    fpr = sum(sig & ~istrue) / (k - numtrue);
    fdr_obs = sum(sig & ~istrue) / max(sum(sig), 1);
    fprintf('%-14s%7d%8.2f%9.4f%10.3f\n', names{i}, sum(sig), tpr, fpr, fdr_obs);
end

% Uncorrected finds the most true effects but a large fraction of its
% discoveries are false. Bonferroni is nearly error-free but misses most
% of the signal. BH-FDR keeps the false fraction of discoveries near 5%
% while retaining far more sensitivity.

%% Visualizing the BH step-up threshold
% Plot sorted p-values against rank, draw the line q*r/m, and reject
% everything up to the last crossing. The threshold adapts to the amount
% of signal in the map.

psort = sort(p);
ranks = 1:k;
bh_line = alpha * ranks / k;
r_star = find(psort <= bh_line, 1, 'last');

figure('Color', 'w');
nshow = 1200;                          % zoom to the smallest p-values
plot(ranks(1:nshow), psort(1:nshow), '.', 'Color', [.4 .4 .4]); hold on
plot(ranks(1:nshow), bh_line(1:nshow), 'r-', 'LineWidth', 2);
yline(alpha / k, ':b', 'Bonferroni \alpha/m');
xline(r_star, '--g', sprintf('last crossing, r = %d', r_star));
xlabel('rank r'); ylabel('p-value');
title('Benjamini-Hochberg step-up (zoom: smallest 1,200 p-values)');
legend({'sorted p-values', 'BH line q\cdotr/m'}, 'Location', 'northwest');

fprintf('BH threshold: reject the %d tests with p <= %.5f\n', r_star, psort(r_star));
fprintf('Bonferroni threshold: p <= %.2e\n', alpha / k);

%% How stable is the control? Repeating the experiment
% FDR control is a statement about expectations across experiments; a
% single study's observed FDR can be higher or lower. Replicate 20 times.

niter = 20;
tpr = zeros(niter, 3); fdr_obs = zeros(niter, 3);

for i = 1:niter
    dat_i = mu + randn(n, k);
    [~, p_i] = ttest(dat_i);
    thr = [alpha, alpha / k, bh_threshold(p_i, alpha)];
    for j = 1:3
        sig = p_i <= thr(j);
        tpr(i, j) = sum(sig & istrue) / numtrue;
        fdr_obs(i, j) = sum(sig & ~istrue) / max(sum(sig), 1);
    end
end

figure('Color', 'w');
subplot(1, 2, 1)
boxplot(tpr, 'Labels', names); title('Sensitivity (TPR)');
subplot(1, 2, 2)
boxplot(fdr_obs, 'Labels', names); hold on
yline(alpha, 'r--', 'nominal q = .05');
title('Observed FDR');

% BH keeps the observed FDR scattered around (slightly below) 5% with far
% higher sensitivity than Bonferroni.

%% Part 4 — Permutation testing and the max-statistic (FWER) threshold
% Permutation tests build the null distribution from the data itself. For
% a one-sample group test, under the null each subject's effect image is
% symmetric around zero, so we sign-flip whole subject images — preserving
% spatial correlation automatically. Tracking the maximum |t| across all
% voxels in each permutation gives an exact FWER threshold: on smooth
% data it adapts to the effective number of independent tests, where
% Bonferroni cannot. (This is the idea behind FSL randomise and PALM.)

rng(11);
n = 25;  k = 4000;
smooth_fwhm = 7;                        % spatial smoothing kernel (voxels)

% smooth, unit-variance spatial noise; true signal d = 0.8 in a 200-voxel region
istrue = false(1, k); istrue(1801:2000) = true;
mu = 0.8 * istrue;

noise = smoothdata(randn(n, k), 2, 'gaussian', smooth_fwhm);
noise = noise ./ std(noise(:));
dat = mu + noise;

tstat = @(X) mean(X) ./ (std(X) ./ sqrt(size(X, 1)));
t_obs = tstat(dat);

nperm = 500;
maxt = zeros(nperm, 1);
for i = 1:nperm
    signs = sign(rand(n, 1) - .5);
    maxt(i) = max(abs(tstat(signs .* dat)));
end

t_perm = quantile(maxt, 0.95);                 % FWER .05 threshold
t_bonf = tinv(1 - 0.025 / k, n - 1);           % Bonferroni (two-tailed)

fprintf('Permutation max-|t| threshold (FWER .05): t > %.2f\n', t_perm);
fprintf('Bonferroni threshold                    : t > %.2f\n', t_bonf);
fprintf('True voxels detected: permutation %d, Bonferroni %d (of %d)\n', ...
    sum(abs(t_obs) > t_perm & istrue), sum(abs(t_obs) > t_bonf & istrue), sum(istrue));

figure('Color', 'w');
subplot(1, 2, 1)
histogram(maxt, 30); hold on
xline(t_perm, 'r-', sprintf('perm = %.2f', t_perm));
xline(t_bonf, 'b--', sprintf('Bonf = %.2f', t_bonf));
xlabel('max |t| per permutation'); ylabel('count');
title('Permutation distribution of the max statistic');

subplot(1, 2, 2)
plot(t_obs, 'Color', [.3 .3 .3]); hold on
yline(t_perm, 'r-', 'permutation'); yline(t_bonf, 'b--', 'Bonferroni');
xregion(1801, 2000, 'FaceColor', 'y');   % true signal region (R2023a+)
xlabel('voxel'); ylabel('t value');
title('Observed t map (smooth data)');

% Because the noise is smooth, the permutation threshold is lower than
% Bonferroni's — extra detections with exact FWER control.

%% Wrapping up
% * Thousands of uncorrected tests guarantee false positives: ~5% of a
%   null map lights up at p < .05.
% * Bonferroni/FWER control makes any false positive unlikely, at a steep
%   cost in power — and it ignores spatial smoothness.
% * BH-FDR controls the expected fraction of false discoveries, adapts to
%   the signal, and balances Type I and Type II errors far better for
%   brain mapping.
% * Permutation max-statistic thresholds control FWER exactly with almost
%   no assumptions.
% See Chapter 29 for how these choices interact with power and sample size.

%% Local functions

function pt = bh_threshold(p, q)
% Benjamini-Hochberg step-up threshold: largest p(r) with p(r) <= (r/m)*q.
% Returns 0 if nothing survives. (CANlab Core's FDR.m is equivalent.)
m = numel(p);
psort = sort(p(:))';
r = find(psort <= (1:m) / m * q, 1, 'last');
if isempty(r), pt = 0; else, pt = psort(r); end
end
