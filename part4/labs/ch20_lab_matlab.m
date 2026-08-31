%% Chapter 20 Lab: Contrasts and Inference with the GLM (MATLAB)
% This lab accompanies Chapter 20, "Contrasts and Inference with the GLM".
% You will fit a three-condition GLM to a simulated voxel, construct and
% test t-contrasts (single conditions, pairwise differences, averages),
% see why contrast scaling and covariate centering change interpretation
% but not inference, run F-contrasts for joint tests, and test a
% conjunction ("A AND B") hypothesis.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples), including
% "Construction and inference with a simple design".
%
% Runtime: under a minute. All data are simulated.

%% 1. A three-condition design
% Imagine three event types: A = famous faces, B = non-famous faces,
% C = houses. We present 42 one-second events (14 per condition) in a
% 400-second run, randomly interleaved. onsets2fmridesign convolves each
% condition's indicator function with the canonical HRF and adds the
% intercept as the LAST column.

TR = 2;                    % repetition time (s)
n_scans = 200;             % number of volumes
run_len = n_scans * TR;    % 400 s

rng(2024);                 % reproducible design and noise
all_onsets = (10:9:384)';  % 42 events, 9 s apart
shuffled = all_onsets(randperm(length(all_onsets)));

ons = {};
ons{1} = sort(shuffled(1:14));     % Condition A onsets (s)
ons{2} = sort(shuffled(15:28));    % Condition B onsets (s)
ons{3} = sort(shuffled(29:42));    % Condition C onsets (s)

X = onsets2fmridesign(ons, TR, run_len);
names = {'A' 'B' 'C' 'Intercept'};

disp(size(X))              % 200 x 4: [A, B, C, intercept]

plotDesign(ons, [], TR);   % onsets and convolved regressors

create_figure('design matrix image');
imagesc(X); colormap gray;
set(gca, 'XTick', 1:4, 'XTickLabel', names);
xlabel('Regressor'); ylabel('Time (TRs)'); title('Design matrix X');

%% 2. Simulate a voxel and fit the GLM
% Ground truth: the voxel responds most to famous faces, somewhat to
% non-famous faces, weakly to houses.
%   beta_A = 1.0, beta_B = 0.6, beta_C = 0.3, intercept = 100

% onsets2fmridesign scales the HRF to a peak of 1, so each task regressor
% peaks near 1. With an intercept of 100, a beta of 1.0 is about a 1%
% signal change, and the noise SD below is 0.15% of baseline.
beta_true = [1.0 0.6 0.3 100]';
sigma_noise = 0.15;

y = X * beta_true + sigma_noise .* randn(n_scans, 1);

% OLS fit, saving the pieces every contrast test will reuse
beta_hat = (X' * X) \ (X' * y);

r       = y - X * beta_hat;         % residuals
[n, p]  = size(X);
dfe     = n - p;                    % error degrees of freedom
sigma2  = (r' * r) / dfe;           % error variance estimate
XtX_inv = inv(X' * X);

disp(table(beta_true, beta_hat, 'RowNames', names))

%% 3. t-contrasts: one question, one number
% A contrast is a linear combination of betas, c'*beta. Its t-statistic:
%   t = c'*beta_hat / sqrt(sigma2 * c' * (X'X)^{-1} * c),  df = dfe
%
% Rules of thumb:
% - Comparisons between conditions: weights sum to ZERO (e.g., [1 -1 0]).
% - Tests against the (unmodeled) baseline need not sum to zero
%   (e.g., [1 0 0] or [1/3 1/3 1/3]) - each beta is 0 under the null.

C = [ 1    0    0    0;      % A vs. baseline
      0    0    1    0;      % C vs. baseline
      1   -1    0    0;      % A - B (famous - nonfamous)
      0.5  0.5 -1    0;      % faces - houses: (A+B)/2 - C
      1/3  1/3  1/3  0 ]';   % task average vs. baseline
con_names = {'A vs baseline'; 'C vs baseline'; 'A - B'; ...
    'faces - houses'; 'task avg vs baseline'};

con_val = C' * beta_hat;                          % contrast estimates
se_con  = sqrt(sigma2 * diag(C' * XtX_inv * C));  % contrast SEs
t_con   = con_val ./ se_con;
p_con   = 2 * (1 - tcdf(abs(t_con), dfe));
true_val = C' * beta_true;

disp(table(true_val, con_val, se_con, t_con, p_con, 'RowNames', con_names))

% Note that A - B has a smaller t than A vs. baseline even though the SEs
% are similar: its true effect (0.4) is simply smaller. How precisely a
% design estimates a PARTICULAR contrast depends on the timing and
% correlation structure of the regressors (design efficiency, Chapter 27).
% The contrast values computed at every voxel become the contrast images
% (COPEs) used in group analysis.

%% 4. Scaling contrast weights does not change inference
% Multiplying c by a constant scales the estimate AND its SE by the same
% factor, so t and P are untouched.

for scale = [1 4]
    c = scale * [1 -1 0 0]';
    val = c' * beta_hat;
    t   = val / sqrt(sigma2 * c' * XtX_inv * c);
    fprintf('c = %d*[1 -1 0]: value = %+6.3f, t = %5.3f\n', scale, val, t);
end

% Exceptions worth caring about:
% 1) Units: for percent-signal-change reporting, scale weights so positive
%    weights sum to +1 and negative weights to -1 (e.g., [.5 .5 -1], not
%    [1 1 -2]) so the contrast is a difference between condition MEANS.
% 2) Group analysis: contrast scale must be identical across participants;
%    with missing runs, re-normalize each subject's weights to sum to +/-1.

%% 5. Centering a covariate changes only the intercept's meaning
% Add a continuous covariate w (mean ~3), raw vs. mean-centered. The task
% betas and covariate slope are identical; only the intercept changes:
% "expected signal at w = 0" becomes "expected signal at average w".

w   = 3 + randn(n_scans, 1);
y_w = y + 0.5 * w;                 % data with a true covariate effect

Xw_raw = [X(:, 1:3) w            ones(n_scans, 1)];
Xw_cen = [X(:, 1:3) w - mean(w)  ones(n_scans, 1)];

b_raw = (Xw_raw' * Xw_raw) \ (Xw_raw' * y_w);
b_cen = (Xw_cen' * Xw_cen) \ (Xw_cen' * y_w);

disp(table(b_raw, b_cen, ...
    'RowNames', {'A' 'B' 'C' 'w (slope)' 'Intercept'}))

%% 6. F-contrasts: testing several effects jointly
% Stack q contrasts as ROWS of a matrix L. The F-statistic
%   F = (L*b)' * inv(L * (X'X)^{-1} * L') * (L*b) / (q * sigma2)
% tests whether the set is JOINTLY nonzero (unsigned, non-directional),
% with (q, dfe) degrees of freedom. It equals the classic full-vs-reduced
% model comparison: how much EXTRA variance do the tested effects explain?

% 1) Omnibus: does the task explain any variance (A, B, C jointly)?
L = [1 0 0 0; 0 1 0 0; 0 0 1 0];
[F_any, p_any] = f_contrast(L, beta_hat, XtX_inv, sigma2, dfe);
fprintf('Any task effect:          F(%d,%d) = %8.2f, p = %g\n', ...
    size(L, 1), dfe, F_any, p_any);

% ... identical to comparing full and reduced (intercept-only) models
X0   = ones(n_scans, 1);
r0   = y - X0 * ((X0' * X0) \ (X0' * y));
F_rss = ((r0' * r0 - r' * r) / 3) / ((r' * r) / dfe);
fprintf('Full vs. reduced (RSS):   F(3,%d) = %8.2f  <- identical\n', dfe, F_rss);

% 2) Do the conditions differ from each other at all? (ANOVA omnibus)
L = [1 -1 0 0; 0 1 -1 0];
[F_diff, p_diff] = f_contrast(L, beta_hat, XtX_inv, sigma2, dfe);
fprintf('Any condition difference: F(%d,%d) = %8.2f, p = %g\n', ...
    size(L, 1), dfe, F_diff, p_diff);

% 3) A one-row F-contrast is the square of the t-contrast (sign is lost)
c = [1 -1 0 0]';
[F_ab, ~] = f_contrast(c', beta_hat, XtX_inv, sigma2, dfe);
t_ab = (c' * beta_hat) / sqrt(sigma2 * c' * XtX_inv * c);
fprintf('F(A-B) = %5.3f = t(A-B)^2 = %5.3f\n', F_ab, t_ab ^ 2);

% Typical F-contrast uses: a joint test of nuisance covariates (e.g., six
% motion parameters), an omnibus test of a factorial model, or condition
% differences across multiple HRF basis functions - whose betas are
% incommensurate ("apples and oranges") and cannot simply be averaged.

%% 7. Conjunction: "A AND B", not "A OR B"
% An F-test asks whether ANY tested combination is nonzero (a logical OR).
% Claims like "this region responds to faces whether or not they are
% famous" are AND claims: both A - C > 0 and B - C > 0. The valid
% minimum-statistic rule: declare a conjunction only if EVERY component
% contrast exceeds the threshold on its own, i.e., min(t) > t_crit.

t_crit = tinv(0.95, dfe);           % one-sided, alpha = .05
c_AC = [1 0 -1 0]'; c_BC = [0 1 -1 0]';

% Voxel 1 responds to both face types; voxel 2 only to famous faces
beta_vox = [0.9 0.8 0.1 100; 0.9 0.1 0.1 100]';
vox_names = {'Voxel 1 (A and B)', 'Voxel 2 (A only)'};

for v = 1:2
    y_v = X * beta_vox(:, v) + sigma_noise .* randn(n_scans, 1);
    b   = (X' * X) \ (X' * y_v);
    r_v = y_v - X * b;
    s2  = (r_v' * r_v) / dfe;

    t_AC = (c_AC' * b) / sqrt(s2 * c_AC' * XtX_inv * c_AC);
    t_BC = (c_BC' * b) / sqrt(s2 * c_BC' * XtX_inv * c_BC);
    is_conj = min(t_AC, t_BC) > t_crit;

    fprintf('%s: t(A-C) = %5.2f, t(B-C) = %5.2f, conjunction = %d\n', ...
        vox_names{v}, t_AC, t_BC, is_conj);
end

% Voxel 2 fails the conjunction even though its F-test over the same two
% contrasts would be highly significant - F and conjunction answer
% different logical questions (OR vs. AND).

%% Explore on your own
% 1. Test [1 -0.5 -0.5 0] ("famous vs. the average of everything else").
%    What is its true value in this simulation?
% 2. Build a 2 x 2 factorial with four conditions (betas 1.0, 0.8, 0.5,
%    0.7) and test main effects and interaction with [1 1 -1 -1],
%    [1 -1 1 -1], and [1 -1 -1 1] (plus a 0 for the intercept).
% 3. In the conjunction demo, raise voxel 2's beta_B toward 0.8. Around
%    what effect size does the conjunction succeed reliably?

%% Local functions

function [F, p] = f_contrast(L, beta_hat, XtX_inv, sigma2, dfe)
% F-test of the joint null L*beta = 0. L has one contrast per ROW.
q  = size(L, 1);
Lb = L * beta_hat;
F  = Lb' * ((L * XtX_inv * L') \ Lb) / (q * sigma2);
p  = 1 - fcdf(F, q, dfe);
end
