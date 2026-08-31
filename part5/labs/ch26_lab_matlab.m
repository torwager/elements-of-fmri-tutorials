%% Chapter 26 Lab: Experiments, Observation, and Causality (MATLAB)
% This lab accompanies Chapter 26, "Experiments, Observation, and
% Causality". Causal structure is one of the few things in statistics you
% can *see* directly, because in simulation you know the ground truth.
% You will:
%
%   1. Recreate the hormone replacement therapy (HRT) story: the same
%      simulated population analyzed observationally vs. as a randomized trial
%   2. Simulate a confounder and show that covariate adjustment removes bias
%   3. Simulate a collider and show that adjustment (or selection) CREATES bias
%   4. Preview mediation analysis on a task -> brain -> behavior chain
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox
% (fitlm). Section 5 optionally uses the CANlab Mediation Toolbox:
%   https://github.com/canlab/MediationToolbox
% Simulation design adapted from CANlab teaching code (github.com/canlab;
% FMRI_simulations covariates_in_regression.m) and CANlab Mediation
% Toolbox tutorials.
%
% Runtime: under a minute. All data are simulated.

%% 1. Observation vs. experiment: the HRT story in numbers
% Observational studies found HRT associated with LOWER heart disease risk;
% a randomized trial showed it CAUSES a higher risk. The culprit:
% healthy-user bias. Healthier people were more likely to elect HRT, and
% being healthy independently protects against heart disease.
%
% We simulate one population where treatment is truly HARMFUL (+0.3 on a
% disease-risk score), then analyze it two ways.

rng(26);
n = 5000;
b_true = 0.3;                        % TRUE causal effect of treatment: harmful

H = randn(n, 1);                     % latent "general health" (unmeasured!)

% Observational study: healthier people opt in to treatment
p_take = 1 ./ (1 + exp(-1.5 * H));   % P(taking HRT) rises with health
T_obs  = double(rand(n, 1) < p_take);
Y_obs  = b_true * T_obs - 1.0 * H + randn(n, 1);   % health lowers disease risk

% Randomized trial: coin-flip assignment, same causal model
T_rct = double(rand(n, 1) < 0.5);
Y_rct = b_true * T_rct - 1.0 * H + randn(n, 1);

obs_effect = mean(Y_obs(T_obs == 1)) - mean(Y_obs(T_obs == 0));
rct_effect = mean(Y_rct(T_rct == 1)) - mean(Y_rct(T_rct == 0));

fprintf('True causal effect of treatment:                 %+.2f (harmful)\n', b_true);
fprintf('Observational estimate (self-selected):          %+.2f  <-- looks protective!\n', obs_effect);
fprintf('Randomized trial estimate:                       %+.2f  <-- correct\n', rct_effect);

%%
% The observational study gets the SIGN wrong: the treated group is
% healthier to begin with. Randomization makes treatment independent of
% health -- and of every confounder we did not think to measure. Note that
% we never used H in the analysis: that is the magic of randomization.

%% 2. A confounder: adjusting for it removes bias
% The confounder C drives both the exposure X and the outcome Y:
%
%   X = a*C + u,      Y = b*X + c*C + e
%
% with a = 0.7, b = +0.4 (the true effect), c = -1.2. Regressing Y on X
% alone gives, in expectation, b + a*c/(a^2 + 1) ~= -0.16: biased right
% past zero to the wrong sign. Adjusting for C recovers b. We repeat the
% experiment 500 times to see the whole sampling distribution.

n = 500; n_reps = 500;
a = 0.7; b = 0.4; c = -1.2;          % C->X, X->Y (true), C->Y paths

est = zeros(n_reps, 2);              % columns: naive, adjusted
for i = 1:n_reps
    C = randn(n, 1);
    X = a * C + randn(n, 1);
    Y = b * X + c * C + randn(n, 1);

    m1 = fitlm(X, Y);                % naive: omit the confounder
    m2 = fitlm([X C], Y);            % adjusted for the confounder
    est(i, 1) = m1.Coefficients.Estimate(2);
    est(i, 2) = m2.Coefficients.Estimate(2);
end

fprintf('True effect b = %.2f\n', b);
fprintf('Naive    (omit C):    mean estimate = %+.3f (theory: %+.3f)\n', ...
    mean(est(:, 1)), b + a * c / (a^2 + 1));
fprintf('Adjusted (include C): mean estimate = %+.3f\n', mean(est(:, 2)));

figure('Color', 'w');
histogram(est(:, 1), 30); hold on;
histogram(est(:, 2), 30);
xline(b, 'k--', 'LineWidth', 2);
xline(0, 'Color', [.5 .5 .5]);
xlabel('Estimated effect of X on Y'); ylabel('Count (500 replications)');
title('Confounder: adjustment removes the bias');
legend({'Naive (omit C)', 'Adjusted for C', 'True effect'});

%%
% Every naive replication is biased -- the whole histogram sits on the
% wrong side of zero -- while adjusted estimates cluster around the truth.
% This is the case where "controlling for" a covariate is exactly right:
% the covariate is a COMMON CAUSE of exposure and outcome. In fMRI, this
% is the logic behind motion regressors, physiological covariates, and
% matching groups on age and sex.

%% 3. A collider: adjusting for it CREATES bias
% Now flip the arrows. X and Y are truly unrelated, but both feed a
% downstream composite S = X + Y + noise (think: an "inclusion score"
% summarizing performance and data quality). S is a COLLIDER: arrows
% collide into it. Conditioning on a collider opens a spurious path
% between its causes.

est_col = zeros(n_reps, 2);          % columns: naive, collider-adjusted
for i = 1:n_reps
    X = randn(n, 1);
    Y = randn(n, 1);                 % true effect of X on Y = 0
    S = X + Y + randn(n, 1);         % collider: caused by BOTH X and Y

    m1 = fitlm(X, Y);                % correct model
    m2 = fitlm([X S], Y);            % "controls for" S -> biased!
    est_col(i, 1) = m1.Coefficients.Estimate(2);
    est_col(i, 2) = m2.Coefficients.Estimate(2);
end

fprintf('True effect of X on Y = 0\n');
fprintf('Naive (no adjustment):   mean estimate = %+.3f\n', mean(est_col(:, 1)));
fprintf('Adjusted for collider S: mean estimate = %+.3f  <-- spurious!\n', mean(est_col(:, 2)));

figure('Color', 'w');
histogram(est_col(:, 1), 30); hold on;
histogram(est_col(:, 2), 30);
xline(0, 'k--', 'LineWidth', 2);
xlabel('Estimated effect of X on Y'); ylabel('Count (500 replications)');
title('Collider: adjustment CREATES the bias');
legend({'Naive (correct)', 'Adjusted for collider S', 'True effect (0)'});

%% 3b. Selection is conditioning too
% You do not need to enter the collider in the model: analyzing only
% observations with high S (keeping only "clean, high-performance" trials
% or participants) conditions on it just the same.

X = randn(2000, 1);
Y = randn(2000, 1);                  % truly unrelated
S = X + Y + randn(2000, 1);
keep = S > median(S);                % post-hoc selection on the collider

r_all = corr(X, Y);
r_sel = corr(X(keep), Y(keep));

figure('Color', 'w');
subplot(1, 2, 1);
scatter(X, Y, 6, 'filled', 'MarkerFaceAlpha', 0.3);
title(sprintf('All observations: r = %+.2f', r_all));
xlabel('X'); ylabel('Y'); axis square
subplot(1, 2, 2);
scatter(X(keep), Y(keep), 6, 'filled', 'MarkerFaceAlpha', 0.3, ...
    'MarkerFaceColor', [0.7 0.2 0.2]);
title(sprintf('Selected on S > median: r = %+.2f', r_sel));
xlabel('X'); axis square

%%
% Among selected observations, X and Y are negatively correlated even
% though they are causally unrelated: to make it past the threshold, an
% observation low on X must be high on Y. This is why post-hoc exclusion
% of trials or subjects based on performance, missing behavioral data, or
% head motion needs scrutiny.
%
% RULE OF THUMB: adjust for common CAUSES (confounders); never adjust for
% common EFFECTS (colliders). The regression cannot tell the difference --
% only your causal model can.

%% 4. Mediation preview: task -> brain -> behavior
% In task fMRI the task X is randomized, but brain activity M and behavior
% Y are only observed. The brain sits in the middle of the causal chain --
% a MEDIATOR:
%
%   M = a*X + e_m,     Y = b*M + c'*X + e_y
%
%   a  : effect of task on brain (experimentally secured by randomization)
%   b  : brain-behavior relationship, controlling for the task
%   c' : direct effect of task on behavior, bypassing this brain measure
%   c  = c' + a*b : total effect;  a*b is the INDIRECT (mediated) effect
%
% We simulate a painful-heat experiment: randomized stimulus intensity X,
% a brain response M, and reported pain Y.

n = 200;
a_true = 0.6; b_true = 0.5; cprime_true = 0.2;

X = [zeros(n/2, 1); ones(n/2, 1)];
X = X(randperm(n));                              % randomized low/high intensity
M = a_true * X + randn(n, 1);                    % brain response (observed)
Y = b_true * M + cprime_true * X + randn(n, 1);  % reported pain

m_a  = fitlm(X, M);                  % a-path:      M ~ X
m_c  = fitlm(X, Y);                  % total (c):   Y ~ X
m_y  = fitlm([X M], Y);              % c' and b:    Y ~ X + M

a_hat      = m_a.Coefficients.Estimate(2);
c_hat      = m_c.Coefficients.Estimate(2);
cprime_hat = m_y.Coefficients.Estimate(2);
b_hat      = m_y.Coefficients.Estimate(3);

fprintf('a  (task -> brain):       %+.3f  (true %.2f)\n', a_hat, a_true);
fprintf('b  (brain -> behavior|X): %+.3f  (true %.2f)\n', b_hat, b_true);
fprintf('c  (total effect):        %+.3f  (true %.2f)\n', c_hat, cprime_true + a_true * b_true);
fprintf('c'' (direct effect):       %+.3f  (true %.2f)\n', cprime_hat, cprime_true);
fprintf('a*b (indirect effect):    %+.3f = c - c'' = %+.3f\n', ...
    a_hat * b_hat, c_hat - cprime_hat);

%% 4b. Bootstrap test of the indirect effect
% Because a*b is a product of estimates, its sampling distribution is
% skewed; the standard test uses the bootstrap: resample participants with
% replacement and read the confidence interval off the a*b distribution.

n_boot = 2000;
ab_boot = zeros(n_boot, 1);
for i = 1:n_boot
    idx = randi(n, n, 1);            % resample participants
    Xb = X(idx); Mb = M(idx); Yb = Y(idx);
    mb_a = fitlm(Xb, Mb);
    mb_y = fitlm([Xb Mb], Yb);
    ab_boot(i) = mb_a.Coefficients.Estimate(2) * mb_y.Coefficients.Estimate(3);
end

ci = prctile(ab_boot, [2.5 97.5]);
p_boot = 2 * min(mean(ab_boot <= 0), mean(ab_boot >= 0));

figure('Color', 'w');
histogram(ab_boot, 40); hold on;
xline(0, 'Color', [.5 .5 .5]);
xline(a_true * b_true, 'k--', 'LineWidth', 2);
xline(ci(1), 'r:', 'LineWidth', 2); xline(ci(2), 'r:', 'LineWidth', 2);
xlabel('Bootstrap indirect effect (a*b)'); ylabel('Count');
title(sprintf('Indirect effect: 95%% CI [%.2f, %.2f], p = %.4f', ci(1), ci(2), p_boot));

%% 5. The same analysis with the CANlab Mediation Toolbox (optional)
% The CANlab Mediation Toolbox wraps this whole analysis -- paths, bootstrap
% tests, and path diagrams -- in a single call, and extends it voxel-wise to
% whole-brain "mediation effect parametric mapping".
% Adapted from CANlab mediation tutorials (canlab.github.io, mediation_1_basics).

if exist('mediation', 'file')
    [paths, stats] = mediation(X, Y, M, 'boot', 'plots', 'verbose', 'bootsamples', 5000);
else
    disp('CANlab Mediation Toolbox not on path; skipping.');
    disp('Get it at https://github.com/canlab/MediationToolbox');
end

%% Wrap-up
% * Randomization makes the IV independent of ALL confounders, even
%   unmeasured ones: the observational and RCT analyses of the same
%   population gave opposite signs.
% * Confounders (common causes): adjust for them, and bias goes away.
% * Colliders (common effects): adjust for them -- or select observations
%   on them -- and bias appears out of nowhere.
% * Mediators carry the effect: adjusting converts a total effect into a
%   direct effect. That changes the question; it does not fix bias.
% * The brain in task fMRI is a mediator: task -> brain is experimental;
%   brain -> behavior needs converging causal evidence.
%
% Next, Chapter 27 applies these principles to designing task fMRI
% experiments.
