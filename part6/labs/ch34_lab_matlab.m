%% Chapter 34 Lab: Structural Equation and Path Models (MATLAB)
% This lab accompanies Chapter 34, "Structural Equation and Path Models".
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part6/ch34-structural-equation-and-path-models
% You will simulate a brain-as-mediator dataset, estimate the mediation
% paths a, b, c', and c with the CANlab Mediation Toolbox, bootstrap the
% indirect effect a*b, see how an unmodeled confounder of the M-Y
% relationship creates spurious mediation, and fit a small three-node
% path model by regression equations.
%
% Requirements: CanlabCore and the MediationToolbox on your MATLAB path,
% plus the Statistics and Machine Learning Toolbox.
%   https://github.com/canlab/CanlabCore
%   https://github.com/canlab/MediationToolbox
% Code adapted from CANlab mediation tutorials
% (canlab.github.io walkthroughs and github.com/canlab/MediationToolbox).
%
% Runtime: about a minute. All data are simulated.

%% 1. Simulate a brain-as-mediator dataset
% We mimic the chapter's running example: a stressor manipulation (X),
% ACC activity as the mediator (M), and stressor-evoked heart rate
% increase as the outcome (Y). One observation per participant
% (single-level mediation), n = 200.
%
% The generative model IS the mediation model, so we know ground truth:
%   M = a*X + noise,   Y = c'*X + b*M + noise
% with true paths a = 0.6, b = 0.5, c' = 0.2.
% Implied total effect: c = c' + a*b = 0.50.

rng(42);                                     % fix the random seed for reproducibility
n = 200;                                     % n = participants (one observation each)
a_true = 0.6; b_true = 0.5; cp_true = 0.2;   % true paths: a (X->M), b (M->Y|X), c' (direct)

X = randn(n, 1);                             % exposure (stressor intensity)
M = a_true .* X + randn(n, 1);               % mediator (ACC activity)
Y = cp_true .* X + b_true .* M + randn(n, 1); % outcome (heart rate)

create_figure('data', 1, 3);
subplot(1, 3, 1); plot(X, M, 'ko'); refline; xlabel('X'); ylabel('M'); title('Path a');
subplot(1, 3, 2); plot(M, Y, 'ko'); refline; xlabel('M'); ylabel('Y'); title('Path b (marginal)');
subplot(1, 3, 3); plot(X, Y, 'ko'); refline; xlabel('X'); ylabel('Y'); title('Total effect c');

%% 2. Estimate the paths with mediation.m, with bootstrap inference
% mediation.m estimates all paths from the two structural regressions
% (M ~ X and Y ~ X + M) and, with the 'boot' option, performs
% bias-corrected, accelerated bootstrap inference on each path --
% including the indirect effect a*b. The Sobel (normal-theory) test is
% overconservative because the product a*b has a skewed sampling
% distribution; the bootstrap respects the skew.
%
% Columns of the paths output: [a  b  c'  c  a*b]

[paths, stats] = mediation(X, Y, M, 'boot', 'verbose', 'bootsamples', 10000, ...
    'names', {'Stressor' 'Heart rate' 'ACC'});

fprintf('\nPath estimates [a b c'' c ab]:\n');
disp(array2table(stats.mean, 'VariableNames', {'a' 'b' 'c_prime' 'c' 'ab'}));

fprintf('Indirect effect a*b = %3.3f, bootstrap p = %3.4f\n', ...
    stats.mean(5), stats.p(5));

% Check the exact decomposition: c - c' = a*b
fprintf('Decomposition check: c - c'' = %3.4f vs. a*b = %3.4f\n', ...
    stats.mean(4) - stats.mean(3), stats.mean(5));

%% 3. Visualize the mediation model
% mediation_path_diagram draws the standard three-variable path diagram
% with coefficients and significance stars. mediation.m with 'plots'
% also produces histograms of the bootstrapped path coefficients --
% look at the a*b histogram and note its skew.

mediation_path_diagram(stats);

[paths2, stats2] = mediation(X, Y, M, 'boot', 'plots', 'verbose', ...
    'names', {'Stressor' 'Heart rate' 'ACC'});  %#ok<ASGLU>

%% 4. An unmodeled confounder creates spurious mediation
% The causal reading of a*b requires NO unmodeled confounding of the
% M-Y relationship. Randomizing X protects paths a and c, but M is only
% observed. Here we simulate a world with NO true mediation (b = 0):
% a lurking variable U (e.g., arousal) drives both M and Y.

U  = randn(n, 1);
M2 = 0.6 .* X + 0.7 .* U + randn(n, 1);              % a = 0.6, confounded
Y2 = 0.3 .* X + 0 .* M2 + 0.7 .* U + randn(n, 1);    % true b = 0 !

% Naive analysis, ignoring U: finds a significant indirect effect
disp('--- Naive mediation (confounder U ignored) ---')
[paths_naive, stats_naive] = mediation(X, Y2, M2, 'boot', 'verbose', ...
    'bootsamples', 10000, 'names', {'Stressor' 'Heart rate' 'ACC'});
fprintf('NAIVE:    a*b = %3.3f, p = %3.4f  <- spurious!\n', ...
    stats_naive.mean(5), stats_naive.p(5));

% Correct analysis: control for U in all regressions with 'covs'
disp('--- Mediation adjusting for the confounder U ---')
[paths_adj, stats_adj] = mediation(X, Y2, M2, 'boot', 'verbose', ...
    'bootsamples', 10000, 'covs', U, 'names', {'Stressor' 'Heart rate' 'ACC'});
fprintf('ADJUSTED: a*b = %3.3f, p = %3.4f  <- no mediation (correct)\n', ...
    stats_adj.mean(5), stats_adj.p(5));

% In real data U is usually unmeasured -- this is why fMRI mediation is
% best framed as pathway discovery rather than causal proof.

%% 5. A three-node path model fit by regression equations
% The recursive SEM of Figure 34.1: ROI 1 drives ROIs 2 and 3 (b12, b13),
% and ROI 2 drives ROI 3 (b23). In matrix form y = B*y + zeta. Because
% the model is recursive (no loops), each structural equation can be
% estimated by OLS:
%   ROI2 ~ ROI1         -> b12
%   ROI3 ~ ROI1 + ROI2  -> b13, b23
% Note the mediation model above is exactly this graph with
% ROI1 = X, ROI2 = M, ROI3 = Y.

T = 300;                                     % T = time points in the simulated series
b12_true = 0.8; b13_true = 0.4; b23_true = 0.5;  % true path coefficients for the three edges

z  = randn(T, 3);                        % independent errors, unit variance
y1 = z(:, 1);
y2 = b12_true .* y1 + z(:, 2);
y3 = b13_true .* y1 + b23_true .* y2 + z(:, 3);

bb12 = glmfit(y1, y2);                   % [intercept; b12]
bb3  = glmfit([y1 y2], y3);              % [intercept; b13; b23]
b12_hat = bb12(2); b13_hat = bb3(2); b23_hat = bb3(3);

fprintf('b12 = %3.3f (true %3.1f), b13 = %3.3f (true %3.1f), b23 = %3.3f (true %3.1f)\n', ...
    b12_hat, b12_true, b13_hat, b13_true, b23_hat, b23_true);

%% 6. Model-implied vs. observed covariance
% Rebuild B and form the model-implied covariance
%   Sigma(theta) = inv(I - B) * Psi * inv(I - B)'
% SEM software estimates parameters by minimizing the discrepancy
% between this matrix and the sample covariance. With 3 variables and
% 6 free parameters (3 paths + 3 error variances) the model is
% just-identified, so it reproduces the sample covariance (nearly)
% exactly -- overall fit can only be tested for over-identified models.

B = [0        0        0;
     b12_hat  0        0;
     b13_hat  b23_hat  0];

Psi = diag([var(y1), var(y2 - b12_hat .* y1), ...
            var(y3 - b13_hat .* y1 - b23_hat .* y2)]);

Sigma_model = (eye(3) - B) \ Psi / (eye(3) - B)';
Sigma_obs   = cov([y1 y2 y3]);

disp('Observed covariance:');      disp(Sigma_obs);
disp('Model-implied covariance:'); disp(Sigma_model);
fprintf('Max abs difference: %3.4f\n', max(abs(Sigma_obs(:) - Sigma_model(:))));

%% 7. Test an individual edge by nested model comparison
% Does the direct ROI1 -> ROI3 connection (b13) improve the model beyond
% the indirect route through ROI2? Dropping the edge gives a nested
% model; for a single OLS path this comparison reduces to the familiar
% partial F-test (equivalently, the t-test on b13).

[b_full, dev_full, stats_full] = glmfit([y1 y2], y3); %#ok<ASGLU>
sse_full = sum(stats_full.resid .^ 2);

[b_red, dev_red, stats_red] = glmfit(y2, y3);         %#ok<ASGLU>
sse_red = sum(stats_red.resid .^ 2);

F = (sse_red - sse_full) / (sse_full / (T - 3));
p = 1 - fcdf(F, 1, T - 3);
fprintf('Edge b13: F(1, %d) = %3.2f, p = %3.2e\n', T - 3, F, p);

%% 8. Where to go next
% - mediation.m supports multilevel (within-person) mediation: pass X,
%   Y, and M as cell arrays with one vector per subject, and add
%   second-level moderators with 'L2M' for moderated mediation.
% - mediation_brain.m runs this analysis at every voxel (Mediation
%   Effect Parametric Mapping), e.g.:
%     results = mediation_brain(X, Y, imgs, 'names', names, ...
%         'mask', mask, 'boot', 'pvals', 5, 'bootsamples', 10000);
% - See the CANlab walkthroughs at canlab.github.io/walkthroughs for
%   brain mediation examples with real data.
