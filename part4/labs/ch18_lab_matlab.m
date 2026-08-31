%% Chapter 18 Lab: The General Linear Model (MATLAB)
% This lab accompanies Chapter 18, "The General Linear Model and Foundations
% of Analysis". You will build a design matrix from event onsets, fit a GLM
% to a simulated voxel time series with ordinary least squares (OLS),
% examine residuals and model fit, and test a simple contrast.
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part4/ch18-the-general-linear-model-and
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples).
%
% Runtime: under a minute. All data are simulated.

%% 1. Define the experimental design: event onsets
% We simulate a 6-minute run (180 volumes, TR = 2 s) with two event types,
% A and B (e.g., famous and non-famous faces). Onsets are in seconds.

TR = 2;                 % repetition time (s)
n_scans = 180;          % number of volumes
run_len = n_scans * TR; % run length in seconds

ons = {};
ons{1} = [10 60 110 160 210 260 310]';   % Condition A onsets (s)
ons{2} = [35 85 135 185 235 285 335]';   % Condition B onsets (s)

%% 2. Build the design matrix by HRF convolution
% onsets2fmridesign builds indicator functions at high temporal resolution,
% convolves them with a canonical HRF (SPM's double-gamma by default), and
% samples the result at the TR. The intercept is added as the LAST column.

X = onsets2fmridesign(ons, TR, run_len);

disp(size(X))     % 180 x 3: [Condition A, Condition B, intercept]

%% 3. Visualize the design
% plotDesign shows the onsets ("stick functions") and the convolved
% regressors. Viewing X as a grayscale image is the classic "design matrix"
% view: one row per volume, one column per regressor.

plotDesign(ons, [], TR);

create_figure('design matrix image');
imagesc(X); colormap gray;
set(gca, 'XTick', 1:3, 'XTickLabel', {'A' 'B' 'Intercept'});
xlabel('Regressor'); ylabel('Time (TRs)');
title('Design matrix X');

%% 4. Simulate a voxel time series with known ground truth
% We create data from the model y = X*beta + noise, so we know the right
% answer: beta_A = 0.8, beta_B = 0.4, intercept = 100 (arbitrary units).
% Noise is IID Gaussian for now -- Chapter 19 adds autocorrelation.

rng(9);                               % seed for reproducible noise
beta_true = [0.8 0.4 100]';           % true effects: [A; B; intercept]
sigma_noise = 0.5;                    % noise SD; modest so this 7-event demo recovers betas clearly

y = X * beta_true + sigma_noise .* randn(n_scans, 1);

create_figure('voxel time series');
plot((0:n_scans-1) * TR, y, 'k-');
xlabel('Time (s)'); ylabel('BOLD signal (a.u.)');
title('Simulated voxel time series');

%% 5. Fit the GLM with ordinary least squares
% The OLS solution in matrix form:  beta_hat = (X'X)^{-1} X'y
% (We use the backslash operator, which is numerically preferable to
% explicitly inverting X'X.)

beta_hat = (X' * X) \ (X' * y);

disp(table(beta_true, beta_hat, ...
    'RowNames', {'Condition A', 'Condition B', 'Intercept'}))

% The estimates should be close to, but not exactly, the true values --
% the difference reflects sampling error from the noise.

%% 6. Residuals and model fit
% Residuals r = y - X*beta_hat are what the model cannot explain. Their
% variance (scaled by the error degrees of freedom, dfe = n - p) estimates
% sigma^2, and R^2 summarizes the proportion of variance explained.

fits = X * beta_hat;                  % fitted values
r    = y - fits;                      % residuals

[n, p] = size(X);
dfe    = n - p;                       % error degrees of freedom
sigma2 = (r' * r) / dfe;              % estimate of error variance

SStot = sum((y - mean(y)) .^ 2);
R2    = 1 - (r' * r) / SStot;

fprintf('dfe = %d, sigma2_hat = %3.2f (true = %3.2f), R^2 = %3.2f\n', ...
    dfe, sigma2, sigma_noise ^ 2, R2);

create_figure('fits and residuals', 2, 1);
subplot(2, 1, 1);
plot((0:n-1) * TR, y, 'k-'); hold on;
plot((0:n-1) * TR, fits, 'r-', 'LineWidth', 2);
legend({'Observed' 'Fitted'}); ylabel('Signal');
title('Observed and fitted time series');
subplot(2, 1, 2);
plot((0:n-1) * TR, r, 'b-');
xlabel('Time (s)'); ylabel('Residual');
title('Residuals (should look like noise around zero)');

%% 7. Standard errors, t-values, and P values
% Var(beta_hat) = sigma^2 * (X'X)^{-1}. The standard error of each beta is
% the square root of the corresponding diagonal element. t = beta / SE is
% compared with a Student's t distribution with dfe degrees of freedom.

XtX_inv = inv(X' * X);
se      = sqrt(sigma2 * diag(XtX_inv));
t_vals  = beta_hat ./ se;
p_vals  = 2 * (1 - tcdf(abs(t_vals), dfe));   % two-tailed

disp(table(beta_hat, se, t_vals, p_vals, ...
    'RowNames', {'Condition A', 'Condition B', 'Intercept'}))

%% 8. A simple contrast: A - B
% Contrasts are linear combinations of betas, c'*beta. The contrast
% [1 -1 0] tests whether Condition A evokes a larger response than B.
% Its t-value uses the same machinery:
%   t = c'*beta_hat / sqrt(sigma2 * c' * (X'X)^{-1} * c)

c = [1 -1 0]';                        % contrast weights: +1 for A, -1 for B, 0 for intercept

con_val = c' * beta_hat;                        % contrast value (0.8 - 0.4 = 0.4 expected)
se_con  = sqrt(sigma2 * c' * XtX_inv * c);      % contrast standard error
t_con   = con_val / se_con;
p_con   = 2 * (1 - tcdf(abs(t_con), dfe));

fprintf('Contrast A - B: value = %3.2f, t(%d) = %3.2f, p = %3.4f\n', ...
    con_val, dfe, t_con, p_con);

%% 9. Cross-check with glmfit
% MATLAB's glmfit (Statistics and Machine Learning Toolbox) should agree
% with our matrix-form solution. glmfit adds its own intercept as the FIRST
% column, so we pass only the task regressors.

[b_glmfit, ~, glm_stats] = glmfit(X(:, 1:2), y);

disp('glmfit betas [intercept; A; B]:'); disp(b_glmfit)
disp('glmfit t-values:'); disp(glm_stats.t)

% The betas, t-values, and P values match our hand computation -- the GLM
% is the same model whether you write the algebra yourself or call a
% packaged routine. CANlab's fmri_data.regress method applies exactly this
% model to every voxel of a 4-D image object at once.

%% Explore on your own
% 1. Move the B onsets closer to the A onsets (e.g., 3 s after each A
%    event). What happens to the correlation between regressors, the
%    standard errors, and the t-values? (This is design efficiency --
%    Chapter 27.)
% 2. Add a linear drift to the simulated data (y = y + 0.02*(1:n)') without
%    modeling it. How do the estimates change? Then add a drift regressor
%    to X and refit.
% 3. Replace the IID noise with AR(1) noise (e = filter(1, [1 -0.5], ...
%    randn(n,1))) and compare the empirical false positive rate of OLS
%    over many simulations. This motivates prewhitening (Chapter 19).
