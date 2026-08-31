%% Chapter 36 Lab: Granger Causal Models (MATLAB)
% This lab accompanies Chapter 36, "Granger Causal Models". You will
% simulate a VAR(1) system with a genuine directed influence and recover it
% with Granger tests; then build the critical counterexample: two regions
% with perfectly symmetric neural coupling but different hemodynamic
% response latencies, where Granger analysis of the BOLD signals reports a
% confident -- and spurious -- directed influence. Finally, you will apply
% a deconvolution remedy.
%
% Requirements: SPM12 on your MATLAB path (for spm_hrf) and the Statistics
% and Machine Learning Toolbox (for fcdf).
%   https://www.fil.ion.ucl.ac.uk/spm/
%
% Runtime: under a minute. All data are simulated.

%% 1. Simulate a VAR(1) system with a true directed influence
% A first-order vector autoregressive model generates each new sample of
% the bivariate state Z(t) = [X(t) Y(t)] from the previous sample plus
% white noise:  Z(t) = A * Z(t-1) + noise.
% Read the coupling matrix A row by row: X depends only on its own past
% (0.5), while Y depends on its own past (0.5) AND on X's past (0.4). So
% the ground truth is a one-way influence: X drives Y at lag 1, with no
% reverse influence. Think of each time step as one TR.

rng(7);
n = 400;                              % time points (think: TRs)
A = [0.5 0.0; ...                     % X(t) <- 0.5*X(t-1)
     0.4 0.5];                        % Y(t) <- 0.4*X(t-1) + 0.5*Y(t-1)

Z = zeros(n, 2);
for t = 2:n
    Z(t, :) = (A * Z(t-1, :)')' + randn(1, 2);
end
X = Z(:, 1);
Y = Z(:, 2);

figure;
plot(1:120, X(1:120), 1:120, Y(1:120), 'LineWidth', 1);
legend('X (driver)', 'Y (receiver)');
xlabel('time (samples)'); ylabel('signal');
title('VAR(1) simulation: X drives Y at lag 1');

%% 2. Test Granger causality in both directions
% A Granger test is a nested-model comparison. To ask "does X Granger
% cause Y?" we compare:
%   Restricted: Y(t) ~ Y(t-1), ..., Y(t-p)          (Y's own history only)
%   Full:       Y(t) ~ Y(t-1..t-p), X(t-1..t-p)     (plus X's history)
% If the full model fits significantly better (F-test on the reduction in
% residual sum of squares), X Granger causes Y. The helper granger_F at
% the bottom of this script implements the test for any lag order p.

[F_xy, p_xy] = granger_F(Y, X, 1);    % X -> Y (true influence)
[F_yx, p_yx] = granger_F(X, Y, 1);    % Y -> X (no influence)

fprintf('X -> Y:  F = %8.2f,  p = %.3g   (true influence)\n', F_xy, p_xy);
fprintf('Y -> X:  F = %8.4f,  p = %.3g   (no influence)\n',  F_yx, p_yx);

% The test recovers the truth emphatically: a huge F for X -> Y and an F
% near zero for Y -> X. Temporal precedence works perfectly here because
% the neural lag is visible at the sampling rate and both series are
% measured the same way -- no differential measurement delay. Hold that
% thought for Section 4.

%% 3. Geweke's directed influence measures
% Geweke quantified directed influence with the log ratio of residual
% variances between the restricted and full models:
%   F_{X->Y} = ln( RSS_restricted / RSS_full )
% If X's history helps predict Y, the full model's residual variance
% shrinks and the measure is > 0; otherwise it is ~0. The difference
% F_{X->Y} - F_{Y->X} is commonly used to infer the NET direction of
% influence -- a convention that will matter below.

Fg_xy = geweke_F(Y, X, 1);
Fg_yx = geweke_F(X, Y, 1);
fprintf('Geweke F_(X->Y) = %.4f\n', Fg_xy);
fprintf('Geweke F_(Y->X) = %.4f\n', Fg_yx);
fprintf('Net influence (X->Y positive): %+.4f\n', Fg_xy - Fg_yx);

%% 4. The hemodynamic confound: equal coupling, unequal HRF lags
% fMRI measures BOLD, not neural activity, and the hemodynamic response
% function varies across brain regions. Now the neural truth is perfectly
% SYMMETRIC coupling (0.3 in each direction; no net directed influence),
% but region 1's HRF peaks early (~4 s) and region 2's peaks late (~7 s)
% -- a difference well within the physiological range.

TR = 1; n2 = 1000; burn = 50;
A_sym = [0.4 0.3; ...                 % region 1 <- itself + region 2, equally
         0.3 0.4];                    % region 2 <- itself + region 1, equally

Zn = zeros(n2 + burn, 2);
for t = 2:(n2 + burn)
    Zn(t, :) = (A_sym * Zn(t-1, :)')' + randn(1, 2);
end
Zn = Zn(burn+1:end, :);               % discard initial transient
neu1 = Zn(:, 1);
neu2 = Zn(:, 2);

% Two plausible regional HRFs (SPM double-gamma; p(1) = delay of response)
h_fast = spm_hrf(TR, [4 16 1 1 6 0 32]);  h_fast = h_fast / max(h_fast);
h_slow = spm_hrf(TR, [7 16 1 1 6 0 32]);  h_slow = h_slow / max(h_slow);

figure;
plot(0:numel(h_fast)-1, h_fast, 0:numel(h_slow)-1, h_slow, 'LineWidth', 1);
legend('region 1 HRF (peak ~4 s)', 'region 2 HRF (peak ~7 s)');
xlabel('time (s)'); ylabel('response (a.u.)');
title('Two plausible regional HRFs');

% Pass each region's neural series through ITS OWN HRF, plus noise
b1 = conv(neu1, h_fast); bold1 = b1(1:n2) + 0.05 * randn(n2, 1);
b2 = conv(neu2, h_slow); bold2 = b2(1:n2) + 0.05 * randn(n2, 1);

figure;
plot(100:199, zscore(bold1(101:200)), 100:199, zscore(bold2(101:200)), ...
    'LineWidth', 1);
legend('BOLD region 1 (fast HRF)', 'BOLD region 2 (slow HRF)');
xlabel('time (s)'); ylabel('z-scored signal');
title('Region 1 leads region 2 -- for purely vascular reasons');

% Granger analysis at both levels: the (unobservable) neural series and
% the (observable) BOLD series.
fprintf('\n%-18s  %10s  %10s\n', '', '1 -> 2', '2 -> 1');
fprintf('%-18s  %10.4f  %10.4f   (symmetric, as designed)\n', ...
    'neural (truth)', geweke_F(neu2, neu1, 1), geweke_F(neu1, neu2, 1));
fprintf('%-18s  %10.4f  %10.4f   (asymmetric: SPURIOUS)\n', ...
    'BOLD (confounded)', geweke_F(bold2, bold1, 1), geweke_F(bold1, bold2, 1));

% At the neural level the two directed influences are nearly identical.
% At the BOLD level the fast-HRF region appears to drive the slow-HRF
% region several times more strongly than the reverse, and both tests are
% "highly significant". Anyone using the net-influence difference would
% confidently -- and wrongly -- conclude region 1 drives region 2. Only
% the measurement changed, not the neural dynamics: temporal precedence
% in BOLD can be hemodynamic rather than neuronal.

%% 5. A deconvolution remedy
% One proposed fix: deconvolve each region's HRF first, reconstructing a
% neural-like signal, and run the Granger analysis on that. Convolution
% is linear, b = H*z with H a Toeplitz matrix built from the HRF, so we
% invert it with ridge-regularized least squares:
%   z_hat = (H'H + lambda*I) \ (H'b)
% Here we grant ourselves a luxury real analyses never have: the TRUE
% regional HRFs. In practice the HRF must be estimated from the data.

lam = 0.05;
neu1_hat = deconvolve(bold1, h_fast, lam);
neu2_hat = deconvolve(bold2, h_slow, lam);

fprintf('\n%-18s  %10s  %10s\n', '', '1 -> 2', '2 -> 1');
fprintf('%-18s  %10.4f  %10.4f\n', 'neural (truth)', ...
    geweke_F(neu2, neu1, 1), geweke_F(neu1, neu2, 1));
fprintf('%-18s  %10.4f  %10.4f\n', 'BOLD (confounded)', ...
    geweke_F(bold2, bold1, 1), geweke_F(bold1, bold2, 1));
fprintf('%-18s  %10.4f  %10.4f\n', 'deconvolved BOLD', ...
    geweke_F(neu2_hat, neu1_hat, 1), geweke_F(neu1_hat, neu2_hat, 1));

% Deconvolution removes the spurious asymmetry -- but the recovered
% influence magnitudes fall far below the neural truth. Even with the
% true HRFs, regularized deconvolution smooths away much of the lag-1
% coupling: the remedy trades a false direction for a large loss of
% sensitivity. Shorter TRs help on both fronts.

%% 6. Takeaways
% * Granger causality formalizes temporal precedence: X Granger causes Y
%   if X's past improves prediction of Y beyond Y's own past.
% * It is exploratory -- no a priori structural model, unlike SEM/DCM --
%   and its claim is about prediction, not mechanism.
% * On neural-scale signals it recovers the true direction beautifully.
% * On BOLD, regional HRF latency differences of 1-2 s masquerade as
%   directed neural influence between symmetrically coupled regions.
% * Deconvolution and faster sampling are partial remedies: they remove
%   the spurious direction at a cost in sensitivity.
% * VAR models assume stationarity: analyze stationary segments, or use
%   models that account for transitions between states.

%% Local functions

function [F, p] = granger_F(target, driver, lagp)
% F-test of "driver Granger-causes target" at lag order lagp.
n = numel(target);
T = (lagp + 1):n;
Xr = ones(numel(T), 1);
Xf = ones(numel(T), 1);
for j = 1:lagp
    Xr = [Xr target(T - j)];                    %#ok<AGROW>
    Xf = [Xf target(T - j) driver(T - j)];      %#ok<AGROW>
end
y = target(T);
rr = y - Xr * (Xr \ y);                         % restricted residuals
rf = y - Xf * (Xf \ y);                         % full residuals
dfe = numel(T) - size(Xf, 2);
F = ((rr' * rr - rf' * rf) / lagp) / ((rf' * rf) / dfe);
p = 1 - fcdf(F, lagp, dfe);
end

function Fg = geweke_F(target, driver, lagp)
% Geweke directed influence: ln(RSS_restricted / RSS_full).
n = numel(target);
T = (lagp + 1):n;
Xr = ones(numel(T), 1);
Xf = ones(numel(T), 1);
for j = 1:lagp
    Xr = [Xr target(T - j)];                    %#ok<AGROW>
    Xf = [Xf target(T - j) driver(T - j)];      %#ok<AGROW>
end
y = target(T);
rr = y - Xr * (Xr \ y);
rf = y - Xf * (Xf \ y);
Fg = log((rr' * rr) / (rf' * rf));
end

function z_hat = deconvolve(bold, h, lam)
% Ridge-regularized deconvolution of an HRF from a time series.
n = numel(bold);
col = zeros(n, 1);
col(1:numel(h)) = h;
H = toeplitz(col, [col(1) zeros(1, n - 1)]);
z_hat = (H' * H + lam * eye(n)) \ (H' * bold);
end
