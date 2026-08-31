%% Chapter 33 Lab: Dynamic Connectivity (MATLAB)
% This lab accompanies Chapter 33, "Dynamic Connectivity". You will:
%   1. Compute sliding-window correlations on data with NO true dynamics and
%      watch convincing "dynamics" appear out of pure sampling variability
%   2. Explore the window-length bias-variance tradeoff
%   3. Track a true regime-switching correlation with a sliding window versus
%      an exponentially weighted (DCC-flavored) estimator
%   4. Test apparent dynamics against a stationary null using
%      phase-randomized surrogates
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox.
% Optional but recommended: Lindquist's Dynamic Correlation Toolbox
%   https://github.com/canlab/Lindquist_Dynamic_Correlation
% which provides sliding_window.m, tapered_sliding_window.m, and DCC.m — a
% full GARCH-based Dynamic Conditional Correlation estimator (DCC.m is also
% distributed with CanlabCore). Code below is adapted from that toolbox's
% Example.m and sliding_window.m.
%
% Runtime: about a minute. All data are simulated.

%% 1. Spurious dynamics from a purely static process
% We simulate two time series (600 "TRs") whose true correlation is a
% CONSTANT r = 0.4 -- a static-correlation null. Any variability we see in
% the windowed correlation is sampling noise, not brain dynamics.

rng(33);
T = 600;                 % number of time points (TRs)
r_true = 0.4;            % true, constant correlation
w = 30;                  % sliding-window length (TRs)

Sigma = [1 r_true; r_true 1];
dat_static = mvnrnd([0 0], Sigma, T);      % correlation NEVER changes

% With the Dynamic Correlation Toolbox on your path you could instead call:
%   Ct = sliding_window(dat_static, w);    % p x p x T array
rho_sw = slidecorr(dat_static, w);

figure;
subplot(2,1,1);
plot(dat_static(:,1)); hold on; plot(dat_static(:,2) + 5);
ylabel('Signal (offset)');
title(sprintf('Two stationary series, constant true r = %.1f', r_true));
subplot(2,1,2);
plot(rho_sw, 'r'); yline(r_true, 'k--');
xlabel('Time (TRs)'); ylabel(sprintf('Windowed r (w = %d)', w));

fprintf('Windowed r ranges from %.2f to %.2f with NO true dynamics.\n', ...
    min(rho_sw), max(rho_sw));

% Interpretation: a correlation estimated from only 30 samples has a large
% standard error (in Fisher-z units, about 1/sqrt(w-3) ~ 0.19), so the
% windowed estimate sweeps across a wide range. Cluster matrices like these
% and you will find "brain states" in pure noise. Real fMRI is worse: noise
% autocorrelation and nuisance regression further reduce the effective
% degrees of freedom in each window.

%% 2. The window-length dilemma (bias-variance tradeoff)
% Longer windows average more samples, so the estimates stabilize -- but
% (Section 3) they also smear out any real changes. Windowing acts as a
% low-pass filter with cutoff ~ 1/w, which is why a common rule of thumb
% sets w to the reciprocal of the lowest frequency in the signal.

windows = [15 30 60 120];
figure;
for i = 1:numel(windows)
    wi = windows(i);
    r = slidecorr(dat_static, wi);
    subplot(numel(windows), 1, i);
    plot(r, 'r'); yline(r_true, 'k--'); ylim([-0.5 1]);
    ylabel(sprintf('w = %d', wi));
    fprintf('w = %3d: SD of windowed r = %.3f\n', wi, std(r, 'omitnan'));
end
xlabel('Time (TRs)');
sgtitle('Same static data, four window lengths');

% The SD of the (spurious) fluctuations shrinks roughly like 1/sqrt(w).
% Before concluding "always use the longest window", see what long windows
% do to REAL dynamics next.

%% 3. Tracking a true regime switch: sliding window vs. DCC-flavored EWMA
% Now the ground truth really changes: the correlation alternates between
% +0.7 and -0.3 every 150 TRs (four "brain states" in sequence). We compare
% two sliding windows with an exponentially weighted moving average (EWMA)
% correlation -- a fixed-parameter cousin of Dynamic Conditional
% Correlation. The forgetting factor lambda plays the role of a Kalman
% gain / learning rate: recent samples get weight (1 - lambda) and older
% evidence decays geometrically. Full DCC additionally removes time-varying
% variance with a GARCH step and estimates its smoothing parameters from
% the data by maximum likelihood -- no window length to choose. With the
% toolbox on your path, try:  Ct = DCC(dat_dyn);

r_t = 0.7 * ones(T, 1);
r_t(mod(floor((0:T-1)' / 150), 2) == 1) = -0.3;   % true correlation path

z = randn(T, 2);
dat_dyn = [z(:,1), r_t .* z(:,1) + sqrt(1 - r_t.^2) .* z(:,2)];

rho30 = slidecorr(dat_dyn, 30);
rho90 = slidecorr(dat_dyn, 90);
rho_ew = ewmacorr(dat_dyn, 0.94);

figure;
stairs(r_t, 'k', 'LineWidth', 2); hold on;
plot(rho30); plot(rho90); plot(rho_ew);
legend({'true r', 'sliding w=30', 'sliding w=90', 'EWMA \lambda=0.94'}, ...
    'Location', 'southwest');
xlabel('Time (TRs)'); ylabel('Correlation');
title('Regime-switching correlation: three estimators');

valid = 120:T;   % skip burn-in where the long window is undefined
ests = {rho30, rho90, rho_ew};
names = {'sliding w=30', 'sliding w=90', 'EWMA lam=0.94'};
for i = 1:3
    e = ests{i};
    fprintf('%-15s RMSE vs true r = %.3f\n', names{i}, ...
        sqrt(mean((e(valid) - r_t(valid)).^2, 'omitnan')));
end

% Interpretation: w = 30 follows the transitions but is noisy within each
% regime; w = 90 is smooth but smears every transition across ~90 TRs
% (bias). The EWMA discounts old samples smoothly and typically achieves a
% better noise/adaptation compromise -- and DCC tunes that compromise from
% the data instead of asking you to pick it.

%% 4. Is it really dynamic? Test against a stationary null
% Statistic: the SD of the windowed correlation over time. Null: phase-
% randomized surrogates. Rotating every Fourier coefficient by a random
% phase -- using the SAME rotation for both series at each frequency --
% preserves each power spectrum (autocorrelation) and the cross-spectrum
% (static correlation), but the surrogate is stationary by construction:
% any true correlation dynamics are destroyed.

nsur = 200;
[obs_s, p_s] = stationaritytest(dat_static, w, nsur);
[obs_d, p_d] = stationaritytest(dat_dyn, w, nsur);

fprintf('Static-null data:      SD of windowed r = %.3f, p = %.3f\n', obs_s, p_s);
fprintf('Regime-switching data: SD of windowed r = %.3f, p = %.3f\n', obs_d, p_d);

% The static data should be non-significant (its variability is exactly
% what a stationary process produces), while the regime-switching data
% should reject decisively. This is the discipline to carry into real
% analyses: never interpret raw sliding-window variability -- or states
% clustered from it -- without a stationary null, careful denoising, and
% checks on motion and arousal.

%% Local functions

function rho = slidecorr(dat, w)
% Sliding-window correlation between the two columns of dat.
% The window ends at time t (samples t-w+1 ... t), mirroring
% sliding_window.m in Lindquist's Dynamic Correlation Toolbox.
T = size(dat, 1);
rho = NaN(T, 1);
for t = w:T
    c = corr(dat(t-w+1:t, :));
    rho(t) = c(1, 2);
end
end

function rho = ewmacorr(dat, lam)
% Exponentially weighted moving correlation (DCC-flavored, fixed lambda).
% Full DCC (see DCC.m in the Dynamic Correlation Toolbox / CanlabCore)
% also fits GARCH variances and estimates lambda-like parameters by
% maximum likelihood.
T = size(dat, 1);
Q = cov(dat);                        % initialize from the static covariance
rho = zeros(T, 1);
for t = 1:T
    zt = dat(t, :)';
    Q = lam * Q + (1 - lam) * (zt * zt');
    rho(t) = Q(1, 2) / sqrt(Q(1, 1) * Q(2, 2));
end
end

function [obs, p] = stationaritytest(dat, w, nsur)
% Compare SD of windowed correlation to a phase-randomized stationary null.
T = size(dat, 1);
rho = slidecorr(dat, w);
obs = std(rho, 'omitnan');
F = fft(dat);
half = 2:ceil(T/2);
nullsd = zeros(nsur, 1);
for s = 1:nsur
    ph = zeros(T, 1);
    ph(half) = 2 * pi * rand(numel(half), 1);  % same rotation, both series
    ph(T - half + 2) = -ph(half);              % conjugate symmetry
    sur = real(ifft(F .* exp(1i * ph)));
    rs = slidecorr(sur, w);
    nullsd(s) = std(rs, 'omitnan');
end
p = (1 + sum(nullsd >= obs)) / (nsur + 1);
end
