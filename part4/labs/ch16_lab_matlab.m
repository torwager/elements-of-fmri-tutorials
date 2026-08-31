%% Chapter 16 Lab: Artifacts and Noise in fMRI (MATLAB)
% This lab accompanies Chapter 16, "Artifacts and Noise in fMRI". You will
% build an fMRI-like time series from its noise ingredients -- slow drift,
% AR(1) autocorrelated noise, transient spikes, and a heartbeat sampled at
% the TR -- then diagnose each one: view the noise in the frequency domain,
% demonstrate temporal aliasing, detect transient outliers with
% RMSSD/DVARS and Mahalanobis distance, and measure how autocorrelation
% inflates naive OLS false positive rates.
%
% Requirements: CanlabCore on your MATLAB path (noise_arp), plus the
% Statistics and Machine Learning Toolbox (tinv, chi2inv, mahal).
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples; see scnlab_outlier_id.m for the full QC workflow).
%
% Runtime: under a minute. All data are simulated.

%% 1. Simulate a noisy voxel from known components
% An 8-minute run: 240 volumes at TR = 2 s. On top of a baseline of 100 we
% add slow drift (sum of slow cosines), AR(1) noise with phi = 0.5
% (each error inherits half of the previous deviation:
% e_t = phi*e_{t-1} + eta_t), a heartbeat simulated at high resolution and
% sampled once per TR, and three transient spikes.

rng(7);
TR = 2; n = 240;
t = (0:n-1)' * TR;                       % volume acquisition times (s)

% Slow drift: sum of slow cosines (1/f-like low-frequency power)
drift = 6*cos(2*pi*t/400) + 3*cos(2*pi*t/180);

% AR(1) noise, phi = 0.5 (CANlab noise generator)
phi = 0.5;
ar1 = noise_arp(n, phi);

% Heartbeat at ~60 bpm with beat-to-beat jitter, built at 50 Hz
dt = 0.02;
t_hi = (0:dt:(n*TR - dt))';
beats = cumsum(normrnd(1, 0.05, round(n*TR*1.3), 1));
beats = beats(beats < t_hi(end));
hb_hi = zeros(size(t_hi));
for i = 1:length(beats)
    hb_hi = hb_hi + exp(-0.5 * ((t_hi - beats(i)) ./ 0.06) .^ 2);
end
hb_tr = hb_hi(round(t / dt) + 1);        % what the scanner sees: 1/TR

% Transient spikes: sudden shifts in 3 volumes
spikes = zeros(n, 1);
spike_vols = [61 62 151];                % MATLAB 1-based volume indices
spikes(spike_vols) = [18 -12 15];

% Composite voxel time series
y = 100 + drift + ar1 + 1.5*hb_tr + spikes;

create_figure('noise components', 5, 1);
comps = {drift, ar1, 1.5*hb_tr, spikes, y};
names = {'Slow drift', 'AR(1) noise (\phi = 0.5)', ...
    'Heartbeat sampled at TR', 'Transient spikes', 'Composite y'};
for i = 1:5
    subplot(5, 1, i);
    plot(t, comps{i}, 'k-', 'LineWidth', 1);
    ylabel(names{i});
end
xlabel('Time (s)');

% Note the heartbeat panel: the true signal beats once per second, but
% sampled every 2 s it looks like slow, irregular wander. That is temporal
% aliasing -- Section 3 shows why.

%% 2. The frequency domain: power spectral density
% The Fourier transform re-expresses the time series as sine waves, each
% with a magnitude and phase. Plotting power against frequency (the PSD)
% shows where the noise lives. Drift dominates the lowest frequencies.

nfft = n;
f = (0:floor(nfft/2)-1)' / (nfft * TR);          % frequency axis (Hz)
P = abs(fft(y - mean(y))) .^ 2;
P = P(1:floor(nfft/2));

create_figure('psd');
semilogy(f(2:end), P(2:end), 'k-', 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('Power (log scale)');
title('PSD of the composite voxel: drift dominates low frequencies');

%% 3. Temporal aliasing: the Nyquist limit
% The Nyquist theorem: a sampling rate Fs represents frequencies only up
% to Fs/2 = 1/(2*TR). Higher frequencies are reflected around the limit
% and reappear below it. Classic demo: a 10 Hz sine sampled at 12 Hz
% (Nyquist = 6 Hz) masquerades as a 12 - 10 = 2 Hz oscillation.

fs_hi = 1000; fs_lo = 12; dur = 5; f_sig = 10;
tt_hi = (0:1/fs_hi:dur)';
tt_lo = (0:1/fs_lo:dur)';
x_hi = sin(2*pi*f_sig*tt_hi);
x_lo = sin(2*pi*f_sig*tt_lo);

create_figure('aliasing', 2, 1);
subplot(2, 1, 1);
plot(tt_hi, x_hi, 'k-', 'LineWidth', 0.5); hold on;
plot(tt_lo, x_lo, 'm.-', 'LineWidth', 1.5);
set(gca, 'XLim', [0 2]);
xlabel('Time (s)');
title('A 10 Hz signal sampled at 12 Hz appears to oscillate at 2 Hz');
legend({'10 Hz signal' 'sampled at 12 Hz'});

% Frequency domain: compute simple periodograms via fft
subplot(2, 1, 2);
fh = (0:floor(length(x_hi)/2)-1)' * fs_hi / length(x_hi);
Ph = abs(fft(x_hi)).^2; Ph = Ph(1:floor(length(x_hi)/2)) ./ max(Ph);
fl = (0:floor(length(x_lo)/2)-1)' * fs_lo / length(x_lo);
Pl = abs(fft(x_lo)).^2; Pl = Pl(1:floor(length(x_lo)/2)) ./ max(Pl);
plot(fh, Ph, 'k-', 'LineWidth', 1); hold on;
plot(fl, Pl, 'm-', 'LineWidth', 1.5);
plot([fs_lo/2 fs_lo/2], [0 1], 'k:');
set(gca, 'XLim', [0 12]);
xlabel('Frequency (Hz)'); ylabel('Normalized power');
legend({'original (10 Hz)' 'aliased (2 Hz)' 'Nyquist limit (6 Hz)'});

% For fMRI at TR = 2 s the Nyquist limit is 1/(2*2) = 0.25 Hz, so our
% ~1 Hz heartbeat aliases into the same low-frequency band as the task --
% and the aliased pattern shifts with small beat-timing changes, so it
% cannot be cleanly filtered or modeled with simple covariates.

%% 4. Multi-voxel outlier detection: RMSSD/DVARS and Mahalanobis distance
% Real QC works on image volumes. We simulate a small "brain" of 200
% voxels sharing the drift and heartbeat (with voxel-specific weights)
% plus voxel-specific AR(1) noise, and corrupt the same three volumes with
% spatially widespread intensity shifts.

V = 200;
w_drift = 0.2 + 1.3*rand(1, V);          % voxel-specific loadings
w_hb = rand(1, V);
noisemat = zeros(n, V);
for v = 1:V
    noisemat(:, v) = noise_arp(n, phi);
end
Y = 100 + drift*w_drift + (1.5*hb_tr)*w_hb + noisemat;
Y(spike_vols, :) = Y(spike_vols, :) + normrnd(0, 6, length(spike_vols), V);

% --- RMSSD (DVARS): root-mean-square image-to-image change --------------
% BOLD changes slowly; artifacts are immediate. Large successive
% differences flag transient outliers. (See also scnlab_outlier_id.m and
% fmri_data.outliers in CanlabCore for the full real-data workflow.)

rmssd = sqrt(mean(diff(Y, 1, 1) .^ 2, 2));      % length n-1
z_rmssd = (rmssd - mean(rmssd)) ./ std(rmssd);
flag_rmssd = find(z_rmssd > 3) + 1;             % diff t-1 -> t flags t

create_figure('rmssd');
plot(2:n, rmssd, 'k-', 'LineWidth', 0.75); hold on;
thr = mean(rmssd) + 3*std(rmssd);
plot([1 n], [thr thr], 'k--');
plot(flag_rmssd, rmssd(flag_rmssd - 1), 'ro', 'MarkerFaceColor', 'r');
xlabel('Time (images)'); ylabel('RMSSD (DVARS)');
title('Image-to-image change flags the corrupted volumes');

disp('True corrupted volumes:'); disp(spike_vols);
disp('Flagged by RMSSD > 3 SD:'); disp(flag_rmssd');

% Note: a spike at volume t produces large differences both into and out
% of the artifact, so the volume after a spike may be flagged as well.
% Adjacent flagged volumes are usually all treated as outliers.

%% 5. Mahalanobis distance: multivariate outliers
% Reduce the images to their first 5 principal components and compute each
% volume's Mahalanobis distance from the multivariate cloud. Under
% approximate normality, squared distances follow a chi-square(5)
% distribution, giving a principled threshold.
%
% With real data and CANlab objects this is one line:
%   [ds, ~, ~, ~, wh_outlier_corr] = mahal(fmri_dat_object, 'noplot', 'corr');
% (see canlab_simple_find_outliers.m in CANlab_help_examples)

Yc = Y - mean(Y);
[U, S, ~] = svd(Yc, 'econ');
k = 5;
scores = U(:, 1:k) * S(1:k, 1:k);        % n x k PC scores per volume

md2 = mahal(scores, scores);             % squared Mahalanobis distance
chi2_thresh = chi2inv(0.999, k);         % expect ~0.1% false flags
flag_mahal = find(md2 > chi2_thresh);

create_figure('mahalanobis');
plot(1:n, md2, 'k-', 'LineWidth', 0.75); hold on;
plot([1 n], [chi2_thresh chi2_thresh], 'k--');
plot(flag_mahal, md2(flag_mahal), 'ro', 'MarkerFaceColor', 'r');
xlabel('Time (images)'); ylabel('Mahalanobis d^2');
title('Multivariate distance from the cloud of images');

disp('Flagged by Mahalanobis:'); disp(flag_mahal');

%% 6. Handle outliers with spike regressors (not deletion)
% Add one indicator column per flagged volume to the design matrix. This
% removes each bad image's influence while preserving the timing and noise
% structure. Deleting time points instead (1) changes the noise properties,
% (2) forces re-timing of the design, and (3) can bias sampling when
% outliers are task-correlated.

flagged = union(flag_rmssd, flag_mahal);
spike_regs = zeros(n, length(flagged));
for i = 1:length(flagged)
    spike_regs(flagged(i), i) = 1;
end
fprintf('Spike regressor matrix: %d x %d -> append to design matrix X\n', ...
    size(spike_regs, 1), size(spike_regs, 2));

%% 7. Autocorrelation: empirical vs. theoretical ACF
% For AR(1) noise the autocorrelation function decays geometrically:
% corr(e_t, e_{t-k}) = phi^k. Check a long sample against theory.

e_long = noise_arp(5000, phi);
max_lag = 15;
acf = ones(max_lag + 1, 1);
for k_lag = 1:max_lag
    c = corrcoef(e_long(1:end-k_lag), e_long(1+k_lag:end));
    acf(k_lag + 1) = c(1, 2);
end

create_figure('acf');
stem(0:max_lag, acf, 'k'); hold on;
plot(0:max_lag, phi .^ (0:max_lag), 'r.-', 'LineWidth', 1);
xlabel('Lag (TRs)'); ylabel('Correlation');
title('ACF of AR(1) noise, \phi = 0.5');
legend({'empirical' 'theoretical \phi^k'});

%% 8. Autocorrelation inflates naive OLS false positives
% Null simulation: fit a blocked task regressor (40-s on/off) to pure
% AR(1) noise -- no signal at all -- 2000 times, and count how often naive
% OLS declares p < .05. Valid inference would give ~5%.

nsim = 2000;
X = [double(sin(2*pi*t/80) > 0), ones(n, 1)];
p = size(X, 2);
tcrit = tinv(0.975, n - p);
cinv = inv(X'*X);

nfp = 0;
for i = 1:nsim
    e = noise_arp(n, phi);               % null data: noise only
    b = X \ e;
    r = e - X*b;
    se = sqrt((r'*r) / (n - p) * cinv(1, 1));
    nfp = nfp + (abs(b(1)/se) > tcrit);
end
fpr_naive = nfp / nsim;

%% 9. Prewhitening restores valid inference
% Transform y and X so the errors become independent. For known phi the
% exact AR(1) whitening matrix is lower-bidiagonal: row t computes
% y_t - phi*y_{t-1}, and row 1 is scaled by sqrt(1 - phi^2). Real software
% estimates phi from the residuals and iterates (Chapters 18-19); see
% fit_gls.m in CanlabCore for a full feasible-GLS implementation.

W = eye(n) - phi * diag(ones(n-1, 1), -1);
W(1, 1) = sqrt(1 - phi^2);
Xw = W * X;
cinvw = inv(Xw'*Xw);

nfp = 0;
for i = 1:nsim
    e = noise_arp(n, phi);
    ew = W * e;
    b = Xw \ ew;
    r = ew - Xw*b;
    se = sqrt((r'*r) / (n - p) * cinvw(1, 1));
    nfp = nfp + (abs(b(1)/se) > tcrit);
end
fpr_white = nfp / nsim;

fprintf('False positive rate, naive OLS on AR(1) noise: %.3f\n', fpr_naive);
fprintf('False positive rate after prewhitening:        %.3f\n', fpr_white);
fprintf('Nominal rate: 0.050\n');

% Naive OLS produces several times too many false positives -- every one a
% "significant activation" in pure noise. Prewhitening restores the
% nominal 5% rate. This is why fMRI software estimates the noise
% autocorrelation and applies generalized least squares at the first level.
