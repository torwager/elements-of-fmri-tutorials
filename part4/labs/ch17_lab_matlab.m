%% Chapter 17 Lab: Image Preprocessing (MATLAB)
% This lab accompanies Chapter 17, "Image Preprocessing". Real preprocessing
% runs inside packages like SPM and fMRIPrep -- here we simulate the
% *concepts* behind four key steps so you can see what each one fixes and
% what it costs: (1) head motion and nuisance regression, (2) slice-timing
% offsets, (3) smoothing-kernel tradeoffs, and (4) high-pass filter design
% for drift removal.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab), including
% CANlab_help_examples: linear_filtering_a_timeseries.m.
%
% Runtime: under a minute. All data are simulated.
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part4/ch17-image-preprocessing

%% 1. Head motion: realignment parameters as nuisance regressors
% Realignment estimates six rigid-body parameters per volume: translations
% in x, y, z (mm) and rotations (radians). We simulate a 200-volume run in
% which the head drifts slowly (random walk) plus one sudden jerk, and a
% voxel whose signal mixes "neural" fluctuations, motion-coupled artifact,
% and noise.

rng(17);                              % seed for reproducibility
n = 200; TR = 1;                      % n = volumes (frames); TR = repetition time (s)

mp = cumsum(0.02 * randn(n, 6));      % 6 motion parameters (random walk)
mp(:, 4:6) = mp(:, 4:6) * 0.005;      % rotations are small (radians)
mp(120:124, 1:3) = mp(120:124, 1:3) + 0.8;    % 0.8 mm jerk at frame 120
mp(120:124, 4:6) = mp(120:124, 4:6) + 0.004;  % small rotation component

neural = noise_arp(n, [.5 .1]);       % AR(2) "neural" fluctuations
w = [2 1.5 1 300 200 100]';           % motion coupling weights
motion_artifact = mp * w;
y = neural + motion_artifact + 0.5 * randn(n, 1);

create_figure('motion parameters', 2, 1);
subplot(2, 1, 1); plot(mp(:, 1:3)); ylabel('Translation (mm)');
legend({'x' 'y' 'z'}); title('Simulated realignment parameters');
subplot(2, 1, 2); plot(rad2deg(mp(:, 4:6))); ylabel('Rotation (deg)');
xlabel('Frame'); legend({'pitch' 'roll' 'yaw'});

%% 1a. Variance explained by motion, and spike regression
% Regress the voxel on the six parameters (plus intercept): the R^2 is the
% fraction of variance the realignment parameters account for. Then compute
% framewise displacement (FWD): the sum of absolute frame-to-frame changes,
% rotations converted to mm as arc length on a 50-mm-radius sphere. Frames
% over threshold get one spike regressor each ("spike regression"), which
% censors them while preserving temporal structure.

X = [mp intercept_model(n)];
r = y - X * (X \ y);
R2 = 1 - var(r) / var(y);
fprintf('Variance explained by 6 motion parameters: R^2 = %.2f\n', R2);

fwd = sum(abs(diff(mp(:, 1:3))), 2) + 50 * sum(abs(diff(mp(:, 4:6))), 2);
fwd = [0; fwd];                       % first frame has no predecessor
thresh = 0.5;                         % censoring threshold (mm); 0.2-0.5 typical
bad = find(fwd > thresh);
fprintf('Mean FWD = %.3f mm; %d frames exceed %.1f mm\n', ...
    mean(fwd), length(bad), thresh);

spikes = zeros(n, length(bad));       % one indicator column per bad frame
for i = 1:length(bad), spikes(bad(i), i) = 1; end

X2 = [mp spikes intercept_model(n)];
r2 = y - X2 * (X2 \ y);
fprintf('R^2 motion only: %.3f;  motion + spikes: %.3f\n', ...
    R2, 1 - var(r2) / var(y));

create_figure('FWD');
plot(fwd, 'k'); hold on;
plot_horizontal_line(thresh);
plot(bad, fwd(bad), 'rv');
xlabel('Frame'); ylabel('FWD (mm)'); title('Framewise displacement');

% Interpretation: realignment parameters soak up most motion-related
% variance; spike regressors handle the jerk. Prevention plus conservative
% censoring beats aggressive exclusion -- motion correlates with age, BMI,
% and clinical status, so excluding many participants biases the sample.

%% 2. Slice-timing offsets and interpolation
% With TR = 2 s and 20 slices acquired interleaved (odd slices first), two
% adjacent slices can be sampled almost a second apart. We build a "true"
% continuous BOLD response with SPM's canonical HRF, sample it at each
% slice's acquisition times, then correct by interpolating to a reference
% slice (the one acquired mid-TR).

TR = 2; n_vol = 100; n_slices = 20;   % TR (s); volumes; slices per volume
dt = 0.05;                            % high-resolution time step (s)
t_hi = (0:dt:(n_vol * TR - dt))';

hrf = spm_hrf(dt);                    % canonical HRF sampled at dt
onsets = 6:24:(n_vol * TR - 20);      % one event onset every 24 s
stim = zeros(size(t_hi));
stim(round(onsets / dt) + 1) = 1;
bold = conv(stim, hrf); bold = bold(1:length(t_hi));

order = [1:2:n_slices 2:2:n_slices];  % interleaved acquisition order
slice_offset = zeros(1, n_slices);
slice_offset(order) = (0:n_slices - 1) * (TR / n_slices);

vol_times = (0:n_vol - 1)' * TR;
picks = [1 11 20];                    % bottom, middle, top slice
labels = {'slice 1 (first)', 'slice 11 (middle)', 'slice 20 (last)'};
sampled = zeros(n_vol, 3);
for i = 1:3
    sampled(:, i) = interp1(t_hi, bold, vol_times + slice_offset(picks(i)));
    fprintf('%s: acquired %.2f s into each TR\n', ...
        labels{i}, slice_offset(picks(i)));
end

ref = 11;                             % reference slice (acquired mid-TR)
ref_offset = slice_offset(ref);
corrected = zeros(n_vol, 3);
for i = 1:3
    corrected(:, i) = interp1(vol_times + slice_offset(picks(i)), ...
        sampled(:, i), vol_times + ref_offset, 'linear', 'extrap');
end

create_figure('slice timing', 2, 1);
subplot(2, 1, 1);
plot(t_hi, bold, 'k-'); hold on; plot(vol_times, sampled, 'o-');
xlim([0 60]); legend([{'true BOLD'} labels]);
title('Uncorrected: same response, different slice times');
subplot(2, 1, 2);
plot(t_hi, bold, 'k-'); hold on; plot(vol_times, corrected, 'o-');
xlim([0 60]); xlabel('Time (s)');
title('After slice-timing correction to the middle slice');

% Interpretation: the last slice is sampled ~1.9 s later than the first --
% enough to bias amplitude estimates if regressors assume reference timing.
% With TR <= 1 s, many pipelines skip this step and shift task regressors
% per slice or use flexible basis sets instead.

%% 3. Smoothing kernel tradeoffs: the matched filter
% Smoothing convolves the image with a Gaussian kernel described by its
% FWHM (= 2.355 * sigma). Matched-filter theory: detection is optimal when
% the kernel matches the signal's spatial extent. We simulate a 2-D slice
% (2-mm voxels) with a small (~7 mm) and a large (~28 mm) activation blob
% plus noise, and measure peak z at each blob center vs. kernel FWHM.

vox_mm = 2; sz = 100;                 % voxel size (mm); grid size (voxels)
[xx, yy] = meshgrid(1:sz, 1:sz);
gauss_blob = @(cx, cy, sig) exp(-((xx - cx).^2 + (yy - cy).^2) / (2 * sig^2));
truth = gauss_blob(28, 50, 1.5) + gauss_blob(70, 50, 6);   % sigma in voxels

fwhms_mm = 0:2:16;                    % smoothing kernels to test (mm)
n_sims = 30;                          % noise realizations per kernel
z_small = zeros(size(fwhms_mm)); z_large = zeros(size(fwhms_mm));
for k = 1:length(fwhms_mm)
    sig_vox = fwhms_mm(k) / 2.355 / vox_mm;
    ps = zeros(n_sims, 1); pl = ps; sds = ps;
    for s = 1:n_sims
        img = truth + randn(sz);
        if sig_vox > 0
            kw = ceil(4 * sig_vox);
            [kx, ky] = meshgrid(-kw:kw, -kw:kw);
            kern = exp(-(kx.^2 + ky.^2) / (2 * sig_vox^2));
            kern = kern / sum(kern(:));
            sm = conv2(img, kern, 'same');
        else
            sm = img;
        end
        ps(s) = sm(50, 28); pl(s) = sm(50, 70);
        bg = sm(5:20, 5:20); sds(s) = std(bg(:));
    end
    z_small(k) = mean(ps) / mean(sds);
    z_large(k) = mean(pl) / mean(sds);
end

create_figure('matched filter');
plot(fwhms_mm, z_small, 'o-'); hold on; plot(fwhms_mm, z_large, 's-');
plot_vertical_line(6); plot_vertical_line(8);
xlabel('Smoothing kernel FWHM (mm)'); ylabel('Peak z at blob center');
legend({'small blob (~7 mm)', 'large blob (~28 mm)'});
title('Detection vs. smoothing: the matched-filter tradeoff');

[~, i1] = max(z_small); [~, i2] = max(z_large);
fprintf('Best FWHM: small blob %d mm, large blob %d mm\n', ...
    fwhms_mm(i1), fwhms_mm(i2));

% Interpretation: each activation is best detected with a kernel near its
% own scale. The conventional 6-8 mm kernel is a compromise -- good for
% cortex-scale signals, costly for fine structure (brainstem, MVPA
% patterns). Smoothing also satisfies Gaussian random field theory
% assumptions for multiple-comparisons correction.

%% 4. High-pass filter design and drift removal
% High-pass filtering removes slow drift by projecting out low-frequency
% discrete-cosine (DCT) regressors KH via the residual-forming matrix
% S = I - KH*pinv(KH). Filtering IS nuisance regression: do not filter and
% then separately regress out covariates -- combine DCT columns, motion
% parameters, spikes, and intercepts into ONE matrix and remove them in a
% single step, applying the same operation to the task regressors.
% Adapted from CANlab_help_examples: linear_filtering_a_timeseries.m

TR = 2; n = 300; hpf = 128;           % TR (s); n = volumes; hpf = high-pass cutoff (s)

[S, KL, KH] = use_spm_filter(TR, n, 'none', 'specify', hpf);
fprintf('%d DCT regressors span frequencies below 1/%d = %.4f Hz\n', ...
    size(KH, 2), hpf, 1 / hpf);

% Combine filter + nuisance covariates into one residual-forming matrix
get_hat = @(X) X * pinv(X);
get_S   = @(H) eye(size(H)) - H;

N = randn(n, 6);                      % stand-in nuisance covs (e.g. motion)
I3 = intercept_model([n/3 n/3 n/3]);  % run intercepts (3 runs)
S_all = get_S(get_hat([KH N I3]));

% Simulated run: 30-s alternating blocks + slow drift + AR(2) noise
task = repmat([ones(15, 1); zeros(15, 1)], 10, 1);      % in TRs
tsec = (0:n - 1)' * TR;
drift = 3 * sin(2 * pi * tsec / 400) + 0.004 * (tsec - mean(tsec));
y_obs = task + drift + noise_arp(n, [.7 .3]) + I3 * randn(3, 1);

y_filt = S_all * y_obs;
fprintf('corr(task, observed) = %.2f;  corr(task, filtered) = %.2f\n', ...
    corr(task, y_obs), corr(task, y_filt));

create_figure('filtering', 2, 1);
subplot(2, 1, 1); plot(tsec, y_obs, 'k'); hold on; plot(tsec, drift, 'r');
title('Observed = 30-s task blocks + drift (red) + AR noise');
subplot(2, 1, 2); plot(tsec, y_filt, 'k'); hold on;
plot(tsec, 2 * (task - mean(task)), 'g');
xlabel('Time (s)');
title('Filtered (128-s cutoff): task blocks (green) survive');

%% 4a. When the filter eats your task
% The filter removes EVERYTHING below the cutoff, including task signal.
% Measure the fraction of task variance removed for alternating block
% designs of increasing length under the 128-s default.

block_lengths = 10:2:90;              % seconds per block
frac_removed = zeros(size(block_lengths));
for k = 1:length(block_lengths)
    bl_tr = round(block_lengths(k) / TR);
    x = repmat([ones(bl_tr, 1); zeros(bl_tr, 1)], ceil(n / (2 * bl_tr)), 1);
    x = x(1:n);
    xc = x - mean(x);
    frac_removed(k) = 1 - var(S * x) / var(xc);
end

create_figure('task overlap');
plot(block_lengths, 100 * frac_removed, 'o-'); hold on;
plot_vertical_line(64);
xlabel('Block length (s), alternating design');
ylabel('% of task variance removed by 128-s filter');
title('A default filter can silently delete your effect');

fprintf('32-s blocks: %.0f%% of task variance removed\n', ...
    100 * frac_removed(block_lengths == 32));
fprintf('64-s blocks: %.0f%% of task variance removed\n', ...
    100 * frac_removed(block_lengths == 64));

% Interpretation: short blocks sail through, but as the design period
% approaches the 128-s cutoff the filter removes most task variance -- with
% 64-s blocks (128-s period) the "correction" deletes the experiment.
% Always check the overlap between design frequencies and filter before
% analyzing, and apply the same filter to task regressors and covariates.

%% Wrap-up
% - Realignment parameters capture head motion; FWD summarizes them for
%   censoring; spike regression censors without destroying temporal
%   structure.
% - Slice-timing offsets shift sampled responses by up to nearly one TR;
%   interpolation or regressor-shifting realigns them.
% - Smoothing is a matched filter: the best kernel depends on the signal
%   extent, so 6-8 mm is a compromise.
% - High-pass filtering is nuisance regression with DCT cosines; design
%   your cutoff around your task frequencies, not the other way around.
%
% Continue to the Chapter 18 lab, where these cleaned-up time series meet
% the General Linear Model.
