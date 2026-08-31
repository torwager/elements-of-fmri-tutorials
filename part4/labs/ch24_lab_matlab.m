%% Chapter 24 Lab: A Mini-Multiverse of Analysis Pipelines (MATLAB)
% This lab accompanies Chapter 24, "Analysis Pipelines: Variations and
% Variability". You will analyze ONE simulated fMRI study through a
% factorial grid of 36 defensible pipeline variants -- high-pass filter
% cutoff x spatial smoothing x outlier handling x autocorrelation
% correction -- and visualize how much the group result "vibrates" across
% pipelines, ending with a specification-curve plot.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path (used for
% onsets2fmridesign; everything else is base MATLAB + Statistics Toolbox).
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples).
%
% Runtime: under a minute. All data are simulated.

%% 1. Simulate one study
% We simulate 24 subjects, each with a 160-volume run (TR = 2 s) and a
% strip of 40 voxels (3 mm each). A blocked task regressor (4-s events
% every 24 s) is convolved with a canonical HRF using onsets2fmridesign.
% The TRUE effect is a Gaussian bump centered on voxel 21, with subject
% amplitudes drawn from a population distribution: real, but moderate.
% Nuisance processes: slow drift, AR(1) noise, and -- for 8 "high-motion"
% subjects -- large global intensity spikes.

TR = 2; n_t = 160; n_v = 40; vox_mm = 3; n_sub = 24;
run_len = n_t * TR;

ons = {};
ons{1} = (12:24:run_len - 24)';            % event onsets (s)
X0 = onsets2fmridesign(ons, TR, run_len);  % CANlab: convolve with SPM HRF
task = X0(:, 1) ./ max(X0(:, 1));          % unit-amplitude task regressor

peak_vox = 21;                                       % true peak voxel
profile = exp(-((1:n_v) - peak_vox).^2 / (2 * 4^2)); % spatial effect profile

rng(24);
sub_amp = 0.25 + 0.15 .* randn(n_sub, 1);  % true subject effect amplitudes
high_motion = false(n_sub, 1);
high_motion(randperm(n_sub, 8)) = true;

phi_true = 0.45;                           % AR(1) autocorrelation
tt = (0:n_t - 1)' ./ (n_t - 1);

data = zeros(n_sub, n_t, n_v);
for s = 1:n_sub
    % AR(1) noise, independent across voxels
    noise = filter(1, [1 -phi_true], randn(n_t, n_v));

    % slow scanner drift: two subject-specific low-frequency cosines
    drift = (1.5 + 3.5 * rand) * cos(2 * pi * 0.9 * tt + 2 * pi * rand) + ...
            (1.5 + 3.5 * rand) * cos(2 * pi * 1.7 * tt + 2 * pi * rand);
    drift = drift * (1 + 0.15 .* randn(1, n_v));

    Y = 100 + sub_amp(s) .* task * profile + drift + noise;

    if high_motion(s)   % motion-like global spikes on 8 random frames
        frames = 5 + randperm(n_t - 10, 8);
        Y(frames, :) = Y(frames, :) + ...
            (8 + 8 .* rand(8, 1)) .* sign(randn(8, 1)) * ones(1, n_v);
    end
    data(s, :, :) = Y;
end

fprintf('True population effect: mean %.2f, sd %.2f\n', ...
    mean(sub_amp), std(sub_amp));

%% 2. Look at the raw data
% One clean subject and one high-motion subject. Note that the slow drift
% is as large as the task effect, and the spikes are much larger.

clean = find(~high_motion, 1); spiky = find(high_motion, 1);
create_figure('example subjects', 2, 1);
for k = 1:2
    s = [clean spiky] * [2 - k; k - 1];  % clean then spiky
    subplot(2, 1, k);
    plot((0:n_t-1) * TR, mean(data(s, :, peak_vox-2:peak_vox+2), 3), 'k-'); hold on;
    plot((0:n_t-1) * TR, 100 + 3 * task, 'r-');
    ylabel('signal'); title(sprintf('Subject %d', s));
    legend({'signal near peak voxel', 'task regressor (scaled)'});
end
xlabel('time (s)');

%% 3. The pipeline grid
% Four analytic choices, each with defensible options:
%   High-pass cutoff:     none / 128 s / 64 s   (DCT cosine regressors)
%   Smoothing FWHM:       0 / 4 / 8 mm          (Gaussian, across voxels)
%   Outlier handling:     none / spike regressors (robust global-signal z)
%   Autocorrelation:      OLS / AR(1) prewhitening
% 3 x 3 x 2 x 2 = 36 pipelines. Each runs the same summary-statistics
% analysis: first-level GLM per subject, the task beta at the a-priori
% peak voxel carried to a group one-sample t-test.

hp_opts = [Inf 128 64];       % Inf = no high-pass filtering
fwhm_opts = [0 4 8];
despike_opts = [false true];
ar1_opts = [false true];

results = table('Size', [0 9], ...
    'VariableTypes', {'double','double','logical','logical','double', ...
                      'double','double','double','logical'}, ...
    'VariableNames', {'hp','fwhm','despike','ar1','mean_beta', ...
                      'sem','t','p','significant'});

for hp = hp_opts
    for fwhm = fwhm_opts
        for despike = despike_opts
            for ar1 = ar1_opts
                betas = zeros(n_sub, 1);
                for s = 1:n_sub
                    Y = squeeze(data(s, :, :));
                    betas(s) = first_level_beta(Y, task, TR, hp, fwhm, ...
                        vox_mm, despike, ar1, peak_vox);
                end
                [~, p, ~, st] = ttest(betas);
                results(end+1, :) = {hp, fwhm, despike, ar1, ...
                    mean(betas), std(betas) / sqrt(n_sub), ...
                    st.tstat, p, p < 0.05}; %#ok<SAGROW>
            end
        end
    end
end

fprintf('%d variants; %d significant at p < .05 (%.0f%%)\n', ...
    height(results), sum(results.significant), ...
    100 * mean(results.significant));
fprintf('group t range: %.2f to %.2f\n', min(results.t), max(results.t));
disp(head(results, 8))

%% 4. The vibration of effects
% The distribution of the group t-statistic across all 36 pipelines --
% same data, same hypothesis, different defensible analysis choices.

t_crit = tinv(0.975, n_sub - 1);

create_figure('vibration of effects');
histogram(results.t, 14, 'FaceColor', [0.27 0.51 0.71]);
hold on; plot([t_crit t_crit], ylim, 'r--', 'LineWidth', 1.5);
xlabel('group t-statistic (same data, different pipelines)');
ylabel('number of pipeline variants');
title('Vibration of effects across 36 analysis pipelines');
legend({'pipelines', sprintf('p = .05 threshold (t = %.2f)', t_crit)});

%% 5. Specification curve
% Top: the 36 group estimates ordered from smallest to largest, with 95%
% CIs, colored by significance. Bottom: which option each variant used.
% Reading a bottom row horizontally shows whether an option systematically
% pushes the estimate up or down.

[~, order] = sort(results.mean_beta);
res = results(order, :);
xx = 1:height(res);
sig = res.significant;

create_figure('specification curve', 2, 1);
subplot(2, 1, 1);
errorbar(xx, res.mean_beta, 1.96 * res.sem, 'LineStyle', 'none', ...
    'Color', [0.8 0.8 0.8]); hold on;
scatter(xx(sig), res.mean_beta(sig), 30, 'r', 'filled');
scatter(xx(~sig), res.mean_beta(~sig), 30, [0.5 0.5 0.5], 'filled');
plot(xlim, [0 0], 'k-');
plot(xlim, mean(sub_amp) * [1 1], 'g:', 'LineWidth', 1.5);
ylabel('group mean beta, peak voxel (95% CI)');
title('Specification curve: one dataset, 36 pipelines');
legend({'95% CI', 'p < .05', 'p >= .05', 'zero', 'true effect'});

spec_labels = {'HP: none', 'HP: 128 s', 'HP: 64 s', 'FWHM: 0 mm', ...
    'FWHM: 4 mm', 'FWHM: 8 mm', 'outliers: none', 'outliers: spike regs', ...
    'autocorr: OLS', 'autocorr: AR(1)'};
spec_masks = [res.hp == Inf, res.hp == 128, res.hp == 64, ...
    res.fwhm == 0, res.fwhm == 4, res.fwhm == 8, ...
    ~res.despike, res.despike, ~res.ar1, res.ar1];

subplot(2, 1, 2); hold on;
n_rows = numel(spec_labels);
for i = 1:n_rows
    m = spec_masks(:, i);
    scatter(xx(m), repmat(n_rows - i + 1, sum(m), 1), 14, 'ks', 'filled');
end
set(gca, 'YTick', 1:n_rows, 'YTickLabel', fliplr(spec_labels));
xlabel('pipeline variants, ordered by estimated effect');
ylim([0.3 n_rows + 0.7]);

%% 6. Which choices matter?
% Average the group t within each level of each factor. Typically: drift
% handling and outlier handling matter a lot, smoothing matters, and AR(1)
% prewhitening barely moves the GROUP result -- in a summary-statistics
% analysis it mainly affects first-level standard errors, not the betas
% carried forward.

disp(groupsummary(results, 'hp', 'mean', {'t', 'significant'}))
disp(groupsummary(results, 'fwhm', 'mean', {'t', 'significant'}))
disp(groupsummary(results, 'despike', 'mean', {'t', 'significant'}))
disp(groupsummary(results, 'ar1', 'mean', {'t', 'significant'}))

%% 7. Wrap-up
% One dataset, one real effect -- and dozens of defensible analyses that
% disagree about "significance". Carp (2012) counted 34,560 plausible
% pipelines; NARPS (Botvinik-Nezer et al., 2020) found only 37% of 70
% teams declared a well-established effect significant on the same data.
% Protect yourself: (1) fix the pipeline BEFORE seeing results
% (preregister, or use a standardized pipeline such as fMRIPrep on
% BIDS-formatted data); (2) if you explore variants, report the multiverse,
% not the best branch; (3) share code, data, and unthresholded maps.
%
% Exercises:
% (a) Set sub_amp = 0.6 + 0.15*randn(n_sub,1) (strong effect). What
%     happens to the fraction of significant variants, and why?
% (b) Add hp = 32 to hp_opts -- aggressive filtering that begins to remove
%     task frequencies. What happens?
% (c) Test voxel 4 instead of the peak voxel (no true effect) and see how
%     often you can find a "significant" pipeline if you pick the best one.

%% Local functions

function b = first_level_beta(Y, task, TR, hp, fwhm, vox_mm, despike, ar1, test_vox)
% One subject through one pipeline variant -> task beta at test_vox.
n_t = size(Y, 1);

if fwhm > 0   % spatial smoothing across the voxel dimension
    sigma_vox = fwhm / 2.355 / vox_mm;
    w = max(3, round(5 * sigma_vox) + 1);
    Y = smoothdata(Y, 2, 'gaussian', w);
end

X = [task ones(n_t, 1)];
if isfinite(hp)   % SPM-style DCT drift regressors for periods > hp
    k = 1:floor(2 * n_t * TR / hp);
    X = [X cos(pi / n_t .* ((0:n_t-1)' + 0.5) * k)];
end
if despike        % one spike regressor per flagged frame
    frames = find_spike_frames(Y);
    S = zeros(n_t, numel(frames));
    for i = 1:numel(frames), S(frames(i), i) = 1; end
    X = [X S];
end

B = X \ Y;        % OLS fit, all voxels at once
if ar1            % AR(1) prewhitening (Cochrane-Orcutt style)
    R = Y - X * B;
    phi = mean(sum(R(2:end, :) .* R(1:end-1, :)) ./ sum(R(1:end-1, :).^2));
    B = (X(2:end, :) - phi * X(1:end-1, :)) \ (Y(2:end, :) - phi * Y(1:end-1, :));
end

b = B(1, test_vox);
end

function frames = find_spike_frames(Y)
% Flag frames with abrupt global intensity jumps (robust z of the
% frame-to-frame difference of the global mean signal).
g = mean(Y, 2);
dg = abs([0; diff(g)]);
mad_dg = median(abs(dg - median(dg))) * 1.4826 + 1e-12;
frames = find((dg - median(dg)) ./ mad_dg > 5);
end
