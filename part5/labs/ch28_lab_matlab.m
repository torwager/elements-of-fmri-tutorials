%% Chapter 28 Lab: Resting-State and Ecological Designs (MATLAB)
% This lab accompanies Chapter 28, "Resting-State and Ecological Designs".
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part5/ch28-resting-state-and-ecological-designs
% You will simulate resting-state multi-ROI time series and compute
% functional connectivity, simulate a shared naturalistic ("movie")
% stimulus and compute inter-subject correlation (ISC), and see how motion
% spikes and vigilance drift inflate connectivity estimates.
%
% Requirements: base MATLAB (Statistics and Machine Learning Toolbox for
% corr/ttest2). All data are simulated; comments point to the CanlabCore
% functions you would use with real data.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples; see canlab_help_prep_for_connectivity).
%
% Runtime: under a minute.

%% 1. Simulate resting-state ROI time series with known network structure
% We simulate a 10-minute resting scan (300 volumes, TR = 2 s) for 8 ROIs
% organized into two networks. Each ROI mixes a slow network-level signal
% (weight w) with slow ROI-specific fluctuations (weight 1-w). Resting
% BOLD is dominated by low frequencies (< ~0.1 Hz), so we low-pass white
% noise with a smoothing kernel.

rng(28);                              % seed for reproducibility
TR    = 2;                            % repetition time (s)
n_t   = 300;                          % 10-minute "scan"
n_roi = 8;                            % regions, 4 per network
network   = [1 1 1 1 2 2 2 2];        % network membership for each ROI
roi_names = arrayfun(@(i) sprintf('N%d-roi%d', network(i), ...
    mod(i - 1, 4) + 1), 1:n_roi, 'UniformOutput', false);

kern   = exp(-0.5 * ((-4:4)' ./ 2) .^ 2);       % Gaussian smoothing kernel
kern   = kern ./ sum(kern);
smoothz = @(z) zscore(conv(z, kern, 'same'));   % slow, unit-variance noise

% One slow shared signal per network
net_sig = [smoothz(randn(n_t, 1)) smoothz(randn(n_t, 1))];

w = 0.7;                              % network signal weight
Y_rest = zeros(n_t, n_roi);
for i = 1:n_roi
    Y_rest(:, i) = w * net_sig(:, network(i)) + ...
                   (1 - w) * smoothz(randn(n_t, 1));
end
Y_rest = zscore(Y_rest);

t = (0:n_t - 1) * TR;
figure('Name', 'Resting ROI time series');
hold on;
colors = [repmat([0 0.35 0.75], 4, 1); repmat([0.8 0.2 0.2], 4, 1)];
for i = 1:n_roi
    plot(t, Y_rest(:, i) + 3 * (i - 1), 'Color', colors(i, :));
end
xlabel('Time (s)'); yticks(3 * (0:n_roi - 1)); yticklabels(roi_names);
title('Simulated resting-state ROIs (blue = network 1, red = network 2)');

% With real data you would instead load and denoise a 4-D image, then
% average within atlas regions:
%   obj = fmri_data(fname);
%   obj_denoised = obj.denoise_timeseries_pipeline(TR, 128, mvmtfname);
%   r = extract_roi_averages(obj_denoised, load_atlas('canlab2018_2mm'));

%% 2. Functional connectivity: the ROI x ROI correlation matrix
% Functional connectivity asks which regions' spontaneous fluctuations
% rise and fall together. The standard summary is the pairwise correlation
% matrix, which recovers the two-network block structure we built in.

FC_rest = corr(Y_rest);

same_net = network' == network;                   % same-network pairs
offdiag  = ~eye(n_roi);
mean_offdiag = @(FC, mask) mean(FC(mask & offdiag));

figure('Name', 'Resting FC');
imagesc(FC_rest, [-1 1]); axis square; colorbar;
colormap(flipud(interp1([0 0.5 1], [0.7 0.1 0.1; 1 1 1; 0.1 0.1 0.7], ...
    linspace(0, 1, 64))));
xticks(1:n_roi); xticklabels(roi_names); xtickangle(45);
yticks(1:n_roi); yticklabels(roi_names);
title('Resting-state functional connectivity');

fprintf('Mean within-network r:  %3.2f\n', mean_offdiag(FC_rest, same_net));
fprintf('Mean between-network r: %3.2f\n', mean_offdiag(FC_rest, ~same_net));

% Strong within-network and much weaker between-network r. Clustering
% and ICA (Chapter 30) find networks by discovering this structure in
% matrices with thousands of regions instead of 8.

%% 3. Naturalistic designs: inter-subject correlation (ISC)
% Now the key move of naturalistic imaging: give every subject the SAME
% rich stimulus. We simulate one stimulus-driven ROI in 10 subjects:
%   Movie: shared stimulus-driven signal + idiosyncratic fluctuations
%   Rest:  idiosyncratic fluctuations only
% Leave-one-out ISC correlates each subject with the mean of the others.
% No event onsets and no HRF model -- the other brains ARE the model.

n_sub = 10;                           % number of subjects
movie_sig = smoothz(randn(n_t, 1));   % shared stimulus-driven time course
a = 0.6;                              % stimulus-driven signal fraction

Y_movie = zeros(n_t, n_sub); Y_solo = zeros(n_t, n_sub);
for s = 1:n_sub
    Y_movie(:, s) = a * movie_sig + (1 - a) * smoothz(randn(n_t, 1));
    Y_solo(:, s)  = smoothz(randn(n_t, 1));      % nothing shared
end

isc_loo = @(Y) arrayfun(@(s) corr(Y(:, s), ...
    mean(Y(:, setdiff(1:n_sub, s)), 2)), 1:n_sub);

isc_movie = isc_loo(Y_movie);
isc_rest  = isc_loo(Y_solo);

figure('Name', 'ISC: movie vs rest');
bar([mean(isc_movie) mean(isc_rest)], 'FaceAlpha', 0.6); hold on;
plot(1 + 0.1 * randn(1, n_sub), isc_movie, 'ko');
plot(2 + 0.1 * randn(1, n_sub), isc_rest, 'ko');
xticklabels({'Movie', 'Rest'}); ylabel('Leave-one-out ISC (r)');
title('Inter-subject correlation');

fprintf('Mean ISC: movie = %3.2f, rest = %3.2f\n', ...
    mean(isc_movie), mean(isc_rest));

% ISC is high wherever -- and only wherever -- activity is driven by the
% shared stimulus. It misses idiosyncratic responses, even strongly
% stimulus-driven ones, because they are not shared across people.

%% 4. Confound 1: motion spikes inflate connectivity
% Head motion changes the signal in many voxels at once. A shared artifact
% adds a common component to every ROI pair, inflating correlations
% everywhere -- including between networks, where true r is ~0. We add 8
% spikes (< 3% of volumes) to the resting data and re-estimate FC, then
% "scrub" (drop) the spike volumes, as flagged in practice by framewise
% displacement and outlier detection.

n_spikes  = 8;                        % motion spikes (< 3% of volumes)
spike_idx = randperm(n_t, n_spikes);
Y_motion  = Y_rest;
for k = 1:n_spikes
    amp = 4 + randn;                  % large global deviation
    Y_motion(spike_idx(k), :) = Y_motion(spike_idx(k), :) + ...
        amp * (0.7 + 0.6 * rand(1, n_roi));
end

FC_motion = corr(Y_motion);
keep = setdiff(1:n_t, spike_idx);
FC_scrubbed = corr(Y_motion(keep, :));

fprintf('\n%-14s %-18s %-18s\n', 'Data', 'Within-network r', 'Between-network r');
FCs = {FC_rest, FC_motion, FC_scrubbed};
labels = {'Clean', 'With spikes', 'Scrubbed'};
for k = 1:3
    fprintf('%-14s %-18.2f %-18.2f\n', labels{k}, ...
        mean_offdiag(FCs{k}, same_net), mean_offdiag(FCs{k}, ~same_net));
end

% A handful of bad volumes noticeably inflates between-network r. This is
% why connectivity pipelines build spike/scrubbing regressors and motion
% covariates -- see fmri_data.denoise_timeseries_pipeline and
% canlab_connectivity_preproc, which implement these steps for real data.
% If one group moves more than another, motion creates spurious group
% differences in "connectivity".

%% 5. Confound 2: vigilance drift masquerades as a group difference
% Roughly half of resting participants are asleep after 10 minutes, and
% drowsiness produces large, slow, widespread fluctuations. We simulate
% two groups with IDENTICAL network structure; the "drowsy" group also has
% a global arousal signal added to every ROI.

n_per_group = 12;                     % subjects per group
sim_subject = @(global_amp) sim_rest_subject(global_amp, n_t, n_roi, ...
    network, w, smoothz);

fc_alert  = zeros(n_per_group, 1);
fc_drowsy = zeros(n_per_group, 1);
for s = 1:n_per_group
    fc_alert(s)  = mean_fc(sim_subject(0),   offdiag);
    fc_drowsy(s) = mean_fc(sim_subject(0.8), offdiag);
end

[~, p, ~, st] = ttest2(fc_drowsy, fc_alert);
fprintf('\nMean FC, alert group:  %3.2f\n', mean(fc_alert));
fprintf('Mean FC, drowsy group: %3.2f\n', mean(fc_drowsy));
fprintf('Spurious group difference: t(%d) = %3.1f, p = %3.2g\n', ...
    st.df, st.tstat, p);

% Identical neural networks, yet the drowsy group shows higher average
% connectivity -- an artifact of arousal. Practical lessons: monitor
% wakefulness (eye tracking, self-report), check group differences in
% sleepiness and motion, match acquisition conditions -- or use an
% engaging naturalistic paradigm, which constrains mental state, keeps
% participants awake, and yields more reliable, behavior-predictive
% connectivity than rest.

%% Explore on your own
% 1. Scan length and reliability: re-run Section 1 with n_t for 3, 6, 9,
%    13, and 27 minutes and correlate each estimated FC with a very long
%    "ground truth" run (n_t = 5000). Gains should diminish after ~9-13
%    minutes of data.
% 2. ISC and response timing: shift movie_sig by a random 0-3 volume lag
%    per subject before mixing. ISC drops even though every subject is
%    stimulus-driven. What does this imply for regions with variable
%    hemodynamics?
% 3. Scrubbing tradeoff: raise n_spikes to 60. How much data can you drop
%    before the scrubbed FC estimate itself becomes unstable?
% 4. Global signal regression: regress the mean ROI signal out of each
%    drowsy subject's data before computing FC. How much of the spurious
%    group difference disappears?

%% Local functions

function Y = sim_rest_subject(global_amp, n_t, n_roi, network, w, smoothz)
% One subject's 8-ROI resting data; global_amp scales an arousal-related
% signal added to ALL ROIs (vigilance fluctuation).
nets = [smoothz(randn(n_t, 1)) smoothz(randn(n_t, 1))];
arousal = smoothz(randn(n_t, 1));
Y = zeros(n_t, n_roi);
for i = 1:n_roi
    Y(:, i) = w * nets(:, network(i)) + (1 - w) * smoothz(randn(n_t, 1)) ...
        + global_amp * arousal;
end
Y = zscore(Y);
end

function m = mean_fc(Y, offdiag)
% Mean off-diagonal correlation across all ROI pairs.
FC = corr(Y);
m = mean(FC(offdiag));
end
