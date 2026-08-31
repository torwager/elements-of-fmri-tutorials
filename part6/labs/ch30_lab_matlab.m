%% Chapter 30 Lab: Introduction to Brain Connectivity (MATLAB)
% This lab accompanies Chapter 30, "Introduction to Brain Connectivity".
% Using simulated data with known ground truth, you will build the core
% objects of functional connectivity analysis:
%
% 1. Multi-ROI time series with community (network) structure
% 2. A parcellated functional connectivity (correlation) matrix
% 3. A seed correlation map over a simulated voxel grid
% 4. Full vs. partial correlation in a chain network
% 5. A shared nuisance signal that inflates FC, and the regression fix
% 6. Split-half reliability of FC estimates as a function of scan length
%
% Sections 1-6 use base MATLAB + Statistics and Machine Learning Toolbox
% (for corr/partialcorr). Section 7 sketches the real-data workflow with
% CANlab tools (CanlabCore + SPM12 on your path) and is left commented.
% Code adapted from CANlab tutorials (github.com/canlab, canlab.github.io).
%
% Runtime: under a minute. All data are simulated.

%% 1. Simulate multi-ROI time series with community structure
% We create 12 ROIs organized into 3 networks of 4 ROIs each. Each network
% has one slow latent signal (think of it as coherent spontaneous activity),
% and each ROI's time series is that latent signal plus ROI-specific noise.
% This is the generative idea behind "functional networks": regions in the
% same network share variance over time.

rng(30);                     % reproducible
n_t     = 240;               % volumes (e.g., 8 min at TR = 2 s)
n_net   = 3;                 % number of networks
roi_per = 4;                 % ROIs per network
n_roi   = n_net * roi_per;

latent = randn(n_t, n_net);                       % one latent signal per network
for k = 1:6
    latent = conv2(latent, ones(5, 1) / 5, 'same');   % smooth: slow fluctuations
end
latent = zscore(latent);

net_assign = kron((1:n_net)', ones(roi_per, 1));  % network label for each ROI
Y = 0.8 * kron(latent, ones(1, roi_per)) + randn(n_t, n_roi);

TR = 2; t = (0:n_t - 1)' * TR;

create_figure_or_figure('ROI time series');
plot(t, Y(:, 1) + 6, 'b-'); hold on;              % two ROIs in network 1
plot(t, Y(:, 2) + 3, 'b--');
plot(t, Y(:, 5), 'r-');                           % one ROI in network 2
xlabel('Time (s)'); ylabel('Signal (offset for display)');
legend({'ROI 1 (net 1)', 'ROI 2 (net 1)', 'ROI 5 (net 2)'});
title('ROIs in the same network share slow fluctuations');

%% 2. The parcellated functional connectivity matrix
% Correlating every ROI with every other ROI yields the node-by-node FC
% matrix -- a miniature "functional connectome". Because ROIs are ordered
% by network, the community structure appears as bright diagonal blocks.

R = corr(Y);

create_figure_or_figure('FC matrix');
imagesc(R, [-1 1]); axis square; colorbar;
title('Functional connectivity (Pearson r)');
xlabel('ROI'); ylabel('ROI');

off = ~eye(n_roi);
within  = net_assign == net_assign';              % same-network pairs
fprintf('Mean within-network r:  %.3f\n', mean(R(within & off)));
fprintf('Mean between-network r: %.3f\n', mean(R(~within)));

%% 3. A seed correlation map
% Seed-based connectivity correlates one region's time series with every
% voxel. We simulate a 48 x 48 "slice" of voxels: three circular blobs are
% assigned to the three networks, and background voxels contain only noise.
% The seed is the latent signal of network 1 -- as if extracted by averaging
% a seed ROI. Voxels in network-1 blobs light up in the map.

nx = 48;
[cx, cy] = meshgrid(1:nx, 1:nx);
centers = [12 12; 34 14; 22 36];                  % blob centers (one per network)
vox_net = zeros(nx);                              % 0 = background
for k = 1:n_net
    d2 = (cx - centers(k, 1)).^2 + (cy - centers(k, 2)).^2;
    vox_net(d2 <= 7^2) = k;
end

vox_ts = randn(n_t, nx * nx);                     % voxel noise
for k = 1:n_net
    idx = find(vox_net(:) == k);
    vox_ts(:, idx) = vox_ts(:, idx) + 0.8 * latent(:, k);
end

seed = latent(:, 1);                              % seed time series (network 1)
seed_map = reshape(corr(seed, vox_ts), nx, nx);

create_figure_or_figure('seed map');
subplot(1, 2, 1); imagesc(vox_net); axis square;
title('Ground truth: voxel network labels'); colorbar;
subplot(1, 2, 2); imagesc(seed_map, [-1 1]); axis square;
title('Seed correlation map (seed = network 1)'); colorbar;

% Only the blob(s) sharing the seed's latent signal show high correlations.
% In a real analysis, per-subject maps like this enter a group-level test.

%% 4. Full vs. partial correlation: the third-variable problem
% Build a chain: A -> B -> C (C receives B's signal, B receives A's).
% A and C have NO direct connection, yet their full correlation is far from
% zero, because B links them. Partial correlation, controlling for B,
% removes the indirect path.

a = zscore(conv(randn(n_t, 1), ones(5, 1) / 5, 'same'));
b = zscore(0.9 * a + 0.5 * randn(n_t, 1));
c = zscore(0.9 * b + 0.5 * randn(n_t, 1));
ABC = [a b c];

R_full = corr(ABC);
R_part = partialcorr(ABC);        % Statistics and Machine Learning Toolbox

disp('Full correlation (A, B, C):');    disp(round(R_full, 2));
disp('Partial correlation (A, B, C):'); disp(round(R_part, 2));

% The A-C entry drops from a substantial value to near zero under partial
% correlation. With many nodes, partial correlations become unstable, and
% regularized inverse-covariance ("graphical lasso") estimators are used.

%% 5. A shared nuisance signal inflates FC -- and regression repairs it
% Motion, respiration, and cardiac signals add common variance to time
% series across the whole brain. We add one slow nuisance signal to every
% ROI and compare FC before contamination, after contamination, and after
% nuisance regression.

g = zscore(conv(randn(n_t, 1), ones(20, 1) / 20, 'same'));   % nuisance signal
Y_bad = Y + 1.2 * g * ones(1, n_roi);                        % hits every ROI

% Residualize each ROI's time series on [nuisance, intercept]
Xn = [g, ones(n_t, 1)];
Y_clean = Y_bad - Xn * (Xn \ Y_bad);

R_bad   = corr(Y_bad);
R_clean = corr(Y_clean);

create_figure_or_figure('nuisance');
subplot(1, 3, 1); imagesc(R, [-1 1]);       axis square; title('True FC');
subplot(1, 3, 2); imagesc(R_bad, [-1 1]);   axis square; title('+ shared nuisance');
subplot(1, 3, 3); imagesc(R_clean, [-1 1]); axis square; title('After regression');

fprintf('Mean between-network r: true %.3f | contaminated %.3f | cleaned %.3f\n', ...
    mean(R(~within)), mean(R_bad(~within)), mean(R_clean(~within)));

% The contaminated matrix is inflated everywhere -- including between
% networks that are truly unconnected. Nuisance regression restores
% estimates close to the truth. On real data the nuisance matrix contains
% motion parameters, ventricle/white-matter signals, and drift terms.

%% 6. Reliability: FC estimates stabilize with scan duration
% Correlations estimated from short scans are noisy. We simulate two long
% independent "sessions" from the same generative model (same true FC) and
% ask: how similar are the edge estimates across sessions as a function of
% the number of volumes used? This split-session correlation of vectorized
% FC matrices is a simple reliability index.

n_long = 1200;
make_session = @(seed_val) simulate_session(seed_val, n_long, n_net, roi_per);
Y1 = make_session(101);
Y2 = make_session(202);

lens = [60 120 240 480 1200];
rel = zeros(size(lens));
mask_ut = triu(true(n_roi), 1);                   % unique edges
for i = 1:numel(lens)
    v1 = corr(Y1(1:lens(i), :));  v1 = v1(mask_ut);
    v2 = corr(Y2(1:lens(i), :));  v2 = v2(mask_ut);
    rel(i) = corr(v1, v2);
end

create_figure_or_figure('reliability');
plot(lens * TR / 60, rel, 'o-', 'LineWidth', 2);
xlabel('Scan length (minutes)'); ylabel('Between-session edge similarity (r)');
title('FC reliability grows with scan duration'); ylim([0 1]); grid on;

% Short scans yield unstable connectomes; reliability climbs steeply with
% duration. This is one reason fingerprinting works so well with long
% high-quality scans (e.g., HCP), and why many groups now collect 20+ min
% of resting-state data per person.

%% 7. Real data with CANlab tools (recipe -- requires CanlabCore, SPM12, data)
% The pipeline below mirrors what you just did, on a real 4-D image.
% Adapted from the CANlab connectivity tutorials:
%   canlab.github.io -> "Prepare fMRI data for connectivity analyses"
%   canlab.github.io -> brainomics_connectivity_demo
%
% % Load one subject's preprocessed 4-D data into an fmri_data object
% dat = fmri_data('sub-01_preprocessed_bold.nii.gz');
%
% % Connectivity preprocessing: removes nuisance covariates (dat.covariates),
% % ventricle & white-matter signals, bandpass filters, and windsorizes.
% TR = 2;
% preprocessed_dat = canlab_connectivity_preproc(dat, 'vw', ...
%     'windsorize', 3, 'bpf', [.008 .25], TR);
%
% % Define a seed from an atlas and extract its average time series
% atlas_obj = load_atlas('canlab2018_2mm');
% m1 = select_atlas_subset(atlas_obj, {'Ctx_1_R'});     % right M1
% m1_timeseries = extract_data(m1, preprocessed_dat);
%
% % Seed map: regress every voxel on the seed time series
% preprocessed_dat.X = m1_timeseries;
% out = regress(preprocessed_dat);        % out.t(1) is the seed-map t image
%
% % Parcellated FC matrix: extract all atlas regions and correlate
% roi_ts = extract_data(atlas_obj, preprocessed_dat);
% R = corr(roi_ts);

%% Local functions

function Y = simulate_session(seed_val, n_t, n_net, roi_per)
% Simulate one session of networked ROI time series (same model as Sec. 1).
rng(seed_val);
latent = randn(n_t, n_net);
for k = 1:6
    latent = conv2(latent, ones(5, 1) / 5, 'same');
end
latent = zscore(latent);
Y = 0.8 * kron(latent, ones(1, roi_per)) + randn(n_t, n_net * roi_per);
end

function create_figure_or_figure(name)
% Use CANlab's create_figure if available; otherwise plain figure.
if exist('create_figure', 'file')
    create_figure(name);
else
    figure('Name', name, 'Color', 'w');
end
end
