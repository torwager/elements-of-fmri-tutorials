%% Lab 6 — Building a Brain Map from Scratch (MATLAB)
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part2/ch06-inferences-about-mind-brain-and
%[text] In this lab you will run the complete statistical brain-mapping loop in
%[text] miniature, on simulated data where *you* planted the true signal:
%[text] simulate a two-condition experiment at many "voxels", build a
%[text] single-subject map from trial-level data, combine subjects into a
%[text] group-level map, threshold it with and without multiple comparisons
%[text] correction, count true and false positives against the known ground
%[text] truth, watch effect sizes inflate in selected voxels, and compare
%[text] voxel-wise mapping against region-of-interest (ROI) tests.
%[text] The final (optional) section runs the same voxel-wise t-test on a real
%[text] dataset using CANlab tools.
%[text]
%[text] Requirements: base MATLAB + Statistics and Machine Learning Toolbox.
%[text] The final section additionally requires CanlabCore + SPM12 on your path.

if isempty(ver('stats'))
    warning('The Statistics Toolbox is required for ttest() below.');
end

%% 1. The generative model: plant the true signal
%[text] The brain mapping framework assumes observed data are a mixture of
%[text] true signal and noise (Figure 6.4 in the book). In a simulation we
%[text] choose the truth: a 40 x 40 grid of voxels, with an [A - B] effect of
%[text] 1 signal unit in two circular regions and exactly zero elsewhere.

rng(6);                                 % fix the random seed for reproducible results
nx = 40; ny = 40;                       % nx, ny = grid size: a 40 x 40 voxel "slice" (1,600 tests)
[xx, yy] = meshgrid(1:nx, 1:ny);

% True signal: effect of 1 inside two circular regions, 0 elsewhere
truth = ((xx-12).^2 + (yy-12).^2 < 25) | ((xx-28).^2 + (yy-25).^2 < 25);
true_effect = 1.0 * truth;

figure;
imagesc(true_effect); axis image off; colormap hot; colorbar;
title('True signal (unobservable in real data)');
fprintf('Voxels with true signal: %d of %d\n', sum(truth(:)), numel(truth));

%% 2. A single-subject map from trial-level data
%[text] One participant performs many trials of condition A (task) and
%[text] condition B (control). At every voxel: baseline 100 units, plus the
%[text] true effect on A trials only, plus trial-to-trial noise. A two-sample
%[text] t-test at each voxel gives the single-subject map (top panel of
%[text] Figure 6.2). "Massively univariate" means the same simple test is
%[text] repeated independently at every voxel.

n_trials = 40;                          % n_trials = trials per condition (A and B each)
noise_sd = 2.0;                         % noise_sd = SD of trial-to-trial noise, in signal units

% Trial-level responses at every voxel [ny x nx x n_trials]
trials_A = 100 + repmat(true_effect, 1, 1, n_trials) + noise_sd * randn(ny, nx, n_trials);
trials_B = 100 + noise_sd * randn(ny, nx, n_trials);

% Two-sample t-test at every voxel: is mean(A) different from mean(B)?
[~, p_single, ~, stat_single] = ttest2(trials_A, trials_B, 'dim', 3);
t_single = stat_single.tstat;

diff_image = mean(trials_A, 3) - mean(trials_B, 3);   % [A - B] difference image

figure;
subplot(1, 2, 1); imagesc(diff_image, [-2 2]); axis image off;
title('Unthresholded difference image'); colorbar;
subplot(1, 2, 2); imagesc(t_single, [-5 5]); axis image off;
title('Single-subject t map'); colorbar;
colormap(gca, 'parula');

%[text] The two true regions are visible — but so are convincing-looking
%[text] blobs that are pure noise. Hypothesis testing exists to tell us which
%[text] blobs are unlikely to be chance.

%% 3. A group-level map: population inference
%[text] Single-subject maps cannot support claims about brains in general.
%[text] The standard group recipe (bottom panel of Figure 6.2): compute one
%[text] [A - B] difference image per participant, then run a one-sample
%[text] t-test across participants at every voxel — treating participants as
%[text] a random effect.

n_sub = 24;                             % n_sub = number of participants
subj_sd = 1.0;                          % subj_sd = SD of between-person variability + noise

% One difference image per participant = true effect + noise [ny x nx x n_sub]
diff_imgs = repmat(true_effect, 1, 1, n_sub) + subj_sd * randn(ny, nx, n_sub);

% One-sample t-test across participants, separately at every voxel
[~, p_map, ~, stat_group] = ttest(diff_imgs, 0, 'dim', 3);
t_map = stat_group.tstat;

figure;
imagesc(t_map, [-8 8]); axis image off; colorbar;
title(sprintf('Group t map (n = %d, df = %d)', n_sub, n_sub - 1));

%% 4. Threshold the map: multiple comparisons
%[text] With 1,600 tests, p < .05 uncorrected *expects* about 73 false
%[text] positives from the 1,462 null voxels alone. We compare three choices:
%[text] uncorrected (per-voxel error rate), Bonferroni (family-wise error
%[text] rate: chance of ANY false positive), and Benjamini-Hochberg FDR
%[text] (expected proportion of significant voxels that are false positives).

alpha = 0.05;                           % alpha = acceptable false positive rate per test
n_vox = nx * ny;                        % n_vox = number of tests in the family (1,600)

sig_unc  = p_map < alpha;               % no correction
sig_bonf = p_map < alpha / n_vox;       % Bonferroni (FWER)

% Benjamini-Hochberg FDR: largest p-value rank under the BH line
p_sorted = sort(p_map(:));
below = p_sorted <= alpha * (1:n_vox)' / n_vox;
if any(below)
    p_crit = p_sorted(find(below, 1, 'last'));
else
    p_crit = 0;
end
sig_fdr = p_map <= p_crit & p_crit > 0;

figure;
titles = {'Unthresholded t map', 'p < .05 uncorrected', ...
          'FDR q < .05', 'Bonferroni FWER < .05'};
masks  = {true(ny, nx), sig_unc, sig_fdr, sig_bonf};
for i = 1:4
    subplot(1, 4, i);
    imagesc(t_map .* masks{i}, [0 8]); axis image off;
    title(titles{i});
end
colormap hot;

%[text] Score each map against the known truth — never possible in real data.
%[text] True positives: significant voxels inside true regions. False
%[text] positives: significant voxels outside them.

fprintf('Ground truth: %d active voxels, %d null voxels\n\n', ...
    sum(truth(:)), sum(~truth(:)));
names = {'Uncorrected p < .05', 'FDR q < .05', 'Bonferroni'};
masks = {sig_unc, sig_fdr, sig_bonf};
for i = 1:3
    m = masks{i};
    tp = sum(m(:) & truth(:));  fp = sum(m(:) & ~truth(:));
    fprintf('%-22s significant: %4d   true pos: %3d   false pos: %3d   sensitivity: %.2f\n', ...
        names{i}, sum(m(:)), tp, fp, tp / sum(truth(:)));
end

%[text] The tradeoff of Figure 6.4 in numbers: uncorrected maps detect nearly
%[text] all true voxels but are contaminated with false positives; Bonferroni
%[text] eliminates false positives but misses real signal at region edges;
%[text] FDR sits in between. No threshold reveals the truth exactly.

%% 5. Significant voxels overestimate their own effects
%[text] Voxels survive thresholding partly because noise happened to favor
%[text] the hypothesis there — so effect sizes estimated in selected voxels
%[text] are biased upward (the "winner's curse"). Every truly active voxel
%[text] has a planted effect of exactly 1.

est_effect = mean(diff_imgs, 3);        % estimated effect per voxel

in_truth    = est_effect(truth);                % all truly active voxels
in_selected = est_effect(truth & sig_bonf);     % ...that also passed Bonferroni
in_falsepos = est_effect(~truth & sig_unc);     % false positives (uncorrected)

fprintf('True effect size (planted):                    1.00\n');
fprintf('Mean estimate, ALL truly active voxels:        %.2f\n', mean(in_truth));
fprintf('Mean estimate, SELECTED (significant) voxels:  %.2f\n', mean(in_selected));
fprintf('Mean |estimate|, false-positive voxels:        %.2f  (true effect: 0.00)\n', ...
    mean(abs(in_falsepos)));

figure;
histogram(in_truth, 20, 'FaceAlpha', 0.6); hold on;
histogram(in_selected, 20, 'FaceAlpha', 0.6);
xline(1.0, '--k', 'true effect = 1');
xlabel('estimated [A - B] effect'); ylabel('voxel count');
legend({'all truly active voxels', 'significant (selected) voxels'});
title('Selection inflates effect estimates');

%[text] Significance tells you WHERE non-zero effects likely exist, not HOW
%[text] LARGE they are.

%% 6. The other end of the continuum: an a priori ROI
%[text] With precise prior knowledge we could average signal in one
%[text] predefined region — a single test, no correction needed (Figure 6.5).
%[text] Maximum power and unbiased estimates... IF the hypothesis is right.
%[text] Compare a correctly placed ROI with one displaced onto
%[text] mostly-null territory.

roi_correct   = (xx-12).^2 + (yy-12).^2 < 25;   % exactly on true region 1
roi_misplaced = (xx-20).^2 + (yy-4).^2  < 25;   % nearby, but wrong

rois = {roi_correct, roi_misplaced};
roi_names = {'Correct a priori ROI', 'Misplaced ROI       '};
for i = 1:2
    dmat = reshape(diff_imgs, [], n_sub);       % voxels x subjects
    roi_avg = mean(dmat(rois{i}(:), :), 1)';    % one value per subject
    [~, p_roi, ~, s] = ttest(roi_avg);
    fprintf('%s  mean effect = %5.2f   t(%d) = %5.2f   p = %.3g\n', ...
        roi_names{i}, mean(roi_avg), n_sub - 1, s.tstat, p_roi);
end

figure;
imagesc(true_effect); axis image off; colormap hot; hold on;
contour(roi_correct, [0.5 0.5], 'g', 'LineWidth', 2);
contour(roi_misplaced, [0.5 0.5], 'c', 'LineWidth', 2);
title('True signal with ROIs: correct (green) vs. misplaced (cyan)');

%[text] The correct ROI gives one decisive, unbiased test. The misplaced ROI,
%[text] sampling mostly null voxels, misses the effect entirely. Prior spatial
%[text] precision buys power — but a wrong prior makes a real effect
%[text] invisible, and the advantage is legitimate only if the ROI was truly
%[text] chosen before seeing the data.

%% 7. (Optional) The same loop on real data with CANlab tools
%[text] The identical logic — one test per voxel, threshold, visualize,
%[text] summarize regions — applied to a real dataset: contrast images for
%[text] 30 participants, [reappraise negative vs. look negative], from
%[text] Wager et al. (2008), Neuron. Requires CanlabCore + SPM12 on your path.
%[text] Adapted from the CANlab voxel-wise t-test walkthrough
%[text] (github.com/canlab/CANlab_help_examples).

if ~isempty(which('load_image_set'))

    img_obj = load_image_set('emotionreg');     % load 30 contrast images
    t = ttest(img_obj);                         % group t-test at every voxel
    t = threshold(t, .05, 'fdr', 'k', 10);      % FDR q < .05, extent >= 10 voxels

    montage(t);  drawnow, snapnow;              % show results on slices

    r = region(t);                              % contiguous significant regions
    table(r);                                   % labeled table of results

else
    disp('CanlabCore not found - skipping the real-data example.');
    disp('See https://canlab.github.io/setup/ for installation.');
end

%[text] Note that CANlab tools perform two-sided tests: hot colors show
%[text] relative activations, cool colors relative deactivations — exactly
%[text] the [Task - Control] contrast logic you just simulated.

%% Wrap-up
%[text] You have run the entire brain-mapping loop with known ground truth:
%[text] maps are statistical constructions, not pictures of activity; group
%[text] maps test generalization across people; thresholds trade error types;
%[text] selected voxels overestimate their effects; and a priori ROIs are
%[text] powerful only when prior knowledge is accurate and honestly a priori.
%[text] Try changing n_sub, the true effect size, or the noise SD, and watch
%[text] how each map responds.
