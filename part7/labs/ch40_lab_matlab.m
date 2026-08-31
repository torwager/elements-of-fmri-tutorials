%% Chapter 40 Lab (MATLAB): Applying a Fixed Brain Signature
%[text] In this lab you play the role of a researcher who has *downloaded* a
%[text] pretrained, population-level pattern (a "signature") and wants to apply it
%[text] to new subjects' data. You will: (1) define a fixed signature weight map;
%[text] (2) simulate condition images for new test subjects; (3) compute pattern
%[text] responses (dot products); (4) test paired forced-choice accuracy and
%[text] ROC / sensitivity / specificity; and (5) see how scanner gain and offset
%[text] distort some similarity metrics but not others.
%[text]
%[text] Requirements: CanlabCore (github.com/canlab/CanlabCore) + SPM12 on your
%[text] path for `roc_plot` and the optional real-data section at the end.
%[text] Sections 1-6 otherwise use base MATLAB + Statistics Toolbox only.
%[text]
%[text] Companion to the Chapter 40 page of "Elements of fMRI Analysis -
%[text] Interactive Tutorials". Parts adapted from CANlab tutorial
%[text] canlab_help_9_apply_a_multivariate_pattern_of_interest.m (github.com/canlab).

%% 1. Define a fixed signature weight map
% The signature is a FIXED set of weights, w. Pretend it was trained on prior
% studies: we never re-fit it to the test data below. Here it lives on a
% 40 x 40 "brain slice" with two positive regions and one negative region.

side = 40;
[xx, yy] = meshgrid(1:side, 1:side);
blob = @(cx, cy, sd) exp(-((xx - cx).^2 + (yy - cy).^2) ./ (2 * sd^2));

W = 1.0 * blob(12, 14, 4) + 0.8 * blob(28, 26, 5) - 0.7 * blob(30, 10, 4);
w = W(:);                                  % vectorized weights, one per "voxel"
n_vox = numel(w);

figure;
imagesc(W, [-1 1]); axis image off; colormap(gca, canlab_colormap_or_default);
colorbar; title('Fixed signature weight map w');

%% 2. Simulate condition images for new test subjects
% 30 new subjects, three conditions per subject:
%   - pain:   expresses the signature strongly (amplitude ~ N(1.0, 0.4))
%   - sound:  an arousing but nonpainful control; weak expression (~ N(0.3, 0.2))
%   - warmth: neutral control; no expression
% Each subject also has a baseline (global) offset shared by all conditions,
% and independent voxel noise -- like real condition images.

rng(2026);
n_sub = 30;

amp_pain  = 1.0 + 0.4 * randn(n_sub, 1);
amp_sound = 0.3 + 0.2 * randn(n_sub, 1);
amp_warm  = zeros(n_sub, 1);

offset   = 0.3 * randn(n_sub, 1);          % subject baseline differences
noise_sd = 4;

make_images = @(amp) amp * w' + offset * ones(1, n_vox) + noise_sd * randn(n_sub, n_vox);
img_pain  = make_images(amp_pain);         % subjects x voxels
img_sound = make_images(amp_sound);
img_warm  = make_images(amp_warm);

figure;
subplot(1, 2, 1); imagesc(W, [-1 1]); axis image off; title('Signature');
subplot(1, 2, 2); imagesc(reshape(img_pain(1, :), side, side), [-8 8]);
axis image off; title('Subject 1, pain image (noisy)');

%% 3. Pattern responses: the dot product
% The pattern response for each image is r = w' * x, a weighted average of
% image values with weights fixed by the signature. With fmri_data objects,
% the same operation is:
%   pexp = apply_mask(data_obj, signature_obj, 'pattern_expression', 'ignore_missing');

resp_pain  = img_pain  * w;
resp_sound = img_sound * w;
resp_warm  = img_warm  * w;

figure; hold on;
plot([1 2 3], [resp_pain resp_sound resp_warm]', '-', 'Color', [.7 .7 .7]);
plot(1 * ones(n_sub, 1), resp_pain,  'ro', 'MarkerFaceColor', [1 .6 .6]);
plot(2 * ones(n_sub, 1), resp_sound, 'o', 'Color', [.9 .6 0]);
plot(3 * ones(n_sub, 1), resp_warm,  'bo', 'MarkerFaceColor', [.6 .6 1]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'pain' 'sound' 'warmth'}, 'XLim', [.5 3.5]);
ylabel('Pattern response  w''x'); title('Signature responses by condition');

fprintf('Mean responses: pain = %3.1f, sound = %3.1f, warmth = %3.1f\n', ...
    mean(resp_pain), mean(resp_sound), mean(resp_warm));

%% 4. Paired forced-choice accuracy
% Within each subject, which condition has the larger response? Baseline
% (offset) differences between subjects cancel in this paired comparison.
% binotest (CanlabCore) gives a binomial test against chance (0.5).

acc_pain_warm  = mean(resp_pain > resp_warm);
acc_pain_sound = mean(resp_pain > resp_sound);

fprintf('Forced-choice accuracy, pain vs warmth: %3.0f%%\n', 100 * acc_pain_warm);
fprintf('Forced-choice accuracy, pain vs sound:  %3.0f%%\n', 100 * acc_pain_sound);

RES = binotest(double(resp_pain > resp_warm), 0.5);   % CanlabCore binomial test
fprintf('Binomial test vs chance (pain vs warmth): p = %3.6f\n', RES.p_val);

%% 5. ROC curves, sensitivity, and specificity
% Single-interval classification: one threshold for everyone. Then the paired
% forced-choice version ('twochoice'), which assumes the positive and null
% observations are entered in the same subject order.

figure;
ROC = roc_plot([resp_pain; resp_warm], [true(n_sub, 1); false(n_sub, 1)], ...
    'color', 'b', 'plothistograms');
fprintf('Single-interval: AUC = %3.2f, sens = %3.2f, spec = %3.2f\n', ...
    ROC.AUC, ROC.sensitivity, ROC.specificity);

figure;
ROC2 = roc_plot([resp_pain; resp_warm], [true(n_sub, 1); false(n_sub, 1)], ...
    'twochoice', 'color', 'r');
fprintf('Forced-choice:   accuracy = %3.2f\n', ROC2.accuracy);

% Forced-choice accuracy is usually higher: between-subject baseline variance
% counts against a single threshold, but cancels within person.

%% 6. Metric choice and calibration across scanners
% Suppose the same subjects were scanned on "Scanner B", which multiplies
% signal by a gain and adds a uniform offset. Compare three pattern-similarity
% metrics:
%   dot product  -- sensitive to gain AND offset
%   cosine       -- invariant to gain, sensitive to offset
%   correlation  -- invariant to gain and uniform offset
% (canlab_pattern_similarity implements these for fmri_data-style matrices.)

gain = 1.8; shift = 5;
img_pain_B = gain * img_pain + shift;

metrics = @(X) deal( ...
    X * w, ...                                                   % dot product
    (X * w) ./ (vecnorm(X, 2, 2) * norm(w)), ...                 % cosine
    corr(X', w));                                                % correlation

[dotA, cosA, corA] = metrics(img_pain);
[dotB, cosB, corB] = metrics(img_pain_B);

figure;
subplot(1, 3, 1); boxplot([dotA dotB], 'Labels', {'Scanner A' 'Scanner B'});
title('Dot product');
subplot(1, 3, 2); boxplot([cosA cosB], 'Labels', {'Scanner A' 'Scanner B'});
title('Cosine similarity');
subplot(1, 3, 3); boxplot([corA corB], 'Labels', {'Scanner A' 'Scanner B'});
title('Correlation');

% Within-scanner forced-choice comparisons are unaffected by gain/offset,
% because both conditions are transformed the same way:
img_warm_B = gain * img_warm + shift;
fprintf('Forced-choice accuracy, Scanner A: %3.2f\n', mean(resp_pain > resp_warm));
fprintf('Forced-choice accuracy, Scanner B: %3.2f\n', ...
    mean(img_pain_B * w > img_warm_B * w));

%% 7. Optional: real data with CANlab pattern-expression idioms
% Requires CanlabCore, SPM12, and Neuroimaging_Pattern_Masks on your path.
% Adapted from canlab_help_9_apply_a_multivariate_pattern_of_interest.m.
% We apply a whole-brain pain-predictive PLS signature (Kragel et al. 2018,
% Nature Neuroscience) to 30 subjects' [regulate - look] contrast images
% (Wager et al. 2008, Neuron).

if ~isempty(which('load_image_set'))

    test_data = load_image_set('emotionreg');       % 30 subjects' contrast images
    [pats, patnames] = load_image_set('pain_cog_emo');
    sig = get_wh_image(pats, 8);                    % whole-brain pain pattern

    % Pattern responses with three metrics. apply_mask resamples the pattern
    % to the data space and treats 0s in data images as missing values.
    pexp     = apply_mask(test_data, sig, 'pattern_expression', 'ignore_missing');
    pexp_cos = apply_mask(test_data, sig, 'pattern_expression', 'ignore_missing', 'cosine_similarity');
    pexp_r   = apply_mask(test_data, sig, 'pattern_expression', 'ignore_missing', 'correlation');

    figure;
    barplot_columns([scale(pexp) scale(pexp_cos) scale(pexp_r)], 'nofigure', ...
        'names', {'dot product (z)' 'cosine (z)' 'correlation (z)'});
    ylabel('Pattern response (z-scored for comparison)');
    title('Pain-signature responses: emotion regulation contrast');

    disp('Correlations among the three metrics across subjects:');
    disp(corr([pexp pexp_cos pexp_r]));

else
    disp('CanlabCore not found on path -- skipping the real-data section.');
end

%% Helper: colormap fallback
function cm = canlab_colormap_or_default
% Use a diverging colormap; fall back to parula-based if none available.
try
    cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
catch
    cm = [linspace(0, 1, 128)' linspace(0, 1, 128)' ones(128, 1); ...
          ones(128, 1) linspace(1, 0, 128)' linspace(1, 0, 128)'];
end
end
