%% Chapter 37 Lab: From Maps to Models -- Population Codes and Decoding (MATLAB)
% This lab accompanies Chapter 37, "Multivariate Brain Analysis: From Maps
% to Models". You will simulate multivoxel patterns in which NO single
% voxel discriminates two conditions but the pattern does, compare a mass
% univariate map with a cross-validated multivariate decoder, and see why
% decoder weights are not a localization map (and how the Haufe transform
% recovers the encoding pattern).
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part7/ch37-multivariate-brain-analysis-from-maps
%
% Requirements: Statistics and Machine Learning Toolbox (for fitcsvm).
% The optional final section uses CanlabCore (https://github.com/canlab).
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples, canlab_help_7_multivariate_prediction_basics.m).
%
% Runtime: under a minute. All data in the core sections are simulated.

%% 1. Two voxels, no single-voxel information
% Two conditions (A and B) shift two voxels in OPPOSITE directions by a
% tiny amount (+/- 0.1 against noise SD 1: per-voxel effect size d = 0.2).
% Crucially, the noise is strongly SHARED between the voxels (r = 0.98),
% as if both ride on a common global fluctuation.

rng(7);                                   % seed, for reproducibility
n = 100;                                  % trials per condition

mu_A = [-0.1  0.1];                       % condition A mean
mu_B = [ 0.1 -0.1];                       % condition B mean
Sigma = [1 .98; .98 1];                   % shared noise, correlation 0.98

% Demean the noise within each condition so the simulated effect equals
% its nominal value exactly (the noise carries no condition information).
EA = mvnrnd([0 0], Sigma, n); EA = EA - mean(EA);
EB = mvnrnd([0 0], Sigma, n); EB = EB - mean(EB);
XA = mu_A + EA;
XB = mu_B + EB;

figure('Color', 'w');
plot(XA(:,1), XA(:,2), 'bo', 'MarkerFaceColor', [.6 .6 1]); hold on;
plot(XB(:,1), XB(:,2), 'o', 'Color', [.9 .5 0], 'MarkerFaceColor', [1 .8 .5]);
xlabel('Voxel 1 activity'); ylabel('Voxel 2 activity');
legend({'Condition A' 'Condition B'});
title('Marginals overlap; the joint pattern separates');
axis equal

%% 2. Quantify: univariate tests vs. the two-voxel pattern
% Each voxel alone: t ~ 1, nonsignificant. The pattern contrast
% (voxel 1 - voxel 2) cancels the shared noise and is hugely significant.

[~, p1, ~, st1] = ttest2(XA(:,1), XB(:,1));
[~, p2, ~, st2] = ttest2(XA(:,2), XB(:,2));
fprintf('Voxel 1 alone:     t = %5.2f, p = %.3f\n', st1.tstat, p1);
fprintf('Voxel 2 alone:     t = %5.2f, p = %.3f\n', st2.tstat, p2);

dA = XA(:,1) - XA(:,2);                   % pattern score, condition A
dB = XB(:,1) - XB(:,2);                   % pattern score, condition B
[~, pd, ~, std_] = ttest2(dA, dB);
fprintf('Pattern (v1 - v2): t = %5.2f, p = %.2e\n', std_.tstat, pd);

% Cross-validated classification: one voxel vs. both voxels
X2 = [XA; XB];
y2 = [ones(n,1); -ones(n,1)];

% Linear SVMs (default box constraint C = 1), 5-fold cross-validation
cv1 = crossval(fitcsvm(X2(:,1), y2, 'KernelFunction', 'linear'), 'KFold', 5);
cvb = crossval(fitcsvm(X2,      y2, 'KernelFunction', 'linear'), 'KFold', 5);
fprintf('CV accuracy, voxel 1 alone: %3.0f%%\n', 100 * (1 - kfoldLoss(cv1)));
fprintf('CV accuracy, both voxels:   %3.0f%%\n', 100 * (1 - kfoldLoss(cvb)));

% Nothing was added to the data -- the information was there all along,
% encoded jointly rather than locally. A miniature population code.

%% 3. Scale up: 200 trials x 120 voxels, one global noise source
% Sixty "signal" voxels respond slightly more to condition A (+/- 0.1);
% the other sixty carry NO task signal. Every voxel rides on a large
% global noise source (SD 1, shared across the whole "brain") plus smaller
% independent noise (SD 0.35) -- like arousal/respiration/scanner effects.

n_tr = 200; V = 120; n_signal = 60;       % trials, voxels, signal voxels
y = repmat([1; -1], n_tr/2, 1);           % +1 = condition A, -1 = B

a = 0.1;                                  % per-voxel signal amplitude
signal = zeros(n_tr, V);
signal(:, 1:n_signal) = repmat(a * y, 1, n_signal);

% Global noise: one value per trial, added to ALL voxels. It is
% condition-independent by design, so we demean it within each condition
% (any chance imbalance would be shared noise masquerading as signal).
g = randn(n_tr, 1);
g(y == 1)  = g(y == 1)  - mean(g(y == 1));
g(y == -1) = g(y == -1) - mean(g(y == -1));

X = signal + repmat(g, 1, V) ...
    + 0.35 * randn(n_tr, V);              % 0.35 = independent voxel noise SD

fprintf('Data: %d trials x %d voxels; signal voxels = 1..%d\n', n_tr, V, n_signal);

%% 4. The mass univariate map comes up empty
% A two-sample t-test at each voxel, Bonferroni-corrected. Each voxel's
% tiny effect (d ~ 0.19) is buried in the shared noise: a few voxels flirt
% with uncorrected significance, but nothing survives correction.

[~, p_map, ~, st] = ttest2(X(y == 1, :), X(y == -1, :));
t_map = st.tstat;

t_bonf = tinv(1 - 0.025 / V, n_tr - 2);   % Bonferroni threshold (two-tailed)
fprintf('Max |t| across voxels:   %4.2f\n', max(abs(t_map)));
fprintf('Bonferroni t threshold:  %4.2f\n', t_bonf);
fprintf('Voxels at p<.05 uncorrected: %d;  surviving Bonferroni: %d\n', ...
    sum(p_map < 0.05), sum(abs(t_map) > t_bonf));

figure('Color', 'w');
imagesc(reshape(t_map, 10, 12)', [-max(abs(t_map)) max(abs(t_map))]);
colormap(redbluecmap_safe); colorbar;
title('Univariate t map (nothing survives correction)');
set(gca, 'XTick', [], 'YTick', []);

%% 5. The multivariate decoder succeeds
% Reverse the equation: all 120 voxels are PREDICTORS of the condition
% label. Train a linear SVM, evaluate with 5-fold cross-validation, and
% compare against the best single voxel under the same test.

mdl = fitcsvm(X, y, 'KernelFunction', 'linear');   % linear SVM (box constraint C = 1)
cvm = crossval(mdl, 'KFold', 5);          % 5-fold cross-validation
acc_pattern = 1 - kfoldLoss(cvm);
fprintf('Whole-pattern decoder accuracy: %5.1f%%\n', 100 * acc_pattern);

acc_single = zeros(V, 1);
for v = 1:V
    cvv = crossval(fitcsvm(X(:, v), y, 'KernelFunction', 'linear'), 'KFold', 5);
    acc_single(v) = 1 - kfoldLoss(cvv);
end
[best_acc, best_vox] = max(acc_single);
fprintf('Best single-voxel accuracy:     %5.1f%% (voxel %d)\n', 100 * best_acc, best_vox);
fprintf('Mean single-voxel accuracy:     %5.1f%%\n', 100 * mean(acc_single));

figure('Color', 'w');
histogram(100 * acc_single, 20, 'FaceColor', [.7 .8 .95]); hold on;
xline(100 * acc_pattern, 'r-', 'LineWidth', 3);
xline(50, 'k:');
xlabel('Cross-validated accuracy (%)'); ylabel('Number of voxels');
legend({'Single voxels', 'Whole pattern', 'Chance'});
title('No voxel decodes; the pattern does');

%% 6. Weights are not localization: the Haufe transform
% Ground truth: only voxels 1-60 carry signal. But the decoder assigns
% comparably large NEGATIVE weights to the no-signal voxels -- it uses
% them as a reference to estimate and subtract the global noise
% (suppressor variables). The Haufe (2014) transform,
%     a  proportional to  Cov(X) * w
% asks the forward question -- how does each voxel covary with the decoder
% score? -- and recovers the true encoding pattern.

w = mdl.Beta;                             % decoder weights (readout recipe)
A_haufe = cov(X) * w;                     % forward / encoding pattern

figure('Color', 'w');
maps = {a * [ones(1, n_signal) zeros(1, V - n_signal)] * 2, w', A_haufe'};
names = {'True encoding (condition difference)', 'Decoder weights w', ...
    'Haufe transform: cov(X) * w'};
for i = 1:3
    subplot(1, 3, i);
    m = maps{i};
    imagesc(reshape(m, 10, 12)', [-max(abs(m)) max(abs(m))]);
    colormap(redbluecmap_safe); colorbar;
    title(names{i});
    set(gca, 'XTick', [], 'YTick', []);
end

fprintf('Mean weight,  signal voxels (1-60):     %+.3f\n', mean(w(1:n_signal)));
fprintf('Mean weight,  noise-only voxels (61+):  %+.3f\n', mean(w(n_signal+1:end)));
fprintf('Mean |Haufe|, signal voxels:            %.3f\n', mean(abs(A_haufe(1:n_signal))));
fprintf('Mean |Haufe|, noise-only voxels:        %.3f\n', mean(abs(A_haufe(n_signal+1:end))));

% Reading the weight map as localization would wrongly conclude that the
% no-signal half of the "brain" encodes the task (with opposite sign!).
% Backward (decoding) models answer "how can the state be read out?";
% forward (encoding) models answer "where is the signal expressed?".

%% 7. OPTIONAL: real-data prediction with CANlab tools
% With CanlabCore + SPM on your path, the same logic runs on real images
% in a few lines using the fmri_data.predict method. The example below
% (adapted from canlab_help_7_multivariate_prediction_basics.m) predicts
% pain ratings from brain images across 33 participants using LASSO-PCR,
% with 5-fold cross-validation holding out all images from a subject
% together. NOTE: downloads a large (~500 MB) dataset from figshare --
% uncomment to run.
%
% fmri_data_file = which('bmrk3_6levels_pain_dataset.mat');
% if isempty(fmri_data_file)
%     fmri_data_file = websave('bmrk3_6levels_pain_dataset.mat', ...
%         'https://ndownloader.figshare.com/files/12708989');
% end
% load(fmri_data_file);                   % loads image_obj (fmri_data)
%
% subject_id = image_obj.additional_info.subject_id;
% holdout_set = zeros(size(subject_id));
% C = cvpartition(length(unique(subject_id)), 'KFold', 5);
% for i = 1:5
%     holdout_set(ismember(subject_id, find(test(C, i)))) = i;
% end
%
% [cverr, stats, optout] = predict(image_obj, ...
%     'algorithm_name', 'cv_lassopcr', 'nfolds', holdout_set);
%
% create_figure('pred vs obs');
% plot(stats.yfit, stats.Y, 'o'); refline;
% xlabel('Predicted pain'); ylabel('Observed pain');
% % stats.weight_obj holds the predictive weight map -- remember: a
% % readout recipe, not an activation map.

%% Explore on your own
% 1. Turn off the shared noise (replace g with independent noise per
%    voxel). Rerun: the univariate map and the decoder should now tell the
%    same story. Why?
% 2. Vary the signal amplitude a from 0.05 to 0.5 and plot max |t| and
%    decoder accuracy against it. Which analysis "wakes up" first?
% 3. Shrink the training set to 40 trials. How stable is the weight map
%    across folds, compared with the Haufe pattern? (Chapters 38-39.)
% 4. Swap the SVM for logistic regression (fitclinear with 'Learner',
%    'logistic'). Does the weight-map story change?

%% Helper: a red-blue colormap without extra toolboxes
function cmap = redbluecmap_safe
% Simple diverging blue-white-red colormap (64 levels).
n2 = 32;
up = linspace(0, 1, n2)'; down = linspace(1, 0, n2)';
cmap = [ [up; ones(n2,1)], [up; down], [ones(n2,1); down] ];
end
