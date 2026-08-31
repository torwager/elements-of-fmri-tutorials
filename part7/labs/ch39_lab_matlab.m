%% Chapter 39 Lab: Training and Testing Predictive Models (MATLAB)
% This lab accompanies Chapter 39, "Training and Testing Predictive Models".
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part7/ch39-training-and-testing-predictive-models
% Using simulated data where the ground truth is known, you will:
%
% 1. See how randomly splitting IMAGES (instead of SUBJECTS) across
%    cross-validation folds leaks subject-level information and produces
%    above-chance accuracy from pure noise
% 2. Fix the leak with grouped (whole-subject) holdout folds and visualize
%    performance with an ROC curve
% 3. Run nested cross-validation to choose a hyperparameter (number of PCR
%    components) without biasing the accuracy estimate
% 4. Watch effect sizes shrink from in-sample fits to an independent,
%    distribution-shifted test cohort
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples: canlab_SVM_on_unpaired_data, hyp_opt_and_mlpcr_demo,
% generalizability_example).
%
% Runtime: a few minutes. All data are simulated -- no downloads.

%% 1. Simulate a dataset with multiple images per subject
% 30 subjects, 8 images each (e.g., 8 single-trial maps per person). Each
% subject has a stable idiosyncratic "fingerprint" pattern present in all
% of their images. The class label is assigned PER SUBJECT, completely at
% random -- the brain data carry ZERO information about it, so the honest
% classification accuracy is 50%.

rng(39);                                  % reproducible

n_sub  = 30;                              % subjects
n_img  = 8;                               % images per subject
n_vox  = 500;                             % simulated voxels
n_obs  = n_sub * n_img;                   % 240 images total

subject_id  = repelem((1:n_sub)', n_img); % subject ID for each image

fingerprint = 1.5 .* randn(n_sub, n_vox); % stable per-subject pattern
Xdat = fingerprint(subject_id, :) + randn(n_obs, n_vox);

y_subject = sign(randn(n_sub, 1));        % random +1/-1 label per subject
y = y_subject(subject_id);                % each image inherits its subject's label

% Pack into an fmri_data object (barebones simulated object: .dat is
% voxels x images, .Y is the outcome)
obj     = fmri_data;
obj.dat = Xdat';                          % voxels x images
obj.Y   = y;

fprintf('Simulated %d images from %d subjects; class balance: %d vs %d subjects\n', ...
    n_obs, n_sub, sum(y_subject > 0), sum(y_subject < 0));

%% 2. The WRONG way: random stratified 5-fold cross-validation
% predict() stratifies folds on the outcome but knows nothing about
% subjects: images from the same person land in both training and test
% sets. The SVM can memorize each subject's fingerprint from their training
% images and "look up" the label for that subject's test images.

[cverr, stats] = predict(obj, 'algorithm_name', 'cv_svm', ...
    'nfolds', 5, 'error_type', 'mcr');

acc_random = 1 - cverr;
fprintf('Random-image-split accuracy: %.3f  <-- far above chance, from pure noise!\n', ...
    acc_random);

%% 3. The RIGHT way: grouped folds (hold out whole subjects)
% Pass a custom vector of holdout-set IDs to 'nfolds': one integer per
% image, with all of a subject's images sharing the same fold ID. Test
% subjects are then never seen during training. (For real studies, see also
% stratified_holdout_set and xval_stratified_holdout_leave_whole_subject_out.)

fold_of_subject = repmat(1:5, 1, n_sub / 5);   % subject -> fold (1..5)
wh_folds = fold_of_subject(subject_id)';       % one fold ID per image

[cverr_g, stats_g] = predict(obj, 'algorithm_name', 'cv_svm', ...
    'nfolds', wh_folds, 'error_type', 'mcr');

acc_grouped = 1 - cverr_g;
fprintf('Grouped (whole-subject) accuracy: %.3f  <-- honest: ~chance\n', acc_grouped);

% The gap between the two numbers is pure subject-level leakage. The same
% principle applies to twins/families (keep them in the same fold) and to
% adjacent time points in autocorrelated series (same fold, or use a
% temporal buffer).

%% 4. Visualize classifier performance with an ROC curve
% roc_plot() takes continuous scores (cross-validated distance from the SVM
% hyperplane) and true class labels, plots the ROC curve, and reports
% sensitivity, specificity, PPV, and AUC. With grouped folds and null
% labels, the curve should hug the diagonal (AUC ~ 0.5).

create_figure('ROC: grouped cross-validation');
ROC = roc_plot(stats_g.dist_from_hyperplane_xval, obj.Y > 0, 'color', 'b');

% Effect size from continuous scores: Cohen's d on the distance from the
% hyperplane is finer-grained than thresholded accuracy (see Chapter 39.5).
scores = stats_g.dist_from_hyperplane_xval;
d_null = cohens_d(scores(obj.Y > 0), scores(obj.Y < 0));
fprintf('Cross-validated effect size (Cohen''s d) for null labels: %.2f\n', d_null);

%% 5. Nested cross-validation for a hyperparameter
% Now a dataset with a REAL continuous outcome: 25 subjects x 4 images,
% and a true multivariate pattern. We use principal component regression
% (cv_pcr) and must choose the number of components -- a hyperparameter.
%
% Picking the k with the best cross-validated performance and reporting
% that same number is biased (the winner was selected on the test folds).
% Nested CV separates the jobs: an INNER grouped-CV loop (within the
% training data only) selects k; the OUTER loop evaluates the whole
% procedure on untouched subjects.

n_sub2 = 25;  n_img2 = 4;  n_obs2 = n_sub2 * n_img2;   % subjects, images per subject, total images
subject_id2 = repelem((1:n_sub2)', n_img2);       % subject ID for each image

w_true = randn(n_vox, 1) ./ sqrt(n_vox);          % true predictive pattern
X2 = randn(n_obs2, n_vox);
subj_intercept = 2 .* randn(n_sub2, 1);           % subject random intercepts
Y2 = 1.5 .* (X2 * w_true) + subj_intercept(subject_id2) + 2 .* randn(n_obs2, 1);   % signal + subject offsets + noise

obj2     = fmri_data;
obj2.dat = X2';
obj2.Y   = Y2;

fold_of_subject2 = repmat(1:5, 1, n_sub2 / 5);    % subject -> outer fold (1..5)
wh_folds2 = fold_of_subject2(subject_id2)';       % grouped outer folds

ks = [1 2 5 10 20];                               % candidate numbers of PCR components (the hyperparameter)

% --- Non-nested (biased) approach first, for comparison:
cverr_all = zeros(size(ks));
for i = 1:numel(ks)
    cverr_all(i) = predict(obj2, 'algorithm_name', 'cv_pcr', ...
        'numcomponents', ks(i), 'nfolds', wh_folds2, ...
        'error_type', 'mse', 'verbose', 0);
end
fprintf('Non-nested: best-of-grid MSE = %.2f at k = %d (optimistic)\n', ...
    min(cverr_all), ks(cverr_all == min(cverr_all)));

% --- Nested CV: outer loop = evaluation, inner loop = selection
yfit_nested = zeros(n_obs2, 1);

for f = 1:5
    te = (wh_folds2 == f);                        % outer test fold
    tr = ~te;
    train_obj = get_wh_image(obj2, find(tr));

    % Inner grouped folds among TRAINING subjects only, relabeled 1..4
    inner_folds = wh_folds2(tr);
    inner_folds(inner_folds > f) = inner_folds(inner_folds > f) - 1;

    cverr_inner = zeros(size(ks));
    for i = 1:numel(ks)                           % inner loop: selection
        cverr_inner(i) = predict(train_obj, 'algorithm_name', 'cv_pcr', ...
            'numcomponents', ks(i), 'nfolds', inner_folds, ...
            'error_type', 'mse', 'verbose', 0);
    end
    [~, best] = min(cverr_inner);

    % Refit on ALL outer-training data with the chosen k.
    % 'nfolds', 1 trains on all observations (no CV) and returns the
    % trained weights in stats.other_output: {vox_weights, intercept}.
    [~, s] = predict(train_obj, 'algorithm_name', 'cv_pcr', ...
        'numcomponents', ks(best), 'nfolds', 1, 'verbose', 0);
    w_hat = s.other_output{1};
    b0    = s.other_output{2};

    % Apply the frozen model to the untouched outer test fold
    yfit_nested(te) = obj2.dat(:, te)' * w_hat + b0;

    fprintf('Outer fold %d: inner loop chose k = %d\n', f, ks(best));
end

mse_nested = mean((yfit_nested - obj2.Y) .^ 2);
fprintf('Nested-CV MSE = %.2f (honest estimate for the full pipeline)\n', mse_nested);

% For serious hyperparameter searches over larger spaces, Bayesian
% optimization (bayesopt in the Statistics and Machine Learning Toolbox)
% is an efficient alternative to grid search -- see the CANlab
% hyp_opt_and_mlpcr_demo tutorial. The same nesting rule applies.

%% 6. Effect-size shrinkage from training to independent test
% Finally, a classification problem with a TRUE (modest) signal. We train
% an SVM once on a training cohort, freeze it, and compare Cohen's d
% computed from the SVM scores:
%   (a) in-sample (same images used for training)      -> inflated
%   (b) cross-validated within the training cohort     -> honest
%   (c) on an independent cohort whose pattern only    -> smaller still:
%       partially overlaps (e.g., a new site/population)   generalization gap

w_class = randn(n_vox, 1) ./ sqrt(n_vox);         % true discriminative pattern
n_per_class = 40;                                 % images per class in each cohort

make_cohort = @(pattern) [randn(n_per_class, n_vox) + 1.0 .* repmat(pattern', n_per_class, 1); ...
                          randn(n_per_class, n_vox) - 1.0 .* repmat(pattern', n_per_class, 1)];

y_class = [ones(n_per_class, 1); -ones(n_per_class, 1)];

Xtrain = make_cohort(w_class);                    % training cohort

% Shifted cohort: the true pattern only ~60% preserved
u = randn(n_vox, 1) ./ sqrt(n_vox);               % random direction for the pattern shift
w_shifted = 0.6 .* w_class + 0.8 .* u;            % ~60% overlap with the true pattern
Xshift = make_cohort(w_shifted);

train_obj    = fmri_data;  train_obj.dat = Xtrain';  train_obj.Y = y_class;

% (a/b) Cross-validated performance in the training cohort
[cverr_tr, stats_tr] = predict(train_obj, 'algorithm_name', 'cv_svm', ...
    'nfolds', 5, 'error_type', 'mcr');

% Train once on ALL training data ('nfolds', 1) and freeze the weights
[~, s_full] = predict(train_obj, 'algorithm_name', 'cv_svm', ...
    'nfolds', 1, 'verbose', 0);
w_svm = s_full.weight_obj.dat;                    % voxel weights

% Scores = pattern expression (dot product of weights with each image)
scores_insample = Xtrain * w_svm;                 % (a) in-sample
scores_cv       = stats_tr.dist_from_hyperplane_xval;  % (b) cross-validated
scores_shift    = Xshift * w_svm;                 % (c) independent shifted cohort

d_in    = cohens_d(scores_insample(y_class > 0), scores_insample(y_class < 0));
d_cv    = cohens_d(scores_cv(y_class > 0),       scores_cv(y_class < 0));
d_shift = cohens_d(scores_shift(y_class > 0),    scores_shift(y_class < 0));

fprintf('\nEffect size (Cohen''s d) of SVM scores:\n');
fprintf('  In-sample (train = test):      d = %.2f   <-- inflated\n', d_in);
fprintf('  Cross-validated:               d = %.2f   <-- honest, same population\n', d_cv);
fprintf('  Independent shifted cohort:    d = %.2f   <-- generalization gap\n', d_shift);

create_figure('Effect-size shrinkage');
bar([d_in d_cv d_shift], 'FaceColor', [.4 .55 .75]);
set(gca, 'XTickLabel', {'In-sample' 'Cross-validated' 'Shifted cohort'});
ylabel('Cohen''s d (SVM scores)');
title('Apparent effect size shrinks as evaluation gets more honest');

% Mapping how performance degrades across increasingly different samples
% (populations, tasks, sites) establishes a model's boundary conditions --
% see test_generalizability() in CanlabCore and the CANlab
% generalizability_example for a full bootstrap-based analysis.

%% 7. Take-home points
% * Errors in training and test data must be independent: hold out whole
%   subjects (and whole families, dyads, or runs), not random images.
% * Every data-dependent operation -- preprocessing, feature selection,
%   hyperparameters -- belongs INSIDE the cross-validation loop; nested CV
%   keeps selection and evaluation on independent observations.
% * Cross-validation evaluates a procedure; the strongest test of the
%   final (frozen) model is prospective application to independent samples,
%   where effect sizes typically shrink as samples get more different.
% * Prefer metrics without hidden flexibility: balanced accuracy and AUC
%   for classifiers, RMSE and out-of-sample R^2 for regression, and
%   effect sizes (Cohen's d from continuous scores) for benchmarking.

%% Local function: Cohen's d with pooled standard deviation

function d = cohens_d(a, b)
% Difference in means divided by the pooled standard deviation
na = numel(a);  nb = numel(b);
sp = sqrt(((na - 1) .* var(a) + (nb - 1) .* var(b)) ./ (na + nb - 2));
d  = (mean(a) - mean(b)) ./ sp;
end
