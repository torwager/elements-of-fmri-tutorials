%% Chapter 38 Lab — Machine Learning Principles and Algorithms (MATLAB)
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part7/ch38-machine-learning-principles-and-algorithms
%[text] This lab makes the core principles of machine learning concrete with small
%[text] simulations: the U-shaped test error curve and the bias–variance tradeoff,
%[text] ridge vs. lasso regularization paths, cross-validation done wrong vs. right
%[text] (feature-selection leakage), and ROC curves for evaluating classifiers.
%[text]
%[text] Requirements: MATLAB with the Statistics and Machine Learning Toolbox.
%[text] The final ROC section uses CANlab Core tools (`create_figure`, `roc_plot`);
%[text] add CanlabCore to your path, or substitute `figure` and `perfcurve`.
%[text] Everything runs on simulated data in well under a minute.

%% 1. Overfitting and the U-shaped test error curve
%[text] We simulate y = sin(2x) + noise, fit polynomials of increasing degree on a
%[text] small training set, and evaluate on an independent test set. Training error
%[text] only falls as complexity grows; test error falls, bottoms out, then rises.

rng(2026);                                    % seed for reproducibility
f_true = @(x) sin(2*x);                       % true function generating the data
sigma = 0.4;                                  % noise SD

n_train = 30; n_test = 400;                   % small training set, large test set
x_train = 4 * rand(n_train, 1);
y_train = f_true(x_train) + sigma * randn(n_train, 1);
x_test  = 4 * rand(n_test, 1);
y_test  = f_true(x_test) + sigma * randn(n_test, 1);

degrees = 1:12;                               % complexity levels to sweep
[train_mse, test_mse] = deal(zeros(size(degrees)));
fits = cell(1, numel(degrees));

for d = degrees
    p = polyfit(x_train, y_train, d);          % fit uses training data ONLY
    train_mse(d) = mean((y_train - polyval(p, x_train)).^2);
    test_mse(d)  = mean((y_test  - polyval(p, x_test)).^2);
    fits{d} = p;
end

xx = linspace(0, 4, 300)';
figure('Color', 'w');
subplot(1, 2, 1); hold on;
plot(x_train, y_train, 'ko', 'MarkerSize', 4);
plot(xx, f_true(xx), 'g-', 'LineWidth', 2);
plot(xx, polyval(fits{1}, xx), 'b--');
plot(xx, polyval(fits{3}, xx), 'r-');
plot(xx, polyval(fits{12}, xx), 'm:');
ylim([-2 2]); xlabel('x'); ylabel('y');
legend({'training data' 'true function' 'degree 1' 'degree 3' 'degree 12'});
title('Fits of increasing complexity');

subplot(1, 2, 2); hold on;
plot(degrees, train_mse, 'b.-');
plot(degrees, test_mse, 'r.-');
yline(sigma^2, ':', 'irreducible noise \sigma^2');
xlabel('Polynomial degree (complexity)'); ylabel('MSE');
legend({'Training error' 'Test error'});
title('Test error is U-shaped');

[~, best] = min(test_mse);
fprintf('Best test-set degree: %d\n', best);

%[text] Degree 1 underfits (high bias); degree 12 overfits (high variance). No
%[text] model beats the irreducible noise floor sigma^2 on test data.

%% 2. The bias–variance tradeoff, visualized
%[text] Bias and variance describe behavior across repeated samples. We redraw the
%[text] training set 50 times, refit each model, and overlay the fitted curves:
%[text] a simple model is consistently wrong (high bias, low variance); a complex
%[text] model is wildly different every time (low bias, high variance).

n_sims = 50;                                  % number of re-drawn training sets
xx = linspace(0.1, 3.9, 200)';                % grid for evaluating fitted curves
show_degrees = [1 3 12];                      % underfit, about right, overfit

figure('Color', 'w');
for i = 1:numel(show_degrees)
    d = show_degrees(i);
    curves = zeros(n_sims, numel(xx));
    for s = 1:n_sims
        xs = 4 * rand(n_train, 1);
        ys = f_true(xs) + sigma * randn(n_train, 1);
        curves(s, :) = polyval(polyfit(xs, ys, d), xx)';
    end
    subplot(1, 3, i); hold on;
    plot(xx, curves', '-', 'Color', [.27 .51 .71 .15]);
    plot(xx, f_true(xx), 'g-', 'LineWidth', 2);
    plot(xx, mean(curves)', 'r--', 'LineWidth', 2);
    bias2 = mean((mean(curves)' - f_true(xx)).^2);
    varc  = mean(var(curves));
    ylim([-2 2]); xlabel('x');
    title(sprintf('degree %d: bias^2=%.3f, var=%.3f', d, bias2, varc));
end
%[text] Expected test error ≈ bias^2 + variance + sigma^2 — exactly why the test
%[text] error curve in Section 1 is U-shaped.

%% 3. Regularization: ridge and lasso coefficient paths
%[text] We simulate 40 correlated features, only 5 of which truly predict y, and
%[text] trace every coefficient as the regularization strength lambda varies.
%[text] Ridge (L2 penalty, sum of squared weights) shrinks all coefficients
%[text] smoothly but never to exactly zero; lasso (L1 penalty, sum of absolute
%[text] weights) drives unimportant coefficients exactly to zero (sparsity).

n = 80; p = 40; p_true = 5;                    % n = observations, p = features, p_true = truly predictive
latent = randn(n, 1);                          % shared factor -> correlated features
X = 0.6 * latent + randn(n, p);
X = zscore(X);
w_true = zeros(p, 1); w_true(1:p_true) = [3 -2.5 2 -1.5 1]';
y = X * w_true + randn(n, 1);

lambdas = logspace(-2, 3, 60);                 % regularization strengths, log-spaced

% Ridge path (ridge() returns standardized-scale coefficients when flag = 0 is
% omitted; use flag = 1 for coefficients on the standardized X we built)
ridge_coefs = ridge(y, X, lambdas, 1);         % p x numel(lambdas)

% Lasso path
[lasso_coefs, FitInfo] = lasso(X, y, 'Lambda', lambdas / n);

figure('Color', 'w');
subplot(1, 2, 1); hold on;
for j = 1:p
    if w_true(j) ~= 0, c = [.86 .08 .24]; else, c = [.83 .83 .83]; end
    plot(log10(lambdas), ridge_coefs(j, :), '-', 'Color', c);
end
yline(0, 'k-');
xlabel('log10(\lambda)'); ylabel('coefficient');
title('Ridge (L2): smooth shrinkage, nothing exactly 0');

subplot(1, 2, 2); hold on;
for j = 1:p
    if w_true(j) ~= 0, c = [.86 .08 .24]; else, c = [.83 .83 .83]; end
    plot(log10(FitInfo.Lambda * n), lasso_coefs(j, :), '-', 'Color', c);
end
yline(0, 'k-');
xlabel('log10(\lambda)');
title('Lasso (L1): sparse — coefficients hit exactly 0');

fprintf('Lasso: number of nonzero coefficients along the path:\n');
disp(unique(sum(abs(lasso_coefs) > 1e-8))');

%[text] Red paths = the 5 true features; gray = noise. With correlated features,
%[text] lasso picks one arbitrary representative per correlated set — on brain
%[text] voxels this yields scattered 'speckles', motivating elastic net
%[text] (lasso(..., 'Alpha', 0.5)), structured penalties, and LASSO-PCR.

%% 4. Cross-validation done WRONG vs. right
%[text] Cross-validation is only valid when every data-dependent choice — feature
%[text] selection included — happens inside the training folds. Here the features
%[text] are PURE NOISE (true accuracy = 50%). Selecting the most outcome-correlated
%[text] features on ALL the data before CV leaks test information and produces
%[text] far-above-chance 'accuracy'; re-selecting inside each fold is honest.

n = 50; p = 2000; k_feats = 20; n_datasets = 20;   % n = participants, p = noise features,
                                                   % k_feats = features kept, n_datasets = repeats
[acc_wrong, acc_right] = deal(zeros(n_datasets, 1));

for s = 1:n_datasets
    X = randn(n, p);                           % pure noise features
    y = [ones(n/2, 1); -ones(n/2, 1)];

    % --- WRONG: select features using ALL data, then cross-validate
    r = corr(X, y);
    [~, sel_all] = maxk(abs(r), k_feats);
    cv = cvpartition(n, 'KFold', 10);
    correct = 0;
    for i = 1:cv.NumTestSets
        tr = training(cv, i); te = test(cv, i);
        w = X(tr, sel_all) \ y(tr);            % least-squares linear classifier
        correct = correct + sum(sign(X(te, sel_all) * w) == y(te));
    end
    acc_wrong(s) = correct / n;

    % --- RIGHT: feature selection re-done inside every training fold
    cv = cvpartition(n, 'KFold', 10);
    correct = 0;
    for i = 1:cv.NumTestSets
        tr = training(cv, i); te = test(cv, i);
        r = corr(X(tr, :), y(tr));
        [~, sel] = maxk(abs(r), k_feats);
        w = X(tr, sel) \ y(tr);
        correct = correct + sum(sign(X(te, sel) * w) == y(te));
    end
    acc_right(s) = correct / n;
end

fprintf('WRONG (select on all data, then CV): %.1f%% (range %.0f-%.0f%%)\n', ...
    100 * mean(acc_wrong), 100 * min(acc_wrong), 100 * max(acc_wrong));
fprintf('RIGHT (select inside each fold):     %.1f%% (range %.0f-%.0f%%)\n', ...
    100 * mean(acc_right), 100 * min(acc_right), 100 * max(acc_right));
fprintf('True accuracy of any model here:     50.0%% (pure noise)\n');
fprintf('Average optimism from leakage:       %+.1f%%\n', ...
    100 * mean(acc_wrong - acc_right));

figure('Color', 'w'); hold on;
histogram(acc_wrong, 0.3:0.04:1, 'FaceColor', [.86 .08 .24]);
histogram(acc_right, 0.3:0.04:1, 'FaceColor', [.27 .51 .71]);
xline(0.5, 'k:', 'chance');
xlabel('10-fold CV accuracy'); ylabel('count of simulated datasets');
legend({'WRONG: selection before CV' 'RIGHT: selection inside CV'});
title('Feature-selection leakage on pure-noise data');

%[text] The gap — often 25–40 percentage points on pure noise — is the optimism
%[text] bought by leakage. With real data the same bias silently inflates accuracy
%[text] on top of any true signal. Cross-validation must wrap the ENTIRE pipeline.

%% 5. ROC curves and AUC
%[text] Accuracy at one threshold hides the sensitivity/specificity tradeoff. The
%[text] ROC curve sweeps the decision threshold over a continuous score; the area
%[text] under it (AUC) summarizes separability (0.5 = uninformative, 1 = perfect).
%[text] Adapted from the CANlab simple_ROC_demo
%[text] (github.com/canlab/ComputationalFoundations). Requires CanlabCore for
%[text] create_figure and roc_plot; or use perfcurve from the Statistics Toolbox.

n_per = 200;                                  % observations per class
labels = [zeros(n_per, 1); ones(n_per, 1)];   % 0 = class A, 1 = class B

scores_good = [randn(n_per, 1); 1.5 + randn(n_per, 1)];   % informative, d' = 1.5
scores_null = [randn(n_per, 1); 0.0 + randn(n_per, 1)];   % uninformative, d' = 0

% Informative model: bowed ROC, AUC well above 0.5
create_figure('ROC_informative');
ROC_good = roc_plot(scores_good, logical(labels), 'color', 'r', 'plothistograms');
title('ROC: informative model (d'' = 1.5)');

% Uninformative model: ROC hugs the diagonal, AUC near 0.5
create_figure('ROC_uninformative');
ROC_null = roc_plot(scores_null, logical(labels), 'color', 'k', 'plothistograms');
title('ROC: uninformative model (d'' = 0)');

fprintf('AUC, informative model:   %.2f\n', ROC_good.AUC);
fprintf('AUC, uninformative model: %.2f\n', ROC_null.AUC);

%[text] The uninformative ROC hugs the diagonal: every gain in sensitivity is paid
%[text] one-for-one in false positives. AUC is threshold-free and base-rate
%[text] insensitive — a standard headline metric for decoding models — though for
%[text] clinical translation the operating point matters as much as the area
%[text] (Chapter 41).

%% Wrap-up
%[text] * Test error is U-shaped in complexity; training error is an ever-improving
%[text]   illusion. Manage the bias–variance tradeoff.
%[text] * Regularization is a continuous complexity dial: ridge shrinks smoothly,
%[text]   lasso sparsifies, elastic net blends the two — trading a little bias for
%[text]   a large variance reduction.
%[text] * Cross-validation is only as honest as its weakest step: choices made with
%[text]   all the data leak and inflate accuracy, even on pure noise.
%[text] * ROC/AUC characterize a classifier across all thresholds.
%[text]
%[text] Chapter 39 builds on this: cross-validation schemes, nested CV for
%[text] hyperparameter tuning, and unbiased performance estimation.
