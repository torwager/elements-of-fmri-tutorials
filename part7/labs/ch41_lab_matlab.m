%% Chapter 41 Lab - Biomarkers and Translational Neuroscience (MATLAB)
% Effect size <-> accuracy <-> AUC <-> NNT conversions, the winner's curse,
% and diagnostic testing with base rates.
%
% Companion to the chapter page "41. Biomarkers and Translational
% Neuroscience". Mirrors the Python lab notebook.
%
% Adapted from CANlab FMRI_simulations (github.com/canlab):
%   effect_size_formulas.m, effect_size_vs_classification_accuracy.m,
%   Cohens_d_to_NNT.mlx, effect_size_inflation_example_sim.m,
%   diagnostic_testing.m
%
% Requirements: MATLAB with Statistics and Machine Learning Toolbox.
% (CanlabCore is not required for this lab.)

%% Part 1 - Effect size converters
% Under a normal, equal-variance model, classification accuracy, AUC, and
% NNT are all functions of the standardized effect size d (plus the
% control event rate, for NNT).

d2r        = @(d) d ./ sqrt(d.^2 + 4);           % d -> point-biserial r
r2d        = @(r) 2 .* r ./ sqrt(1 - r.^2);      % r -> d
d2acc_forc = @(d) normcdf(d ./ sqrt(2));         % forced-choice acc = AUC
d2acc_sing = @(d) normcdf(d ./ 2);               % single-interval accuracy
acc2d_sing = @(a) 2 .* norminv(a);               % accuracy -> required d
d2nnt      = @(d, cer) 1 ./ (normcdf(d + norminv(cer)) - cer); % Furukawa
auc2nnt    = @(auc) 1 ./ (2 .* auc - 1);         % Kraemer & Kupfer

% A reference table for common effect sizes
d = [0.2 0.5 0.8 1.5 2.32]';    % last row: NPS response to evoked pain
T = table(d, d2r(d), d2acc_sing(d), d2acc_forc(d), d2nnt(d, 0.5), ...
    'VariableNames', {'d' 'r' 'acc_single' 'AUC' 'NNT_cer50'});
disp(T)

fprintf('d required for 90%% accuracy, single interval: %3.2f\n', acc2d_sing(0.9));
fprintf('d required for 90%% accuracy, forced choice:   %3.2f\n', sqrt(2) * norminv(0.9));

%% Plot: accuracy and AUC as a function of d
% Forced choice (comparing one draw from each group) is easier than
% single-interval classification (deciding about one draw alone), so the
% same d yields higher accuracy.

d = 0:0.05:3;

figure('Color', 'w'); hold on
h1 = plot(d, d2acc_forc(d), 'LineWidth', 3);
h2 = plot(d, d2acc_sing(d), 'LineWidth', 3);
h3 = plot([0 3], [.5 .5], 'k:', 'LineWidth', 2);

% crosshairs at d = 0.8 (a "large" effect in conventional terms)
plot([.8 .8], [.5 d2acc_sing(.8)], 'k--')
plot([0 .8], [d2acc_sing(.8) d2acc_sing(.8)], 'k--')

xlabel('Effect size (d)'); ylabel('Classification accuracy')
set(gca, 'YLim', [.4 1], 'FontSize', 14)
legend([h1 h2 h3], {'Forced choice (= AUC)' 'Single interval' 'Chance'}, ...
    'Location', 'southeast')
title('Accuracy by effect size')

%% Part 2 - Verify the formulas by direct simulation
% Simulate controls ~ N(0,1) and patients ~ N(d,1) and measure accuracy
% empirically.

rng(1)
d_true = 0.8;
n = 100000;
x1 = randn(n, 1);            % controls
x2 = d_true + randn(n, 1);   % patients

% Forced choice: pick the higher of one patient and one control
emp_forced = mean(x2 > x1);

% Single interval: classify one observation against the optimal cutpoint,
% which is midway between the means when variances are equal
cut = d_true / 2;
emp_single = mean([x2 > cut; x1 <= cut]);

fprintf('Forced choice:  theory %5.3f, simulated %5.3f\n', d2acc_forc(d_true), emp_forced);
fprintf('Single interval: theory %5.3f, simulated %5.3f\n', d2acc_sing(d_true), emp_single);

%% Part 3 - From d to number needed to treat (NNT)
% NNT = 1 / (EER - CER): how many patients must be treated for one
% additional "responder" relative to control. Converting a continuous d to
% NNT requires a response threshold, expressed via the control event rate.

% Verify Furukawa's formula by simulation: treatment ~ N(d,1),
% control ~ N(0,1), response = value above threshold thr
d_true = 0.5;
thr = 0;                                 % response threshold on the scale
cer = 1 - normcdf(thr);                  % control event rate = 0.5
x_ctrl  = randn(n, 1);
x_treat = d_true + randn(n, 1);
eer_emp = mean(x_treat > thr);
cer_emp = mean(x_ctrl > thr);
nnt_emp = 1 / (eer_emp - cer_emp);

fprintf('NNT at d = %.1f, CER = %.2f: formula %.1f, simulated %.1f\n', ...
    d_true, cer, d2nnt(d_true, cer), nnt_emp);
fprintf('Kraemer & Kupfer (threshold-free): NNT = %.1f\n', auc2nnt(d2acc_forc(d_true)));

%% NNT depends strongly on where the response threshold is set
% For a fixed d, a stringent response criterion (low CER) implies a much
% larger NNT than a lenient one. Clinical claims must state the criterion.

d_vals = [0.2 0.5 0.8 1.2];
cer_vals = 0.02:0.01:0.90;

figure('Color', 'w'); hold on
for i = 1:length(d_vals)
    plot(cer_vals, d2nnt(d_vals(i), cer_vals), 'LineWidth', 3);
end
set(gca, 'YScale', 'log', 'FontSize', 14)
xlabel('Control event rate (CER)'); ylabel('NNT (log scale)')
legend(cellstr(num2str(d_vals', 'd = %.1f')), 'Location', 'northeast')
title('NNT as a function of d and response criterion')

%% Part 4 - The winner's curse: inflation of discovery effect sizes
% When effects are selected because they passed a significance threshold
% (out of many tests), the selected effect sizes are inflated. Simulate a
% "discovery" study with many voxels, keep the significant ones, and
% re-measure the same voxels in a "replication" study.

rng(42)
n_vox  = 5000;                 % independent tests (e.g., voxels)
n_sub  = 20;                   % subjects per study
d_true = 0.5;                  % same true effect everywhere

% Discovery study: one-sample t-test per voxel
Y_disc = d_true + randn(n_sub, n_vox);           % subjects x voxels
[~, p, ~, stats] = ttest(Y_disc);
d_disc = stats.tstat ./ sqrt(n_sub);             % observed d per voxel

alpha_v = 0.001;                                 % voxelwise threshold
sig = p < alpha_v & d_disc > 0;
fprintf('%d of %d voxels significant at p < %.3f\n', sum(sig), n_vox, alpha_v);

% Replication study: same voxels, new subjects, no selection
Y_rep = d_true + randn(n_sub, n_vox);
d_rep = mean(Y_rep) ./ std(Y_rep);

fprintf('True d = %.2f\n', d_true);
fprintf('Mean observed d, significant voxels (discovery):  %.2f\n', mean(d_disc(sig)));
fprintf('Mean observed d, same voxels (replication):       %.2f\n', mean(d_rep(sig)));

figure('Color', 'w'); hold on
histogram(d_disc(sig), 30, 'FaceColor', [.8 .3 .3], 'FaceAlpha', .5);
histogram(d_rep(sig), 30, 'FaceColor', [.3 .3 .8], 'FaceAlpha', .5);
plot([d_true d_true], get(gca, 'YLim'), 'k--', 'LineWidth', 2);
xlabel('Observed effect size (d)'); ylabel('Number of voxels')
legend({'Discovery (selected)' 'Replication (same voxels)' 'True d'})
title('Winner''s curse: selected effects shrink on replication')

%% Expected inflation grows with the number of tests and shrinks with n
% Expected maximum |Z| across v null tests, converted to d for several
% sample sizes: what pure selection can produce even with NO true effect.
% Adapted from expected_effect_size_inflation_curve.m

n_tests = round(logspace(0, 5, 30));
n_vals  = [10 20 30 100];

nsim = 2000;
d_max = zeros(length(n_tests), length(n_vals));
for v = 1:length(n_tests)
    z_max = max(randn(nsim, n_tests(v)), [], 2);   % max across tests
    p_max = 1 - normcdf(mean(z_max));              % expected max as p-value
    for i = 1:length(n_vals)
        t_equiv = tinv(1 - p_max, n_vals(i) - 1);
        d_max(v, i) = t_equiv ./ sqrt(n_vals(i));
    end
end

figure('Color', 'w');
plot(log(n_tests), d_max, 'LineWidth', 3);
legend(cellstr(num2str(n_vals', 'n = %d')), 'Location', 'northwest')
xlabel('Log number of tests performed'); ylabel('Expected max d under the null')
set(gca, 'FontSize', 14)
title('Effect size inflation from selection alone')

%% Part 5 - Diagnostic testing: base rates and PPV
% PPV is strongly driven by specificity and prevalence; sensitivity has
% less impact. Adapted from diagnostic_testing.m

calc_ppv = @(sens, spec, prev) sens .* prev ./ ...
    (sens .* prev + (1 - spec) .* (1 - prev));
calc_npv = @(sens, spec, prev) spec .* (1 - prev) ./ ...
    ((1 - sens) .* prev + spec .* (1 - prev));

% Even a "98/98" test is wrong 2 times out of 3 at 1% prevalence:
fprintf('sens 98%%, spec 98%%, prev 10%%:  PPV = %.2f\n', calc_ppv(.98, .98, .10));
fprintf('sens 98%%, spec 98%%, prev  1%%:  PPV = %.2f\n', calc_ppv(.98, .98, .01));
fprintf('sens 98%%, spec 99.9%%, prev 0.1%%: PPV = %.2f\n', calc_ppv(.98, .999, .001));

%% PPV curves
prev_lines = [.5 .2 .1 .05 .01];

figure('Color', 'w');

% PPV by sensitivity (specificity fixed)
subplot(1, 3, 1); hold on
sens_vals = .8:.005:1;
for j = 1:length(prev_lines)
    plot(sens_vals, calc_ppv(sens_vals, .98, prev_lines(j)), 'LineWidth', 3);
end
xlabel('Sensitivity'); ylabel('PPV'); title('Specificity = 98%')
legend(cellstr(num2str(100 * prev_lines', 'Prev %2.0f%%')), 'Location', 'southeast')

% PPV by specificity (sensitivity fixed)
subplot(1, 3, 2); hold on
spec_vals = .8:.005:1;
for j = 1:length(prev_lines)
    plot(spec_vals, calc_ppv(.98, spec_vals, prev_lines(j)), 'LineWidth', 3);
end
xlabel('Specificity'); ylabel('PPV'); title('Sensitivity = 98%')

% PPV by prevalence and specificity (sensitivity fixed at 90%)
subplot(1, 3, 3)
prev_vals = .01:.01:.5;
spec_vals = .8:.005:1;
[S, P] = meshgrid(spec_vals, prev_vals);
contourf(S, P, calc_ppv(.90, S, P), 20)
xlabel('Specificity'); ylabel('Prevalence'); title('PPV, sensitivity = 90%')
colorbar

%% A realistic brain-biomarker scenario
% A chronic-pain classifier performs at ~90% sensitivity / 80% specificity
% against healthy controls. Deployed where 20% of patients have the
% condition, what does a positive test mean?

sens = .90; spec = .80; prev = .20;
fprintf('Realistic biomarker: PPV = %.2f, NPV = %.2f\n', ...
    calc_ppv(sens, spec, prev), calc_npv(sens, spec, prev));

% Improving specificity to 95% helps more than improving sensitivity:
fprintf('  ...with spec -> 95%%:  PPV = %.2f\n', calc_ppv(.90, .95, prev));
fprintf('  ...with sens -> 99%%:  PPV = %.2f\n', calc_ppv(.99, .80, prev));

%% Wrap-up
% Take-homes:
% 1. d, r, accuracy, AUC, and NNT are interconvertible under normal
%    assumptions - report the clinically meaningful ones.
% 2. A "large" d of 0.8 gives only ~66% single-interval accuracy;
%    individual-level decisions need d in the 2+ range.
% 3. Selected (significant) effects are inflated; only prospective tests
%    on independent data give unbiased performance estimates.
% 4. PPV depends on the base rate: specificity and prevalence, more than
%    sensitivity, determine what a positive test means.
