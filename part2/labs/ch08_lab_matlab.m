%% Chapter 8 Lab - Valid and Invalid Inferences (MATLAB)
%[text] # Chapter 8 Lab - Valid and Invalid Inferences
%[text] This lab reproduces several of the "statistical lies" simulations that accompany
%[text] the book: cases where standard brain-mapping analyses mislead. You generate the
%[text] data, so you know the ground truth - and can watch thresholding, voxel selection,
%[text] and circular analysis distort it.
%[text]
%[text] Sections 1-5 use base MATLAB + Statistics Toolbox only. The optional final
%[text] section runs whole-brain versions with CANlab tools (CanlabCore + SPM12 on path).
%[text]
%[text] *Simulations adapted from the CANlab FMRI\_simulations repository*
%[text] *(github.com/canlab), Principles\_of\_fMRI\_statistical\_lies.*

%% Section 1: The winner's curse - thresholding inflates effect sizes
%[text] Every voxel has the SAME true effect, Cohen's d = 0.5, measured in N = 30
%[text] participants. Thresholding keeps only voxels whose *estimated* d cleared a bar,
%[text] so the survivors are the ones where noise pushed the estimate up.

rng(1);
N = 30; n_vox = 20000; d_true = 0.5;

data  = randn(N, n_vox) + d_true;            % every voxel truly active at d = 0.5
d_hat = mean(data) ./ std(data);             % estimated effect size per voxel
t     = d_hat .* sqrt(N);                    % one-sample t-statistic
p     = 2 * tcdf(-abs(t), N - 1);

for alpha = [.001 .0001]
    sig = p < alpha;
    fprintf('p < %.4f -> %5d of %d voxels significant | mean estimated d among them: %.3f (truth: %.1f)\n', ...
        alpha, sum(sig), n_vox, mean(d_hat(sig)), d_true);
end

% Visualize the selection
d_min = tinv(1 - .001/2, N - 1) / sqrt(N);   % minimum estimated d for p < .001
figure('Color', 'w'); hold on
histogram(d_hat, 80, 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
histogram(d_hat(p < .001), 80, 'FaceColor', [.8 .2 .3], 'EdgeColor', 'none');
xline(d_true, 'k-', 'LineWidth', 2); xline(d_min, 'r--', 'LineWidth', 1.5);
xlabel('Estimated Cohen''s d'); ylabel('Number of voxels');
title('Thresholding reports only the lucky half of the distribution');
legend({'All voxels (true d = 0.5)', 'Significant, p < .001', 'True d', 'Selection threshold'});

%[text] The significant voxels report d near 0.75 - a ~50% overstatement - and
%[text] stricter correction makes the inflation WORSE. Note also that most truly
%[text] active voxels were missed entirely (the power problem of Section 8.3).

%% Section 2: Null brain-behavior correlations that look spectacular
%[text] N = 20 participants, 10,000 voxels of pure noise: the true brain-behavior
%[text] correlation is exactly zero everywhere. We search the "brain" for the voxel
%[text] most correlated with behavior - the "voodoo correlations" scenario.

rng(7);
N = 20; n_vox = 10000;

behavior = randn(N, 1);
brain    = randn(N, n_vox);                  % NULL data: no true correlation anywhere

r = corr(behavior, brain);                   % 1 x n_vox correlations
[r_best, best] = max(abs(r));
r_best = r(best);
t_best = r_best * sqrt((N - 2) / (1 - r_best^2));
p_best = 2 * tcdf(-abs(t_best), N - 2);

fprintf('Best voxel of %d: r = %.2f, p = %.5f (uncorrected)\n', n_vox, r_best, p_best);
fprintf('Voxels with |r| > 0.5: %d;  |r| > 0.6: %d\n', sum(abs(r) > .5), sum(abs(r) > .6));

figure('Color', 'w');
plot(behavior, brain(:, best), 'o', 'Color', [0 0 .5], 'MarkerFaceColor', 'b'); lsline
text(min(behavior), max(brain(:, best)), sprintf('r = %.2f', r_best), 'FontSize', 14);
xlabel('Behavior'); ylabel('Brain activity (best voxel)');
title('A completely null brain, after searching 10,000 voxels');

%[text] How bad is it? Repeat at several sample sizes and track the maximum null |r|.

N_vals = [10 20 40 100]; n_iter = 20;
rmax = zeros(n_iter, numel(N_vals));
for j = 1:numel(N_vals)
    for i = 1:n_iter
        beh = randn(N_vals(j), 1);
        rmax(i, j) = max(abs(corr(beh, randn(N_vals(j), n_vox))));
    end
end

figure('Color', 'w');
errorbar(N_vals, mean(rmax), std(rmax), 'o-', 'Color', [.8 .2 .3], 'LineWidth', 2);
xlabel('Sample size (N)'); ylabel('Max |r| across brain'); ylim([0 1]);
title('Maximum null correlation across 10,000 voxels');

%[text] With N = 10-20, the best null voxel routinely shows |r| > 0.7 - the range of
%[text] many celebrated small-study brain-behavior correlations. Valid alternatives:
%[text] estimate correlations in independent data, or test predictive models on
%[text] held-out participants.

%% Section 3: Artifacts masquerade as signal
%[text] A t-test lights up any CONSISTENT effect, neural or not. This 2-D "slice"
%[text] contains a genuine neural region (d = 0.5) and a stronger artifact rim around
%[text] the "ventricles" (d = 1.2), like pulsation, inflow, or motion effects.

rng(3);
sz = 64; N = 20;
[xx, yy] = meshgrid(1:sz, 1:sz);

neural   = double((xx - 20).^2 + (yy - 18).^2 < 36) * 0.5;          % circle in "cortex"
in_vent  = ((xx - 32).^2 / 64  + (yy - 38).^2 / 25) < 1;
in_rim   = ((xx - 32).^2 / 100 + (yy - 38).^2 / 49) < 1;
artifact = double(in_rim & ~in_vent) * 1.2;                         % ventricle rim

truth = neural + artifact;
data  = randn(N, sz * sz) + repmat(truth(:)', N, 1);

d_hat = mean(data) ./ std(data);
t     = d_hat .* sqrt(N);
p     = 2 * tcdf(-abs(t), N - 1);
t_map = reshape(t .* (p < .001), sz, sz);                           % thresholded map

figure('Color', 'w');
subplot(1, 3, 1); imagesc(neural);   axis image off; title('True NEURAL signal (d = 0.5)');
subplot(1, 3, 2); imagesc(artifact); axis image off; title('True ARTIFACT (d = 1.2)');
subplot(1, 3, 3); imagesc(t_map);    axis image off; title('Thresholded t map, p < .001');
colormap hot; colorbar

%[text] The map shows two "activations," and the artifact is the bigger one. Nothing
%[text] in the statistics distinguishes them. Defenses are physiological, not
%[text] statistical: artifact-aware preprocessing, nuisance modeling, and skepticism
%[text] about blobs near ventricles, edges, and large vessels.

%% Section 4: A false double dissociation from circular ROI selection
%[text] Tasks A and B activate the brain IDENTICALLY: the same 500 active voxels,
%[text] the same true effect (d = 0.4). With N = 20 and p < .001, power is poor, so
%[text] each map catches a different noise-driven subset. We then select a "Task A
%[text] region" from the A map, a "Task B region" from the B map, and extract each
%[text] region's response to BOTH tasks - the classic circular analysis.
%[text] (Adapted from lie9\_false\_double\_dissociation.m.)

rng(11);
N = 20; n_vox = 5000; n_active = 500; d_true = 0.4;

truth = zeros(1, n_vox); truth(1:n_active) = d_true;   % same true map for BOTH tasks
taskA = randn(N, n_vox) + truth;
taskB = randn(N, n_vox) + truth;

pfun = @(x) 2 * tcdf(-abs(mean(x) ./ (std(x) ./ sqrt(size(x, 1)))), size(x, 1) - 1);
roiA = pfun(taskA) < .001;                             % circular ROI selection
roiB = pfun(taskB) < .001;
fprintf('Task A map: %d sig. voxels | Task B map: %d | overlap: %d\n', ...
    sum(roiA), sum(roiB), sum(roiA & roiB));

means = [mean(mean(taskA(:, roiA))) mean(mean(taskB(:, roiA))); ...
         mean(mean(taskA(:, roiB))) mean(mean(taskB(:, roiB)))];

figure('Color', 'w');
b = bar(means); b(1).FaceColor = [.3 .5 .7]; b(2).FaceColor = [.9 .6 .2];
set(gca, 'XTickLabel', {'ROI from Task A map', 'ROI from Task B map'});
ylabel('Mean activity'); yline(d_true, 'k:', 'true effect (both tasks)');
legend({'Task A', 'Task B'});
title('A "double dissociation" - from identical true effects');

%[text] A beautiful crossover interaction, manufactured from noise around identical
%[text] effects. The antidote is INDEPENDENCE: select ROIs in half the participants,
%[text] test the dissociation in the other half.

half = N / 2;
roiA1 = pfun(taskA(1:half, :)) < .01;                  % select on subjects 1-10 only
roiB1 = pfun(taskB(1:half, :)) < .01;

means_ind = [mean(mean(taskA(half+1:end, roiA1))) mean(mean(taskB(half+1:end, roiA1))); ...
             mean(mean(taskA(half+1:end, roiB1))) mean(mean(taskB(half+1:end, roiB1)))];

disp('Independent test (subjects 11-20); rows: ROI from A, ROI from B; cols: Task A, Task B')
disp(means_ind)
fprintf('(True effect for both tasks in active voxels: %.1f)\n', d_true);

%[text] With independent selection and testing, the crossover vanishes: each ROI
%[text] responds about equally to both tasks - the dissociation was never real.
%[text] (Means sit below the true d = 0.4 because selecting on only 10 subjects is
%[text] noisy and admits some null voxels; the independent test is unbiased, not
%[text] magically precise.)

%% Section 5: Box 8.1 in numbers - sensitivity, specificity, and PPV
%[text] Diagnostic use of a test depends on the base rate. Even an excellent test
%[text] collapses for rare conditions.
%[text] (Adapted from diagnostic\_testing.m, github.com/canlab FMRI\_simulations.)

calc_ppv = @(sens, spec, prev) sens .* prev ./ (sens .* prev + (1 - spec) .* (1 - prev));

examples = [.98 .98 .10; .98 .98 .01; .98 .999 .001; .90 .90 .20; .90 .80 .20];
for i = 1:size(examples, 1)
    fprintf('sens=%.3f spec=%.3f prevalence=%5.1f%%  ->  PPV = %.2f\n', ...
        examples(i, 1), examples(i, 2), 100 * examples(i, 3), ...
        calc_ppv(examples(i, 1), examples(i, 2), examples(i, 3)));
end

prev = linspace(.001, .5, 300);
figure('Color', 'w'); hold on
plot(prev, calc_ppv(.98, .98, prev), 'LineWidth', 2);
plot(prev, calc_ppv(.90, .90, prev), 'LineWidth', 2);
plot(prev, calc_ppv(.90, .80, prev), 'LineWidth', 2);
xlabel('Prevalence (base rate)'); ylabel('Positive predictive value'); ylim([0 1]);
legend({'sens 98%, spec 98%', 'sens 90%, spec 90%', 'sens 90%, spec 80%'}, 'Location', 'southeast');
title('PPV depends heavily on the base rate');

%[text] With 98% sensitivity AND specificity, a positive test for a 1%-prevalence
%[text] condition is right only about a third of the time - why brain biomarkers need
%[text] near-perfect specificity to be clinically useful for rare conditions.

%% Optional: whole-brain versions with CANlab tools
%[text] If CanlabCore + SPM12 are on your path, run whole-brain versions of Lies 1
%[text] and 2 using the same code style as the book's companion repository.
%[text] (Adapted from lie1\_high\_threshold.m and lie4\_null\_correlation.m,
%[text] github.com/canlab FMRI\_simulations.)

if exist('fmri_data', 'file') == 2

    % --- Lie 1, whole-brain: true signal vs. what survives a stringent threshold
    N = 20;
    [obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', N);   % simulated 4-D dataset
    mt = mean(true_obj);                                        % true signal image
    t  = ttest(obj, .05, 'unc');                                % observed t map

    orthviews(mt);                                              % where truth lives
    t = threshold(t, .00005, 'unc');                            % stringent threshold
    montage(t);                                                 % sparse, "focal" blobs

    % --- Lie 4, whole-brain: null brain-behavior correlation map
    obj = sim_data(fmri_data, 'n', N, 'null', 'smoothness', 10);
    corr_matrix = obj.additional_info{1}.corr_matrix;           % vector-matrix corr
    r = corr_matrix(obj.Y, obj.dat');

    r2t = @(r, n) r .* sqrt((n - 2) ./ (1 - r.^2));
    t2p = @(t, n) 2 .* (1 - tcdf(abs(t), n - 2));
    r_obj = statistic_image('volInfo', obj.volInfo, 'p', t2p(r2t(r, N), N), 'dat', r, 'type', 'r');
    r_obj = threshold(r_obj, .005, 'unc');
    montage(r_obj);                                             % null data, real-looking blobs
else
    disp('CanlabCore not found - skipping optional whole-brain section.');
end
