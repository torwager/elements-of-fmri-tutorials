%[text] %% Chapter 7 Lab — Forward and Reverse Inference (MATLAB)
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part2/ch07-forward-and-reverse-inference
%[text] Standard brain maps tell us the probability of activation given a psychological
%[text] state — *forward inference*, or sensitivity. What we usually want is the
%[text] probability of the state given activation — *reverse inference*, or positive
%[text] predictive value (PPV). In this lab you will implement Bayes' rule as a
%[text] reverse-inference calculator, verify it with a Monte Carlo simulation, explore
%[text] how PPV depends on specificity and base rate, and build a toy "reverse
%[text] inference map" showing why forward and reverse maps of the same state can look
%[text] very different. An optional final section connects this to real data with
%[text] CANlab's Neurosynth similarity tools.
%[text]
%[text] Core sections use base MATLAB only. The optional last section requires
%[text] CanlabCore + SPM12 on your path. Companion reading: the
%[text] [Chapter 7 tutorial page](https://torwager.github.io/elements-of-fmri-tutorials/book/part2/ch07-forward-and-reverse-inference).

%% 1. Bayes' rule as a reverse-inference calculator
%[text] Forward inference is P(Brain | Psy) — the **sensitivity** of a brain measure to
%[text] a psychological state. Reverse inference is P(Psy | Brain) — the **PPV**.
%[text] Bayes' rule connects them:
%[text]
%[text] $$PPV = \frac{Sens \times BR}{Sens \times BR + (1-Spec)(1-BR)}$$
%[text]
%[text] where BR is the **base rate** of the state and Spec is the **specificity**,
%[text] P(~Brain | ~Psy). Reproduce the chapter's caudate / punishment-motivation example.

% PPV from sensitivity, specificity, and base rate (Bayes' rule)
ppv = @(sens, spec, br) (sens .* br) ./ (sens .* br + (1 - spec) .* (1 - br));

scenario = ["Chapter example"; "Perfect sensitivity"; "Rare state (1% base rate)"; ...
            "Highly specific region"; "Specific AND common state"];
sens = [0.90; 1.00; 0.90; 0.90; 0.90];         % sensitivity P(Brain | Psy) per scenario
spec = [0.80; 0.80; 0.80; 0.99; 0.99];         % specificity P(~Brain | ~Psy) per scenario
base_rate = [0.10; 0.10; 0.01; 0.10; 0.50];    % base rate P(Psy) per scenario

PPV = ppv(sens, spec, base_rate);
T = table(scenario, sens, spec, base_rate, round(PPV, 3), ...
    'VariableNames', {'Scenario', 'Sensitivity', 'Specificity', 'BaseRate', 'PPV'});
disp(T)

%[text] The first three rows reproduce the chapter's numbers: PPV = 0.33 despite 90%
%[text] sensitivity; perfect sensitivity barely helps (0.36); a 1% base rate collapses
%[text] the PPV to 0.04. Only a *highly specific* region and a reasonably *common*
%[text] state make activation strong evidence for the state.

%% 2. Monte Carlo check
%[text] Bayes' rule is exact, but it is worth seeing it emerge from raw counts.
%[text] Simulate 200,000 "moments of mental life": the state is present with
%[text] probability = base rate; the region activates with probability = sensitivity
%[text] (state present) or 1 - specificity (state absent). Then count: among moments
%[text] with activation, how often was the state actually present?

rng(7);                                     % seed = 7 -> reproducible simulated data
n = 200000;                                 % number of simulated "moments of mental life"
s = 0.90; sp = 0.80; br = 0.10;             % sensitivity, specificity, base rate (chapter example)

psy = rand(n, 1) < br;                      % is the state present?
p_act = s .* psy + (1 - sp) .* ~psy;        % P(activation) this moment
brain = rand(n, 1) < p_act;                 % does the region activate?

counts = [sum(~psy & ~brain), sum(~psy & brain); ...
          sum( psy & ~brain), sum( psy & brain)];
disp(array2table(counts, 'VariableNames', {'Brain_No', 'Brain_Yes'}, ...
    'RowNames', {'Psy_No', 'Psy_Yes'}))

fprintf('Empirical P(Psy | Brain) = %.3f\n', mean(psy(brain)));
fprintf('Analytic  PPV            = %.3f\n', ppv(s, sp, br));

%[text] The empirical fraction matches Bayes' rule. The crosstab shows *why* PPV is
%[text] low: the state is rare, so a modest false-alarm rate applied to the huge
%[text] number of state-absent moments generates far more false alarms than hits.

%% 3. How selective must a region be?
%[text] Insist on a confident reverse inference — PPV >= 0.9. How specific must
%[text] activation be? Sweep specificity at several base rates (sensitivity = 0.9).

spec_grid = linspace(0.5, 0.999, 300);      % specificity grid, chance (0.5) to near-perfect
base_rates = [0.5 0.25 0.10 0.01];          % base rates P(Psy), common to rare

figure('Color', 'w'); hold on;
for b = base_rates
    plot(spec_grid, ppv(0.90, spec_grid, b), 'LineWidth', 2, ...
        'DisplayName', sprintf('base rate = %g', b));
end
yline(0.9, '--', 'PPV = 0.9', 'Color', [.5 .5 .5]);
xlabel('Specificity, P(~Brain | ~Psy)'); ylabel('PPV = P(Psy | Brain)');
title('PPV vs. specificity (sensitivity fixed at 0.90)');
legend('Location', 'west'); hold off;

% Required specificity to reach a target PPV, solving the PPV equation:
required_spec = @(sens, br, p) 1 - sens .* br .* (1 - p) ./ (p .* (1 - br));

br_grid = [0.50 0.25 0.10 0.05 0.01 0.001]';   % base rates, common (50%) to very rare (0.1%)
Req = table(br_grid, required_spec(0.9, br_grid, 0.80), ...
    required_spec(0.9, br_grid, 0.90), required_spec(0.9, br_grid, 0.95), ...
    'VariableNames', {'BaseRate', 'SpecFor_PPV80', 'SpecFor_PPV90', 'SpecFor_PPV95'});
disp(Req)

%[text] At a 50% base rate, ~89% specificity suffices for PPV = 0.9. At 10% you
%[text] already need 99%; at 1% you need ~99.9% — the region would have to activate
%[text] almost never in the absence of the state. For rare states, no single
%[text] anatomical region is plausibly that selective, which is why the field has
%[text] turned to multivariate patterns explicitly tested for specificity.

%% 4. Base rates rule: the screening problem
%[text] The same arithmetic governs medical screening. Mammography: ~90% sensitivity
%[text] and 80–99% specificity, yet with a ~1.5% ten-year base rate (age 40–50) the
%[text] PPV of a positive screen is only ~6% at 80% specificity, ~40% at 98%.

br_grid = logspace(-3, 0, 400);   % base rates from 0.1% to 100% (log-spaced)
br_mammo = 0.0147;                % 10-year breast cancer base rate, women age 40-50 (ACS 2017-2018)

figure('Color', 'w');
semilogx(br_grid, ppv(0.90, 0.80, br_grid), 'LineWidth', 2); hold on;
semilogx(br_grid, ppv(0.90, 0.98, br_grid), 'LineWidth', 2);
plot(br_mammo, ppv(0.90, 0.80, br_mammo), 'o', 'MarkerSize', 8, 'LineWidth', 2);
plot(br_mammo, ppv(0.90, 0.98, br_mammo), 'o', 'MarkerSize', 8, 'LineWidth', 2);
xline(br_mammo, ':', '10-yr base rate, age 40-50', 'Color', [.5 .5 .5]);
xlabel('Base rate P(Psy)  [log scale]'); ylabel('PPV');
title('PPV vs. base rate (sensitivity = 0.90)');
legend({'specificity = 0.80 (US mammography)', 'specificity = 0.98 (Denmark)'}, ...
    'Location', 'northwest'); hold off;

fprintf('PPV at base rate 1.47%%: spec 0.80 -> %.3f, spec 0.98 -> %.3f\n', ...
    ppv(0.9, 0.80, br_mammo), ppv(0.9, 0.98, br_mammo));

%% 5. A toy reverse-inference map (the Neurosynth idea)
%[text] A **forward-inference map** shows P(activation | term); a **reverse-inference
%[text] (association) map** asks how diagnostic activation is of the term, accounting
%[text] for how often each voxel activates across all *other* kinds of studies.
%[text]
%[text] Simulate a small 2D "brain" with two active regions:
%[text] - **Region S** (selective): activates strongly for the target state, rarely otherwise.
%[text] - **Region G** (general): activates strongly for the target state AND many
%[text]   other states — like "salience network" regions that appear in a large
%[text]   fraction of all fMRI studies.
%[text] From each voxel's sensitivity and false-alarm rate we compute the likelihood
%[text] ratio LR = P(act | Psy) / P(act | ~Psy) and the posterior P(Psy | act).

[xx, yy] = meshgrid(1:64, 1:48);
gauss2d = @(cy, cx, w) exp(-((xx - cx).^2 + (yy - cy).^2) ./ (2 * w^2));

% P(act | Psy): BOTH regions respond strongly to the target state
sens_map = 0.05 + 0.85 * (gauss2d(16, 18, 5) + gauss2d(30, 46, 6));
sens_map = min(max(sens_map, 0), 0.95);

% P(act | ~Psy): only Region G also responds to many OTHER states
fa_map = 0.05 + 0.65 * gauss2d(30, 46, 6);
fa_map = min(max(fa_map, 0), 0.90);

br = 0.10;                                           % prior probability P(Psy) of the state
lr_map = sens_map ./ fa_map;                         % likelihood ratio
posterior_map = ppv(sens_map, 1 - fa_map, br);       % P(Psy | act)

figure('Color', 'w', 'Position', [100 100 1200 320]);
titles = {'Forward map: P(act | Psy)', 'Likelihood ratio', ...
          sprintf('Reverse map: P(Psy | act), base rate = %g', br)};
maps = {sens_map, lr_map, posterior_map};
for i = 1:3
    subplot(1, 3, i);
    imagesc(maps{i}); axis image off; colorbar;
    title(titles{i});
end
colormap hot;

fprintf('Region S (selective): sens = %.2f, P(act|~Psy) = %.2f, P(Psy|act) = %.2f\n', ...
    sens_map(16, 18), fa_map(16, 18), posterior_map(16, 18));
fprintf('Region G (general):   sens = %.2f, P(act|~Psy) = %.2f, P(Psy|act) = %.2f\n', ...
    sens_map(30, 46), fa_map(30, 46), posterior_map(30, 46));

%[text] In the forward map, regions S and G look almost identical — a standard group
%[text] analysis would report both. But only the selective region S carries real
%[text] evidence about the state: activation there raises P(Psy) from the 10% prior
%[text] to ~67%, while region G leaves it near 13%. The entire difference lies in
%[text] P(act | ~Psy), which standard brain maps never measure.

%% 6. Optional: meta-analytic reverse inference with CANlab tools
%[text] With CanlabCore + SPM12 installed, you can compare a real brain pattern
%[text] against the Neurosynth meta-analytic database — which terms' maps does it
%[text] resemble most and least? High similarity to one construct and low similarity
%[text] to alternatives is evidence for specificity, the ingredient reverse
%[text] inference needs. Uncomment to run.

% % Adapted from CANlab tutorials:
% % github.com/canlab/CANlab_help_examples (neurosynth_topic_similarity_and_wedge_plot.m)
% % Requires the Neurosynth Feature Set 1 file on your MATLAB path.
%
% test_dat = load_image_set('pain_pdm');           % example: pain-predictive pattern
% image_obj = test_dat.get_wh_image(1);
%
% [image_by_feature_correlations, top_feature_tables] = ...
%     neurosynth_feature_labels(image_obj, 'images_are_replicates', false, 'noverbose');
%
% % Wedge plot of the most and least similar Neurosynth terms
% r_to_plot = [top_feature_tables{1}.testr_high; top_feature_tables{1}.testr_low];
% textlabels = [top_feature_tables{1}.words_high(:)' top_feature_tables{1}.words_low(:)'];
% create_figure('wedge_plot');
% tor_wedge_plot(r_to_plot, textlabels, 'outer_circle_radius', .3, ...
%     'colors', {[1 .7 0] [.4 0 .8]}, 'bicolor', 'nofigure');

%[text] For a pain-predictive pattern, the top-ranked terms are "pain",
%[text] "stimulation", "noxious", "somatosensory" — and the *least* similar terms
%[text] ("object", "recognition", "semantic") are just as informative: they show the
%[text] pattern is not simply tracking general visual or cognitive engagement.
%[text]
%[text] **Wrap-up:** PPV depends strongly on specificity and base rate, weakly on
%[text] sensitivity. Valid reverse inference requires measuring P(act | ~Psy) across
%[text] many confusable alternatives — the logic behind Neurosynth association maps
%[text] and specificity-tested multivariate signatures.
