%% Chapter 25 Lab: Neuroimaging Meta-Analysis (MATLAB)
% This lab accompanies Chapter 25, "Neuroimaging Meta-Analysis". You will
% simulate a literature of studies, pool effect sizes with fixed- and
% random-effects models (forest plot), quantify heterogeneity, diagnose
% publication bias with a funnel plot, and build a coordinate-based
% MKDA-style kernel density map with a Monte Carlo FWER threshold.
%
% Sections 1-6 use base MATLAB only, on simulated data. Section 7 shows how
% the same logic runs on real coordinate databases with the CANlab MKDA
% toolbox (github.com/canlab/Canlab_MKDA_MetaAnalysis).
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples, canlab_meta_analysis_walkthrough1/2).
%
% Runtime: about a minute. All data are simulated.

%% 1. Simulate a literature of studies
% Each study i measures the same effect in n_i participants. True effects
% vary across studies (tau = 0.2) around a grand mean mu = 0.3 -- different
% scanners, task variants, populations. Observed effect = true effect +
% sampling noise with variance v_i = 1/n_i (one-sample design).

rng(7);
k   = 25;
n   = randi([10 80], k, 1);              % per-study sample sizes
mu_true  = 0.3;                          % grand mean effect
tau_true = 0.2;                          % between-study SD

theta = mu_true + tau_true * randn(k, 1);    % true study effects
v     = 1 ./ n;                              % sampling variances
d     = theta + sqrt(v) .* randn(k, 1);      % observed effect sizes
se    = sqrt(v);                             % standard errors

disp(table((1:k)', n, d, se, 'VariableNames', {'study' 'n' 'effect' 'se'}))

%% 2. Fixed- vs random-effects pooling
% Fixed effects (FE): one common true effect, weights w = 1/v.
% Random effects (RE): true effects vary; weights w* = 1/(v + tau2), with
% tau2 estimated by DerSimonian-Laird from Cochran's Q. Only RE conclusions
% generalize beyond the studies in hand (see Chapter 21 for the same
% distinction in group analysis).

w_fe  = 1 ./ v;
fe    = sum(w_fe .* d) / sum(w_fe);
fe_se = 1 / sqrt(sum(w_fe));

Q     = sum(w_fe .* (d - fe).^2);                    % Cochran's Q
c     = sum(w_fe) - sum(w_fe.^2) / sum(w_fe);
tau2  = max(0, (Q - (k - 1)) / c);                   % DerSimonian-Laird
I2    = max(0, (Q - (k - 1)) / Q);                   % % between-study variation

w_re  = 1 ./ (v + tau2);
re    = sum(w_re .* d) / sum(w_re);
re_se = 1 / sqrt(sum(w_re));

fprintf('True grand mean: %.3f (tau = %.2f)\n', mu_true, tau_true);
fprintf('Fixed effects:   %.3f +/- %.3f\n', fe, fe_se);
fprintf('Random effects:  %.3f +/- %.3f\n', re, re_se);
fprintf('Q = %.1f (df = %d), tau2 = %.3f, I2 = %.0f%%\n', Q, k - 1, tau2, 100 * I2);

% The RE standard error is larger -- that is honesty, not weakness: it
% carries the real between-study variability that FE ignores.

%% 3. Forest plot
% One row per study (whiskers = 95% CI), pooled estimates at the bottom.

[~, order] = sort(d);
figure('Color', 'w'); hold on;
for row = 1:k
    i = order(row);
    plot(d(i) + [-1.96 1.96] * se(i), [row row], '-', 'Color', [.6 .6 .6]);
    plot(d(i), row, 's', 'MarkerFaceColor', [.27 .51 .71], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 4 + 6 * w_re(i) / max(w_re));
end
plot(fe + [-1.96 1.96] * fe_se, [-1 -1], '-', 'Color', [.2 .2 .2], 'LineWidth', 4);
plot(re + [-1.96 1.96] * re_se, [-2.5 -2.5], '-', 'Color', [.7 .1 .1], 'LineWidth', 4);
xline(0, 'k-'); xline(mu_true, 'r--');
text(fe, -1, '  fixed effects', 'FontSize', 8);
text(re, -2.5, '  random effects', 'FontSize', 8, 'Color', [.7 .1 .1]);
set(gca, 'YTick', 1:k, 'YTickLabel', ...
    arrayfun(@(i) sprintf('Study %02d (n=%d)', i, n(i)), order, 'UniformOutput', false), ...
    'FontSize', 7, 'YLim', [-3.5 k + 1]);
xlabel('Standardized effect size'); title('Forest plot: 25 simulated studies');

%% 4. Publication bias and the funnel plot
% Real literatures are filtered: we simulate 80 attempted studies but a
% study is "published" only if p < .05 (or with 15% luck otherwise). The
% funnel plot (effect vs standard error, precise studies on top) shows the
% missing lower-left corner, and the pooled estimate is inflated.

k_all   = 80;
n_all   = randi([10 100], k_all, 1);
v_all   = 1 ./ n_all;
th_all  = mu_true + tau_true * randn(k_all, 1);
d_all   = th_all + sqrt(v_all) .* randn(k_all, 1);
se_all  = sqrt(v_all);

z_all   = abs(d_all ./ se_all);                      % each study's own test
p_all   = erfc(z_all / sqrt(2));                     % two-tailed normal p
pub     = p_all < 0.05 | rand(k_all, 1) < 0.15;      % publication filter

% Random-effects pooling of published studies only (same formulas as above,
% packaged as the local function pool_re at the end of this file)
re_all = pool_re(d_all, v_all);
re_pub = pool_re(d_all(pub), v_all(pub));

fprintf('Published %d of %d studies\n', sum(pub), k_all);
fprintf('RE estimate, all studies:      %.3f\n', re_all);
fprintf('RE estimate, published only:   %.3f  (true = %.3f)\n', re_pub, mu_true);

figure('Color', 'w'); hold on;
h_pub = scatter(d_all(pub),  se_all(pub), 22, [.27 .51 .71], 'filled');
h_fd  = scatter(d_all(~pub), se_all(~pub), 22, [.7 .7 .7]);   % the file drawer
se_grid = linspace(0.001, max(se_all) * 1.05, 50);
plot(mu_true - 1.96 * se_grid, se_grid, 'k--');               % pseudo-95% funnel
plot(mu_true + 1.96 * se_grid, se_grid, 'k--');
h_true = xline(mu_true, 'k-');
h_re   = xline(re_pub, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse');                                % precise on top
xlabel('Effect size'); ylabel('Standard error');
legend([h_pub h_fd h_true h_re], ...
    {'published', 'file drawer', 'true mean', 'RE (published)'}, ...
    'Location', 'northwest', 'FontSize', 8);
title('Funnel plot: publication filter removes the lower-left corner');

%% 5. Coordinate-based meta-analysis: simulate peak coordinates
% Most studies publish only peak coordinates. We simulate 21 studies on a
% 2D "axial slice" (MNI-like mm coordinates, elliptical brain mask):
%   - 20 honest studies: 85% detect a true region T at (-40, 22), reporting
%     1-3 peaks near it (scatter SD 7 mm), plus 1-5 uniform noise peaks
%   - 1 rogue study: 40 peaks tightly clustered at (45, -60) -- the
%     fixed-effects trap (like one study run at a very liberal threshold)

[xx, yy]  = meshgrid(-90:2:90, -126:2:90);           % 2 mm grid
in_brain  = (xx / 72).^2 + ((yy + 18) / 95).^2 <= 1; % elliptical mask
gx = xx(in_brain); gy = yy(in_brain);                % brain grid points

focus_T = [-40 22];                                  % true region
focus_R = [45 -60];                                  % rogue study's region

n_studies = 21;
peaks_by_study = cell(n_studies, 1);
n_subj = zeros(n_studies, 1);

for i = 1:20                                         % honest studies
    n_subj(i) = randi([10 60]);
    pk = sample_in_brain(randi([1 5]));              % noise peaks
    if rand < 0.85                                   % detects region T
        m_true = randi([1 3]);
        pk = [pk; focus_T + 7 * randn(m_true, 2)];   %#ok<AGROW>
    end
    peaks_by_study{i} = pk;
end
n_subj(21) = 12;                                     % rogue study
peaks_by_study{21} = focus_R + 4 * randn(40, 2);

all_peaks = cat(1, peaks_by_study{:});
fprintf('%d studies, %d peaks total (rogue study alone: 40)\n', ...
    n_studies, size(all_peaks, 1));

figure('Color', 'w'); hold on; axis equal;
plot(gx, gy, '.', 'Color', [.93 .93 .93], 'MarkerSize', 2);
h_hon = plot(all_peaks(1:end-40, 1), all_peaks(1:end-40, 2), 'o', ...
    'MarkerSize', 3, 'Color', [.27 .51 .71]);
h_rog = plot(all_peaks(end-39:end, 1), all_peaks(end-39:end, 2), 'o', ...
    'MarkerSize', 3, 'Color', [.7 .1 .1]);
h_T   = plot(focus_T(1), focus_T(2), 'ko', 'MarkerSize', 16);
xlabel('x (mm)'); ylabel('y (mm)');
title('Reported peak coordinates, 21 studies');
legend([h_hon h_rog h_T], {'honest studies', 'rogue study (40 peaks)', ...
    'true region T'}, 'Location', 'southeast', 'FontSize', 8);

%% 6. KDA vs MKDA maps, and a Monte Carlo FWER threshold
% Both use a spherical kernel, radius r = 10 mm.
% KDA:  count of ALL peaks within r mm (peak = unit of analysis).
% MKDA: binary indicator map per STUDY (1 within r mm of any of its peaks),
%       averaged with sqrt(N) study weights -> weighted proportion of
%       studies activating near each voxel. One peak-rich study can no
%       longer dominate.

r = 10;

% KDA: peak density
kda = zeros(size(gx));
for p = 1:size(all_peaks, 1)
    kda = kda + double((gx - all_peaks(p, 1)).^2 + (gy - all_peaks(p, 2)).^2 <= r^2);
end

% MKDA: study indicator maps and weighted proportion
maps = zeros(n_studies, numel(gx));
for i = 1:n_studies
    maps(i, :) = study_indicator(peaks_by_study{i}, gx, gy, r);
end
w = sqrt(n_subj); w = w / sum(w);
mkda = (w' * maps)';

figure('Color', 'w');
subplot(1, 2, 1);
imagesc_brain(kda, in_brain, xx, yy); title('KDA: peak density');
subplot(1, 2, 2);
imagesc_brain(mkda, in_brain, xx, yy); title('MKDA: weighted prop. of studies');

% The KDA hotspot is the rogue region; in MKDA it nearly vanishes while the
% true region T (activated by most studies) dominates. This is why original
% peak-level KDA/ALE "should not be used" (Chapter 25).

% Monte Carlo null: peaks randomly located in the brain mask, holding each
% study's number of peaks fixed. (Real MKDA randomizes each study's
% contiguous blobs, preserving within-study peak clustering.) Save the MAX
% weighted proportion per iteration -> FWER-controlling threshold.
n_iter = 200;                            % use >= 10,000 for a real analysis
max_stat = zeros(n_iter, 1);
n_peaks = cellfun(@(p) size(p, 1), peaks_by_study);

for it = 1:n_iter
    null_maps = zeros(n_studies, numel(gx));
    for i = 1:n_studies
        null_maps(i, :) = study_indicator(sample_in_brain(n_peaks(i)), gx, gy, r);
    end
    max_stat(it) = max(w' * null_maps);
end

sorted_stat = sort(max_stat);                    % base-MATLAB 95th percentile
thresh = sorted_stat(ceil(0.95 * n_iter));
fprintf('FWER threshold: weighted proportion >= %.3f\n', thresh);
fprintf('Significant grid points: %d (max MKDA = %.2f at region T)\n', ...
    sum(mkda >= thresh), max(mkda));

figure('Color', 'w');
mkda_thr = mkda; mkda_thr(mkda < thresh) = 0;
imagesc_brain(mkda_thr, in_brain, xx, yy);
title(sprintf('MKDA thresholded at FWER p < .05 (>= %.2f)', thresh));

%% 7. Real MKDA with the CANlab toolbox
% The pipeline above is a 2D miniature of what the CANlab MKDA toolbox does
% in 3D with a curated coordinate database. The code below is the real
% workflow, adapted from the CANlab meta-analysis walkthroughs
% (CANlab_help_examples/canlab_meta_analysis_walkthrough1.m). It is shown,
% not run: it needs Canlab_MKDA_MetaAnalysis + CanlabCore + SPM12 on your
% path and a coordinate database file.
%
% The database is a tab-delimited text file, one row per reported peak,
% with special columns x, y, z (MNI or T88), study, Contrast (a unique
% index per independent contrast map), Subjects (sample size), FixedRandom
% ('Fixed' coordinates are automatically downweighted), and CoordSys.
%
%   % --- Locate the database and make an analysis directory -------------
%   dbfilename  = 'Agency_meta_analysis_database.txt';
%   dbname      = which(dbfilename);          % read_database reads 'dbname'
%   analysisdir = fullfile(pwd, 'meta_analysis_example');
%   mkdir(analysisdir); cd(analysisdir)
%
%   % --- Read coordinates and set up study-level maps --------------------
%   clear DB
%   read_database;                  % script: builds DB from 'dbname'
%   DB = Meta_Setup(DB, 10);        % 10 mm kernel radius; sqrt(N) weights
%
%   % --- Indicator maps + Monte Carlo + thresholded results ---------------
%   Meta_Activation_FWE('all', DB, 500, 'nocontrasts', 'noverbose');
%   % 'all' = setup + mc + results. Use >= 5000-10000 iterations for a
%   % publishable analysis; 'nocontrasts' skips the interactive prompt for
%   % contrasts across task types.
%
%   % --- Visualize with CANlab object-oriented tools ----------------------
%   img = fmri_data('Activation_FWE_all.img', 'noverbose');
%   o2  = montage(img);                          % slice montage
%   o2  = addpoints(o2, DB.xyz, 'MarkerFaceColor', [.5 0 0], ...
%                   'Marker', 'o', 'MarkerSize', 4);   % overlay the peaks
%   r   = region(img);                           % region-class object
%   [rpos, rneg] = table(r);                     % labeled results table
%   surface(img, 'cutaway', 'ycut_mm', -30);     % surface rendering
%
% Meta_Activation_FWE builds one indicator map per contrast, weights each by
% sqrt(N) (downweighting fixed-effects studies), randomizes each map's
% contiguous blobs within gray matter on every Monte Carlo iteration, and
% thresholds with the max-statistic distribution -- exactly the logic you
% built by hand in Section 6, in 3D.

%% Local functions

function [re, re_se] = pool_re(d, v)
% Random-effects (DerSimonian-Laird) pooled estimate.
k  = numel(d);
w  = 1 ./ v;
fe = sum(w .* d) / sum(w);
Q  = sum(w .* (d - fe).^2);
c  = sum(w) - sum(w.^2) / sum(w);
tau2 = max(0, (Q - (k - 1)) / c);
ws = 1 ./ (v + tau2);
re = sum(ws .* d) / sum(ws);
re_se = 1 / sqrt(sum(ws));
end

function pts = sample_in_brain(m)
% m uniform random locations within the elliptical 2D brain mask.
pts = zeros(0, 2);
while size(pts, 1) < m
    cand = [rand(2 * m, 1) * 180 - 90, rand(2 * m, 1) * 216 - 126];
    ok = (cand(:, 1) / 72).^2 + ((cand(:, 2) + 18) / 95).^2 <= 1;
    pts = [pts; cand(ok, :)]; %#ok<AGROW>
end
pts = pts(1:m, :);
end

function ind = study_indicator(peaks, gx, gy, r)
% Binary map over brain grid points: 1 within r mm of any of this study's peaks.
mind2 = inf(size(gx));
for p = 1:size(peaks, 1)
    mind2 = min(mind2, (gx - peaks(p, 1)).^2 + (gy - peaks(p, 2)).^2);
end
ind = double(mind2 <= r^2)';
end

function imagesc_brain(vals, in_brain, xx, yy)
% Display a brain-masked vector of values as an image in mm coordinates.
img = nan(size(xx));
img(in_brain) = vals;
imagesc(xx(1, [1 end]), yy([1 end], 1), img, 'AlphaData', ~isnan(img));
set(gca, 'YDir', 'normal'); axis image; colormap hot; colorbar;
xlabel('x (mm)'); ylabel('y (mm)');
end
