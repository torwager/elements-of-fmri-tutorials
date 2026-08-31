%% Lab 31 - Multivariate Decomposition: PCA and ICA
% In this lab you will build fMRI-like data where the ground truth is known --
% two "networks," each a spatial map paired with a time course, linearly mixed
% with noise -- and then ask PCA and ICA to recover them. Because the two maps
% overlap, the true sources are NOT orthogonal: PCA finds the right subspace
% but blends the sources, while spatial ICA's independence criterion unmixes
% them. You will also use scree/permutation plots to choose dimensionality,
% implement FastICA from scratch in a few lines, and finish with a miniature
% group ICA plus dual regression on a simulated multi-subject dataset.
%
% Requirements: CanlabCore + SPM12 on your MATLAB path
% (onsets2fmridesign, create_figure), plus the Statistics Toolbox (corr).
% Adapted from CANlab tutorials (github.com/canlab).
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part6/ch31-multivariate-decomposition-pca-and-ica
% The Python lab notebook mirrors this script.

%% 1. Build ground-truth sources and mix them
% The decomposition model is X = A * S + E, where X is time x voxels, the
% columns of A (T x k) are component time courses, and the rows of S (k x V)
% are spatial maps. Our two maps share a 60-voxel overlap zone, so they are
% not orthogonal -- but, being sparse on/off blocks, they remain close to
% statistically independent. Overlap costs orthogonality, not independence.

rng(7);                                        % seed for reproducibility
T = 200; V = 360; TR = 2;                      % T = time points, V = voxels, TR = repetition time (s)

% --- Spatial maps: two overlapping blocks of voxels ---
s1 = zeros(1, V); s1(1:150)  = 1;              % network 1: voxels 1-150
s2 = zeros(1, V); s2(91:255) = 1;              % network 2: voxels 91-255 (overlap!)
S_true = [s1; s2];                             % k x V

% --- Time courses: event onsets convolved with the canonical (SPM) HRF ---
ons1 = [20 100 180 260 340]';                  % network 1 event onsets (sec)
ons2 = [60 140 220 300 380]';                  % network 2 event onsets (sec)
X1 = onsets2fmridesign({ons1}, TR, T * TR);    % last column is the intercept
X2 = onsets2fmridesign({ons2}, TR, T * TR);
A_true = [X1(:, 1) X2(:, 1)];                  % T x k

% --- Mix and add noise ---
X = A_true * S_true + 0.3 * randn(T, V);       % mix + Gaussian noise (SD = 0.3)

create_figure('sources', 2, 2);
subplot(2, 2, 1); plot(S_true', 'LineWidth', 2);
title('True spatial maps (overlapping)'); xlabel('Voxel'); legend({'map 1', 'map 2'});
subplot(2, 2, 2); plot(A_true, 'LineWidth', 1.5);
title('True time courses'); xlabel('Time (TR)');
subplot(2, 2, 3); imagesc(X); colormap gray;
title('Mixed data X (time x voxels)'); xlabel('Voxel'); ylabel('Time (TR)');
subplot(2, 2, 4); plot(X(:, 120));
title('One voxel in the overlap zone (voxel 120)'); xlabel('Time (TR)');

% How overlapping are the sources? The property that matters for PCA is
% orthogonality, so measure the cosine of the angle between the vectors
% (0 = orthogonal, 1 = identical).
cos_maps = (s1 * s2') / (norm(s1) * norm(s2))
cos_timecourses = (A_true(:, 1)' * A_true(:, 2)) / (norm(A_true(:, 1)) * norm(A_true(:, 2)))

% The maps overlap substantially (cosine ~0.4): the sources are far from
% orthogonal. PCA's spatial components are forced to be mutually orthogonal,
% so it CANNOT return these maps as they are.

%% 2. PCA via the SVD: scree plot, permutation null, eigenimages
% PCA of the mean-centered data comes from the singular value decomposition
% Xc = U * S * V'. The columns of V are eigenimages (spatial modes), the
% columns of U their time courses, and the squared singular values (divided
% by their sum) give the proportion of variance each component explains.
% A permutation null -- shuffling each voxel's time series independently --
% shows how large components get by chance.

Xc = X - mean(X);                              % mean-center each voxel
[U, S_sv, V_sv] = svd(Xc, 'econ');
sv = diag(S_sv);
varexp = 100 * sv .^ 2 / sum(sv .^ 2);

nperm = 20; null10 = zeros(nperm, 10);         % nperm = number of permutations (quick null; use more for publication)
for p = 1:nperm
    Xp = zeros(size(Xc));
    for j = 1:V, Xp(:, j) = Xc(randperm(T), j); end
    svp = svd(Xp);
    vep = 100 * svp .^ 2 / sum(svp .^ 2);
    null10(p, :) = vep(1:10)';
end

create_figure('scree', 1, 2);
subplot(1, 2, 1);
plot(varexp(1:10), 'ko-', 'LineWidth', 2); hold on;
plot(max(null10), 'r--', 'LineWidth', 2);
title('Scree plot'); xlabel('Component'); ylabel('% variance explained');
legend({'data', 'permutation max (null)'});
subplot(1, 2, 2);
plot(cumsum(varexp(1:10)), 'ko-', 'LineWidth', 2);
title('Cumulative % variance'); xlabel('Number of components');

% Two components tower above a flat noise floor, with an elbow at component 3.
% Only the first two beat the permutation maximum: keep k = 2.

% Did PCA recover the sources? The subspace, yes -- the sources, no:
pc_maps = V_sv(:, 1:2)';                       % eigenimages (rows)
pc_time = U(:, 1:2) .* sv(1:2)';               % their time courses

create_figure('pca_maps', 1, 2);
subplot(1, 2, 1); plot(S_true', 'LineWidth', 2);
title('True spatial maps'); xlabel('Voxel');
subplot(1, 2, 2); plot(pc_maps', 'LineWidth', 2);
title('PCA eigenimages (orthogonal, variance-ranked)'); xlabel('Voxel');

disp('|corr| between true maps (rows) and PC maps (columns):');
disp(abs(corr(S_true', pc_maps')));

% Each PC correlates moderately with BOTH true maps: PC 1 is a variance-
% weighted blend spanning both networks, and PC 2 is a difference component --
% the orthogonal remainder. This is not a bug; it is what PCA is for.

%% 3. Rank-k reconstruction: PCA nails the subspace
% Truncating the SVD after k components gives the best rank-k approximation
% of the data. The error plummets through k = 2, then flattens: beyond rank 2
% the components only reproduce noise. This is why PCA is the standard
% pre-reduction step before ICA.

recon_err = zeros(1, 9);
for k = 0:8
    Xk = U(:, 1:k) * S_sv(1:k, 1:k) * V_sv(:, 1:k)';
    recon_err(k + 1) = norm(Xc - Xk, 'fro') / norm(Xc, 'fro');
end
create_figure('reconstruction');
plot(0:8, recon_err, 'ko-', 'LineWidth', 2);
xlabel('Rank k'); ylabel('Relative reconstruction error');
title('Truncated SVD: error vs. rank');

%% 4. Spatial ICA, implemented from scratch (FastICA)
% ICA models X = A * S with the rows of S statistically independent -- a
% stronger requirement than uncorrelated, but one that (unlike orthogonality)
% our overlapping sparse maps nearly satisfy. The maps are strongly
% non-Gaussian (0/1 weights), and mixtures of non-Gaussian sources are more
% Gaussian than the sources themselves -- so maximizing the non-Gaussianity of
% the recovered maps points the axes back at the true sources.
%
% We implement symmetric FastICA (tanh contrast) in a few lines -- see the
% local function fastica_symm at the end of this script. For spatial ICA the
% voxels are the samples: we remove each image's spatial mean, whiten via the
% SVD (the PCA pre-reduction step), then iterate the fixed-point update.
% In practice you would use GIFT (icatb_fastICA), FSL MELODIC, or CanlabCore's
% ica() method on an fmri_data object -- same model, industrial-strength code.

k = 2;                                         % k = number of components to estimate
[ic_maps, W_unmix] = fastica_symm(Xc, k);      % k x V independent spatial maps
ic_time = Xc * pinv(ic_maps);                  % T x k time courses (least squares)

% Match ICA components to true sources (order and sign are arbitrary in ICA)
cc = corr(S_true', ic_maps');                  % true x IC, signed
[~, comp_order] = max(abs(cc), [], 2);
comp_signs = sign(cc(sub2ind(size(cc), (1:k)', comp_order)));
ic_maps_m = ic_maps(comp_order, :) .* comp_signs;
ic_time_m = ic_time(:, comp_order) .* comp_signs';

create_figure('ica_vs_truth', 2, 2);
subplot(2, 2, 1); plot(S_true', 'LineWidth', 2);
title('True spatial maps'); xlabel('Voxel');
subplot(2, 2, 2); plot(ic_maps_m', 'LineWidth', 2);
title('ICA maps (matched for sign/order)'); xlabel('Voxel');
subplot(2, 2, 3); plot(A_true, 'LineWidth', 1.5);
title('True time courses'); xlabel('Time (TR)');
subplot(2, 2, 4); plot(ic_time_m, 'LineWidth', 1.5);
title('ICA time courses'); xlabel('Time (TR)');

disp('|corr| between true maps (rows) and ICA maps (columns):');
disp(abs(corr(S_true', ic_maps_m')));
disp('|corr| between true and ICA time courses:');
disp(abs(corr(A_true, ic_time_m)));

% Each ICA component correlates near 1.0 with exactly one source: same
% two-dimensional subspace as PCA, but the independence criterion rotated the
% axes onto the sources. Note that we had to match order and flip signs by
% hand -- ICA components come with arbitrary sign, scale, and order, and are
% not ranked by variance.

%% 5. A miniature group analysis: temporal concatenation + dual regression
% Real studies run ICA on a group. The standard recipe: temporally
% concatenate all subjects' data into one tall (N*T) x V matrix, run one
% spatial ICA to get group maps G, then use dual regression per subject:
%   Stage 1 (spatial regression):   A_i = Y_i * pinv(G)   -> time courses
%   Stage 2 (temporal regression):  S_i = pinv(A_i) * Y_i -> subject maps
% Our six simulated subjects share the two networks but differ in ways dual
% regression should detect: each expresses network 2 with a different
% amplitude, and network 2's borders shift slightly across subjects.

rng(11);                                       % fresh seed for the group study
n_subj = 6;                                    % number of simulated subjects
amp2_true = linspace(0.5, 2, n_subj);          % network-2 amplitude per subject

subj_data = cell(1, n_subj); subj_maps_true = cell(1, n_subj);
for i = 1:n_subj
    s2_i = zeros(1, V);
    s2_i((91 + 3 * (i - 1)):(255 + 3 * (i - 1))) = 1;   % shifted network 2
    S_i = [s1; s2_i];

    o1 = 20 + 80 * (0:4)' + TR * randi([-2 2], 5, 1);   % jittered onsets
    o2 = 60 + 80 * (0:4)' + TR * randi([-2 2], 5, 1);
    Xi1 = onsets2fmridesign({o1}, TR, T * TR);
    Xi2 = onsets2fmridesign({o2}, TR, T * TR);
    A_i = [Xi1(:, 1) amp2_true(i) * Xi2(:, 1)];

    subj_data{i} = A_i * S_i + 0.3 * randn(T, V);
    subj_maps_true{i} = S_i;
end

% --- Group spatial ICA on the temporally concatenated data ---
Y_group = cat(1, subj_data{:});                % (N*T) x V
Y_group = Y_group - mean(Y_group);
G = fastica_symm(Y_group, k);                  % group maps, k x V

% Match group ICs to the canonical sources for readable output
cc = corr(S_true', G');
[~, comp_order] = max(abs(cc), [], 2);
comp_signs = sign(cc(sub2ind(size(cc), (1:k)', comp_order)));
G = G(comp_order, :) .* comp_signs;

create_figure('group_maps');
plot(G', 'LineWidth', 2); title('Group ICA spatial maps'); xlabel('Voxel');
legend({'group IC 1', 'group IC 2'});

% --- Dual regression for each subject ---
amp2_hat = zeros(1, n_subj); S_hat = cell(1, n_subj);
for i = 1:n_subj
    Yc = subj_data{i} - mean(subj_data{i});
    A_i = Yc * pinv(G);                        % stage 1: T x k time courses
    S_i = pinv(A_i) * Yc;                      % stage 2: k x V subject maps
    amp2_hat(i) = std(A_i(:, 2));              % expression amplitude, network 2
    S_hat{i} = S_i;
end

% Does each recovered map match its OWN subject's truth best?
map_corr = zeros(n_subj);
for i = 1:n_subj
    for j = 1:n_subj
        map_corr(i, j) = corr(S_hat{i}(2, :)', subj_maps_true{j}(2, :)');
    end
end

create_figure('dual_regression', 1, 2);
subplot(1, 2, 1); plot(amp2_true, amp2_hat, 'ko', 'MarkerFaceColor', 'k');
xlabel('True network-2 amplitude'); ylabel('Recovered (std of stage-1 time course)');
title('Amplitude recovery');
subplot(1, 2, 2); imagesc(map_corr); colorbar;
xlabel('True subject'); ylabel('Recovered subject');
title('Stage-2 map vs true map correlation');

fprintf('Correlation(true amplitude, recovered amplitude): %.2f\n', ...
    corr(amp2_true', amp2_hat'));
fprintf('Mean map correlation with own true map:     %.2f\n', mean(diag(map_corr)));
fprintf('Mean map correlation with others'' true maps: %.2f\n', ...
    mean(map_corr(~eye(n_subj))));

% Dual regression recovers both kinds of individual differences from a single
% group decomposition: stage-1 time courses track each subject's expression
% amplitude almost perfectly, and stage-2 maps match their own subject's truth
% better than other subjects' (the correlation matrix is diagonally dominant).
% In a real study these subject-level outputs feed group inference: voxel-wise
% t-tests on the stage-2 maps, and GLMs relating stage-1 time courses to
% tasks, behavior, or clinical status.

%% 6. Notes and further directions
% - With CanlabCore, spatial ICA of a real dataset is one line on an
%   fmri_data object:  icadat = ica(obj, 20);  The mixing and separation
%   matrices are stored in icadat.additional_info{1:2}, and the ica() help
%   shows the dual-regression back-reconstruction:
%       B = pinv(icadat.dat) * obj.dat;        % stage 1
%       S_hat = pinv(B') * obj.dat';           % stage 2
% - Dual regression works against ANY set of spatial maps, including labeled
%   template networks from large consortium studies -- see the CANlab dual
%   regression walkthrough (canlab_dual_regression_example) in
%   CANlab_help_examples, which uses Neurosynth topic maps as seeds.
% - Production group ICA: GIFT (icatb_fastICA and friends) or FSL MELODIC.
%
% Explore further: shrink the overlap to zero AND give the networks different
% amplitudes -- with orthogonal maps and unequal variances, PCA's axes snap
% onto the true sources. Increase the noise or the subject-to-subject map
% shifts and watch group ICA and dual regression degrade gracefully.

%% Local function: symmetric FastICA with a tanh contrast
function [S_ica, W] = fastica_symm(Xdat, k)
% Spatial ICA of a (time x voxels) matrix: returns k spatially independent
% maps (k x voxels) and the unmixing matrix W applied to the whitened data.
% Steps: remove each image's spatial mean, whiten via SVD (PCA reduction to
% k dimensions), then iterate the FastICA fixed-point update with symmetric
% decorrelation (Hyvarinen's algorithm, tanh nonlinearity).

nvox = size(Xdat, 2);
Xr = Xdat - mean(Xdat, 2);                     % remove spatial mean per image
[~, ~, Vr] = svd(Xr, 'econ');
Z = sqrt(nvox) * Vr(:, 1:k)';                  % k x V whitened spatial signals

W = orth(randn(k));                            % random orthogonal start
for iter = 1:200                               % 200 = max fixed-point iterations
    G = tanh(W * Z);
    Wnew = G * Z' / nvox - diag(mean(1 - G .^ 2, 2)) * W;
    [Evec, Eval] = eig(Wnew * Wnew');          % symmetric decorrelation:
    Wnew = Evec * diag(1 ./ sqrt(diag(Eval))) * Evec' * Wnew;  % (WW')^-1/2 W
    if max(abs(abs(diag(Wnew * W')) - 1)) < 1e-9   % convergence tolerance
        W = Wnew; break
    end
    W = Wnew;
end
S_ica = W * Z;                                 % k x V independent maps
end
