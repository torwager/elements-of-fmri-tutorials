%% Chapter 13 Lab: Fundamental MRI Physics (MATLAB)
% This lab accompanies Chapter 13, "Fundamental MRI Physics". You will plot
% T1 recovery and T2 decay curves for different tissues, use the MR signal
% equation to see how TR and TE create tissue contrast, render a digital
% head phantom as proton-density-, T1-, and T2-weighted images, and explore
% k-space with the 2D Fourier transform.
%
% Requirements: base MATLAB only (fft2/ifft2 are built in). CanlabCore is
% not required for this lab.
% Portions adapted from the CANlab tutorial "Lab 1: T1 decay and basic
% MATLAB" (canlab.github.io, github.com/canlab).
%
% Runtime: seconds. All data are simulated.
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part3/ch13-fundamental-mri-physics

%% 1. T1 recovery and T2 decay: exponential curves
% After an RF pulse, longitudinal magnetization recovers as
%   Mz(t)  = M0 * (1 - exp(-t/T1))
% and transverse magnetization decays as
%   Mxy(t) = M0 * exp(-t/T2).
% The famous "63%" rule falls out of the base e: at t = T1 the recovery
% curve reaches 1 - 1/e = 0.63, and at t = T2 the decay curve has fallen
% to 1/e = 0.37. Check the constant:

1 / exp(1)     % approximately 0.3679

% Define the relaxation functions as anonymous functions of time and a
% tissue-specific constant:

t1relax = @(t, T1) 1 - exp(-t ./ T1);   % longitudinal recovery
t2decay = @(t, T2) exp(-t ./ T2);       % transverse decay

% Illustrative tissue constants (ms) and relative proton densities:
%                 T1      T2     rho
% Gray matter    1000     100    0.85
% White matter    600      80    0.70
% CSF            3000    2000    1.00

T1s = [1000 600 3000];   % T1 per tissue (ms): [gray, white, CSF]
T2s = [100 80 2000];     % T2 per tissue (ms): [gray, white, CSF]
rhos = [0.85 0.70 1.00]; % relative proton density: [gray, white, CSF]
names = {'Gray matter', 'White matter', 'CSF'};
colors = [1 .5 0; 1 0 0; 0 0 1];        % orange, red, blue

%% 2. Plot the relaxation curves for each tissue
% Each curve crosses the dashed line exactly at its own time constant.
% At any single readout time the tissues sit at different heights --
% that vertical separation IS tissue contrast.

t = 1:4000;                             % time after excitation (ms)

figure('Name', 'Relaxation curves');

subplot(1, 2, 1); hold on;
for i = 1:3
    plot(t, t1relax(t, T1s(i)), 'Color', colors(i, :), 'LineWidth', 3);
    plot(T1s(i), 1 - 1/exp(1), 'o', 'Color', colors(i, :), ...
        'MarkerFaceColor', colors(i, :), 'HandleVisibility', 'off');
end
plot([0 4000], [1 - 1/exp(1), 1 - 1/exp(1)], 'k--');
xlabel('Time (ms)'); ylabel('M_z / M_0');
title('T1 (longitudinal) relaxation');
legend([names, {'63% recovered'}], 'Location', 'southeast');

subplot(1, 2, 2); hold on;
for i = 1:3
    plot(t, t2decay(t, T2s(i)), 'Color', colors(i, :), 'LineWidth', 3);
    plot(T2s(i), 1/exp(1), 'o', 'Color', colors(i, :), ...
        'MarkerFaceColor', colors(i, :), 'HandleVisibility', 'off');
end
plot([0 600], [1/exp(1), 1/exp(1)], 'k--');
xlim([0 600]);
xlabel('Time (ms)'); ylabel('M_{xy} / M_0');
title('T2 (transverse) relaxation');
legend([names, {'decayed to 37%'}]);

%% 3. The signal equation: how TR and TE create contrast
% A pulse sequence excites the spins every TR ms and reads out TE ms after
% each excitation. The measured signal from a tissue is approximately
%   S = rho * (1 - exp(-TR/T1)) * exp(-TE/T2)
% The first factor is longitudinal recovery (T1 term, set by TR); the
% second is remaining transverse signal (T2 term, set by TE).
%
% Classic regimes:
%   long TR + short TE          -> proton-density image
%   intermediate TR + short TE  -> T1-weighted image
%   long TR + intermediate TE   -> T2-weighted image

signal = @(TR, TE, T1, T2, rho) rho .* (1 - exp(-TR ./ T1)) .* exp(-TE ./ T2);

TRs = linspace(50, 6000, 400);          % TR values to sweep (ms)
TEs = linspace(1, 400, 400);            % TE values to sweep (ms)

figure('Name', 'Signal equation');

subplot(1, 2, 1); hold on;
for i = 1:3
    plot(TRs, signal(TRs, 10, T1s(i), T2s(i), rhos(i)), ...
        'Color', colors(i, :), 'LineWidth', 3);
end
plot([500 500], [0 1], 'k:');           % a T1-weighted choice
xlabel('TR (ms)'); ylabel('Signal');
title('Signal vs. TR (TE = 10 ms)');
legend([names, {'TR = 500 (T1-weighted)'}], 'Location', 'southeast');

subplot(1, 2, 2); hold on;
for i = 1:3
    plot(TEs, signal(4000, TEs, T1s(i), T2s(i), rhos(i)), ...
        'Color', colors(i, :), 'LineWidth', 3);
end
plot([100 100], [0 1], 'k:');           % a T2-weighted choice
xlabel('TE (ms)'); ylabel('Signal');
title('Signal vs. TE (TR = 4000 ms)');
legend([names, {'TE = 100 (T2-weighted)'}]);

% Left panel: at intermediate TR, white > gray > CSF (T1 contrast); at
% long TR all tissues converge toward their proton densities.
% Right panel: at intermediate TE, CSF >> gray > white (T2 contrast).

%% 4. A digital head phantom in three contrasts
% Build a simple "head" from ellipses -- a gray-matter shell, white-matter
% interior, and two CSF-filled ventricles -- then evaluate the signal
% equation voxelwise at three (TR, TE) settings.

n = 128;                                 % image matrix size (n-by-n voxels)
[x, y] = meshgrid(linspace(-1, 1, n));   % voxel coordinate grids on [-1, 1]
ellipse = @(cx, cy, axr, ayr) ((x - cx) ./ axr).^2 + ((y - cy) ./ ayr).^2 < 1;

labels = zeros(n);                       % 0 = background
labels(ellipse(0, 0, 0.72, 0.92)) = 1;   % gray-matter shell
labels(ellipse(0, 0, 0.52, 0.72)) = 2;   % white-matter interior
labels(ellipse(-0.13, -0.05, 0.09, 0.33)) = 3;   % left ventricle (CSF)
labels(ellipse(+0.13, -0.05, 0.09, 0.33)) = 3;   % right ventricle (CSF)

weighted_image = @(TR, TE) ...
    (labels == 1) .* signal(TR, TE, T1s(1), T2s(1), rhos(1)) + ...
    (labels == 2) .* signal(TR, TE, T1s(2), T2s(2), rhos(2)) + ...
    (labels == 3) .* signal(TR, TE, T1s(3), T2s(3), rhos(3));

settings = {'Proton density (TR=4000, TE=10)', 4000, 10; ...
            'T1-weighted (TR=500, TE=10)',      500, 10; ...
            'T2-weighted (TR=4000, TE=100)',   4000, 100};

figure('Name', 'Contrast weighting');
for s = 1:3
    subplot(1, 3, s);
    imagesc(weighted_image(settings{s, 2}, settings{s, 3}), [0 1]);
    axis image off; title(settings{s, 1}, 'FontSize', 9);
end
colormap gray;

% Proton density: all tissues bright and similar (water content only).
% T1-weighted: white > gray, ventricles dark -- the classic anatomical scan.
% T2-weighted: ventricles glow, white matter dark. Pathology with extra
% water (edema, tumors) lights up here, hence its clinical popularity.
%
% Try it: set TR very short AND TE long -- both factors suppress signal at
% once, giving a dim, useless image. That combination is never used.

%% 5. k-space: images as sums of spatial frequencies
% The scanner measures the 2D Fourier transform of the slice, one spatial
% frequency (kx, ky) at a time. We can go the other way: transform the
% phantom into k-space with fft2 and see what lives where.

img = weighted_image(500, 10);           % the T1-weighted phantom

F = fftshift(fft2(img));                 % image -> k-space (center in middle)
[kx, ky] = meshgrid((1:n) - n/2 - 1);    % k-space coordinate grids
k_radius = sqrt(kx.^2 + ky.^2);          % distance from k-space center

center_mask = k_radius < 10;             % keep radius < 10 (of 64): low spatial frequencies only

figure('Name', 'k-space filtering');
subplot(1, 4, 1); imagesc(img); axis image off; title('Image');
subplot(1, 4, 2); imagesc(log(1 + abs(F))); axis image off;
title('k-space (log magnitude)');
subplot(1, 4, 3);
imagesc(abs(ifft2(ifftshift(F .* center_mask)))); axis image off;
title('Center of k-space only');
subplot(1, 4, 4);
imagesc(abs(ifft2(ifftshift(F .* ~center_mask)))); axis image off;
title('Edges of k-space only');
colormap gray;

% Fraction of the total k-space energy inside the small center disk:
energy_fraction = sum(abs(F(center_mask)).^2, 'all') / sum(abs(F).^2, 'all')

% 1. k-space is brightest at its center: most image energy is in a few
%    low spatial frequencies.
% 2. Center only -> blurry but recognizable brain (gross shape/contrast).
%    This is why motion during central k-space lines is so damaging.
% 3. Edges only -> just outlines (fine detail and tissue boundaries).

%% 6. A single k-space point is a stripe pattern across the whole image
% Each k-space point encodes one 2D sinusoid spanning the entire field of
% view. Its distance from the k-space center sets the stripe frequency and
% its polar angle sets the orientation.

points = [3 0; 0 3; 10 6];               % (kx, ky) examples
ptitles = {'(k_x=3, k_y=0): vertical stripes', ...
           '(k_x=0, k_y=3): same freq, rotated 90\circ', ...
           '(k_x=10, k_y=6): higher freq, oblique'};

figure('Name', 'Single k-space points');
for s = 1:3
    Fp = zeros(n);
    Fp(n/2 + 1 + points(s, 2), n/2 + 1 + points(s, 1)) = 1;
    Fp(n/2 + 1 - points(s, 2), n/2 + 1 - points(s, 1)) = 1;  % conjugate pair
    subplot(1, 3, s);
    imagesc(real(ifft2(ifftshift(Fp)))); axis image off;
    title(ptitles{s}, 'FontSize', 9);
end
colormap gray;

% The first two points are equidistant from the origin, so their stripes
% have the same spacing but orientations 90 degrees apart. A real image is
% a weighted sum of thousands of these gratings; the inverse FFT is the
% bookkeeping that adds them up.

%% 7. Summary
% - T1 recovery and T2 decay are exponentials with tissue-specific time
%   constants; the "63%" rule follows from the base e.
% - S = rho*(1 - exp(-TR/T1))*exp(-TE/T2): TR and TE select which tissue
%   property dominates (PD, T1-weighted, or T2-weighted contrast).
% - Images are acquired in k-space; the center carries gross structure and
%   most of the energy, the periphery carries fine detail, and the inverse
%   Fourier transform assembles the image.
% - T2*-weighted (functional) images use similar timing to T2-weighted
%   ones but add sensitivity to field inhomogeneities from deoxygenated
%   hemoglobin -- the bridge to BOLD physiology in Chapter 14.
