---
title: "17. Image Preprocessing"
subject: "Part 4: Signal Processing and Analysis"
---

# Image Preprocessing

:::{admonition} What you will learn
:class: tip
- The standard fMRI preprocessing pipeline — distortion correction, slice-timing correction, realignment, coregistration, normalization, smoothing, temporal filtering — and what each step fixes and costs
- How EPI distortion is corrected with B0 field maps or reverse-blip ("topup") acquisitions, and why uncorrected distortion breaks coregistration and normalization
- How rigid-body realignment estimates six motion parameters, how framewise displacement (FWD) summarizes them, and how censoring and spike regression handle bad volumes
- Why normalization templates matter (MNI305, MNI152, ICBM2009, IXI555, and friends), and how too much or too little warping flexibility each produces errors
- How to choose smoothing kernels and high-pass filter cutoffs so that you remove artifacts without removing your task signal
:::

:::{admonition} 🖥️ Ways to run this chapter's code
:class: seealso
- **In your browser, no setup:** open the [interactive Python lab](./labs/ch17-lab-python.ipynb) and click the **⏻ power icon** at the top right of the notebook. Run cells top-to-bottom, starting with the first (setup/import) cell.
- **In the cloud:** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part4/labs/ch17-lab-python.ipynb) · [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part4/labs/ch17_lab_matlab.m)
- The code tabs on **this page** are static previews with copy buttons — the labs are where code runs.
:::

:::{div}
:class: run-quick
**Run this code:** [⚡ In-browser lab](./labs/ch17-lab-python.ipynb) · [Colab](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part4/labs/ch17-lab-python.ipynb) · [MATLAB Online](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part4/labs/ch17_lab_matlab.m)
:::

## Overview

Before any statistical analysis, fMRI data pass through a sequence of preprocessing steps with three aims: (1) minimize the influence of artifacts from data acquisition and head movement; (2) transform the data so it better meets statistical assumptions; and (3) for group analysis, standardize the locations of brain regions across individuals. The typical pipeline — reconstruction, distortion correction, slice-timing correction, motion correction (realignment), coregistration of structural and functional images, spatial normalization, smoothing, temporal filtering, and physiological noise correction — is summarized below. A recurring theme runs through every step: each correction removes a mixture of artifact *and* real signal, so each involves a tradeoff that is worth understanding rather than accepting on faith.

:::{figure} images/ch17_fig1_preproc_pipeline.png
:alt: Preprocessing pipeline showing the structural T1 stream (coregister to functional, warp to atlas template) and the functional stream (denoise, distortion correction, slice timing, realignment, normalization, smoothing)
:width: 90%

An example preprocessing pipeline. The structural (T1) image is coregistered to the functional images and warped to an atlas template; the resulting warping parameters are then applied to the functional time series, which has itself been denoised, distortion-corrected, slice-time corrected, and realigned. Smoothing comes last. *(Figure 17.1 from the book.)*
:::

**Reconstruction and artifact checks.** Images are acquired in k-space and reconstructed into image space at the scanner, stacked into 3-D volumes (one per TR), and aggregated into 4-D NIfTI files whose headers record voxel size, TR, TE, and more. This is the moment to check for scanner artifacts and distortions — sometimes with surprising results: one early PET study localized "anticipatory pain activation" to the temporal pole, which turned out to be the jaw of subjects clenching their teeth.

**EPI distortion correction.** Susceptibility artifacts cause signal loss and geometric warping of functional EPI images, especially near air–tissue boundaries (orbitofrontal cortex, temporal poles, brainstem), producing nonlinear mismatches with the anatomical image. Two main corrections exist, both based on extra scans acquired before the functional runs. *Unwarping* collects a B0 field map, estimates voxel displacement from field inhomogeneity, and applies the inverse. *Reverse-blip* methods acquire brief EPI scans with opposite phase-encoding directions (A–P and P–A); distortion is maximal along the phase-encoding axis with opposite sign in the two scans, so algorithms such as FSL's `topup` interpolate to the midpoint. Distortion correction is recommended — without it, functional images do not align well with structural images or atlases — but it is not magic: correction can *over*-correct, introducing new artifacts.

:::{figure} images/ch17_fig2_distortion_correction.png
:alt: EPI image with susceptibility artifacts, estimated distortion field map in the phase encoding direction, and distortion correction error with overcorrected brainstem
:width: 95%

Distortion correction. Left: a functional EPI image with characteristic distortion in orbitofrontal cortex and brainstem. Center: estimated distortion along the phase-encoding direction from a reverse-blip acquisition (red = expansion, blue = compression). Right: correction improves the match to the T1 anatomy but can also over-correct — here functional signal is expanded beyond the brainstem boundary. *(Figure 17.2 from the book.)*
:::

**Slice-timing correction.** Slices within a volume are acquired at different times, but analyses treat each volume as a snapshot. Slice-timing correction interpolates each voxel's time series to the acquisition time of a reference slice (often the middle of the TR). Getting the actual acquisition order right — sequential vs. interleaved — is essential; assuming the wrong order does harm. Interpolation always introduces some error, can interact badly with head motion, and can smear transient "spike" artifacts in time. With modern rapid acquisitions (TR ≤ 1 s), it is increasingly common to skip slice-timing correction entirely and instead shift the *task regressors* per slice or use flexible hemodynamic basis sets. Correction matters more as TR grows, for brief events, and for resting-state connectivity at long TRs, where timing offsets between regions distort correlations.

**Motion correction (realignment).** Head motion is among the most serious artifact sources in fMRI: a voxel's time series is only meaningful if it samples the same brain location at every time point. Realignment registers each volume to a reference (the first or mean image) using a rigid-body transformation with 6 parameters — translations in $x, y, z$ and rotations (roll, pitch, yaw) — estimated by iteratively minimizing squared differences between images. The six parameter time courses are saved and routinely used as nuisance covariates in the GLM. Realignment does not fix everything: spin-history effects persist during and after movement. Additional defenses include head restraints and participant training, prospective motion correction, motion covariates in the GLM, and *censoring* high-motion volumes. Censoring uses **framewise displacement**, the sum of absolute frame-to-frame changes in the six parameters, with rotations converted to millimeters as arc length on a sphere of radius 50 mm (roughly the cortex-to-center distance):

$$
\mathrm{FWD}_t = |\Delta d_{x,t}| + |\Delta d_{y,t}| + |\Delta d_{z,t}| + 50\left(|\Delta \alpha_t| + |\Delta \beta_t| + |\Delta \gamma_t|\right)
$$

with rotations in radians. Volumes exceeding a threshold (0.2–0.5 mm is typical for resting state) are removed by "scrubbing" or — preferably — modeled with *spike regression*, one nuisance regressor per bad image, which preserves the natural temporal correlation structure for inference. Typical mean FWD is about 0.05 mm in healthy young adults, 0.10 mm in the population-based UK Biobank, and 0.2–0.4 mm in children; a mean-FWD cutoff around 0.25 mm excludes only the worst runs in adults. But beware selection bias: motion is heritable, correlated with body-mass index and clinical variables, and can correlate with task states — so aggressive exclusion of images, runs, or participants can mask true effects or manufacture spurious ones. Conservative censoring plus vigorous motion *prevention* (padding, tape for tactile feedback, coaching, mock-scanner training) is the recommended combination.

**Coregistration.** The session's high-resolution structural image must be aligned to the functional images, both for visualization and so that normalization estimated from the detailed T1 can be applied to the functionals. A rigid-body transform is used, but minimizing squared intensity differences is inappropriate across modalities — tissue intensities are ordered W > G > V in functional images and V > G > W in structural images — so coregistration instead maximizes **mutual information**. Coregistration is a common failure point, largely because susceptibility distortions differ across modalities; checking it for every participant is essential quality control.

**Normalization.** Group analysis requires each voxel to correspond to the same brain region in every subject. Normalization registers each subject's anatomy to a standardized atlas space defined by a template — usually a group-average brain in volume or surface space. Linear (affine) registration adds scaling and shearing to the rigid-body parameters but is too coarse on its own; nonlinear algorithms (SPM's unified segmentation, FSL's FNIRT, ANTs) use smooth 3-D basis functions to locally stretch and shrink the image. Surface-based normalization (FreeSurfer) instead treats the cortex as a sheet and aligns folding patterns — consistently outperforming volumetric methods for cortex, and enabling smoothing that does not mix signals across gyri or into non-gray tissue. Flexibility is a double-edged sword: too many unconstrained basis functions overfit local features and can grossly distort overall brain shape (a "local minimum" solution), while too-rigid warps leave individuals poorly aligned. Inter-subject registration is one of the largest sources of error in group analysis, so inspect every normalized brain.

The choice of template matters more than many realize. "MNI space" is not one thing: MNI305, MNI152, ICBM2009 (symmetric and asymmetric variants), IXI555, and others differ in resolution, sharpness, the locations of some structures, and even the origin, shifting the whole brain relative to other templates. Report the exact template you used, and check that anatomical atlases used for localization match it.

:::{figure} images/ch17_fig5_mni_templates.png
:alt: Axial slices of six templates registered to MNI space: MNI305, MNI152, IXI555, Keuken7T, ICBM2009sym, CIT168
:width: 95%

Templates normalized to MNI space are broadly similar but differ in resolution, smoothness, and the locations of some structures. MNI305 and MNI152 were common targets in the 2000s; higher-definition templates such as IXI555 and ICBM2009 have gradually replaced them. *(Figure 17.5 from the book.)*
:::

**Spatial smoothing.** Smoothing convolves the data with a 3-D Gaussian kernel characterized by its full width at half maximum, $\mathrm{FWHM} = 2\sqrt{2\ln 2}\,\sigma \approx 2.355\,\sigma$. It compensates for residual anatomical misalignment after normalization and satisfies the smoothness assumptions of Gaussian random field theory for multiple-comparisons correction. The cost is partial-volume mixing: averaging truly active voxels with null neighbors dilutes signal. In theory, a kernel matched to the spatial extent of the activation is optimal (the matched-filter idea), but in practice one fixed kernel — typically 6 or 8 mm — is applied everywhere, even though signals in cortex are smoother than in the brainstem. MVPA studies often skip smoothing to preserve fine-grained patterns, though modest smoothing can actually help when decoding across participants.

**Temporal filtering.** Each voxel's time series contains slow scanner drift and low-frequency noise; high-pass filtering removes it, either by adding low-frequency cosine covariates to the GLM or by premultiplying data and regressors with a filtering matrix (frequency-domain filtering via the Fourier transform is equivalent). The critical rule: **your task frequencies must not overlap the filter's stopband**. SPM's default 128-s (0.0078 Hz) cutoff will remove nearly all true signal from designs with alternating blocks of about 32 s or longer — producing null results no matter how strong the effect. Resting-state analyses typically band-pass at 0.01–0.08 Hz, which brackets the canonical HRF's peak frequency (~0.03 Hz) and excludes much respiratory (0.08–0.4 Hz) and aliased cardiac noise, though interest in faster BOLD components is growing. Physiological noise correction (e.g., RETROICOR, RVHRCOR) uses recorded respiration and pulse to build nuisance regressors — with the same caveat that closes the loop on this chapter's theme: tasks change physiology, so corrections remove real signal along with artifact, and they help only when the artifact removed outweighs the signal lost.

## Hands-on tutorial

Real preprocessing runs inside packages like SPM, FSL, and fMRIPrep — but every step is easier to reason about once you have simulated it yourself. Here we look at two steps you can fully understand in a few lines of code: using realignment parameters as nuisance regressors, and designing a high-pass filter that removes drift without removing your task. The full labs add slice-timing and smoothing-kernel simulations.

**Step 1 — Motion parameters as nuisance regressors.** We simulate a voxel time series contaminated by head motion, then quantify how much variance the six realignment parameters explain.

:::{note}
The tabs below are **static previews** (with copy buttons) showing the key step in each language. To run and modify this code, use the [interactive in-browser lab](./labs/ch17-lab-python.ipynb) or the Colab / MATLAB Online links above.
:::

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Requires CanlabCore + SPM12 on your MATLAB path
% Adapted from CANlab tutorials (github.com/canlab)
rng(17); n = 200; TR = 1;

mp = cumsum(0.02 * randn(n, 6));          % 6 motion params (random walk)
mp(120:124, :) = mp(120:124, :) + 0.8;    % a sudden head jerk

true_signal = noise_arp(n, [.5 .1]);      % "neural" fluctuations, AR(2)
y = true_signal + mp * [2 1.5 1 3 2 1]' + 0.5 * randn(n, 1);

% Variance explained by motion parameters (R^2 from nuisance regression)
X  = [mp ones(n, 1)];
r  = y - X * (X \ y);
R2 = 1 - var(r) / var(y);
fprintf('Motion explains %.0f%% of voxel variance\n', 100 * R2)

% Framewise displacement: |diffs| summed, rotations (radians) x 50 mm
fwd = sum(abs(diff(mp(:, 1:3))), 2) + 50 * sum(abs(diff(mp(:, 4:6))), 2);
create_figure('FWD'); plot(fwd); plot_horizontal_line(0.5);
xlabel('Frame'); ylabel('FWD (mm)');
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np

rng = np.random.default_rng(17)
n, TR = 200, 1.0

mp = np.cumsum(0.02 * rng.standard_normal((n, 6)), axis=0)  # 6 motion params
mp[120:125] += 0.8                                          # a sudden head jerk

# AR(2)-like "neural" fluctuations plus motion-coupled artifact
true_signal = np.convolve(rng.standard_normal(n), [1, .5, .1])[:n]
y = true_signal + mp @ np.array([2, 1.5, 1, 3, 2, 1]) \
    + 0.5 * rng.standard_normal(n)

# Variance explained by motion parameters (R^2 from nuisance regression)
X = np.column_stack([mp, np.ones(n)])
resid = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
print(f"Motion explains {100 * (1 - resid.var() / y.var()):.0f}% of variance")

# Framewise displacement: |diffs| summed, rotations (radians) x 50 mm
fwd = (np.abs(np.diff(mp[:, :3], axis=0)).sum(1)
       + 50 * np.abs(np.diff(mp[:, 3:], axis=0)).sum(1))
print(f"{(fwd > 0.5).sum()} frames exceed a 0.5 mm censoring threshold")
```
:::
::::

**Step 2 — High-pass filtering without destroying the task.** We build a discrete-cosine high-pass filter (the same construction SPM uses) and apply it as a residual-forming matrix, following the CANlab principle that filtering and nuisance regression should happen in *one* step.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Adapted from CANlab_help_examples: linear_filtering_a_timeseries.m
TR = 2; n = 300; hpf = 128;                       % cutoff in seconds

[S, KL, KH] = use_spm_filter(TR, n, 'none', 'specify', hpf);
% KH: low-frequency cosine regressors; S: residual-forming matrix

task  = repmat([ones(15, 1); zeros(15, 1)], 10, 1);   % 30-s alternating blocks
drift = 4 * cos((1:n)' * 2 * pi / 280);               % slow scanner drift
y_obs = task + drift + noise_arp(n, [.7 .3]);

y_filt = S * y_obs;                                % drift removed
fprintf('corr(task, observed) = %.2f, corr(task, filtered) = %.2f\n', ...
    corr(task, y_obs), corr(task, y_filt))

% Danger check: how much TASK variance does the filter remove?
task_lost = 1 - var(S * task) / var(task - mean(task));
fprintf('Filter removes %.0f%% of task variance\n', 100 * task_lost)
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np
from scipy.linalg import pinv

TR, n, hpf = 2.0, 300, 128.0                      # cutoff in seconds

# Discrete-cosine high-pass basis (SPM's construction)
k = int(np.floor(2 * n * TR / hpf + 1))
t = np.arange(n)
KH = np.column_stack([np.sqrt(2 / n) * np.cos(np.pi * (2 * t + 1) * j / (2 * n))
                      for j in range(1, k)])
S = np.eye(n) - KH @ pinv(KH)                     # residual-forming matrix

task = np.tile(np.r_[np.ones(15), np.zeros(15)], 10)  # 30-s alternating blocks
drift = 4 * np.cos(np.arange(n) * 2 * np.pi / 280)    # slow scanner drift
rng = np.random.default_rng(0)
y_obs = task + drift + rng.standard_normal(n)

y_filt = S @ y_obs                                # drift removed
print(f"corr(task, observed) = {np.corrcoef(task, y_obs)[0,1]:.2f}, "
      f"filtered = {np.corrcoef(task, y_filt)[0,1]:.2f}")

# Danger check: how much TASK variance does the filter remove?
task_c = task - task.mean()
print(f"Filter removes {100 * (1 - (S @ task).var() / task_c.var()):.0f}% "
      "of task variance")
```
:::
::::

With 30-s alternating blocks (60-s period), the 128-s filter leaves the task nearly untouched while removing the drift. Re-run the same code with 64-s blocks (128-s period) and watch the filter remove most of your task variance — the exact failure mode the chapter warns about. The full labs continue with a slice-timing offset demo and a smoothing-kernel tradeoff simulation.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch17-lab-python.ipynb) or download the [MATLAB live script](./labs/ch17_lab_matlab.m), which mirrors it using CANlab tools.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part4/labs/ch17-lab-python.ipynb) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part4/labs/ch17_lab_matlab.m)

## Thought questions

1. Every preprocessing step removes a mixture of artifact and true signal. Pick two steps (e.g., high-pass filtering and physiological noise correction) and describe a concrete experiment in which each would *reduce* your power to detect a true effect. What could you change — in the design or the preprocessing — to shift the balance back?
2. A collaborator proposes excluding all participants with mean FWD above 0.1 mm from a study comparing adolescents with and without ADHD. Drawing on what you know about who moves more and why, explain the selection-bias risks, and propose a defensible alternative motion-handling strategy.
3. Slice-timing correction and motion correction each involve interpolation, and their errors interact — a voxel that moves samples different brain locations at different times. Discuss the arguments for the order of these two steps, and why rapid multiband acquisition changes the calculus (for both task and resting-state analyses).
4. Two labs analyze the same dataset. One normalizes volumetrically to MNI152 with aggressive nonlinear warping; the other uses surface-based normalization with modest smoothing on the cortical sheet. For which brain structures and which scientific questions would you expect their results to diverge most, and why?
5. Your task alternates 45-s blocks of pain and rest, and your software's default high-pass cutoff is 128 s. Walk through what happens to your task effect, how you would detect the problem from the data and design (not just from a null result), and two distinct ways to fix it.

## Quiz yourself

:::{dropdown} **Q1.** What are the three broad aims of fMRI preprocessing?
**Answer:** (1) Minimize the influence of artifacts from acquisition and head movement; (2) check and transform the data to better meet statistical assumptions; and (3) for group analyses, standardize the locations of brain regions across individuals.
:::

:::{dropdown} **Q2.** What are the two main approaches to EPI distortion correction, and what extra data does each require?
**Answer:** Unwarping uses a B0 field map scan to estimate signal displacement from field inhomogeneity and applies the inverse. Reverse-blip methods (e.g., FSL's topup) acquire brief additional EPI scans with opposite phase-encoding directions (A–P and P–A) and interpolate the images to the midpoint of the two opposite distortions.
:::

:::{dropdown} **Q3.** How many parameters does rigid-body realignment estimate, and what are they?
**Answer:** Six: three translations (shifts along x, y, and z) and three rotations (roll, pitch, and yaw). They are estimated by iteratively minimizing the squared difference between each volume and a reference image, and saved for use as nuisance covariates.
:::

:::{dropdown} **Q4.** How is framewise displacement (FWD) computed, and how are rotations handled?
**Answer:** FWD is the sum of the absolute frame-to-frame changes (temporal derivatives) of the six realignment parameters. Rotations are first converted to millimeters by computing arc length on a sphere of radius 50 mm — approximately the distance from the cortex to the center of the brain.
:::

:::{dropdown} **Q5.** Why is minimizing squared intensity differences inappropriate for coregistering structural to functional images, and what is used instead?
**Answer:** Because tissue classes have different intensity orderings in the two modalities — white > gray > ventricles in functional (T2*) images but ventricles > gray > white in structural (T1) images — the same anatomy does not have matching intensities. Coregistration therefore maximizes mutual information between the images instead.
:::

:::{dropdown} **Q6.** Why is spike regression generally preferred over scrubbing for censoring high-motion volumes?
**Answer:** Spike regression adds one nuisance regressor per bad image within the GLM, which removes the volume's influence while preserving the natural temporal correlation structure of the time series — important for valid statistical inference (P values). Scrubbing physically deletes volumes, disrupting temporal structure.
:::

:::{dropdown} **Q7.** Give two reasons fMRI data are spatially smoothed, and the main cost of smoothing.
**Answer:** Smoothing blurs residual anatomical misalignment left over from normalization, and it makes images smooth enough for Gaussian random field theory to give accurate corrected P values. The cost is partial-volume mixing: signal from truly active regions is averaged with null neighbors, reducing effective signal and blurring fine spatial patterns (one reason MVPA studies often skip it).
:::

:::{dropdown} **Q8.** Your design alternates 40-s task and 40-s rest blocks. What happens if you apply SPM's default 128-s high-pass filter, and why?
**Answer:** The task's fundamental frequency (1/80 s = 0.0125 Hz) is close to the filter's 0.0078 Hz cutoff, and designs with alternating blocks of roughly 32 s or longer have substantial power at or below the cutoff — so the filter removes much of the true task-related signal, sharply reducing power (potentially to near zero for longer blocks) regardless of the true effect size.
:::
