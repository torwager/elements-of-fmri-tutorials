---
title: "31. Multivariate Decomposition: PCA and ICA"
subject: "Part 6: Brain Connectivity"
---

# Multivariate Decomposition: PCA and ICA

:::{admonition} What you will learn
:class: tip
- How decomposition methods factor a time $\times$ voxels data matrix into spatial maps and time courses, turning ~100,000 voxels into a handful of latent components
- How PCA and the singular value decomposition (SVD) find orthogonal components ordered by variance explained, and how eigenimages and their time courses are read off $U$, $S$, and $V$
- Why ICA's statistical independence criterion — stronger than uncorrelatedness — can unmix sources that PCA blends together, and why non-Gaussianity is the key ingredient
- How to choose how many components to keep (scree elbows, variance thresholds, permutation tests)
- How group ICA with temporal concatenation and dual regression yields subject-specific time courses and spatial maps for group inference
:::

## Overview

Everything up to this point in the book has been largely *univariate*: a model is fit to one voxel at a time, and the results are stitched into maps. Multivariate decomposition methods flip this around. They operate on many variables at once — typically the full $T \times V$ matrix $X$ of one subject's fMRI data, with $T$ time points as rows and $V$ voxels as columns — and factor it into a small number of components, each pairing a **spatial map** (a pattern over voxels) with a **time course**. The components are *latent variables*: not directly measured, but expressed across many measured voxels. Decomposition methods are workhorses of functional connectivity analysis, and ICA in particular is the most widely used method for identifying resting-state networks. The broader family includes factor analysis, multidimensional scaling, and non-negative matrix factorization — the members differ mainly in the constraints they impose (orthogonality vs. independence), whether they model noise explicitly, and how they normalize the data — but PCA and ICA are the bedrock, and the focus here.

**Principal components analysis** reduces the dimensionality of a set of correlated variables while retaining as much variance as possible. It transforms the original variables into new ones — the principal components — that are uncorrelated and ordered by the variance they explain. Formally, PCA eigen-decomposes the data covariance (or correlation) matrix, $\Sigma = Q \Lambda Q^T$, where the columns of $Q$ are eigenvectors defining linear combinations of voxels and the diagonal of $\Lambda$ holds eigenvalues proportional to each component's share of variance. In practice the components are computed via the singular value decomposition of the mean-centered data:

$$
X = U S V^T
$$

where $U$ ($T \times T$) and $V$ ($V \times V$) are orthonormal and $S$ is diagonal with singular values sorted largest to smallest. In fMRI terms, each column of $V$ is an **eigenimage** — a spatial mode capturing covariance structure across voxels — and the corresponding column of $U$ (scaled by its singular value) is that eigenimage's time course. The data are exactly the sum of rank-one layers, $X = \sum_j s_j \mathbf{u}_j \mathbf{v}_j^T$, and truncating the sum after $k$ terms gives the best possible rank-$k$ reconstruction. The proportion of variance explained by component $j$ is $s_j^2 / \sum_i s_i^2$.

:::{figure} images/ch31_fig1_pca_svd_overview.png
:alt: Panel A shows the SVD X = U S V-transpose with time courses in U and eigenimages in V-transpose; panel B shows the data matrix rebuilt as a sum of rank-one outer products, each pairing a time course with a spatial map
:width: 85%

Overview of PCA via the singular value decomposition. (A) The time $\times$ voxels data matrix $X$ is factored as $U S V^T$: $V^T$ holds the eigenimages (spatial maps), $U$ their associated time courses, and $S$ the singular values. (B) Equivalently, $X$ is a sum of rank-one components $s_j \mathbf{u}_j \mathbf{v}_j^T$, each the outer product of a time course and a spatial map, weighted by its singular value. *(Figure 31.1 from the book.)*
:::

Because components are ordered, later ones can often be discarded with little loss, making PCA a natural data-reduction and exploration tool — for spotting outlier volumes in quality control, for finding distributed patterns related to tasks, and as a preprocessing step for ICA and for the predictive models of Part 7. Common rules for choosing how many components to keep: look for an "elbow" where the scree plot of variance explained decelerates; retain enough components to explain a target proportion of variance (e.g., 90%); keep components with eigenvalues greater than 1 (more than one original variable's worth); or compare against decompositions of permuted (randomized) data and keep components that beat chance.

PCA's guarantees come with a catch: the components must be *orthogonal*, and they chase *variance*. Real signal sources — brain networks, physiological artifacts, scanner drift — have no obligation to be orthogonal to one another or to line up with directions of maximal variance. When two sources overlap or correlate, the first PC typically captures a variance-weighted blend of both and later PCs capture orthogonal remainders: the subspace is right, but the axes within it are rotated away from the true sources. **Independent components analysis** addresses exactly this. It models the data as a linear mixture of latent sources,

$$
X = A\,S
$$

where $X$ is $T \times V$, the mixing matrix $A$ ($T \times k$) holds each component's time course in its columns, and the source matrix $S$ ($k \times V$) holds spatially independent maps in its rows. Neither $A$ nor $S$ is observed; ICA seeks an unmixing matrix $W$ such that $\hat{S} = WX$ recovers the sources — blindly, without knowing the mixing process.

:::{figure} images/ch31_fig2_ica_overview.png
:alt: The data matrix X equals the mixing matrix A times the source matrix S; columns of A are component time courses and rows of S are spatially independent maps
:width: 85%

Overview of spatial ICA. The time $\times$ voxels data $X$ is modeled as a mixing matrix $A$ — whose columns are component time courses — times a source matrix $S$ whose rows are statistically independent spatial maps. *(Figure 31.2 from the book.)*
:::

What makes this possible is a stronger criterion than PCA's: statistical **independence**, $p(s_1, s_2) = p(s_1)\,p(s_2)$, based on higher-order properties of the distributions rather than just second moments. Independence implies uncorrelatedness, but the converse holds only for Gaussian signals — uncorrelated non-Gaussian signals can still carry higher-order dependence. ICA exploits two facts: many real-world signals are non-Gaussian, and by the Central Limit Theorem, *mixing* non-Gaussian signals pushes the mixtures toward Gaussianity. So ICA searches for the unmixing that makes the recovered components as non-Gaussian as possible — via algorithms such as infomax (maximizing entropy) or FastICA and JADE (working on kurtosis and related contrasts). At least $k-1$ of the sources must be non-Gaussian; with multiple Gaussian sources the model is identifiable only up to an orthogonal rotation. The flexibility carries costs: components have arbitrary sign, scale, and order, and are not ranked by importance, so researchers must inspect and label them. In fMRI, the standard variant is **spatial ICA** — maximizing independence of the spatial maps, consistent with the idea that different networks involve distinct sets of voxels — rather than temporal ICA, since tasks and artifacts routinely induce correlated time courses. Before ICA, data are typically mean-centered, whitened, and reduced with PCA. Calling ICA components "networks" is something of a misnomer — nothing guarantees that regions loading on a component are interconnected — but the usage is standard, and ICA-derived resting-state networks replicate remarkably well across subjects and sessions. ICA is also widely used for denoising during preprocessing.

In practice ICA is usually run on a **group** of subjects. The most common approach is *temporal concatenation*: stack each subject's $T \times V$ matrix into an $NT \times V$ matrix $Y$, and fit $Y = AS$, so all subjects share one set of group spatial maps $S$ while each subject's block of $A$ ($A^{(i)}$) provides subject-specific time courses. Subject-specific maps are then recovered by *back-reconstruction*: either invert the model per subject (solve $X^{(i)} = A^{(i)} S^{(i)}$ for $S^{(i)}$), or use the two-stage **dual regression**. In stage 1 (spatial regression), each subject's data are regressed onto the group maps, yielding a time course per component — each map's expression at every time point, controlling for the other maps. In stage 2 (temporal regression), the same subject's data are regressed onto those time courses, yielding subject-specific spatial maps. Group inference follows naturally: t-tests on the subject-level maps identify voxels reliably loading on each component, and the subject-level time courses can serve as dependent variables in GLMs relating components to tasks and behavior. Because ICA solutions can vary across studies (noise sensitivity; sign and scale ambiguity), a popular strategy is to run dual regression against *fixed, labeled* template networks from large-scale resting-state studies — trading within-subject independence for comparability — while newer approaches such as Group Information-guided ICA and the NeuroMark pipeline estimate subject-level components constrained to match high-quality group templates. Popular implementations include MELODIC (FSL) and the GIFT toolbox (MATLAB).

## Hands-on tutorial

The best way to understand what PCA and ICA each can and cannot do is to build data where the truth is known. We simulate two source "networks" — each a spatial map paired with a time course — mix them into a time $\times$ voxels matrix with noise, and then ask each method to recover them. The key manipulation: the two spatial maps *overlap*, so the true sources are **not orthogonal** — though, being sparse blocks, they remain close to statistically independent. The labs then push further: scree plots for dimensionality, and a miniature dual regression on a simulated "group."

**Step 1 — Simulate mixed sources.** Two spatial maps (overlapping blocks of voxels) and two event-driven time courses combine as $X = AS + E$.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Requires CanlabCore + SPM12 on your MATLAB path
% Adapted from CANlab tutorials (github.com/canlab)
rng(7);
T = 200; V = 360; TR = 2;

% Two overlapping spatial maps (non-orthogonal sources)
s1 = zeros(1, V); s1(1:150)  = 1;              % "network" 1: voxels 1-150
s2 = zeros(1, V); s2(91:255) = 1;              % "network" 2: voxels 91-255 (overlap!)
S_true = [s1; s2];

% Two event-related time courses via the canonical HRF
a1 = onsets2fmridesign({[20 100 180 260 340]'}, TR, T*TR);
a2 = onsets2fmridesign({[60 140 220 300 380]'}, TR, T*TR);
A_true = [a1(:,1) a2(:,1)];

X = A_true * S_true + 0.3 * randn(T, V);       % mix + noise
s1 * s2' / (norm(s1) * norm(s2))               % cosine ~0.4: maps not orthogonal
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np
from scipy.stats import gamma

rng = np.random.default_rng(7)
T, V, TR = 200, 360, 2.0

# Two overlapping spatial maps (non-orthogonal sources)
s1 = np.zeros(V); s1[0:150]  = 1              # "network" 1: voxels 1-150
s2 = np.zeros(V); s2[90:255] = 1              # "network" 2: voxels 91-255 (overlap!)
S_true = np.vstack([s1, s2])

# Two event-related time courses via a double-gamma HRF
t = np.arange(0, 32, TR)
hrf = gamma.pdf(t, 6) - gamma.pdf(t, 16) / 6
def timecourse(onsets_s):
    stick = np.zeros(T); stick[(np.array(onsets_s) / TR).astype(int)] = 1
    return np.convolve(stick, hrf)[:T]
A_true = np.column_stack([timecourse([20, 100, 180, 260, 340]),
                          timecourse([60, 140, 220, 300, 380])])

X = A_true @ S_true + 0.3 * rng.standard_normal((T, V))   # mix + noise
print(s1 @ s2 / (np.linalg.norm(s1) * np.linalg.norm(s2)))
# cosine ~0.4: the maps overlap, so they are not orthogonal
```
:::
::::

**Step 2 — PCA mixes, ICA unmixes.** PCA (via SVD) captures the right two-dimensional subspace — the scree plot shows two components towering over the noise floor — but its orthogonal axes are rotated blends of the two correlated sources. FastICA, run on the same two-dimensional reduction, recovers maps that match the true sources almost perfectly.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% PCA via SVD of the mean-centered data
Xc = X - mean(X);
[U, Ssv, Vsv] = svd(Xc, 'econ');
var_explained = 100 * diag(Ssv).^2 / sum(diag(Ssv).^2);
figure; plot(var_explained(1:10), 'ko-'); title('Scree plot');

pc_maps = Vsv(:, 1:2)';                        % eigenimages (rows)
disp('|corr| between true sources and PCA maps:')
disp(abs(corr(S_true', pc_maps')))             % PCs blend the sources

% Spatial ICA: the full lab implements FastICA from scratch in ~15 lines;
% with GIFT installed use icatb_fastICA, or CanlabCore's ica() on fmri_data
```
:::
:::{tab-item} Python
:sync: python

```python
from sklearn.decomposition import FastICA

# PCA via SVD of the mean-centered data
Xc = X - X.mean(axis=0)
U, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
var_explained = 100 * sv**2 / np.sum(sv**2)
print("Variance explained by PCs 1-4:", var_explained[:4].round(1))

pc_maps = Vt[:2]                              # eigenimages (rows)
print("|corr| true sources vs PCA maps:")
print(np.abs(np.corrcoef(S_true, pc_maps)[:2, 2:]).round(2))  # blended

# Spatial ICA: voxels are the samples, so decompose the transpose
ica = FastICA(n_components=2, random_state=0, whiten="unit-variance")
ic_maps = ica.fit_transform(Xc.T).T           # k x V independent spatial maps
ic_time = ica.mixing_                         # T x k component time courses
print("|corr| true sources vs ICA maps:")
print(np.abs(np.corrcoef(S_true, ic_maps)[:2, 2:]).round(2))  # ~1.0
```
:::
::::

You should see PCA map–source correlations well below 1 (each PC correlates with *both* sources), while each ICA map correlates near 1.0 with exactly one source. Same subspace, different axes — the independence criterion is what points the axes at the sources. The full labs complete the arc: reconstructing $X$ from rank-one layers, the scree/permutation view of dimensionality, and a mini dual regression that recovers subject-specific time courses and maps from a concatenated "group" dataset.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch31-lab-python.ipynb) or download the [MATLAB live script](./labs/ch31_lab_matlab.m), which mirrors it using CANlab tools.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch31-lab-python.ipynb) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch31_lab_matlab.m)

## Thought questions

1. PCA components are orthogonal and variance-ranked; ICA components are independent, unordered, and of arbitrary sign and scale. For each property, give a concrete fMRI scenario where it is an asset and one where it is a liability.
2. Spatial ICA assumes networks are *spatially* independent, yet the chapter notes that regions loading on one component need not be interconnected, and overlapping "hub" regions belong to multiple networks. How would violations of spatial independence distort the recovered maps, and what might temporal ICA recover instead in the same data?
3. ICA works because mixtures of non-Gaussian sources are more Gaussian than the sources themselves. What kinds of fMRI signals and artifacts are strongly non-Gaussian, and what would happen to an ICA decomposition if the dominant sources of variance were nearly Gaussian?
4. A lab runs group ICA on its own 30-subject study; another performs dual regression against published template networks from a 1,000-subject consortium. Compare the two strategies for (a) sensitivity to this sample's idiosyncratic networks, (b) replicability across studies, and (c) interpretability of a patient-vs-control difference on "the default mode network."
5. Choosing the model order $k$ changes what ICA finds: small $k$ merges networks, large $k$ splits them into subnetworks. Is there a "true" number of brain networks? Defend a strategy for choosing and reporting $k$ in a study linking network expression to a clinical variable.

## Quiz yourself

:::{dropdown} **Q1.** In a decomposition of a time $\times$ voxels fMRI matrix, what two things does each component pair together?
**Answer:** A spatial map (a pattern of weights over voxels — an eigenimage in PCA, an independent map in ICA) and an associated time course describing how strongly that pattern is expressed at each time point.
:::

:::{dropdown} **Q2.** In the SVD $X = USV^T$ of mean-centered fMRI data, where are the eigenimages, their time courses, and the variance information?
**Answer:** The columns of $V$ (rows of $V^T$) are the eigenimages; the columns of $U$ hold the corresponding time courses; and the diagonal of $S$ holds the singular values, whose squares (divided by their sum) give each component's proportion of variance explained.
:::

:::{dropdown} **Q3.** Name three common rules for deciding how many principal components to retain.
**Answer:** Look for an elbow (deceleration) in the scree plot of variance explained; keep enough components to reach a variance threshold such as 90%; keep components with eigenvalues > 1; or use permutation testing to keep components explaining more variance than in randomized data. (Any three.)
:::

:::{dropdown} **Q4.** Why is statistical independence a stronger requirement than uncorrelatedness, and for what kind of signals are the two equivalent?
**Answer:** Independence requires the joint distribution to factorize, $p(s_1,s_2) = p(s_1)p(s_2)$, constraining all higher-order moments, whereas uncorrelatedness only constrains second moments. The two coincide only for Gaussian signals — uncorrelated non-Gaussian signals can still have higher-order dependence, which is exactly what ICA exploits.
:::

:::{dropdown} **Q5.** In the ICA model $X = AS$ for fMRI, what do $A$ and $S$ contain, and what does ICA actually estimate?
**Answer:** $A$ (the mixing matrix) contains each component's time course in its columns; $S$ (the source matrix) contains spatially independent maps in its rows. Neither is observed: ICA estimates an unmixing matrix $W$ such that $\hat{S} = WX$ approximates the sources, without knowledge of the mixing process.
:::

:::{dropdown} **Q6.** Why can ICA fail if more than one source is Gaussian?
**Answer:** Any orthogonal rotation of jointly Gaussian variables has exactly the same distribution, so with multiple Gaussian sources the unmixing matrix is only identifiable up to an orthogonal transformation — there is no unique independent solution to find.
:::

:::{dropdown} **Q7.** Why is spatial ICA preferred over temporal ICA for fMRI?
**Answer:** Spatial ICA matches the assumption that different brain networks engage largely distinct sets of voxels (independent spatial maps), while temporal independence is less plausible — tasks and artifacts routinely produce correlated time courses across components. Spatial ICA also suits fMRI's many-voxels, fewer-time-points geometry.
:::

:::{dropdown} **Q8.** Describe the two stages of dual regression and what each produces.
**Answer:** Stage 1 (spatial regression): regress each subject's data onto the group spatial maps, yielding a subject-specific time course per component, controlling for the other maps. Stage 2 (temporal regression): regress the same subject's data onto those time courses, yielding subject-specific spatial maps. The time courses and maps then support group-level tests relating components to tasks, behavior, or group status.
:::
