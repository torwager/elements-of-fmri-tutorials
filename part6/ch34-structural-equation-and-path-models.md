---
title: "34. Structural Equation and Path Models"
subject: "Part 6: Brain Connectivity"
---

# Structural Equation and Path Models

:::{admonition} What you will learn
:class: tip
- How structural equation models (SEMs) and path analysis model the covariance structure among variables using systems of regression equations
- How to read a path diagram and estimate path coefficients for a small network of brain regions
- How mediation analysis decomposes a total effect into direct ($c'$) and indirect ($a \times b$) components, and why the bootstrap is the standard inference for $a \times b$
- How brain-as-mediator designs link experimental manipulations, brain activity, and outcomes in a single model
- Why unmodeled confounding — especially of the mediator–outcome relationship — limits causal claims, and how moderation (PPI) analysis differs from mediation
:::

:::{admonition} 🖥️ Ways to run this chapter's code
:class: seealso
- **In your browser, no setup:** open the [interactive Python lab](./labs/ch34-lab-python.ipynb) and click the **⏻ power icon** at the top right of the notebook. Run cells top-to-bottom, starting with the first (setup/import) cell.
- **In the cloud:** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch34-lab-python.ipynb) · [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch34_lab_matlab.m)
- The code tabs on **this page** are static previews with copy buttons — the labs are where code runs.
:::

:::{div}
:class: run-quick
**Run this code:** [⚡ In-browser lab](./labs/ch34-lab-python.ipynb) · [Colab](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch34-lab-python.ipynb) · [MATLAB Online](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch34_lab_matlab.m)
:::

## Overview

Structural equation modeling is a family of techniques for modeling relationships — more formally, the variance–covariance structure — among a set of variables. Where the GLM predicts one outcome from a set of predictors, an SEM specifies a *system* of equations: multiple predictors, multiple outcomes, and possibly latent factors that are not observed directly but are expressed in the measured variables. The emphasis is on comparing models with different structures and finding one that gives a parsimonious, low-error account of the observed correlations. When all variables are observed (no latent factors), the SEM is called **path analysis**, and the two most widely used inferences in this framework are tests of **mediation** — whether the relationship between two variables is transmitted through a third — and **moderation** — whether one variable changes the strength of the relationship between two others.

Applied to fMRI, SEM is a classic tool for *effective connectivity*: a set of brain regions is chosen a priori, along with a hypothesized set of directed connections among them. The strength of each connection is a **path coefficient** — the expected change in activity in one region per unit change in a region that influences it. Writing $y$ for the vector of regional activities at one time point, the path model is

::::{div}
:class: eq-tip
$$
y = B\,y + \zeta
$$
:::{div}
:class: eq-tip-text
y — vector of regional activities at one time point · B — matrix of directed path coefficients (zeros on the diagonal and for absent connections) · ζ — vector of independent errors
:::
::::
:::{div}
:class: eq-where
*where* $y$ *is the vector of regional activities at one time point,* $B$ *the matrix of directed path coefficients (with zeros on the diagonal and wherever no connection is hypothesized), and* $\zeta$ *a vector of independent errors.*
:::

Rearranging so $y$ appears on only one side, $y = (I - B)^{-1}\zeta$, implies a model covariance matrix

::::{div}
:class: eq-tip
$$
\Sigma(\theta) = (I - B)^{-1}\,\Psi\,(I - B)^{-T}
$$
:::{div}
:class: eq-tip-text
Σ(θ) — model-implied covariance · θ — free parameters (paths and error variances) · B — path matrix · Ψ — error covariance · I — identity matrix
:::
::::
:::{div}
:class: eq-where
*where* $\Sigma(\theta)$ *is the covariance matrix implied by the free parameters* $\theta$ *(the paths and error variances),* $B$ *the path matrix,* $\Psi$ *the covariance of the errors* $\zeta$*, and* $I$ *the identity matrix.*
:::

Estimation minimizes the discrepancy between this model-implied covariance and the sample covariance, typically by maximum likelihood. Model fit is assessed with a $\chi^2$ test comparing the two covariance matrices — and, importantly, the logic is reversed from ordinary hypothesis testing: a *nonsignificant* result means the model is adequate (we cannot reject it), not that it is true or best. Specific edges are tested by comparing nested models with and without the connection using likelihood ratio tests.

:::{figure} images/ch34_fig1_three_roi_sem.png
:alt: Three regions of interest with directed paths b12, b13, b23, and the matrix equation y = By + zeta
:width: 70%
:class: book-figure

A simple recursive three-variable SEM. ROI 1 influences ROIs 2 and 3, and ROI 2 influences ROI 3; the path coefficients between the three nodes appear as elements of the matrix $B$. *(Figure 34.1 from the book. © the authors and MIT Press; reproduced with permission — not covered by this site's CC-BY license.)*
:::

Standard SEM has known limitations for fMRI time series: it assumes observations are independent across time (ignoring temporal autocorrelation, which distorts standard errors), it handles dynamics inflexibly, and it does not incorporate experimental inputs. Extensions such as the *unified SEM*, which couples a path model with a vector autoregressive model, and the *extended unified SEM*, which adds experimental inputs, address these gaps. Group analysis is another subtlety: naively concatenating subjects can yield results inconsistent with single-subject analyses, motivating stacked or multilevel SEMs — or simply t-tests across subject-level path estimates, the familiar summary-statistics approach.

**Mediation analysis** asks whether the effect of an exposure $X$ on an outcome $Y$ is at least partially transmitted through an intervening mediator $M$. The full model is a pair of regressions,

::::{div}
:class: eq-tip
$$
M_i = i_1 + a\,X_i + e_{M,i}
$$
$$
Y_i = i_2 + c'\,X_i + b\,M_i + e_{Y,i}
$$
:::{div}
:class: eq-tip-text
Xᵢ, Mᵢ, Yᵢ — exposure, mediator, outcome for participant i · a — X→M path · b — M→Y path controlling for X · c′ — direct effect of X on Y · i₁, i₂ — intercepts · e — residual errors
:::
::::
:::{div}
:class: eq-where
*where* $X_i$*,* $M_i$*, and* $Y_i$ *are the exposure, mediator, and outcome for participant* $i$*;* $a$*,* $b$*, and* $c'$ *the path coefficients;* $i_1$ *and* $i_2$ *the intercepts; and* $e_{M,i}$*,* $e_{Y,i}$ *the residual errors.*
:::

and the reduced model, without the mediator, is $Y_i = i_3 + c\,X_i + e_i$, where $i_3$ is its intercept and $c$ the total effect. Path $a$ is the effect of exposure on mediator, path $b$ the effect of mediator on outcome controlling for exposure, and $c'$ the **direct effect** — what remains of the $X \to Y$ relationship once $M$ is accounted for. The total effect decomposes exactly:

::::{div}
:class: eq-tip
$$
c = c' + a b
$$
:::{div}
:class: eq-tip-text
c — total effect of X on Y · c′ — direct effect (X→Y controlling for M) · ab — indirect effect transmitted through M
:::
::::
:::{div}
:class: eq-where
*where* $c$ *is the total effect of* $X$ *on* $Y$*,* $c'$ *the direct effect, and* $ab$ *the indirect effect transmitted through the mediator* $M$*.*
:::

so the **indirect (mediated) effect** is the product $a \times b$, and testing mediation means testing $H_0\!: ab = 0$ (equivalently, $c - c' = 0$). The classical Sobel test uses a large-sample normal approximation for the standard error of $\widehat{ab}$, but it is overconservative; modern practice uses **bootstrap tests** — resampling cases, re-estimating $ab$ thousands of times, and forming confidence intervals and P values from the bootstrap distribution.

:::{figure} images/ch34_fig2_mediation_moderation.png
:alt: Panel A shows x to m to y mediation triangle; panel B shows m pointing at the x-to-y arrow, moderation
:width: 80%
:class: book-figure

Mediation and moderation. (A) The three-variable path diagram used in mediation analysis: $m$ partially explains the relationship between $x$ and $y$. (B) Moderation: the value (level) of $m$ changes the strength of the $x$–$y$ relationship. *(Figure 34.2 from the book. © the authors and MIT Press; reproduced with permission — not covered by this site's CC-BY license.)*
:::

Mediation is a natural fit for fMRI because it links experiment, brain, and behavior in one model. A canonical **brain-as-mediator** design: participants are randomized to a stressful or innocuous challenge ($X$), task-evoked activity in an anterior cingulate cortex (ACC) region is the mediator ($M$), and stressor-evoked heart rate increase is the outcome ($Y$). Mediation holds when stress shifts ACC activity (path $a$) and ACC activity tracks heart rate within groups (path $b$) strongly enough that the direct effect $c'$ shrinks relative to $c$. Variables can be organized in other ways too — a brain region mediating the effect of another region on behavior, or one region mediating the link between two other regions. With trial-by-trial data, **multilevel mediation** estimates paths within person and lets path strengths vary (and be moderated) across people, and **mediation effect parametric mapping** searches the whole brain for voxels that act as mediators, analogous to a seed analysis.

:::{figure} images/ch34_fig3_brain_mediation.png
:alt: Stressor to ACC to heart rate path diagram above two scatterplots contrasting mediation and no mediation
:width: 85%
:class: book-figure

Brain-based mediation. Left: ACC activity mediates the stressor–heart rate relationship — the $a$ and $b$ effects account for the group difference, so $c - c' > 0$ and the parallel fit lines nearly coincide. Right: significant $a$ and $b$ effects but little mediation, because $ab$ is small relative to the total effect and a large direct effect $c'$ remains. *(Figure 34.3 from the book. © the authors and MIT Press; reproduced with permission — not covered by this site's CC-BY license.)*
:::

These models are transparent and flexible, but their assumptions deserve respect. The causal reading of $a \times b$ requires, among other things, **no unmodeled confounding** of the mediator–outcome relationship: even when $X$ is randomized, $M$ is merely observed, so any third variable that drives both $M$ and $Y$ (arousal, attention, global signal) can manufacture a spurious "indirect effect." These assumptions are easy to state but hard to verify and often violated, which is why mediation in neuroimaging is best framed as *pathway discovery and description* rather than definitive causal inference — you will create exactly this failure mode in the lab.

Finally, **moderation** asks a different question: does the $X$–$Y$ relationship change with the level of a third variable $M$? It is tested with an interaction term in a standard regression, $Y_i = b_0 + b_1 X_i + b_2 M_i + b_3 (X_i \times M_i) + e_i$, where $b_0$ is the intercept, $b_1$ and $b_2$ the main effects of exposure and moderator, and $b_3$ the interaction coefficient — rejecting $H_0\!: b_3 = 0$ establishes moderation. The time-series version is the widely used **psychophysiological interaction (PPI)** analysis: a seed region's time course, a task variable, and their interaction enter a first-level GLM at every voxel, testing where connectivity with the seed depends on task context. Implementations differ on whether to deconvolve the HRF to form the interaction at the "neural" level ([SPM](https://www.fil.ion.ucl.ac.uk/spm/)) or interact the convolved signals directly (FSL) — deconvolution is only exact when the HRF is known and the system is linear — and generalized PPI (gPPI) models all task conditions with separate interaction terms for greater power in complex designs.

## Hands-on tutorial

In this tutorial you will simulate a brain-as-mediator dataset ($X \to M \to Y$ with a direct path), estimate all four paths, and bootstrap the indirect effect — the core mediation workflow. The MATLAB code uses `mediation.m` from the CANlab Mediation Toolbox; the Python code builds the same analysis from regressions so you can see every moving part.

**Step 1 — Simulate mediation data and estimate the paths.** We generate $n = 200$ observations with true paths $a = 0.6$, $b = 0.5$, $c' = 0.2$, then recover them: path $a$ from the regression of $M$ on $X$, paths $b$ and $c'$ from the regression of $Y$ on $X$ and $M$, and the total effect $c$ from the reduced model.

:::{note}
The tabs below are **static previews** (with copy buttons) showing the key step in each language. To run and modify this code, use the [interactive in-browser lab](./labs/ch34-lab-python.ipynb) or the Colab / MATLAB Online links above.
:::

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Requires CanlabCore + MediationToolbox on your MATLAB path
% Adapted from CANlab mediation tutorials (github.com/canlab/MediationToolbox)
rng(42);                               % fix the random seed for reproducibility
n = 200;                               % n = participants (one observation each)
X = randn(n, 1);                       % exposure (e.g., stressor intensity)
M = 0.6 * X + randn(n, 1);             % mediator (e.g., ACC activity), true a = 0.6
Y = 0.2 * X + 0.5 * M + randn(n, 1);   % outcome (e.g., heart rate), b = 0.5, c' = 0.2

% Estimate paths; 'boot' = bootstrap the indirect effect for inference
% (10,000 samples; use at least ~5,000 for stable tails)
[paths, stats] = mediation(X, Y, M, 'boot', 'verbose', 'bootsamples', 10000, ...
    'names', {'Stressor' 'Heart rate' 'ACC'});

% Columns of paths: [a  b  c'  c  a*b]
mediation_path_diagram(stats);         % path diagram with coefficients and stars
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np

rng = np.random.default_rng(42)                 # fix the random seed for reproducibility
n = 200                                         # n = participants (one observation each)
X = rng.standard_normal(n)                      # exposure (e.g., stressor intensity)
M = 0.6 * X + rng.standard_normal(n)            # mediator (e.g., ACC), true a = 0.6
Y = 0.2 * X + 0.5 * M + rng.standard_normal(n)  # outcome, true b = 0.5, c' = 0.2

def ols(y, *preds):                             # OLS betas, intercept first
    D = np.column_stack([np.ones(len(y)), *preds])
    return np.linalg.lstsq(D, y, rcond=None)[0]

a      = ols(M, X)[1]                           # Path a:   M ~ X
cp, b  = ols(Y, X, M)[1:3]                      # Paths c', b:  Y ~ X + M
c      = ols(Y, X)[1]                           # Path c:   total effect, Y ~ X

print(f"a={a:.3f}  b={b:.3f}  c'={cp:.3f}  c={c:.3f}")
print(f"indirect a*b={a * b:.3f}  check c - c'={c - cp:.3f}")   # identical
```
:::
::::

**Example output:**

```text
a=0.519  b=0.531  c'=0.130  c=0.406
indirect a*b=0.275  check c - c'=0.275
```

The estimates land near the true values ($a = 0.6$, $b = 0.5$, $c' = 0.2$), and $\widehat{c} - \widehat{c'} = \widehat{a}\widehat{b}$ exactly — the decomposition is an algebraic identity for OLS.

**Step 2 — Bootstrap inference for the indirect effect.** The sampling distribution of $\widehat{ab}$ is skewed, so we resample cases with replacement, re-estimate $a \times b$ each time, and read the 95% confidence interval off the bootstrap distribution. `mediation.m` does this automatically with the `'boot'` option; in Python we write the loop ourselves.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% mediation.m already bootstrapped a*b ('boot', 10,000 samples).
% stats fields hold [a  b  c'  c  a*b] estimates and bootstrap P values:
fprintf('a*b = %3.3f, bootstrap p = %3.4f\n', stats.mean(5), stats.p(5));

% The printed table also reports the Sobel-style parametric test --
% compare: the bootstrap test is typically less conservative.
```
:::
:::{tab-item} Python
:sync: python

```python
n_boot = 5000                                   # bootstrap samples; ~5,000+ for stable CI tails
boot_ab = np.empty(n_boot)
for i in range(n_boot):
    idx = rng.integers(0, n, n)                 # resample cases with replacement
    a_i = ols(M[idx], X[idx])[1]
    b_i = ols(Y[idx], X[idx], M[idx])[2]
    boot_ab[i] = a_i * b_i

lo, hi = np.percentile(boot_ab, [2.5, 97.5])
p_boot = 2 * min((boot_ab <= 0).mean(), (boot_ab >= 0).mean())
print(f"a*b = {a * b:.3f}, 95% bootstrap CI [{lo:.3f}, {hi:.3f}], p = {p_boot:.4f}")
```
:::
::::

**Example output:**

```text
a*b = 0.275, 95% bootstrap CI [0.176, 0.386], p = 0.0000
```

(None of the 5,000 bootstrap samples crossed zero, so the two-sided P value prints as 0 — report it as $p < 0.001$.)

The interval should exclude zero decisively — the data were built with real mediation. The full labs push further: they show how an **unmodeled confounder** of the $M$–$Y$ relationship produces a significant "indirect effect" when the true $b$ path is zero (and how adjusting for the confounder repairs it), and they fit the three-ROI path model of Figure 34.1 by regression equations, testing an individual edge by model comparison.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch34-lab-python.ipynb) or download the [MATLAB live script](./labs/ch34_lab_matlab.m), which mirrors it using the CANlab Mediation Toolbox.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch34-lab-python.ipynb) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch34_lab_matlab.m)

## Thought questions

1. In a study where the exposure $X$ is randomized but the brain mediator $M$ is only measured, which of the paths $a$, $b$, and $c'$ inherit a causal interpretation from the randomization, and which do not? What kinds of variables could confound the $M$–$Y$ relationship in a typical fMRI experiment, and what would each do to $\widehat{ab}$?
2. Two SEMs with reversed arrows (ROI 1 → ROI 2 versus ROI 2 → ROI 1) can imply exactly the same covariance matrix and therefore fit the data equally well. What does this "equivalent models" problem mean for claims of effective connectivity based on SEM fit, and what could you add to the study design to break the symmetry?
3. Suppose you find significant $a$ and $b$ paths but the indirect effect $ab$ is small relative to the total effect $c$, as in the right panel of Figure 34.3. A colleague argues the brain region is "clearly involved" and should be called a mediator anyway. How would you adjudicate, and what additional evidence (statistical or experimental) would strengthen or weaken the mediation claim?
4. Multilevel mediation estimates paths from trial-by-trial data within each person, while single-level mediation uses one observation per person. What confounds does each approach protect against, and which is vulnerable to individual-difference confounds like age or head motion? When would the two give different answers?
5. SPM's PPI deconvolves the seed time course to form the interaction at the "neural" level, while FSL interacts the HRF-convolved signals directly. Lay out the assumptions each approach makes, and explain why the choice matters more for event-related than for blocked designs.

## Quiz yourself

:::{dropdown} **Q1.** What is the difference between SEM and path analysis?
**Answer:** Both model relationships among variables with systems of equations, but SEM is the superset: it allows latent (unobserved) factors, correlated errors, and non-recursive relations. Path analysis is the special case in which all variables are observed — as in the three-variable mediation model. Latent variables are seldom used in neuroimaging applications.
:::

:::{dropdown} **Q2.** In the mediation model, what do paths $a$, $b$, $c'$, and $c$ each represent?
**Answer:** Path $a$ is the effect of the exposure $X$ on the mediator $M$; path $b$ is the effect of $M$ on the outcome $Y$ controlling for $X$; path $c'$ is the direct effect of $X$ on $Y$ that remains after accounting for $M$; and path $c$ is the total effect of $X$ on $Y$ from the reduced model without the mediator.
:::

:::{dropdown} **Q3.** Write the decomposition of the total effect, and state what the indirect effect is.
**Answer:** $c = c' + ab$. The indirect (mediated) effect is the product $a \times b$ — the portion of the $X \to Y$ effect transmitted through $M$ — and testing mediation means testing $H_0\!: ab = 0$, which is equivalent to testing $c - c' = 0$.
:::

:::{dropdown} **Q4.** Why has the bootstrap replaced the Sobel test for inference on $a \times b$?
**Answer:** The Sobel test relies on a large-sample normal approximation to the standard error of $\widehat{ab}$, but the product of two coefficients has a skewed sampling distribution, making the Sobel test overconservative (P values larger than necessary). Bootstrapping resamples cases, re-estimates $ab$ many times, and builds confidence intervals and P values from the empirical distribution, respecting the skew.
:::

:::{dropdown} **Q5.** In an SEM overall fit test, what does a *nonsignificant* $\chi^2$ mean — and what does it not mean?
**Answer:** The test compares the model-implied covariance matrix with the observed one, so a nonsignificant result means the model's covariance is not significantly different from the data's — the model is adequate and cannot be rejected. It does not mean the model is true or the best possible one; other models (including equivalent models with different arrows) may fit as well.
:::

:::{dropdown} **Q6.** Give the three ways of organizing variables for a brain mediation analysis described in the chapter.
**Answer:** (1) A brain region mediates the effect of a task manipulation on an outcome (e.g., dmPFC mediating threat effects on heart rate); (2) a brain region mediates the relationship between another brain region and an outcome (e.g., regions mediating vlPFC's effect on reappraisal success); (3) a brain region mediates the relationship between two other brain regions, using time-series/trial data (multilevel) or one contrast value per subject (single-level).
:::

:::{dropdown} **Q7.** Even with randomized $X$, why can a significant $a \times b$ fail to reflect true mediation?
**Answer:** Because the mediator is observed rather than manipulated, the $b$ path (the $M$–$Y$ relationship controlling for $X$) can be created by an unmodeled confounder that influences both $M$ and $Y$. Randomization protects the $a$ and $c$ paths but not $b$, so the causal "no unmodeled confounding" assumption must hold for the indirect effect to be interpreted causally — which is why mediation in fMRI is often framed as pathway discovery rather than causal proof.
:::

:::{dropdown} **Q8.** What is a PPI analysis, and what three regressors define it at the first level?
**Answer:** Psychophysiological interaction analysis is a time-series moderation test of whether connectivity with a seed region depends on task context. The first-level GLM includes the seed region's time course (the physiological variable), the task/psychological variable, and their interaction; the interaction term's parameter map is taken to the group level, and a significant interaction indicates task-modulated connectivity.
:::

:::{div}
:class: book-tile
![Cover of Elements of Functional Magnetic Resonance Imaging](../cover-small.jpg)
**The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::

---

[⌂ Back to home](https://torwager.github.io/elements-of-fmri-tutorials/) · [Table of contents](../contents.md) · [How to use this site](../how-to-use.md)
