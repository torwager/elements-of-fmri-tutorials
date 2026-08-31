---
title: "36. Granger Causal Models"
subject: "Part 6: Brain Connectivity"
---

# Granger Causal Models

:::{admonition} What you will learn
:class: tip
- What it means for one time series to "Granger cause" another, and how vector autoregressive (VAR) models formalize temporal precedence as evidence of directed influence
- How Geweke's decomposition splits the total dependence between two series into two directed influences plus an instantaneous component
- How Granger causality differs from SEM and DCM: exploratory rather than confirmatory, with no a priori structural model required
- Why regional differences in hemodynamic response latency can create *spurious* Granger causality in BOLD data — and how deconvolution has been proposed as a remedy
- When Granger causality is, and is not, an informative tool for fMRI connectivity analysis
:::

## Overview

Granger causality mapping uses multivariate time series to infer *directed* connectivity from **temporal precedence**: if knowing the past of region $X$ improves prediction of the present of region $Y$ — beyond what $Y$'s own past already provides — we say that $X$ *Granger causes* $Y$. The idea originated in economics, where Clive Granger proposed it as a pragmatic, testable stand-in for causality. Unlike structural equation models (Chapter 34) and dynamic causal models (Chapter 35), Granger causality does not require you to specify a structural model of which regions connect to which in advance. It simply asks, for any pair (or set) of regions, whether one series' history carries predictive information about the other. That makes it more exploratory than confirmatory — potentially most useful in the earlier stages of scientific inquiry, when you do not yet have strong hypotheses about network structure.

The formal machinery is the **vector autoregressive (VAR) model**. Consider two (possibly multivariate) time series $X_t$ and $Y_t$, assumed stationary. First fit each series on its own past, using models of order $p$:

$$
X_t = \sum_{j=1}^{p} A_j X_{t-j} + \epsilon_t, \qquad
Y_t = \sum_{j=1}^{p} B_j Y_{t-j} + \eta_t
$$

where $\epsilon_t$ and $\eta_t$ are zero-mean white noise with covariance matrices $\Sigma_1$ and $\Sigma_2$, and the coefficient matrices $A_j$ and $B_j$ capture the strength of self-influence at each lag $j$. Next, stack the two series into $Z_t = (X_t^{\top}, Y_t^{\top})^{\top}$ and fit a joint model:

$$
Z_t = \sum_{j=1}^{p} C_j Z_{t-j} + \nu_t, \qquad
\operatorname{Cov}(\nu_t) = \Sigma =
\begin{bmatrix}
\Sigma_{xx} & \Sigma_{xy} \\
\Sigma_{yx} & \Sigma_{yy}
\end{bmatrix}
$$

The joint model lets the current value of each series depend on the past of *both*. If adding the cross-regressive terms significantly improves the fit — that is, if $Y$'s history predicts $X$ after controlling for $X$'s own history, or vice versa — a Granger-causal relationship is inferred. In practice this is a nested-model comparison: an F-test (or likelihood ratio test) of the full model against the restricted one.

Geweke proposed an elegant way to quantify these influences using the innovation covariances. The **total linear dependence** between $X$ and $Y$ decomposes into three parts:

$$
F_{X,Y} = \ln \frac{|\Sigma_1|\,|\Sigma_2|}{|\Sigma|}
        = F_{Y \to X} + F_{X \to Y} + F_{X \cdot Y}
$$

where

$$
F_{Y \to X} = \ln \frac{|\Sigma_1|}{|\Sigma_{xx}|}, \qquad
F_{X \to Y} = \ln \frac{|\Sigma_2|}{|\Sigma_{yy}|}, \qquad
F_{X \cdot Y} = \ln \frac{|\Sigma_{xx}|\,|\Sigma_{yy}|}{|\Sigma|}
$$

$F_{Y \to X}$ exceeds zero when past values of $Y$ improve prediction of the current $X$ (the joint model's innovation variance $|\Sigma_{xx}|$ shrinks below the restricted model's $|\Sigma_1|$), and symmetrically for $F_{X \to Y}$. The third term, $F_{X \cdot Y}$, captures *instantaneous* dependence — shared variance at zero lag that neither history explains. A common summary is the **difference** $F_{X \to Y} - F_{Y \to X}$, used to infer which region's history is the more influential.

Applied to fMRI, Granger causality is often used in a whole-brain, seed-based manner. A **Granger causality map (GCM)** is computed with respect to a single reference region (a seed region or voxel; see Chapter 30): it maps both *sources* of influence on the seed and *targets* of influence from the seed across the entire brain — regions that Granger cause the seed, and regions the seed Granger causes.

The approach has been genuinely controversial for fMRI, and the central critique is hemodynamic. BOLD is an indirect, sluggish measure of neuronal activity, and the hemodynamic response function varies across brain regions and individuals (Chapter 18). If region $A$'s vasculature responds a second faster than region $B$'s, then $A$'s BOLD signal will *temporally precede* $B$'s even when the underlying neural activity is simultaneous — or even when $B$'s neural activity leads. Temporal precedence in the signal may reflect hemodynamic, not neuronal, causes, producing spurious directed influence. Proponents have argued that shorter TRs help; critics have argued that Granger causality mapping should not be applied to BOLD data at all. One proposed middle path is to first *deconvolve* the fMRI time series — removing each region's estimated hemodynamic response to recover latent neural-like signals — and apply Granger analysis to those. Two further practical cautions: the VAR framework assumes **stationarity**, so time series spanning distinct states (e.g., baseline and task blocks) should be partitioned into stationary segments or handled with models that account for state transitions; and standard VAR models are linear, though nonlinear extensions exist (e.g., estimating the AR model in a kernel-defined feature space). Finally, remember the scope of the claim: Granger causality indicates only *whether* activity in one region helps predict activity in another — it makes no statement about *how* that influence occurs.

## Hands-on tutorial

In this tutorial you will build the entire Granger story from scratch on simulated data: first a bivariate VAR(1) system with a genuine directed influence, which Granger tests recover correctly — then a system with perfectly *symmetric* neural coupling whose two regions have different hemodynamic latencies, which produces a confidently wrong directional inference at the BOLD level.

**Step 1 — Simulate a VAR(1) with a true directed influence and test both directions.** Region $X$ drives region $Y$ at lag 1 ($Y_t$ receives $0.4\,X_{t-1}$), with no influence in the reverse direction.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Simulate VAR(1): X -> Y at lag 1, no reverse influence
rng(7);
n = 400;
A = [0.5 0.0; 0.4 0.5];          % X(t) <- 0.5*X(t-1);  Y(t) <- 0.4*X(t-1) + 0.5*Y(t-1)
Z = zeros(n, 2);
for t = 2:n
    Z(t, :) = (A * Z(t-1, :)')' + randn(1, 2);
end
X = Z(:, 1); Y = Z(:, 2);

% Granger test as a nested-model F-test (lag p = 1)
% Restricted: Y(t) ~ Y(t-1);  Full: Y(t) ~ Y(t-1) + X(t-1)
T  = (2:n)';
rr = Y(T) - [ones(size(T)) Y(T-1)] * ([ones(size(T)) Y(T-1)] \ Y(T));
rf = Y(T) - [ones(size(T)) Y(T-1) X(T-1)] * ...
            ([ones(size(T)) Y(T-1) X(T-1)] \ Y(T));
F  = ((rr'*rr - rf'*rf) / 1) / ((rf'*rf) / (numel(T) - 3));
p  = 1 - fcdf(F, 1, numel(T) - 3);
fprintf('X -> Y:  F = %.1f, p = %.2g\n', F, p);   % large F, tiny p
% Swap X and Y in the code above to test Y -> X (F near 0)
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests

rng = np.random.default_rng(7)
n = 400
A = np.array([[0.5, 0.0],      # X(t) <- 0.5 X(t-1)
              [0.4, 0.5]])     # Y(t) <- 0.4 X(t-1) + 0.5 Y(t-1)
Z = np.zeros((n, 2))
for t in range(1, n):
    Z[t] = A @ Z[t-1] + rng.standard_normal(2)
X, Y = Z[:, 0], Z[:, 1]

# grangercausalitytests: does column 2 Granger cause column 1?
res_xy = grangercausalitytests(np.column_stack([Y, X]), maxlag=[1])
res_yx = grangercausalitytests(np.column_stack([X, Y]), maxlag=[1])
F_xy, p_xy, *_ = res_xy[1][0]["ssr_ftest"]   # X -> Y: F ~ 92, p ~ 1e-19
F_yx, p_yx, *_ = res_yx[1][0]["ssr_ftest"]   # Y -> X: F ~ 0,  p ~ 0.99
```
:::
::::

**Step 2 — The HRF confound: equal neural coupling, unequal hemodynamic lags.** Now the neural coupling is perfectly symmetric ($0.3$ in both directions), but region 1's HRF peaks at ~4 s while region 2's peaks at ~7 s. At the BOLD level, Granger analysis confidently reports that the fast-HRF region drives the slow-HRF one — a directed influence that does not exist in the neural dynamics.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Requires SPM12 on the path (for spm_hrf)
TR = 1; n = 1000;
As = [0.4 0.3; 0.3 0.4];                    % symmetric neural coupling
Z = zeros(n, 2);
for t = 2:n
    Z(t, :) = (As * Z(t-1, :)')' + randn(1, 2);
end

h_fast = spm_hrf(TR, [4 16 1 1 6 0 32]);    % HRF peaking early (~4 s)
h_slow = spm_hrf(TR, [7 16 1 1 6 0 32]);    % HRF peaking late  (~7 s)
b1 = conv(Z(:, 1), h_fast); b1 = b1(1:n) + 0.05 * randn(n, 1);
b2 = conv(Z(:, 2), h_slow); b2 = b2(1:n) + 0.05 * randn(n, 1);

% Re-run the nested F-tests from Step 1 on b1, b2:
% "fast -> slow" now shows a much larger F than "slow -> fast",
% even though the neural coupling is exactly symmetric.
```
:::
:::{tab-item} Python
:sync: python

```python
from scipy.stats import gamma

def hrf(t, peak):                            # double-gamma HRF
    h = gamma.pdf(t, peak) - gamma.pdf(t, 16) / 6
    return h / h.max()

t_hrf = np.arange(0, 30, 1.0)               # TR = 1 s
h_fast, h_slow = hrf(t_hrf, 4.0), hrf(t_hrf, 7.0)

As = np.array([[0.4, 0.3], [0.3, 0.4]])     # symmetric neural coupling
Z = np.zeros((1000, 2))
for t in range(1, 1000):
    Z[t] = As @ Z[t-1] + rng.standard_normal(2)

b1 = np.convolve(Z[:, 0], h_fast)[:1000] + 0.05 * rng.standard_normal(1000)
b2 = np.convolve(Z[:, 1], h_slow)[:1000] + 0.05 * rng.standard_normal(1000)

grangercausalitytests(np.column_stack([b2, b1]), maxlag=[1])  # fast -> slow: huge F
grangercausalitytests(np.column_stack([b1, b2]), maxlag=[1])  # slow -> fast: modest F
```
:::
::::

The full labs push both steps further: they compute Geweke's directed-influence measures $F_{X \to Y}$ and $F_{Y \to X}$ from restricted-vs-full residual variances, verify that the *neural* series in Step 2 are symmetric while the *BOLD* series are not, and show that deconvolving each region with its own HRF before testing removes the spurious asymmetry (at a real cost in sensitivity).

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch36-lab-python.ipynb) or download the [MATLAB live script](./labs/ch36_lab_matlab.m), which mirrors it using SPM's HRF tools.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch36-lab-python.ipynb) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch36_lab_matlab.m)

## Thought questions

1. Granger causality and DCM occupy opposite ends of an exploratory–confirmatory spectrum. Sketch a research program on, say, fear learning in which each method would be used at a different stage, and explain what each contributes that the other cannot.
2. The core critique of Granger causality for fMRI is that temporal precedence in BOLD may be hemodynamic rather than neuronal. Suppose you observed that region A Granger causes region B during a task but not during rest, with identical scanning parameters. Does the task-vs-rest contrast rescue the causal interpretation? What confounds survive, and what additional data would you want?
3. Deconvolution promises to remove regional HRF differences before Granger analysis, but it requires estimating each region's HRF and amplifies noise. Under what combinations of TR, scan duration, and expected neural lag would you judge deconvolution-based Granger analysis worth attempting, and when would you conclude the enterprise is hopeless?
4. The VAR framework assumes stationarity, yet most task fMRI deliberately alternates between states. Describe two concrete strategies from the chapter for handling a blocked task design, and discuss what each one sacrifices.
5. "X Granger causes Y" is a statement about prediction, not mechanism. Give an example — from neuroscience or elsewhere — where a variable robustly Granger causes another with no plausible direct influence, and identify the general causal structure (e.g., a common driver with unequal delays) that produces such cases.

## Quiz yourself

:::{dropdown} **Q1.** In one sentence, what does it mean for time series $X$ to Granger cause time series $Y$?
**Answer:** Past values of $X$ improve the prediction of the current value of $Y$, over and above what $Y$'s own past values already predict.
:::

:::{dropdown} **Q2.** What class of statistical model is used to test for Granger causality, and what key assumption about the time series does it make?
**Answer:** Vector autoregressive (VAR) models, in which each series is regressed on lagged values of itself and (in the full model) the other series. The framework assumes the series are stationary — constant mean and covariance over time.
:::

:::{dropdown} **Q3.** How does Granger causality differ from SEM and DCM in what the analyst must specify in advance?
**Answer:** SEM and DCM are confirmatory: they require a priori specification of a structural model (which regions are connected), typically comparing a few candidate models. Granger causality requires no structural model — it asks directly whether one region's history predicts another's — making it more exploratory.
:::

:::{dropdown} **Q4.** In Geweke's framework, the total linear dependence $F_{X,Y}$ decomposes into three terms. Name them.
**Answer:** The directed influence of $X$ on $Y$ ($F_{X \to Y}$), the directed influence of $Y$ on $X$ ($F_{Y \to X}$), and the instantaneous influence $F_{X \cdot Y}$ — dependence at zero lag not explained by either history.
:::

:::{dropdown} **Q5.** What is a Granger causality map (GCM)?
**Answer:** A whole-brain map computed with respect to a single reference (seed) region, showing both regions whose activity Granger causes the seed (sources) and regions the seed Granger causes (targets).
:::

:::{dropdown} **Q6.** Why can two regions with identical, simultaneous neural activity still show a significant directed Granger influence in BOLD data?
**Answer:** Because the hemodynamic response function varies across regions. If one region's vascular response is faster, its BOLD signal temporally precedes the other's, and the Granger test attributes this hemodynamic lag difference to directed neural influence — a spurious result.
:::

:::{dropdown} **Q7.** What two remedies have been proposed for the hemodynamic confound in fMRI Granger analysis?
**Answer:** Sampling faster (shorter TRs), which proponents argue mitigates the problem, and deconvolving each region's estimated HRF from its time series first, so that Granger analysis is applied to reconstructed neural-like signals rather than raw BOLD.
:::

:::{dropdown} **Q8.** Your experiment alternates 30-second rest and task blocks. Why is fitting a single VAR across the whole run problematic, and what should you do instead?
**Answer:** The series is non-stationary — its mean and covariance differ between states — violating the VAR assumption. You should either partition the data into stationary segments (e.g., analyze task and rest separately) or use a model that explicitly accounts for transitions between stationary periods.
:::
