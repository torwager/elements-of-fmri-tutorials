---
title: "29. Statistical Power and Sample Size"
subject: "Part 5: Experimental Design"
---

# Statistical Power and Sample Size

:::{admonition} What you will learn
:class: tip
- What effect sizes are (Cohen's $d$, Pearson's $r$), how they differ from statistical significance, and why they drive power analysis
- Why post hoc effect sizes from significant voxels are systematically inflated — the "winner's curse" of selection bias — and how stricter thresholds make it worse
- How to compute statistical power analytically and by simulation, as a function of sample size, effect size, and the multiple comparisons correction
- How to invert a power analysis to find the minimum detectable effect size for a planned sample
- What the BWAS debate implies for individual-differences research, and how to balance the number of participants against scan time per participant
:::

## Overview

Power analysis is a study planning tool: it is performed *before* data collection, to choose a sample size and design that give a good chance of detecting the effects you care about. In fMRI this is complicated because power depends on many factors — the psychological manipulation, task design, acquisition hardware, preprocessing, and analysis choices — and it differs across brain regions. Fortunately, all of these influences converge on a "final common pathway": the **effect size** in the group analysis. If you can specify the effect size you expect, or the minimum effect you would like to be able to detect, then power and sample size calculations become straightforward and can be done with standard tools (or a few lines of code, as in the labs below).

Effect sizes are unit-free descriptions of how large an effect is, independent of sample size. Cohen's $d$ divides a mean effect by its standard deviation — for a one-sample test, $d = \mu / \sigma$ — and benchmark values of $d$ = 0.2, 0.5, and 0.8 are conventionally called small, medium, and large. Pearson's $r$ plays the same role for brain–behavior correlations, and the two are interconvertible:

$$
r = \frac{d}{\sqrt{d^2 + 4}}, \qquad d = \frac{2r}{\sqrt{1 - r^2}}
$$

Crucially, effect size is distinct from statistical significance. Test statistics and P values depend on sample size: with a large enough $N$, arbitrarily small — and practically meaningless — effects become "significant." Effect sizes, by contrast, speak to practical importance, replicability, and the sample sizes future studies will need. For a one-sample t-test the two are linked by $t = d\sqrt{N}$, so an observed t-statistic can be converted back into an effect size estimate ($\hat{d} = t/\sqrt{N}$) — a trick used throughout the labs.

There is a catch, however: typical neuroimaging studies are poorly suited to *estimating* effect sizes. A mass-univariate analysis tests 50,000–350,000 voxels and then reports effect sizes only for the small subset that survived thresholding. Voxels reach significance when they have true signal *and* noise that happens to favor the hypothesis, so the selected voxels' estimated effects are biased upward — like a mediocre golfer reporting only his 10 best holes. This **selection bias**, or "winner's curse," grows worse with more tests, smaller samples, and smaller true effects — and, counterintuitively, *more stringent* correction for multiple comparisons increases the inflation, because only the luckiest noise survives. The same logic underlies P-hacking, HARKing, and circular ("double-dipping") analyses: selecting, weighting, or sorting results in a way that depends on the test outcome short-circuits inference. Remedies include pre-registration and Registered Reports, and — increasingly practical in neuroimaging — discovery–validation designs that lock in a region, network, or multivariate pattern in one sample and test that single effect in an independent sample, yielding unbiased effect size estimates.

:::{figure} images/ch29_fig1_effect_size_inflation.png
:alt: Simulated true signal with d = 0.5, signal plus noise, and 3D surface plots comparing true effect sizes to inflated post hoc effect sizes in significant voxels
:width: 90%

Selection bias in action. Voxels in blue were assigned a true effect of Cohen's $d = 0.5$, and spatially correlated noise was added for each of 30 simulated participants. After a group t-test thresholded at P < 0.001, every significant voxel's estimated effect size (right, "Post Hoc") exceeds the true value (center, "True") — because voxels whose noise favors the hypothesis are selected, and the rest are discarded. *(Figure 29.1 from the book.)*
:::

**Statistical power** is the probability of rejecting the null hypothesis when it is false — of detecting a true effect if it exists. Power increases with the true effect size and with sample size, and decreases as the significance threshold becomes more stringent. This last point is where multiplicity exacts its price: correcting for multiple comparisons across the brain pushes the effective per-test alpha from 0.05 down to roughly 0.001 (a level that often approximates FDR q < 0.05) or to ~4 × 10⁻⁶ (typical familywise error control with permutation tests). The power curves below quantify the cost. Detecting a brain–behavior correlation of $r = 0.5$ with 80% power requires about N = 28 in a single a priori ROI, N = 53 at P < 0.001, and N = 100 or more with FWER correction. For a small effect of $r = 0.1$, the corresponding numbers are N = 781, N = 1,537, and nearly N = 3,000. Comparing two groups roughly quadruples the total sample needed relative to a one-sample test. As a planning anchor, a reasonable expectation for task effects in individual voxels is $d \approx 0.5$ (based on a reference set of Human Connectome Project task activations): detecting it in a single ROI requires N = 34, and with FWER control approximately N = 121 — while a two-group comparison with FWER control requires roughly N = 466.

:::{figure} images/ch29_fig2_power_curves.png
:alt: Four panels showing power as a function of sample size for correlations at P<.05 and P<.001, and number of subjects needed as a function of effect size with FWER correction
:width: 95%

Power curves as a function of effect size and sample size. (A) Power for a single test at P < 0.05 two-tailed: detecting $r$ = 0.1–0.5 with 80% power requires N = 781, 193, 84, 46, and 28. (B) At P < 0.001 (often approximating FDR correction), the required samples roughly double: N = 1,537, 378, 163, 88, and 53. (C, D) Sample size needed for 80% power with permutation-based FWER correction, for one-sample tests (d) and correlations (r). *(Figure 29.2 from the book.)*
:::

These sobering numbers frame the modern **power failure** debate. Underpowered studies not only miss true effects; their "significant" findings are more likely to be false positives, and the effects that do survive are inflated by the winner's curse — a recipe for irreplicable findings. The problem is most acute for **Brain-Wide Association Studies (BWAS)** linking individual differences in brain measures to behavior. Marek, Tervo-Clemmens and colleagues (2022) analyzed tens of thousands of participants and found that the largest brain–behavior correlations in single voxels or edges were around $r = 0.1$ — requiring hundreds to thousands of participants for adequate power. But the story has a constructive coda: with appropriate multivariate models that aggregate signal across the brain, effects up to about $r = 0.4$ are attainable in the same data, cutting required samples roughly 16-fold and bringing individual-differences prediction within reach of samples in the hundreds. Task-based states and decoding analyses can yield much larger effects still (multivariate patterns can exceed $d = 3$), permitting detection in small samples.

Finally, power is not just about the number of participants: the amount and quality of data per participant matter too. With a fixed scanning budget, the optimal allocation depends on the ratio of between-subject variance ($\sigma_B^2$) to within-subject variance ($\sigma_W^2$). Large $\sigma_B^2$ — true individual differences — favors more participants; large $\sigma_W^2$ — noisy within-person measurements — favors longer scans. For typical group analyses, between-subject variance is rate-limiting, so the rule of thumb is: scan as many participants as possible, with at least ~30 minutes of functional time each. An empirical analysis of a working memory study (below) found that 40 total scan hours were best spent on ~38–40 participants scanned for about an hour each. For studies predicting individual-level outcomes, or in homogeneous-population designs like visual psychophysics, deeper scanning of fewer participants can be preferable. When in doubt, larger samples are better — with only 40 subjects and FWER correction, power was around 10% even for a robust working memory contrast.

:::{figure} images/ch29_fig3_scan_time_allocation.png
:alt: Working memory contrast map and power curves for allocating 40 total scan hours across different numbers of subjects, with a pie chart showing roughly balanced within- and between-subject variance
:width: 90%

Balancing scan time and number of subjects. (A) Contrast map for an N-back working memory task. (B) Power as a function of sample size given a fixed budget of 40 total scan hours, for three thresholds. Within- and between-subject variance were roughly balanced (pie chart); the maximum-power solution is ~38 subjects scanned for about one hour each. *(Figure 29.3 from the book.)*
:::

## Hands-on tutorial

In this tutorial you will build a power calculator from scratch and use it to answer the planning questions above: How many participants do I need? What is the smallest effect I can detect? And how badly will the winner's curse inflate my post hoc effect sizes? The calculations mirror the simulation code behind the book's power figures (github.com/canlab).

**Step 1 — Analytic power curves.** Power for a one-sample t-test comes from the noncentral t distribution: with true effect $d$ and sample size $N$, the t-statistic is distributed with noncentrality $\delta = d\sqrt{N}$, and power is the probability that it exceeds the critical value $t_{crit}$. One convention to fix up front: a *planned* test is two-tailed, while thresholds applied to statistic maps ("$p < .001$") are conventionally directional, so we carry a `tails` argument.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Power for a one-sample t-test (requires Statistics Toolbox)
% Adapted from power_calc.m in the CANlab tools (github.com/canlab)
d = 0.5;  alpha = .05;  tails = 2;  N = 2:200;

tcrit = tinv(1 - alpha/tails, N - 1);    % critical t (tails = 1 for maps)
delta = d .* sqrt(N);                    % noncentrality parameter
pow   = 1 - nctcdf(tcrit, N - 1, delta); % power at each N

ncrit = N(find(pow >= .80, 1));          % N needed for 80% power
fprintf('d = %.1f: N = %d for 80%% power at p < %.2f\n', d, ncrit, alpha)

figure; plot(N, pow, 'LineWidth', 2); hold on
plot([2 200], [.8 .8], 'k--')
xlabel('Sample size (N)'); ylabel('Power')
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def power_one_sample(d, n, alpha=0.05, tails=2):
    """Power of a one-sample t-test with true effect size d.
    tails=2 for a planned test; tails=1 for directional map thresholds."""
    df = n - 1
    t_crit = stats.t.ppf(1 - alpha / tails, df)
    return 1 - stats.nct.cdf(t_crit, df, d * np.sqrt(n))

d, alpha = 0.5, 0.05
N = np.arange(2, 201)
pow_curve = power_one_sample(d, N, alpha)

n_crit = N[np.argmax(pow_curve >= 0.80)]     # N needed for 80% power
print(f"d = {d}: N = {n_crit} for 80% power at p < {alpha}")

plt.plot(N, pow_curve); plt.axhline(0.8, ls="--", c="k")
plt.xlabel("Sample size (N)"); plt.ylabel("Power")
```
:::
::::

Repeating this across effect sizes and alpha levels reproduces the book's power curves — and shows the multiplicity tax directly: the same $d = 0.5$ effect that needs N = 34 at P < 0.05 needs roughly N = 122 under whole-brain FWER correction, and about 460 participants if the question is a patient–control difference. Inverted, the same code reconstructs Figure 29.2C/D from scratch: with FWER correction, N = 30 buys 80% power only for $d \approx 1.17$, and N = 100 only for $d \approx 0.56$.

**Step 2 — The winner's curse.** Simulate many "voxels" that all share the same true effect ($d = 0.5$), threshold the map, and compare post hoc effect sizes in significant voxels with the truth.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Simulate 20,000 voxels, all with true d = 0.5, N = 30 subjects
% Concept from effect_size_inflation_example_sim.m (github.com/canlab)
rng(29); n = 30; nvox = 20000; d_true = 0.5;

dat  = d_true + randn(n, nvox);          % subjects x voxels
[~, ~, ~, st] = ttest(dat);              % group t-test at each voxel
p_dir = 1 - tcdf(st.tstat, n - 1);       % directional p, as in mapping
d_hat = st.tstat ./ sqrt(n);             % observed effect sizes

sig = p_dir < .001;
fprintf('True d = %.2f | mean d-hat: all voxels %.2f, significant only %.2f\n', ...
    d_true, mean(d_hat), mean(d_hat(sig)))
```
:::
:::{tab-item} Python
:sync: python

```python
# Simulate 20,000 voxels, all with true d = 0.5, N = 30 subjects
rng = np.random.default_rng(29)
n, nvox, d_true = 30, 20000, 0.5

dat = d_true + rng.standard_normal((n, nvox))   # subjects x voxels
t, _ = stats.ttest_1samp(dat, 0.0)
p_dir = stats.t.sf(t, n - 1)                    # directional p, as in mapping
d_hat = t / np.sqrt(n)                          # observed effect sizes

sig = p_dir < 0.001
print(f"True d = {d_true} | mean d-hat: all voxels {d_hat.mean():.2f}, "
      f"significant only {d_hat[sig].mean():.2f}")
```
:::
::::

Every voxel has the same true effect, yet the significant subset averages $\hat{d} \approx 0.76$ — about 50% inflation, and the bias worsens at stricter thresholds or smaller N (at $p < 4\times10^{-6}$ with N = 15, the surviving voxels average $\hat{d} \approx 2$, four times the truth). The full labs go further: power for correlations via the Fisher z transform, one- and two-group sample sizes under Bonferroni/FDR/FWER-level alphas, minimum detectable effect sizes for a planned N, and a reanalysis of the univariate-vs-multivariate power gap highlighted by the BWAS debate.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch29-lab-python.ipynb) or download the [MATLAB live script](./labs/ch29_lab_matlab.m), which mirrors it using CANlab-style power utilities.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part5/labs/ch29-lab-python.ipynb)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part5/labs/ch29_lab_matlab.m)

## Thought questions

1. Stricter multiple comparisons correction reduces false positives but *increases* the inflation of post hoc effect sizes among surviving voxels. Explain why both statements are true simultaneously, and what this implies for using published whole-brain results to plan a new study's sample size.
2. A colleague argues that small-sample studies are self-correcting: "If we found it with N = 15, the effect must be huge, so it's the effects worth caring about." Using the concepts of power, confidence interval width, and the winner's curse, construct the counterargument — and identify any circumstances in which your colleague might have a point (e.g., multivariate task decoding with $d > 3$).
3. The Marek et al. (2022) BWAS findings ($r \approx 0.1$ for univariate brain–behavior correlations) triggered calls for consortium-scale samples, yet multivariate models achieve $r \approx 0.4$ in the same data. How should these two facts jointly shape the design of an individual-differences study with a budget for 300 participants? Consider what is gained and lost by moving from voxel-level maps to predictive patterns.
4. You have 60 hours of scanner time for a group study of a cognitive task. Using the concepts of within-subject variance ($\sigma_W^2$) and between-subject variance ($\sigma_B^2$), describe how you would decide between 60 participants × 1 hour and 20 participants × 3 hours — and how your answer would change if the goal were instead to build a within-person predictive model, or to study a rare patient population.
5. Pre-registration and discovery–validation designs both aim to eliminate selection bias, but they make different tradeoffs between analytic flexibility and inferential rigor. For a new fMRI study of emotion regulation in adolescents, sketch a workflow that preserves the ability to make data-driven methodological choices while still delivering an unbiased effect size estimate for the primary hypothesis.

## Quiz yourself

:::{dropdown} **Q1.** What is statistical power?
**Answer:** Power is the probability of rejecting the null hypothesis when it is false — that is, the probability of detecting a true effect if it exists. It increases with larger true effect sizes and larger samples, and decreases with more stringent significance thresholds.
:::

:::{dropdown} **Q2.** How does an effect size like Cohen's $d$ differ from a t-statistic or P value?
**Answer:** Cohen's $d$ is a unit-free, sample-size-independent description of an effect's magnitude (mean effect divided by standard deviation). Test statistics and P values depend on sample size and reflect statistical significance, not magnitude — with a large enough N, even tiny, practically unimportant effects become significant. The two are related by $t = d\sqrt{N}$ for a one-sample test.
:::

:::{dropdown} **Q3.** What is the "winner's curse" in neuroimaging, and what three factors make it worse?
**Answer:** When effect sizes are reported only for voxels that survived thresholding, the estimates are biased upward, because voxels are selected partly for having favorable noise. The bias grows with (1) more tests (voxels) conducted, (2) smaller samples, and (3) smaller true effects. Stricter statistical thresholds also increase the inflation among surviving voxels.
:::

:::{dropdown} **Q4.** Roughly how many participants are needed to detect a brain–behavior correlation of r = 0.5 with 80% power in (a) one a priori ROI at P < 0.05, (b) at P < 0.001, and (c) with whole-brain FWER correction?
**Answer:** Approximately (a) N = 28 for a single ROI, (b) N = 53 at P < 0.001 (a level that often approximates FDR correction), and (c) N = 100 or more for FWER control. Requirements grow dramatically for smaller effects — r = 0.1 needs N = 781 for even a single test.
:::

:::{dropdown} **Q5.** What is the rule of thumb for sample sizes when comparing two groups (e.g., patients vs. controls) rather than testing one group?
**Answer:** You need approximately 4 times the total sample size for a between-group comparison. For example, an effect of d = 0.4 detectable with N = 50 in a one-sample test requires about N = 101 per group (~200 total) for 80% power in a two-group comparison.
:::

:::{dropdown} **Q6.** Why are pre-registration and discovery–validation approaches effective against selection bias?
**Answer:** Both lock down which effect is tested and how, before the confirmatory test is run. Testing a single pre-specified effect (a region, network, or multivariate pattern) in an independent sample eliminates the dependence between selection criteria and test outcome (circularity), so the false positive rate is controlled and the replication-sample effect size is an unbiased estimate of the true effect.
:::

:::{dropdown} **Q7.** What did the Marek et al. (2022) BWAS analyses find, and what is the multivariate counterpoint?
**Answer:** Across large-scale datasets, the largest univariate brain–behavior correlations were only about r = 0.1, implying that hundreds to thousands of participants are needed for reliable voxel- or edge-level individual-differences findings. However, multivariate models aggregating signal across the brain can reach about r = 0.4 in the same data — reducing required sample sizes roughly 16-fold and making prediction feasible with samples in the hundreds.
:::

:::{dropdown} **Q8.** With a fixed total scanning budget, when should you favor more participants versus more scan time per participant?
**Answer:** The balance depends on the ratio of between-subject to within-subject variance. When between-subject variance dominates (typical for group analyses of cognitive and affective tasks), add participants — with at least ~30 minutes of functional time each. When within-subject variance dominates, or for predicting individual-level outcomes, longer or repeated scans per person become more valuable. Empirically, one working memory study found ~40 hours best allocated as ~1 hour each for ~38–40 participants.
:::
