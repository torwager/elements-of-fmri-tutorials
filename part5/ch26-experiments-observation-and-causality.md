---
title: "26. Experiments, Observation, and Causality"
subject: "Part 5: Experimental Design"
---

# Experiments, Observation, and Causality

:::{admonition} What you will learn
:class: tip
- The difference between experimentally manipulated and observed variables, and why it determines what causal claims a study can support
- What internal and external validity mean — and why they tend to trade off
- Why randomization licenses causal inference: it makes the independent variable statistically independent of all confounders, measured and unmeasured
- How to read simple causal diagrams (DAGs) and tell confounders, colliders, and mediators apart — so you know when adjusting for a covariate removes bias and when it *creates* bias
- The major bias families identified in clinical trials (selection, attrition, performance, detection, reporting) and practical strategies for avoiding them in fMRI studies
:::

## Overview

Variables in neuroimaging studies come in two flavors: those we **manipulate** and those we **observe**. In task fMRI, the stimuli and instructions we present — faces vs. houses, high vs. low reward, painful vs. non-painful heat — are independent variables (IVs) under experimental control, usually manipulated within-person so that each participant experiences every level. Brain activity at each voxel, behavior, and physiology are dependent variables (DVs) that we observe but do not control. This distinction maps onto two goals that every study shares. **Internal validity** is the ability to make correct causal inferences — to be confident that the activation we see is caused by the manipulated task attribute and not something else. **External validity** is the ability of results to generalize to new contexts, measures, and populations. Frustratingly, the two trade off: tightly controlled laboratory experiments maximize internal validity but test a narrow slice of behavior in a selected sample, while large observational studies (resting-state correlations with clinical traits, naturalistic movie-watching) sample richer contexts and populations but mix together multiple causal influences. The strongest research programs combine both: experiments to establish that an effect is causal, and observational or large-scale studies to establish its size and reach in the real world.

The single most important tool for internal validity is **randomization**. Without it, associations are open to **confounding** — a spurious component of the X–Y association created by a third variable that influences both. The examples are close to home. If most trials of condition A appear early in the run and condition B late, the A − B contrast is confounded with arousal, fatigue, and practice. If the word "Ouch!" appears on the screen only during painful stimulation, visual and semantic processing are confounded with pain. And if participants *choose* which stimulus they receive on each trial, the trial type becomes entangled with novelty, decision processes, and the fluctuating internal states that drove the choice — the observed activation could be caused by any of them. Randomizing the order and timing of conditions breaks all of these links at once, including links to confounders nobody thought of.

The same logic underlies randomized controlled trials (RCTs) in medicine, and their history holds a cautionary tale. Large observational studies in the 1980s found that women taking hormone replacement therapy (HRT) had roughly *three-fold lower* risk of coronary heart disease, and HRT was widely recommended for prevention; some 15 million American women were taking it by 2001. When a large RCT (the Women's Health Initiative) finally reported in 2002, women randomized to HRT showed a *29% increase* in heart disease. The observational association had the wrong sign because of **healthy-user bias**: the kind of person who elected to take HRT also exercised, adhered to medications, smoked less, and was healthier in ways that lowered heart disease risk on their own. The causal structure is easiest to see as a diagram:

```{mermaid}
flowchart LR
    subgraph OBS["Observational study"]
        H(("General health<br/>latent, unmeasured")) -- "a +" --> T["Taking HRT"]
        H -- "c −" --> D["Heart disease"]
        T -- "b +" --> D
    end
    subgraph RCT["Randomized trial"]
        R[/"Random<br/>assignment"/] --> T2["Taking HRT"]
        T2 -- "b +" --> D2["Heart disease"]
    end
    style H fill:#fef3c7,stroke:#d97706,color:#78350f
    style R fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
```

In the observational world, the raw HRT–disease association blends the true causal path $b$ (harmful) with a spurious path through the latent confounder — healthy people both take HRT (path $a$) and avoid heart disease (path $c$) — and the spurious component dominates, flipping the sign. Randomization severs every arrow *into* the treatment, leaving only $b$. Measured proxies for the confounder (education, medication adherence) can be entered as covariates, and they help exactly to the extent that they capture the latent "healthy user" construct — which they rarely do completely. *(This diagram is redrawn after Figure 26.1 from the book.)*

The diagram language generalizes, and it repays learning because **covariate adjustment is not always benign**. Consider the simple linear system $X = aC + u$ and $Y = bX + cC + e$. Regressing $Y$ on $X$ alone yields a biased estimate of the causal effect:

$$
\mathbb{E}\big[\hat{\beta}_X\big] \;=\; b \;+\; ac\,\frac{\mathrm{Var}(C)}{\mathrm{Var}(X)}
$$

Adjusting for the **confounder** $C$ removes the second term and recovers $b$. But the same operation applied to a different causal structure backfires. A **collider** is a variable caused by *both* X and Y (for example, a composite "data quality" or inclusion score). Conditioning on a collider — whether by entering it as a covariate or by selecting only observations above some threshold — *creates* a spurious X–Y association where none existed. This is why post-hoc selection of trials or participants based on performance, missing data, or head motion deserves real scrutiny. A **mediator** lies on the causal path from X to Y; adjusting for it removes the very effect you may want to measure, converting a *total* effect into a *direct* effect. Same regression, three different structures, three different verdicts:

```{mermaid}
flowchart LR
    subgraph S1["Confounder - adjust"]
        C1((C)) --> X1[X]
        C1 --> Y1[Y]
        X1 -- "?" --> Y1
    end
    subgraph S2["Collider - do NOT adjust"]
        X2[X] --> S((S))
        Y2[Y] --> S
        X2 -. "?" .-> Y2
    end
    subgraph S3["Mediator - depends on question"]
        X3[X] --> M((M))
        M --> Y3[Y]
        X3 -. "direct" .-> Y3
    end
    style C1 fill:#dcfce7,stroke:#22c55e,color:#14532d
    style S fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style M fill:#fef3c7,stroke:#d97706,color:#78350f
```

Clinical trialists have catalogued the ways studies go wrong, and every one of them applies to fMRI. **Selection bias**: the study group differs systematically from the intended population (the healthy-user bias is one form), limiting generalization and — when selection differs across treatment arms — internal validity too. **Attrition bias**: dropout that differs across groups and relates to the outcome, as when control participants who feel no benefit leave the study. **Performance bias**: groups differ in exposures other than the treatment, such as extra follow-up care for the drug arm. **Detection bias**: outcomes are measured differently across groups. **Reporting bias**: what gets analyzed, published, and cited depends on the results ("file-drawer" effects, outcome switching, spin). Well-designed RCTs counter these with strict randomization, blinding of participants and assessors, identical procedures across arms, intent-to-treat analysis, checks that randomization actually balanced the groups (chance failures are common when group sizes fall below ~50), pre-registration of a primary outcome, and public registries. Neuroimaging is increasingly adopting the same practices — and one striking demonstration of why: in a meta-analysis of antidepressant trials, drug–placebo differences shrank toward zero when side-effect-driven *unblinding* was controlled, showing how expectation alone can manufacture a "treatment effect."

For fMRI specifically, three strategies eliminate or manage confounds. **Matching or stratification**: design the study so IVs cannot be explained by known alternatives — stratify trial types across time blocks so fatigue and habituation are orthogonal to task contrasts, and match groups on age, sex, and baseline scores. **Independent manipulation of the suspected confound**: if pain-evoked activation might really reflect general negative emotion, manipulate negative emotion without pain and show the region is insensitive to it. **Statistical control**: enter what you cannot manipulate — head motion regressors, spike artifacts, physiological noise — as covariates, with the DAG lessons above in mind. Finally, keep the big causal picture in view. In a typical task fMRI study, the task is randomized but the brain is only *observed*: the design supports the claim that the task caused the activation, but the path from brain activity to behavior remains correlational, because neither brain activity nor performance was assigned at random. This **brain-as-mediator** framing — randomized X, observed brain M, observed outcome Y — is the natural bridge from experimental design to mediation analysis, and 'triangulating' effects across methods (invasive recordings, EEG/MEG, brain stimulation) is what ultimately turns correlated pathways into causal claims.

```{mermaid}
flowchart LR
    X["Task X<br/>(randomized)"] -- "a" --> M["Brain response M<br/>(observed)"]
    M -- "b" --> Y["Behavior Y<br/>(observed)"]
    X -- "c' (direct)" --> Y
```

## Hands-on tutorial

Causal structure is one of the few things in statistics you can *see* directly in simulation, because you know the ground truth. Below we build the two central demonstrations in a few lines each; the full labs extend them with replications, selection effects, and a mediation preview.

**Step 1 — A confounder: adjusting fixes the bias.** We recreate the HRT story in miniature: $C$ (think "general health") drives both the exposure $X$ and the outcome $Y$. The true causal effect of $X$ on $Y$ is $+0.4$, but the naive regression gets the sign wrong — and adjusting for $C$ recovers the truth.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Adapted from CANlab FMRI_simulations (covariates_in_regression.m)
rng(26);
n = 1000;
C = randn(n, 1);                    % confounder ("general health")
X = 0.7*C + randn(n, 1);            % exposure, influenced by C
Y = 0.4*X - 1.2*C + randn(n, 1);    % true causal effect of X = +0.40

naive    = fitlm(X, Y);             % omits the confounder
adjusted = fitlm([X C], Y);         % adjusts for the confounder

fprintf('True effect: 0.40 | naive: %.2f | adjusted: %.2f\n', ...
    naive.Coefficients.Estimate(2), adjusted.Coefficients.Estimate(2));
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np, statsmodels.api as sm

rng = np.random.default_rng(26)
n = 1000
C = rng.standard_normal(n)                      # confounder ("general health")
X = 0.7 * C + rng.standard_normal(n)            # exposure, influenced by C
Y = 0.4 * X - 1.2 * C + rng.standard_normal(n)  # true causal effect of X = +0.40

naive    = sm.OLS(Y, sm.add_constant(X)).fit()                        # omits C
adjusted = sm.OLS(Y, sm.add_constant(np.column_stack([X, C]))).fit()  # adjusts

print(f"True effect: 0.40 | naive: {naive.params[1]:.2f} | "
      f"adjusted: {adjusted.params[1]:.2f}")
```
:::
::::

**Step 2 — A collider: adjusting *creates* bias.** Now $X$ and $Y$ are truly unrelated, but both feed into a downstream composite $S$ (think of an inclusion score built from performance and data quality). Left alone, the regression correctly finds nothing; "controlling for" $S$ manufactures a strong spurious effect.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
X = randn(n, 1);
Y = randn(n, 1);                    % true causal effect of X = 0
S = X + Y + randn(n, 1);            % collider: caused by BOTH X and Y

naive    = fitlm(X, Y);             % correct: no effect
adjusted = fitlm([X S], Y);         % adjusting for the collider = biased!

fprintf('True effect: 0.00 | naive: %.2f | adjusted: %.2f\n', ...
    naive.Coefficients.Estimate(2), adjusted.Coefficients.Estimate(2));
```
:::
:::{tab-item} Python
:sync: python

```python
X = rng.standard_normal(n)
Y = rng.standard_normal(n)              # true causal effect of X = 0
S = X + Y + rng.standard_normal(n)      # collider: caused by BOTH X and Y

naive    = sm.OLS(Y, sm.add_constant(X)).fit()                        # correct
adjusted = sm.OLS(Y, sm.add_constant(np.column_stack([X, S]))).fit()  # biased!

print(f"True effect: 0.00 | naive: {naive.params[1]:.2f} | "
      f"adjusted: {adjusted.params[1]:.2f}")
```
:::
::::

The naive estimate in Step 2 hovers near zero while the "adjusted" one lands around −0.5 with an impressive t-statistic — a completely spurious finding produced by the adjustment itself. The full labs go on to (1) contrast an observational study with an RCT on the same simulated population, (2) repeat both demonstrations across hundreds of replications, (3) show that *selecting* observations on a collider (as in post-hoc trial exclusion) biases results just like covarying for one, and (4) preview mediation analysis on a simulated task → brain → behavior chain.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch26-lab-python.ipynb) or download the [MATLAB live script](./labs/ch26_lab_matlab.m), which mirrors it and connects to the CANlab Mediation Toolbox.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part5/labs/ch26-lab-python.ipynb) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part5/labs/ch26_lab_matlab.m)

## Thought questions

1. Naturalistic movie-watching designs expose participants to rich, engaging stimuli, but close-ups of faces tend to occur at emotionally charged moments in the plot. Analyze this design in terms of internal and external validity, and propose a hybrid design that adds an experimental manipulation to recover some internal validity without destroying the naturalism.
2. In a pain study, participants move their heads more during high-intensity stimulation. Consider three causal readings of head motion — confounder, mediator (motion-induced signal change caused by the painful event), and part of a collider (trials excluded when motion exceeds a threshold). What does including motion regressors, or excluding high-motion trials, do to your pain contrast under each reading?
3. A study compares brain connectivity in medicated depressed patients vs. healthy controls and finds group differences. Drawing on the HRT example, list at least three variables that could play the role of "general health" here, say which ones can be measured as proxies, and explain what covariate adjustment can and cannot fix when the true confounder is latent.
4. A colleague argues, "Our task was fully randomized, so our finding that amygdala activity drives self-reported fear is causally established." Which link in the task → brain → behavior chain does randomization actually secure, which does it not, and what converging evidence would strengthen the brain → behavior claim?
5. Pre-registration and registries were developed to combat reporting bias in clinical trials. Which specific analytic decisions in a typical fMRI study (first-level modeling, trial exclusion, ROI selection, contrast choice) create analogous opportunities for bias, and which would you commit to in advance?

## Quiz yourself

:::{dropdown} **Q1.** What is the difference between internal and external validity, and how are they related?
**Answer:** Internal validity is the ability to make correct causal inferences about study variables; external validity is the ability of the results to generalize to new contexts, measures, and populations. They tend to trade off: tight laboratory control boosts internal validity but narrows the sample and context, while diverse observational settings boost external validity but admit confounding.
:::

:::{dropdown} **Q2.** Why does randomizing the levels of an independent variable license causal inference?
**Answer:** With random assignment (and adequate sample size), the IV is statistically independent of all other potential causes of the outcome — measured and unmeasured alike. Any systematic difference in the outcome across levels can therefore be attributed to the manipulation rather than to confounders.
:::

:::{dropdown} **Q3.** What is a confound? Give an example from task fMRI.
**Answer:** A confound is a third variable that creates a spurious component of the association between the IV and DV by influencing both (or by being associated with the IV and causing the DV). Example: presenting condition A trials early in a run and condition B trials late confounds the A − B contrast with arousal, fatigue, and practice effects.
:::

:::{dropdown} **Q4.** In the hormone replacement therapy story, why did observational studies and the randomized trial reach opposite conclusions?
**Answer:** Observational studies suffered from healthy-user (selection) bias: women who elected to take HRT were healthier in many ways that independently lowered heart disease risk, making HRT look protective. Randomization removed the influence of the "healthy user" confounder on treatment status, revealing HRT's true causal effect — a 29% *increase* in heart disease.
:::

:::{dropdown} **Q5.** What is a collider, and what happens if you adjust for one?
**Answer:** A collider is a variable causally influenced by both X and Y. Conditioning on it — entering it as a covariate or selecting observations based on it — induces a spurious association between X and Y. Unlike a confounder, a collider should be left alone; adjustment *creates* bias rather than removing it.
:::

:::{dropdown} **Q6.** What is attrition bias, and how does intent-to-treat analysis help?
**Answer:** Attrition bias arises when dropout (or exclusion) differs across groups and is related to the outcome — for example, control participants who feel no benefit leaving the study, inflating the apparent treatment effect. Intent-to-treat analysis includes all participants in the groups to which they were originally randomized, preserving the balance created by randomization.
:::

:::{dropdown} **Q7.** Name the three basic strategies for eliminating confounds in fMRI studies.
**Answer:** (1) Matching or stratification — designing so IVs are orthogonal to known alternative causes (e.g., stratifying trial types across time blocks; matching groups on age and sex). (2) Independent manipulation of a suspected confound to test whether it can produce the effect. (3) Statistical control — entering unavoidable nuisance variables (motion, spikes, physiological noise) as covariates.
:::

:::{dropdown} **Q8.** In a typical task fMRI study, which causal link is experimentally secured and which is not — and why does this make the brain a "mediator"?
**Answer:** Randomizing the task secures the task → brain link: activation differences can be attributed to the manipulation. But brain activity and behavior are both merely observed, so the brain → behavior link is correlational. The brain sits between the randomized cause and the behavioral outcome — a mediator (X → M → Y) whose a-path is experimental but whose b-path requires converging evidence (e.g., stimulation, lesions, multi-method triangulation).
:::
