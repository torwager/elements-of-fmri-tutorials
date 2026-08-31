---
title: "7. Forward and Reverse Inference"
subject: "Part 2: Brain Mapping"
---

# Forward and Reverse Inference

:::{admonition} What you will learn
:class: tip
- The difference between forward inference, $P(\text{Brain} \mid \text{Psy})$, and reverse inference, $P(\text{Psy} \mid \text{Brain})$ — and why standard brain maps only support the former
- How the diagnostic-testing concepts of sensitivity, specificity, base rate, and positive predictive value (PPV) map onto brain mapping
- How Bayes' rule links forward and reverse inference, and why PPV depends strongly on specificity and base rate but only weakly on sensitivity
- Why "affirming the consequent" is a logical fallacy — and how it sneaks into interpretations of brain activation
- How meta-analytic databases (e.g., Neurosynth) and multivariate predictive models make valid reverse inference possible
:::

:::{admonition} 🖥️ Ways to run this chapter's code
:class: seealso
- **In your browser, no setup:** open the [interactive Python lab](./labs/ch07-lab-python.ipynb) and click the **⏻ power icon** at the top right of the notebook. Run cells top-to-bottom, starting with the first (setup/import) cell.
- **In the cloud:** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part2/labs/ch07-lab-python.ipynb) · [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part2/labs/ch07_lab_matlab.m)
- The code tabs on **this page** are static previews with copy buttons — the labs are where code runs.
:::

:::{div}
:class: run-quick
**Run this code:** [⚡ In-browser lab](./labs/ch07-lab-python.ipynb) · [Colab](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part2/labs/ch07-lab-python.ipynb) · [MATLAB Online](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part2/labs/ch07_lab_matlab.m)
:::

## Overview

Brain mapping is well suited to one kind of question and poorly suited to another. When we induce a psychological state — pain, fear, punishment motivation — and ask which regions respond, we are making a **forward inference**: we learn about the probability of brain activity given a psychological state, written $P(\text{Brain} \mid \text{Psy})$. In the language of diagnostic testing, this is the **sensitivity** (or hit rate) of the brain measure. Standard brain maps are precisely this: maps of forward inference. The stronger the statistical effect, the more likely a region "lights up" when the state is present.

Often, though, what we really want is the opposite direction. If the caudate nucleus is active, can we conclude that a person is experiencing punishment motivation? That is a **reverse inference** — the probability of a mental state given brain activity, $P(\text{Psy} \mid \text{Brain})$ — known in the diagnostic literature as the **positive predictive value (PPV)** of the test. If we could make such inferences reliably, we could "decode" mental states from brain activity: infer whether someone is in pain, engaging cognitive effort, or hiding information.

:::{figure} images/ch07_fig1_forward_reverse_inference.png
:alt: Forward inference maps psychological states to brain activity via sensitivity; reverse inference maps brain activity to psychological states via PPV; a two-by-two table relates hits, misses, false alarms, and correct rejections to sensitivity, specificity, PPV, and NPV
:width: 90%

Forward and reverse inference applied to brain mapping, and their relationship with diagnostic testing measures. "Brain" refers to a brain measure being present (e.g., activation above threshold); "Psy" refers to an underlying psychological state being present. Sensitivity is $P(\text{Brain} \mid \text{Psy})$; specificity is $P(\sim\text{Brain} \mid \sim\text{Psy})$; PPV is $P(\text{Psy} \mid \text{Brain})$. *(Figure 7.1 from the book.)*
:::

The catch is that forward and reverse inference are not interchangeable. A region can respond reliably to a task without being informative about that task, because most brain structures respond to many things — there is a many-to-many mapping between psychological constructs and brain regions. The caudate is strongly activated by reward, cognitive control, motor behavior, and more; even 99% sensitivity to punishment decisions would not make caudate activation good evidence of punishment motivation. Treating it as such is the classical logical fallacy of **affirming the consequent**: if all dogs prefer ice cream over fruit, and Mary prefers ice cream, it does not follow that Mary is a dog. $P(\text{Ice Cream} \mid \text{Dog}) = 1$ says nothing about $P(\text{Dog} \mid \text{Ice Cream})$ until we consider $P(\text{Ice Cream} \mid \sim\text{Dog})$.

**Bayes' rule** makes the relationship exact:

$$
P(\text{Psy} \mid \text{Brain}) = \frac{P(\text{Brain} \mid \text{Psy}) \; P(\text{Psy})}{P(\text{Brain})}
$$

Expanding the denominator over the state being present or absent gives the PPV as a function of three quantities — sensitivity, **specificity** ($P(\sim\text{Brain} \mid \sim\text{Psy})$, i.e., how rarely the region activates in the *absence* of the state), and the **base rate** $P(\text{Psy})$ (the proportion of time the state occurs at all):

$$
\text{PPV} = \frac{\text{Sens} \times P(\text{Psy})}{\text{Sens} \times P(\text{Psy}) + (1 - \text{Spec}) \times \left(1 - P(\text{Psy})\right)}
$$

Plugging in numbers is sobering. Suppose the caudate responds to punishment motivation with 90% sensitivity and 80% specificity, and people experience that state 10% of the time (probably generous). Observing caudate activity then implies punishment motivation with probability of only **33%** — nowhere near the 90% sensitivity. Raising sensitivity to a perfect 100% barely helps (PPV = 36%). But drop the base rate to 1% and the PPV collapses to **4%**. PPV is driven by specificity and base rate, not sensitivity — and standard brain maps measure only sensitivity. The same arithmetic bedevils medical screening: mammography has roughly 90% sensitivity and 80–99% specificity, but with a ~1.5% ten-year base rate of breast cancer in 40–50-year-old women, the PPV of a positive screen is only about 6% in the U.S. (at 80% specificity) versus about 40% in Denmark (at 98% specificity).

Valid reverse inference is possible — it just requires evaluating specificity against the alternatives, not merely showing that a region activates. That means testing a brain measure across many confusable states: task conditions, mental constructs, and behaviors that could produce the same activation. Meta-analytic databases such as Neurosynth aggregate results from thousands of studies, allowing brain patterns to be compared against hundreds of psychological terms and supporting empirically grounded reverse-inference maps. Large multi-task datasets and naturalistic paradigms push in the same direction. And because averaging over a whole structure (the hippocampus, say) mixes many circuits with many functions, multivariate pattern analysis — training algorithms to find the *fine-grained features* that best discriminate a state from thousands of alternatives — is a particularly promising route. Properly validated, such brain measures could let us assess pain, emotion, and cognition in people who cannot report them, and reveal which brain systems treatments actually act on.

## Hands-on tutorial

In this tutorial you will implement Bayes' rule as a reverse-inference calculator, reproduce the chapter's caudate example, and see quantitatively how selective a region must be — and how common a mental state must be — before its activation supports a confident reverse inference. The full labs extend this to a simulated likelihood-ratio "reverse inference map" in the spirit of Neurosynth.

**Step 1 — Compute PPV from sensitivity, specificity, and base rate.** We implement the PPV equation and evaluate the chapter's toy example: 90% sensitivity, 80% specificity, 10% base rate.

:::{note}
The tabs below are **static previews** (with copy buttons) showing the key step in each language. To run and modify this code, use the [interactive in-browser lab](./labs/ch07-lab-python.ipynb) or the Colab / MATLAB Online links above.
:::

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% PPV from sensitivity, specificity, and base rate (Bayes' rule)
ppv = @(sens, spec, br) (sens .* br) ./ ...
      (sens .* br + (1 - spec) .* (1 - br));

% Chapter example: caudate and punishment motivation
fprintf('Sens 0.90, Spec 0.80, BR 0.10 -> PPV = %.2f\n', ppv(.90, .80, .10));
fprintf('Sens 1.00 (perfect!)         -> PPV = %.2f\n', ppv(1.0, .80, .10));
fprintf('Base rate 0.01               -> PPV = %.2f\n', ppv(.90, .80, .01));
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np

def ppv(sens, spec, base_rate):
    """PPV from sensitivity, specificity, and base rate (Bayes' rule)."""
    sens, spec, br = np.asarray(sens), np.asarray(spec), np.asarray(base_rate)
    return sens * br / (sens * br + (1 - spec) * (1 - br))

# Chapter example: caudate and punishment motivation
print(f"Sens 0.90, Spec 0.80, BR 0.10 -> PPV = {ppv(.90, .80, .10):.2f}")
print(f"Sens 1.00 (perfect!)          -> PPV = {ppv(1.0, .80, .10):.2f}")
print(f"Base rate 0.01                -> PPV = {ppv(.90, .80, .01):.2f}")
```
:::
::::

You should get PPV = 0.33, 0.36, and 0.04 — matching the chapter. Perfect sensitivity buys almost nothing; a rarer state destroys the inference.

**Step 2 — How selective must a region be?** We sweep specificity from 0.5 to 1 at several base rates and ask: how specific must activation be for a confident (PPV ≥ 0.9) reverse inference?

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
spec = linspace(.5, .999, 200);
base_rates = [.5 .1 .01];

figure; hold on;
for br = base_rates
    plot(spec, ppv(.90, spec, br), 'LineWidth', 2, ...
        'DisplayName', sprintf('base rate = %g', br));
end
yline(.9, '--', 'PPV = 0.9'); xlabel('Specificity'); ylabel('PPV');
title('Reverse inference needs specificity AND a decent base rate');
legend('Location', 'northwest');
```
:::
:::{tab-item} Python
:sync: python

```python
import matplotlib.pyplot as plt

spec = np.linspace(0.5, 0.999, 200)
fig, ax = plt.subplots(figsize=(6, 4))
for br in [0.5, 0.1, 0.01]:
    ax.plot(spec, ppv(0.90, spec, br), lw=2, label=f"base rate = {br}")
ax.axhline(0.9, ls="--", color="gray")
ax.set(xlabel="Specificity", ylabel="PPV = P(Psy | Brain)",
       title="Reverse inference needs specificity AND a decent base rate")
ax.legend()
```
:::
::::

At a 10% base rate, even 99% specificity yields a PPV of about 0.91 — barely clearing the bar — and at a 1% base rate no realistic specificity suffices for a single region. This is why open-ended tests across many alternative states, and multivariate patterns tuned for discrimination, are central to modern reverse inference.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch07-lab-python.ipynb) or download the [MATLAB live script](./labs/ch07_lab_matlab.m), which mirrors it and adds an optional CANlab/Neurosynth similarity example.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part2/labs/ch07-lab-python.ipynb) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part2/labs/ch07_lab_matlab.m)

## Thought questions

1. A widely cited study reports that the insula activates in 30% of all published fMRI studies, spanning pain, emotion, interoception, salience, and language. A newspaper headline announces that a defendant's insula activation during testimony "reveals he felt genuine remorse." Deconstruct this claim using sensitivity, specificity, and base rate — which quantity does each piece of evidence in the study actually speak to?
2. Neurosynth computes both "forward inference" maps ($P(\text{activation} \mid \text{term})$) and "reverse inference" (association) maps that account for how often activation appears across all other terms. Why can these two maps look dramatically different for the same term (e.g., "pain"), and which regions would you expect to survive in the reverse map but not dominate the forward map?
3. The base rate $P(\text{Psy})$ is straightforward for a disease (prevalence), but what does it even mean for a mental state like "punishment motivation" occurring in daily life versus within a specific experimental task context? How does constraining the context (e.g., "given the participant is doing one of these 10 tasks") change the Bayesian calculation, and is that constrained inference still useful?
4. Multivariate predictive models are said to enable valid reverse inference because they are optimized and *tested* for discrimination among alternatives. What would a convincing specificity evaluation for a "pain signature" look like — which alternative conditions must be tested, and why is averaging activity over a whole anatomical region unlikely to achieve the same specificity?
5. Screening programs (e.g., mammography under age 40) are sometimes abandoned not because the test is insensitive but because the PPV is too low. Propose an fMRI-based "biomarker" use case, estimate plausible sensitivity, specificity, and base rate, and argue whether deployment would be justified.

## Quiz yourself

:::{dropdown} **Q1.** Define forward inference and reverse inference as conditional probabilities.
**Answer:** Forward inference is the probability of brain activity given a psychological state, $P(\text{Brain} \mid \text{Psy})$. Reverse inference is the probability of the psychological state given brain activity, $P(\text{Psy} \mid \text{Brain})$.
:::

:::{dropdown} **Q2.** Which diagnostic-testing quantities correspond to forward and reverse inference?
**Answer:** Forward inference corresponds to sensitivity (hit rate, recall); reverse inference corresponds to the positive predictive value (PPV) of the test.
:::

:::{dropdown} **Q3.** What is specificity, and how is it written in probability notation?
**Answer:** Specificity is the probability that the brain measure does *not* respond when the state is absent: $P(\sim\text{Brain} \mid \sim\text{Psy})$, which equals 1 minus the false alarm rate.
:::

:::{dropdown} **Q4.** Write Bayes' rule relating $P(\text{Psy} \mid \text{Brain})$ to $P(\text{Brain} \mid \text{Psy})$.
**Answer:** $P(\text{Psy} \mid \text{Brain}) = P(\text{Brain} \mid \text{Psy}) \, P(\text{Psy}) / P(\text{Brain})$ — sensitivity times the base rate of the state, divided by the overall probability of observing the brain measure across all states.
:::

:::{dropdown} **Q5.** In the chapter's caudate example (90% sensitivity, 80% specificity, 10% base rate), roughly what is the PPV, and what happens when sensitivity rises to 100%?
**Answer:** The PPV is only about 33%. Raising sensitivity to 100% increases it only to about 36% — PPV is driven by specificity and base rate, not sensitivity.
:::

:::{dropdown} **Q6.** What logical fallacy is committed when strong forward inference is used to claim reverse inference, and what is the ice-cream illustration?
**Answer:** Affirming the consequent. Even if every dog prefers ice cream over fruit ($P(\text{Ice Cream} \mid \text{Dog}) = 1$), someone who prefers ice cream is not therefore a dog — the inference ignores $P(\text{Ice Cream} \mid \sim\text{Dog})$, i.e., how often the observation occurs without the state.
:::

:::{dropdown} **Q7.** Why can't averaging activity across an entire structure like the hippocampus establish that its activation implies memory?
**Answer:** The hippocampus (like most structures) contains many circuits participating in many mental functions, so its average activity is not specific to memory. Establishing reverse inference requires showing specificity against alternative processes — for which fine-grained multivariate patterns, tested across many alternatives, are more promising.
:::

:::{dropdown} **Q8.** Name two developments that make open-ended tests of specificity across many mental states feasible.
**Answer:** Large meta-analytic databases aggregating thousands of studies (e.g., Neurosynth), large multi-task datasets, and naturalistic experiments exposing participants to thousands of images and concepts — any two of these.
:::

:::{div}
:class: book-tile
📖 **The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::

---

[⌂ Back to home](https://torwager.github.io/elements-of-fmri-tutorials/) · [Table of contents](../contents.md) · [How to use this site](../how-to-use.md)
