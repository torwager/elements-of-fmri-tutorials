---
title: "3. Imaging and Society"
subject: "Part 1: Motivation"
---

# Imaging and Society

:::{admonition} What you will learn
:class: tip
- What biomarkers and *neuromarkers* are, and the different roles they can play in medicine, law, and beyond
- Why inferences about an individual person are far harder than inferences about group averages (the "group to individual," or G2i, problem)
- How base rates, standardization, cost, and access complicate the real-world use of any biomarker — not just brain-based ones
- Why the decision to *adopt* a neuromarker is a societal question that goes beyond the science of validating it
- How to weigh the harms of premature commercialization against the harms of preserving a flawed status quo
:::

## Overview

Imaging technologies already permeate daily life. Satellite imagery originally developed for military purposes now guides us to restaurants and hiking trails, and image-recognition systems have spread into forensics, agriculture, law, entertainment, marketing, and political science. Neuroimaging is following a similar trajectory: it is increasingly invoked in clinics, courtrooms, and commercial products. That makes it crucial to understand what brain images can — and cannot — tell us about the mind, disease, and human performance. Used with appropriate, disciplined inference, imaging can be a powerful tool for improving the human condition. Used to overreach, it joins a long historical tradition of biological observations pressed into service to "sell" ideas, treatments, and training programs, causing harm and muddying our understanding of which measures are valid and when.

In medicine, imaging is routine, and its value rests on **biomarkers**: biological measures that serve as indicators of disease-relevant physiological processes. Some biomarkers are qualitative — a bright region on a chest CT, read by an expert radiologist, revealing fluid-filled lungs that need immediate draining. Others are quantitative — echocardiography turns moving images of the heart into a numeric estimate of left ventricular ejection fraction, an actionable indicator of heart failure. One of the book's authors experienced both kinds firsthand when a chest CT, ordered just before his wife was to be discharged from an emergency room days after childbirth, revealed a life-threatening condition that physical examination had missed. Imaging succeeds most clearly when pathology is visible to a trained eye and points to a clear course of action, or when a well-validated procedure converts measurements into a clinically meaningful score.

Brain-based biomarkers — variously called **neuromarkers**, signatures, or predictive models — are an expanding frontier. They can index stable traits or transient states; draw on structural or functional measures from MRI, fMRI, PET, EEG, MEG, or fNIRS; and serve many distinct goals: differential diagnosis of neurological disorders, biologically grounded subtyping of patients, surrogate endpoints when symptom measures are unreliable, prognosis of future disease course, and mechanistic targets for treatment development. The U.S. FDA distinguishes more than seven major biomarker categories with different uses. And because fMRI, EEG, and MEG measure the living brain in action, they raise a further possibility: neuromarkers of a person's *mental state* — what they perceive, remember, intend, or feel. The past decade of "decoding" studies has shown that far more about the mind can be read from patterns of brain activity than most researchers once thought possible.

Yet when it comes to using fMRI in the clinic or courtroom, most neuroimaging experts agree that "we're not there yet." That has not stopped the influx: one review found a four-fold increase in the use of neuroscience evidence in the legal system from 2005 to 2015, with neurological evidence appearing in 10–12% of murder trials and roughly 25% of death penalty trials, typically to support claims about competence, self-control, or insanity. Making valid claims about an *individual* requires that brain–behavior relationships be (a) reproducible across individuals and studies, (b) generalizable across scanners, contexts, and populations, and (c) large enough in effect size to license inferences about a single person. This last hurdle is formalized as the **"group to individual" (G2i) problem**: individual inference demands data quality and effect sizes roughly an order of magnitude beyond what group-average inference requires — plus the often-unexamined assumption that this person's brain works like the brains of the comparison group. Consider a defendant whose dorsolateral prefrontal cortex (dlPFC) activity sits two standard deviations below the mean, offered as evidence that a brain abnormality caused his compulsive shoplifting. For that claim to hold, dlPFC activity would need to be strongly tied to delay of gratification (and ideally shoplifting), the *same* region would need to be implicated in this particular brain, and his brain would have to be organized like those of the reference population. None of these steps is trivial.

Other obstacles apply to biomarkers of every kind. For conditions with low base rates, a positive test can be deeply misleading — the probability of actually having the condition may remain low even after a positive result (this statistical reality prompted the 2016 revision of mammogram screening recommendations). Biomarkers can be expensive, burden health-care systems, and be unequally accessible; tests that begin as optional can drift toward becoming obligatory. And before a measure is trusted in high-stakes settings, it must be standardized and validated: DNA evidence became a courtroom gold standard only after a rigorous process of establishing norms, quality-control procedures, and laboratory accreditation. Doing the same for neuroimaging — standardizing acquisition and analysis across a bewildering diversity of pipelines — is a daunting prospect. Importantly, none of these problems is unique to brain imaging; many widely accepted forms of evidence, from clinical diagnoses to polygraphs, fingerprints, and other forensic techniques, are less accurate than their reputations suggest.

Even where the science is solid, the decision to *deploy* a neuromarker in clinical, legal, commercial, military, or educational settings is about more than science. It involves cost-effectiveness, fair access, potential misuse, commercial viability, and unintended consequences. Neuroimaging is scientifically validated for detecting some forms of brain pathology, including Alzheimer's disease — but whether everyone at risk should be scanned depends on cost, on who pays, and on whether the result would change treatment or life decisions. And even when a scan says little about an individual, imaging can yield population-level insight: resting-state connectivity studies suggested that dementia-related neuropathology spreads along anatomically connected networks, supporting the hypothesis that misfolded proteins propagate from neuron to neuron much like prion diseases — a new disease concept opening new treatment avenues. Pain tells a similar story: brain activity measured with fMRI corresponds reliably to pain experience across individuals and studies, yet neuroimaging is not used clinically for pain, and there are legitimate fears it could be misapplied as a "pain lie detector" in legal settings where financial stakes are high.

Meanwhile, companies are forging ahead of the evidence — marketing fMRI lie detection to courts, promising EEG-based diagnosis of ADHD, autism, depression, concussion, and more, and claiming to reveal consumers' and voters' hidden preferences (or even pets' secret emotions). Such claims typically far outstrip even the most optimistic reading of the literature, and they cause real harm. But the opposite error is also harmful: the status quo these technologies would replace is often poor. Symptom-based psychiatric diagnoses are unreliable; courtroom experts opine on handwriting and bite marks with little scientific backing; the assumption that infants (or particular ethnic groups) feel less pain led to decades of documented harm; women and minorities remain systematically under-treated for pain; and tens of thousands of people yearly undergo back surgery without clear spinal pathology, 20–40% of whom end up worse off ("failed back surgery syndrome"). Chronic pain alone costs the U.S. an estimated $600 billion per year — more than cancer and heart disease combined. The critical question is therefore comparative: *does a new technology help us make better decisions than we would make without it?* The path forward is a "middle way" — cautious, cumulative, rigorously principled progress, built on open sharing across laboratories — pursued with urgency, because sitting on the sidelines has costs too.

## Key ideas in pictures

Chapter 3 has no figures in the book; the diagrams below summarize its central arguments.

**The group-to-individual (G2i) inference chain.** Every link must hold before a group-level finding can support a claim about one person — and each link demands stronger evidence than the one before.

```{mermaid}
flowchart TD
    A["Group-level finding<br/>e.g., lower dlPFC activity is associated with<br/>poorer delay of gratification, on average"] --> B["Reproducibility<br/>Does the association replicate<br/>across studies and individuals?"]
    B --> C["Generalizability<br/>Does it hold across scanners, contexts,<br/>and populations — including this person's?"]
    C --> D["Effect size<br/>Is the relationship strong enough<br/>(~10x group-level requirements)<br/>to classify a single person?"]
    D --> E["Individual mapping<br/>Is the same region, with the same function,<br/>implicated in this particular brain?"]
    E --> F["Defensible individual inference<br/>e.g., a dlPFC abnormality plausibly<br/>contributed to this defendant's behavior"]
    style A fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style F fill:#dcfce7,stroke:#22c55e,color:#14532d
```

**The adoption decision is bigger than the science.** Scientific validation is necessary but not sufficient; deploying a neuromarker in the real world raises questions science alone cannot answer.

```{mermaid}
flowchart LR
    subgraph S["Scientific validation"]
        S1["Reliable measurement"]
        S2["Reproducible brain-outcome links"]
        S3["Standardized acquisition & analysis"]
        S4["Accredited laboratories<br/>(the DNA-evidence precedent)"]
    end
    subgraph D["Societal adoption"]
        D1["Cost-effectiveness & who pays"]
        D2["Fair and equal access"]
        D3["Base rates & predictive value"]
        D4["Potential misuse<br/>('pain lie detector', neuromarketing)"]
        D5["Does it change decisions<br/>or treatment at all?"]
    end
    S --> D --> U["Use in clinic, courtroom,<br/>school, or marketplace"]
    style S fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style D fill:#fef3c7,stroke:#f59e0b,color:#713f12
    style U fill:#dcfce7,stroke:#22c55e,color:#14532d
```

**Two ways to cause harm.** The chapter's closing argument is that both overreach and excessive caution have victims; the goal is a "middle way" judged by a single comparative question.

```{mermaid}
flowchart TD
    Q["Does the new brain measure help us make<br/>better decisions than we would otherwise?"]
    Q --> O["Overreach<br/>Commercial claims that outstrip evidence:<br/>fMRI lie detection, EEG diagnosis of<br/>ADHD/autism/depression, neuromarketing"]
    Q --> I["Inaction (status quo bias)<br/>Unreliable symptom-based diagnoses,<br/>unscientific forensic testimony,<br/>under-treated pain, failed back surgeries,<br/>$600B/yr cost of chronic pain"]
    O --> H1["Harm: bad decisions, eroded trust,<br/>obscured understanding of valid measures"]
    I --> H2["Harm: preventable suffering from<br/>decisions made with poor evidence"]
    Q --> M["The middle way<br/>Cautious, cumulative, principled progress —<br/>open sharing, rigorous validation, urgency"]
    style Q fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style O fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style I fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style M fill:#dcfce7,stroke:#22c55e,color:#14532d
```

## Thought questions

1. Courts routinely admit expert testimony on handwriting, bite marks, and polygraph-adjacent techniques with weak scientific foundations, while many scholars argue fMRI evidence should be excluded until it meets far stricter benchmarks. Is it defensible to hold neuroimaging to a *higher* standard than evidence already in use — or does consistency demand that we either raise the bar for all forensic evidence or admit brain imaging under the same permissive rules? Who is harmed by each choice?

2. You are a judge in the shoplifting case described in this chapter: the defense offers a scan showing dlPFC metabolism two standard deviations below the mean, alongside published group studies linking dlPFC to delay of gratification. Walk through the G2i chain and specify concretely what additional evidence — about reproducibility, effect size, and this defendant's own brain — you would require before giving the scan any evidentiary weight. Is there any realistic body of evidence that would satisfy you, and if not, what does that imply about neuroimaging in criminal responsibility determinations generally?

3. A validated brain-based pain neuromarker would address genuine injustices — under-treatment of women and minorities whose pain reports are distrusted, and unnecessary surgeries guided by misleading spine MRIs. But insurers and defendants in tort cases would have strong financial incentives to demand such scans, and a test that begins as optional may become effectively obligatory. Design a policy for deploying a pain neuromarker that captures the benefits while preventing it from becoming a mandatory "pain lie detector." What happens, under your policy, to the claimant whose scan is negative but who insists the pain is real?

4. Amyloid PET and related measures can indicate elevated Alzheimer's risk years before symptoms, but disease-modifying treatment options remain limited. Given the chapter's argument that a scan's value depends on whether it changes decisions, who should decide whether an asymptomatic person is scanned — and should the results be shareable with (or discoverable by) long-term-care insurers, employers, or family members with their own genetic stake? How does the population-level research value of widespread scanning weigh against individual risks of unactionable bad news?

5. Consumer neurotechnology companies currently sell EEG-based diagnosis, neuromarketing insight into voters' "hidden preferences," and similar services with claims that far outstrip the literature — yet these operate largely outside FDA-style regulation. Drawing on the DNA-evidence precedent (standardization, quality control, laboratory accreditation), sketch what a credible regulatory regime for commercial neurotechnology would require. Would your regime have blocked anything of value, and would a scientifically honest company be able to survive under it?

## Quiz yourself

:::{dropdown} **Q1.** What is a biomarker?
**Answer:** A biological measure that serves as an indicator of a disease-relevant physiological process or function — for example, bright fluid-filled regions on a chest CT, or a lab value derived from an image. Biomarkers guide diagnosis, treatment selection, and monitoring.
:::

:::{dropdown} **Q2.** The chapter's emergency-room story illustrates two distinct kinds of image-based biomarkers. What are they?
**Answer:** Qualitative biomarkers read by an expert (a radiologist recognizing fluid in the lungs on a chest CT) and quantitative biomarkers, where a standardized procedure converts image measurements into a clinically actionable numeric score (echocardiography yielding left ventricular ejection fraction, an indicator of heart failure).
:::

:::{dropdown} **Q3.** What is a neuromarker, and what kinds of things can neuromarkers measure?
**Answer:** A neuromarker (also called a brain signature or predictive model) is a biomarker based on neuroimaging — structural or functional measures from MRI, fMRI, PET, EEG, MEG, fNIRS, and related methods. Neuromarkers can index stable traits or transient states, and can serve diagnosis, patient subtyping, prognosis, surrogate endpoints, and mechanistic treatment targets.
:::

:::{dropdown} **Q4.** How has the use of neuroscience evidence in the legal system changed in recent years?
**Answer:** It rose sharply — one review found a four-fold increase from 2005 to 2015 across major felonies, with neurological evidence used in about 10–12% of murder trials and roughly 25% of death penalty trials, typically to support claims about competence, behavioral control, or insanity, and in tort cases to corroborate brain injury and chronic pain.
:::

:::{dropdown} **Q5.** What three measurement properties must a brain–behavior relationship have before it can support claims about an individual person?
**Answer:** It must be (a) reproducible — stably related to the outcome across individuals and studies; (b) generalizable — holding across scanners, contexts, and populations; and (c) sufficiently strong — with effect sizes large enough to permit inference about a single person rather than only a group average.
:::

:::{dropdown} **Q6.** What is the "group to individual" (G2i) problem?
**Answer:** The problem of using an effect established in a group on average to draw conclusions about one individual. Individual-level inference requires data quality and effect sizes roughly an order of magnitude greater than group-level inference, plus assumptions that the individual's brain is organized like those of the comparison population — assumptions that can be unacceptably naive.
:::

:::{dropdown} **Q7.** Why can a positive biomarker test be misleading for conditions with low base rates?
**Answer:** When a condition is rare in the population, even a fairly accurate test yields many false positives relative to true positives, so the probability of actually having the condition remains low after a positive result. This statistical fact motivated the 2016 change in consensus recommendations to reduce mammogram screening for breast cancer.
:::

:::{dropdown} **Q8.** What precedent does DNA evidence set for how neuroimaging biomarkers might become accepted in court?
**Answer:** DNA testing became a courtroom gold standard only after a rigorous process of establishing detailed norms and regulations, quality-control procedures for validating each sample, and accreditation of testing laboratories. Neuroimaging would need the analogous machinery — standardized acquisition and analysis, per-scan quality validation, and accredited labs — which is daunting given the diversity and complexity of imaging pipelines.
:::

:::{div}
:class: book-tile
📖 **The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::

---

[⌂ Back to home](https://torwager.github.io/elements-of-fmri-tutorials/) · [Table of contents](../contents.md) · [How to use this site](../how-to-use.md)
