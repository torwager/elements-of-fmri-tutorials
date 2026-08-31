---
title: "10. Why Imaging Is Not Phrenology"
subject: "Part 2: Brain Mapping"
---

# Why Imaging Is Not Phrenology

:::{admonition} What you will learn
:class: tip
- What phrenology actually claimed, why it counts as pseudoscience, and why "imaging is the new phrenology" is a false equivalence
- How incentives in publishing and media coverage inflate claims about brain findings — and the three hidden inferences (predictive value, necessity, sufficiency) smuggled into "the neural basis of X" headlines
- The three criteria that separate scientific from pseudoscientific claims about brain–mind associations: reproducible evidence, unbiased hypothesis tests, and plausible mechanisms
- Why prospective, falsifiable prediction — testing decoding models on independent data specified in advance — is the strongest form of validation neuroimaging offers
- Concrete evidence that fMRI yields robust, reproducible brain–mind mappings: meta-analytic classifiers, pattern responses in individuals, and validated signatures such as the Neurologic Pain Signature
:::

## Overview

Pseudoscientific "brain reading" has a long history, and neuroimaging has not always kept clear of it. Early functional imaging studies sometimes committed errors of inference and leaps of logic serious enough that critics asked whether the whole enterprise was pseudoscience — some comparing it to phrenology, the 19th-century practice of reading personality from bumps on the skull. The comparison is worth taking seriously, because answering it forces us to articulate what makes any claim about brain–mind associations scientific. The short answer developed in this chapter: neuroimaging signals arise from real, well-characterized neurophysiological mechanisms; task-evoked activity replicates across laboratories; and the field increasingly tests strong, predefined, falsifiable predictions. These are exactly the hallmarks of science that phrenology lacked.

Much of the criticism has been aimed at studies claiming to find an area *devoted to* a complex mental process — a "depression area," "love area," or "reward center." The incentive structure pushes in this direction: journals and funders prize work billed as novel and transformative, and the use of positive words in scientific papers rose nearly 10-fold over 40 years, with paradigm-shift words like "novel," "innovative," and "unprecedented" rising 150-fold. More often, though, it is popular media coverage rather than the papers themselves that inflates the claims. The result is a familiar leap of logic. As Chapter 7 explained, the vast majority of studies establish an association between a mental state and regional activity — *forward inference*, $P(\text{brain} \mid \text{mental state})$ — but do not establish that the activity distinguishes that state from related ones like general arousal or emotion (*reverse inference*), nor that the association is strong or generalizable enough to be meaningful. Concluding that a study found "the neural basis for love" quietly imports three untested inferences: (a) that we can tell someone is in love from their brain (positive predictive value); (b) that without the region one would be incapable of love (necessity); and (c) that stimulating the region would produce love (sufficiency). Testing these would require a long research program of stimulation, lesion, and imaging studies across many alternative mental states — not a single association map.

When people sense extreme, unjustified claims, they often swing to an equally extreme opposing position — and "fMRI is the new phrenology" is a case in point. Phrenology was developed by Franz Joseph Gall and his student Johann Spurzheim, who claimed to identify personality traits from patterns of skull protrusions. Gall was, notably, among the earliest proponents of *localizationism* — the idea that particular mental functions are accomplished by separate, localized brain regions — a view later supported empirically by Broca, Hughlings Jackson, Fritsch and Hitzig, and others, who showed that lesions in particular areas produce distinct patterns of deficits. The tension between localizationist and more holistic, integrative accounts of brain function continues today, and both views have merit in different cases. But Gall's method was another matter. He reportedly conceived the theory after noticing that a classmate with excellent verbal memory had bulging eyes; he went on to catalog twenty-seven faculties — from "Reproductive Instinct" and "Tendency to Murder" to "Wisdom," "Poetic Talent," and "Religious Sentiment" — each tied to a skull protrusion on the basis of case anecdotes. "Tendency to Murder" was localized above the ears partly because an acquaintance with bulges there had left business to become a butcher, and another had become a public executioner.

Phrenology was controversial even in its own era and is now a premier example of **pseudoscience**: claims that wear the costume of scientific investigation but rest on cherry-picked evidence and invalid tests. Equating skull bumps with Broca's findings is a false equivalence for two reasons. First, Gall and Spurzheim's associations were not reproducible, whereas Broca's finding that left inferior frontal lesions impair speech has been replicated many times and corroborated by converging methods, including neuroimaging. Second, phrenology has no plausible mechanism: variation in skull thickness and shape places little or no constraint on the function of the tissue inside. Across every technique used to link biology and mind — lesion studies, invasive electrophysiology, EEG and MEG, calcium imaging, optogenetic and chemogenetic stimulation, PET and fMRI — the same three criteria separate science from pseudoscience: (a) the strength and reproducibility of the evidence, (b) hypothesis testing with statistical methods free from bias, and (c) the existence of plausible mechanisms underlying the association.

Neuroimaging passes these tests. Mechanistically, the links between BOLD signals, regional blood flow, and electrical activity in neuronal populations are firmly established: BOLD fMRI is closely coupled to local field potentials (Chapter 14). Empirically, associations between fMRI activity and stimuli, task demands, and behavior are highly reproducible across studies and laboratories, documented in hundreds of meta-analyses. And predictively, fMRI patterns can be used to *decode* the type or intensity of stimuli and tasks from brain activity alone, with models tested unbiasedly in independent samples. Decoding is the most powerful of these validations because it is falsifiable in the strongest sense: the model predicts an objectively observable state of the world — $P(\text{stimulus type} \mid \text{fMRI activity})$ — in advance, without post hoc adjustment. Making and testing precise, falsifiable predictions is how scientific claims earn trust in every field; Einstein became famous when his 1915 predictions about gravitational lensing were confirmed during the 1919 solar eclipse, and as Milton Friedman put it, "The only relevant test of the validity of a hypothesis is comparison of prediction with experience."

Two examples make this concrete. Neurosynth.org aggregates reported activation peaks from thousands of studies, yielding meta-analytic maps for thousands of psychological terms. Using maps for working memory, emotion, and pain, a Gaussian Naive Bayes classifier was applied — untouched — to 281 images from *independent* individual participants (79 experiencing heat pain, 108 viewing aversive pictures, 94 performing working memory). It picked the correct task with 65–95% sensitivity and 80–98% specificity, far above the ~38% chance level. In the second example, the Neurosynth working-memory map was used as a *pattern of interest*: a weighted average of each participant's data image, with the meta-analytic map supplying the weights, gives a single **pattern response** per person. Although the pattern came from entirely different studies, populations, and scanners, the response was stronger for working memory than rest in 100% of test participants, and stronger for a harder (3-back) than an easier (2-back) task in 95% — effect sizes of d = 1.97 and 1.09, far beyond conventional benchmarks for "large" effects. A stable, reproducible pattern, strong enough to support inferences about individual people.

These are not isolated cases. Reproducible findings span distinct emotion types, social cognition, emotion regulation, memory, and brain disorders from Alzheimer's to PTSD and chronic pain. Trait-like clinical biomarkers are harder, but progress is real: one large study (N = 871) classified autism at 67% accuracy, and another with over 1,000 participants identified four depression "biotypes," classifying patients versus controls at close to 90% accuracy and predicting response to brain-stimulation treatment. Perhaps the clearest exhibit is the **Neurologic Pain Signature (NPS)**, an fMRI-based measure that responds to painful stimulation across diverse populations worldwide and tracks pain intensity strongly enough to make accurate predictions about individuals. Validated in more than 50 unique samples, and tested in 20 cohorts across North America, Europe, and Asia (N = 603), it responded as expected in 95% of participants, with effect sizes roughly 10 times larger than typical correlations between stable traits and resting-state activity or brain structure. None of this means the neuroimaging literature is free of pseudoscientific claims — rigor varies widely — but many associations are robust enough that they are simply not in dispute. In short, they have all the hallmarks of solid science that phrenology did not.

## Key ideas in pictures

:::{figure} images/ch10_fig1_phrenology_map.png
:alt: A 19th-century phrenological map showing mental faculties assigned to regions of the skull
:width: 70%

A phrenological map: putative associations between protrusions on the skull and mental faculties such as "firmness," "benevolence," and "destructiveness." Gall assigned his twenty-seven faculties to skull locations on the basis of case anecdotes — an unreproducible, cherry-picked method with no plausible mechanism, since skull shape places essentially no constraint on the function of the brain tissue beneath. *(Figure 10.1 from the book.)*
:::

**The hidden inferences behind "the neural basis of X."** A forward-inference finding — activity associated with a mental state — is routinely misread as licensing three much stronger claims, each of which requires its own program of evidence.

```{mermaid}
flowchart TD
    A["Published finding (forward inference)<br/>Mental state X is associated with<br/>activity in region R, on average"] --> B["'Region R is the neural basis of X'"]
    B -.->|implicitly assumes| C["Positive predictive value<br/>We can tell someone is in state X<br/>from activity in R<br/>(requires reverse inference: testing R<br/>against arousal, emotion, other states)"]
    B -.->|implicitly assumes| D["Necessity<br/>Without R, state X is impossible<br/>(requires lesion evidence)"]
    B -.->|implicitly assumes| E["Sufficiency<br/>Stimulating R produces state X<br/>(requires stimulation evidence)"]
    style A fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style B fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style C fill:#fef3c7,stroke:#f59e0b,color:#713f12
    style D fill:#fef3c7,stroke:#f59e0b,color:#713f12
    style E fill:#fef3c7,stroke:#f59e0b,color:#713f12
```

**What separates science from pseudoscience — for any brain–mind method.** The same three criteria apply whether the technique is lesion mapping, electrophysiology, optogenetics, or fMRI. Phrenology fails all three; well-conducted neuroimaging meets all three.

```{mermaid}
flowchart LR
    subgraph P["Phrenology"]
        P1["Anecdotal case studies,<br/>cherry-picked evidence"]
        P2["No unbiased hypothesis tests;<br/>claims not reproducible"]
        P3["No plausible mechanism:<br/>skull shape does not constrain<br/>brain function"]
    end
    subgraph N["Neuroimaging done well"]
        N1["Reproducible associations:<br/>hundreds of meta-analyses;<br/>consistent task-evoked activity<br/>across labs"]
        N2["Predefined, falsifiable predictions<br/>tested in independent samples<br/>without post hoc adjustment"]
        N3["Established mechanism:<br/>BOLD is coupled to local field<br/>potentials and blood flow"]
    end
    C["Criteria for scientific claims<br/>(a) strong, reproducible evidence<br/>(b) unbiased hypothesis testing<br/>(c) plausible mechanism"] --> P
    C --> N
    P --> V1["Pseudoscience"]
    N --> V2["Valid brain-mind mapping"]
    style C fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style P fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style N fill:#dcfce7,stroke:#22c55e,color:#14532d
    style V1 fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style V2 fill:#dcfce7,stroke:#22c55e,color:#14532d
```

**Prediction as the strongest validation.** Decoding models close the loop that association studies leave open: they commit to a falsifiable prediction about the observable world before the data are seen.

```{mermaid}
flowchart LR
    A["Train a model on prior data<br/>e.g., Neurosynth meta-analytic maps<br/>for working memory, emotion, pain"] --> B["Fix the model in advance<br/>no tuning to the test sample"]
    B --> C["Predict observable states<br/>in independent participants:<br/>P(stimulus or task | fMRI activity)"]
    C --> D{"Prediction<br/>correct?"}
    D -->|"Yes: e.g., 3-way task classification<br/>65-95% sensitivity, 80-98% specificity<br/>vs ~38% chance; NPS responds in 95%<br/>of participants across 20 cohorts"| E["Claim survives a strong,<br/>unbiased test"]
    D -->|No| F["Claim is falsified —<br/>and the field learns which states<br/>fMRI does not capture"]
    style A fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style E fill:#dcfce7,stroke:#22c55e,color:#14532d
    style F fill:#fef3c7,stroke:#f59e0b,color:#713f12
```

## Thought questions

1. Gall was an early localizationist, and localizationism itself was later vindicated in part by Broca, Fritsch and Hitzig, and others. Since phrenology's *core theoretical commitment* was not entirely wrong, what exactly made it pseudoscience? Use the chapter's three criteria to argue whether a theory can be partly correct and still pseudoscientific — and whether a methodologically rigorous field can still produce wrong theories.

2. Lesion and stimulation studies are the classic route to claims of necessity and sufficiency, and double dissociations are often treated as the gold standard for functional specialization. But suppose region A supports a resource that task X relies on heavily and task Y only weakly, while region B does the reverse: a double dissociation appears without either region being a dedicated "module," and distributed, population-coded representations could produce the same pattern. How should this possibility change the way we interpret a "love area" claim even *after* lesion and stimulation evidence is added to the imaging evidence the chapter says is missing?

3. The Neurosynth working-memory pattern classified harder versus easier tasks in 95% of individuals, yet the chapter warns against concluding that any region is "the working memory area." Reconcile these: what precisely does a high-performing distributed decoding model license us to say about brain–mind mappings, and what does it still leave undetermined about the causal role and functional specialization of any voxel that contributes to it?

4. The chapter attributes much claim inflation to incentives — journals rewarding "transformative" findings and journalists needing to sell stories. Falsifiable prediction is offered as a scientific corrective, but does it fix the incentive problem? Consider: who chooses which predictions get tested and reported, what happens to failed predictions, and whether a decoding model with impressive accuracy can itself become the next over-sold "brain reading" product. Design one concrete norm (for journals, funders, or media) that would push the field toward the chapter's standard.

5. The NPS shows effect sizes roughly 10 times larger than typical correlations between stable traits and resting-state activity or brain structure. What does this contrast suggest about *which kinds* of brain–mind mappings (transient evoked states vs. stable traits and disorders) are currently strong enough for individual-level inference — and how should that shape your expectations for the clinical biomarkers (autism, depression biotypes) described at the end of the chapter?

## Quiz yourself

:::{dropdown} **Q1.** What was phrenology?
**Answer:** A 19th-century practice, developed by Franz Joseph Gall and his student Johann Spurzheim, of identifying personality traits and mental faculties by examining the pattern of bumps on a person's skull. Gall cataloged twenty-seven faculties, each supposedly tied to a particular skull protrusion.
:::

:::{dropdown} **Q2.** What is pseudoscience, and why does phrenology qualify?
**Answer:** Pseudoscience is a set of claims that appear to be based on scientific investigation and hypothesis testing but in fact rest on selective "cherry picking" of evidence or other invalid tests. Phrenology qualifies because its associations came from unreplicated case anecdotes, its claims were not reproducible, and it had no plausible mechanism — skull shape places little or no constraint on the function of the brain tissue inside.
:::

:::{dropdown} **Q3.** A study reports that a brain region activates when people view pictures of their romantic partners. What three implicit inferences are smuggled into the headline "scientists find the brain's love center"?
**Answer:** (a) That we can tell whether someone is in love from their brain activity (positive predictive value — a reverse inference requiring the region to distinguish love from arousal, thinking about others, or emotion in general); (b) that without the region one would be incapable of love (necessity); and (c) that stimulating the region would produce love (sufficiency). None of these follows from a forward-inference association.
:::

:::{dropdown} **Q4.** What three criteria distinguish scientific from pseudoscientific claims about brain–mind associations, regardless of the technique used?
**Answer:** (a) The strength and reproducibility of the evidence; (b) the use of hypothesis testing and statistical methods free from bias; and (c) the existence of plausible mechanisms underlying the association.
:::

:::{dropdown} **Q5.** Why is Broca's finding not "just phrenology with lesions"?
**Answer:** Because it meets the criteria phrenology fails: the finding that left inferior frontal lesions impair speech is highly reproducible and has been confirmed by many investigators, it has been corroborated by converging evidence from other methods including neuroimaging, and it rests on a plausible mechanism — damage to the brain tissue that performs the function, not variation in the overlying skull.
:::

:::{dropdown} **Q6.** What mechanistic evidence grounds the link between fMRI signals and neural activity?
**Answer:** Links between neuroimaging signals, regional blood flow, and electrical activity in neuronal populations are firmly established: BOLD fMRI signals are closely coupled to local field potentials in the brain. This gives neuroimaging the plausible physiological basis that phrenology never had.
:::

:::{dropdown} **Q7.** Why is decoding (predicting stimuli or tasks from brain activity) a stronger validation than a standard activation study?
**Answer:** Because it makes falsifiable predictions about objectively observable variables — P(stimulus type | fMRI activity) — specified in advance and tested in independent samples without post hoc adjustment. If the predictions are precise and the test unbiased, the results give an honest picture of which states fMRI activity does and does not capture, in line with the principle that hypotheses are validated by comparing prediction with experience.
:::

:::{dropdown} **Q8.** What did the two Neurosynth-based examples and the Neurologic Pain Signature demonstrate about the robustness of fMRI brain–mind mappings?
**Answer:** A Naive Bayes classifier built from meta-analytic maps classified working memory, emotion, and pain in 281 independent participants with 65–95% sensitivity and 80–98% specificity (chance ~38%). The Neurosynth working-memory pattern — derived from entirely different studies, populations, and scanners — distinguished task from rest in 100% and harder from easier tasks in 95% of individuals (d = 1.97 and 1.09). The NPS responded to painful stimulation in 95% of 603 participants across 20 international cohorts, with effect sizes about 10 times larger than typical trait–brain correlations — showing that fMRI mappings can be robust, reproducible, sensitive, and specific.
:::
