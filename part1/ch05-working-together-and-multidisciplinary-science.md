---
title: "5. Working Together and Multidisciplinary Science"
subject: "Part 1: Motivation"
---

# Working Together and Multidisciplinary Science

:::{admonition} What you will learn
:class: tip
- Why the brain serves as a "common language" that connects psychology, neuroscience, medicine, statistics, engineering, computer science, and beyond
- What distinct kinds of expertise a state-of-the-art fMRI study draws on, and what each contributes at every stage from design to interpretation
- Why no single discipline sees the whole picture — and why integration across fields is essential for a cumulative science of mind and brain
- What makes multidisciplinary collaboration genuinely hard: unfamiliar concepts and vocabulary, misaligned incentives, and the humility it demands
- What "bridge" scientists are, and why developing proficiency across several fields makes both individual researchers and teams more effective
:::

## Overview

Many fields study aspects of the mind, the brain, and the body — but largely in isolation from one another. Psychologists study the mind and behavior, traditionally with minimal reference to biology. Neuroscientists study the brain from ion channels to systems, but only a fraction of that work concerns mental experience or behavior. Clinical researchers study the treatment and prevention of brain-centered disorders, yet only a small percentage of their studies directly assess brain function. Statistics, engineering, and computer science lead the study of complex computational processes, but have only recently turned toward biological systems. Educators and athletic trainers study learning, lawyers and judges reason about intention, suffering, and capacity — all processes that arise in the brain — yet rarely with any brain measurement or neuroscientific grounding. Each discipline offers a crucial but incomplete window into the same frontier: the study of the mind and brain is, ultimately, the study of us.

The situation recalls the well-known story of the blind people examining an elephant: one feels something long, rubbery, and flexible; another a smooth, firm surface; a third a flat, delicate membrane. The mind and brain are a very large elephant. Their study spans spatial scales from molecules to cells to systems, and time scales from the nanosecond gating of ion channels to brain–mind relationships across a human lifetime; it ranges from elementary reflexes to our most abstract thoughts. Each discipline specializes in a different piece. Some separation is inevitable — scientific knowledge grows exponentially, and mastering even a specialized subfield's literature is hard enough — but connections across fields are vital if we are to form a comprehensive understanding of the mind and brain and use it to make informed, beneficial, and compassionate decisions across every sphere of human life.

Neuroimaging is inherently multidisciplinary in both directions: it *enables* interaction across fields, and it *requires* it. The brain provides a common foundation and language around which practitioners from diverse disciplines can come together. Consider a neuroscientist studying the molecular basis of learning, a pharmacologist interested in antipsychotic drugs, a psychiatrist examining depression, a social psychologist investigating altruism, and a marketing professor studying purchasing decisions. Each is studying the function of the mesolimbic dopamine system — in different ways, with a unique perspective to contribute. The same logic operates within a field: psychological theories of motivation, emotion regulation, and prejudice can proliferate in separate vocabularies while relying on the same core neural systems. By grounding concepts in measurable processes that can be compared and contrasted whatever words describe them, neuroimaging helps build a cumulative, unified science.

At the same time, doing cutting-edge neuroimaging well demands expertise no single field supplies. Collecting high-quality data requires physics and engineering. Analyzing it requires signal processing, statistics, and increasingly machine learning. Interpreting findings — and tailoring analyses to neuroscientific, psychological, or clinical hypotheses — requires expertise in those domains and in neuroanatomy.

To make this concrete, imagine a study of how an antidepressant affects the prefrontal cortex and how those effects relate to symptom trajectories. We know a great deal about the molecular pharmacology of such drugs, but the science of their effects on large-scale brain systems supporting thought, emotion, and decision-making is in its infancy. A "dream team" for this study might span seven fields. The psychiatrist and psychologist define the clinical phenomena, sampling frame, and outcome measures, and — with the statistician and pharmacologist — design a task that isolates depression-related mental processes and efficiently detects drug effects. The psychologist and statistician ensure the design is efficient, well powered, and capable of supporting valid causal inference. The pharmacologist contributes knowledge of the drug's cellular and molecular mechanisms and absorption kinetics, and of drug effects on vasculature and blood gases that could masquerade as neural signals. The neuroscientist brings knowledge of how the drug penetrates the brain and acts on neurons, glia, and neural systems, helping generate hypotheses about which regions should be affected. The MRI physicist or biomedical engineer designs an acquisition protocol that yields high-quality images and minimizes artifacts in the regions that matter most, and understands how vascular and physiological drug effects can alter the fMRI signal independent of neural function. The computer scientist manages the data volume and provides pipelines for processing, analysis, and quality control. During analysis, the statistician scrutinizes the data and the assumptions behind each test; afterward, the neuroanatomist localizes the effects, and the psychiatrist, neuroscientist, and psychologist interpret them against the existing literature.

This portrait does not mean every study needs seven experts — that would be impractical. The real lesson is that the best science builds on foundational knowledge from all of these disciplines, and that neuroimaging researchers should develop key proficiencies across fields. A scientist using fMRI might come from any of these home disciplines but likely has some capability in nearly all of them. Such researchers can run high-quality studies with small, agile teams — and their breadth helps them recognize when, and whom, to ask for advice. A good metaphor for this kind of collaboration is a *confluence*: many great rivers running together, their ideas and techniques intermingling. Scientific confluence benefits science and society well beyond fMRI itself — it is fertile ground for new ideas and techniques, and it builds a multidisciplinary scientific literacy that transfers broadly across fields and jobs.

None of this is easy. Multidisciplinary science demands communication, patience, humility, and resources. Collaborators must learn and discuss unfamiliar concepts in unfamiliar terms, and each must be willing to be a novice in someone else's area of expertise. It is often easier to stay within a narrow specialty, communicating in the specialized language of one's own field; collaboration instead requires caring about problems outside one's defined interests, publishing in venues one's home field may not read, and spending time teaching teammates basics that are old news at home but innovative in an interdisciplinary context. Funding agencies have recognized these obstacles and created mechanisms and support structures to promote cross-field teams. But perhaps the most powerful answer is individual: scientists who cultivate multiple types of expertise — "bridge" scientists — shrink the gulf between the psychologist and the physicist, or the pharmacologist and the statistician, and act as the glue holding larger teams together. A little knowledge goes a long way; knowing a few words of a colleague's scientific language can dramatically improve communication. Offering a route to that kind of multidisciplinary fluency is one of the reasons this book exists.

## Key ideas in pictures

Chapter 5 has no figures in the book; the diagrams below summarize its central ideas.

**The brain as a common language.** Fields with divergent concepts, vocabularies, and methods converge on the same object of study. Neuroimaging grounds their questions in shared, measurable processes — here, five very different researchers all studying the mesolimbic dopamine system.

```{mermaid}
flowchart TD
    N["Neuroscientist<br/>learning"] --> B
    P["Pharmacologist<br/>antipsychotics"] --> B
    Y["Psychiatrist<br/>depression"] --> B
    S["Social psychologist<br/>altruism"] --> B
    M["Marketing professor<br/>purchasing"] --> B
    B["One system:<br/>mesolimbic dopamine,<br/>measured in common"]
    B --> U["Concepts grounded in<br/>measurable processes —<br/>a cumulative science"]
    style B fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style U fill:#dcfce7,stroke:#22c55e,color:#14532d
```

**Who does what in an fMRI study.** The hypothetical antidepressant study's "dream team," mapped onto the stages of a study. Most contributions span stages, and several kinds of expertise overlap — which is exactly why bridge scientists are so valuable.

```{mermaid}
flowchart TD
    subgraph D["1 · Design"]
        direction LR
        D1["Psychiatrist &<br/>psychologist<br/>phenomena · task"]
        D2["Statistician<br/>power · causal<br/>inference"]
        D3["Pharmacologist<br/>mechanisms ·<br/>kinetics"]
        D4["Neuroscientist<br/>regional<br/>hypotheses"]
    end
    subgraph A["2 · Acquisition"]
        A1["MRI physicist / engineer<br/>protocol · artifact control ·<br/>non-neural signal effects"]
    end
    subgraph C["3 · Processing & analysis"]
        direction LR
        C1["Computer scientist<br/>pipelines · QC"]
        C2["Statistician<br/>tests &<br/>assumptions"]
    end
    subgraph I["4 · Interpretation"]
        direction LR
        I1["Neuroanatomist<br/>localization"]
        I2["Psychiatrist · neuroscientist<br/>· psychologist<br/>meaning & literature"]
    end
    D --> A --> C --> I
    style D fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style A fill:#fef3c7,stroke:#f59e0b,color:#713f12
    style C fill:#fce7f3,stroke:#ec4899,color:#831843
    style I fill:#dcfce7,stroke:#22c55e,color:#14532d
```

**From separate streams to confluence.** Collaboration is hard for identifiable reasons, and the chapter identifies concrete forces that overcome them — chief among them, individual scientists who become bridges between fields.

```{mermaid}
flowchart TD
    O["Obstacles<br/>unfamiliar vocabulary ·<br/>incentives favor specialty ·<br/>time cost of teaching"]
    O --> R1["Communication<br/>& humility<br/>being a novice again"]
    O --> R2["Funding built for<br/>cross-field teams"]
    O --> R3["Bridge scientists<br/>the glue of teams"]
    R1 --> F["Confluence<br/>rivers of ideas<br/>running together"]
    R2 --> F
    R3 --> F
    F --> G["New ideas & techniques ·<br/>broad scientific literacy"]
    style O fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style F fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style G fill:#dcfce7,stroke:#22c55e,color:#14532d
```

## Thought questions

1. The chapter argues that the brain can serve as a "common language" that grounds concepts from different fields in shared, measurable processes. But translation can also flatten meaning: does a social psychologist's "altruism" and a marketing professor's "purchasing preference" really become commensurable because both engage the mesolimbic dopamine system? When does grounding concepts in a shared neural substrate genuinely unify theories, and when does it merely create an illusion of agreement?

2. In the antidepressant study example, several failure modes cross disciplinary boundaries — for instance, a drug's effect on vasculature and blood gases could mimic a neural effect, a fact sitting at the intersection of pharmacology, physiology, and MRI physics. Identify another plausible cross-disciplinary failure mode in an fMRI study you might design, and specify which combinations of expertise would be needed to notice it and which single-discipline team would likely miss it.

3. The chapter endorses both large multidisciplinary "dream teams" and small, agile teams of broadly trained researchers. What are the respective advantages and risks of each model — for scientific rigor, for training the next generation, and for who gets credit? Given that a seven-expert team is "highly impractical" for most studies, where should the field invest: in team-building infrastructure or in multidisciplinary education?

4. Academic incentives — hiring, promotion, and prestige — largely reward depth within a home discipline, while this chapter argues science needs bridge scientists who publish in unfamiliar venues and spend time teaching collaborators material that is not novel in their own field. If you were designing a department's tenure criteria or a funding agency's review process, what concrete changes would make bridge science a viable career strategy rather than a sacrifice?

5. The elephant parable implies that each field's picture of the mind and brain is real but partial. Choose a field discussed in the chapter that currently makes consequential decisions with little neuroscientific input — law, education, or athletic training, for example. What is one decision in that field that a mature multidisciplinary brain science could genuinely improve, what could go wrong if neuroimaging were imported without the accompanying expertise, and how do the cautions from the imaging-and-society discussion (Chapter 3) apply?

## Quiz yourself

:::{dropdown} **Q1.** In what sense does the brain provide a "common language" for different scientific fields?
**Answer:** The brain is a shared foundation that many disciplines — psychology, neuroscience, medicine, statistics, engineering, computer science, education, law — are all studying from different angles. Neuroimaging grounds each field's concepts in the same measurable processes, so theories couched in different vocabularies can be compared, contrasted, and integrated.
:::

:::{dropdown} **Q2.** What point does the parable of the blind people and the elephant illustrate in this chapter?
**Answer:** Each discipline perceives only one part of the mind and brain — one feels the trunk, another the side, another the ear — and mistakes its partial view for the whole. The study of mind and brain spans enormous spatial scales (molecules to systems) and time scales (nanoseconds to lifetimes), so every field's window is crucial but incomplete, and integration is required for a comprehensive understanding.
:::

:::{dropdown} **Q3.** Why is neuroimaging described as *inherently* multidisciplinary?
**Answer:** In two ways: it enables cross-field interaction by giving diverse disciplines a common, measurable object of study; and it requires cross-field expertise, because collecting high-quality data takes physics and engineering, analyzing it takes signal processing, statistics, and machine learning, and interpreting it takes neuroscientific, psychological or clinical, and neuroanatomical knowledge.
:::

:::{dropdown} **Q4.** Which seven fields make up the "dream team" in the chapter's hypothetical antidepressant study?
**Answer:** An MRI physicist or biomedical engineer, a psychologist, a statistician, a pharmacologist, a computer scientist, a psychiatrist, and a neuroscientist/neuroanatomist.
:::

:::{dropdown} **Q5.** In the antidepressant example, what does the pharmacologist contribute beyond knowledge of the drug's molecular mechanism?
**Answer:** Knowledge of the drug's absorption kinetics into brain tissue, which shapes the experimental design and analysis, and knowledge of drug effects on brain vasculature and blood gas levels — effects that can produce artifacts in the fMRI signal that mimic or mask true neural effects.
:::

:::{dropdown} **Q6.** Why might an MRI physicist matter for *interpreting* a drug study, not just for collecting good images?
**Answer:** Because drugs can affect vasculature and physiology, and the fMRI signal is vascular in origin, the physicist or biomedical engineer can identify how such effects might alter the measured signal independently of neural function — knowledge essential for distinguishing genuine neural drug effects from physiological artifacts.
:::

:::{dropdown} **Q7.** What is a "bridge" scientist, and why are they valuable?
**Answer:** A researcher who develops expertise across multiple disciplines, shrinking the gulf between fields so that, say, the psychologist and the physicist have a shared vocabulary. Bridge scientists act as the glue holding larger teams together, can run high-quality studies with small agile teams, and know when and whom to ask for specialist advice.
:::

:::{dropdown} **Q8.** What does the chapter identify as the main challenges of multidisciplinary collaboration, and what helps overcome them?
**Answer:** Challenges include learning unfamiliar concepts and terminology, the humility of being a novice in a collaborator's field, incentives favoring narrow specialization, publishing in venues outside one's home field, and time spent teaching teammates basics. Remedies include a spirit of collaboration and patience, dedicated funding mechanisms for multidisciplinary teams, and individual scientists developing multiple types of expertise — even a little knowledge of another field's "language" dramatically improves communication.
:::

:::{div}
:class: book-tile
![Cover of Elements of Functional Magnetic Resonance Imaging](../cover-small.jpg)
**The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::

---

[⌂ Back to home](https://torwager.github.io/elements-of-fmri-tutorials/) · [Table of contents](../contents.md) · [How to use this site](../how-to-use.md)
