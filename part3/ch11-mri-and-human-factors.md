---
title: "11. MRI and Human Factors"
subject: "Part 3: MRI Environment and MRI Signal"
---

# MRI and Human Factors

:::{admonition} What you will learn
:class: tip
- The three main hardware components of an MR scanner — magnet, radiofrequency coils, and gradient coils — and what each contributes to imaging
- Why the static magnetic field is *always on*, and how ferromagnetic, paramagnetic, and diamagnetic materials behave differently around it
- The major safety concerns in the MR environment: projectiles, implants, RF heating (SAR), acoustic noise, and how screening prevents accidents
- How stimuli are delivered and behavior is measured inside the scanner, and why precise timing synchronization matters
- How human factors — head movement, claustrophobia, anxiety, sleepiness, and the unusual scanner environment itself — shape data quality and study design
:::

## Overview

An MR scanner is a large, remarkably versatile piece of hardware built around three main components. The first is a superconducting electromagnet that generates an extremely strong static magnetic field, usually 1.5 to 7 Tesla (T). A standard 1.5T clinical field is about 30,000 times stronger than the Earth's magnetic field, and an "ultra-high field" 7T magnet is 140,000 times stronger. The second component is the radiofrequency (RF) coils, which lie close to the head and both transmit and receive energy at the resonant frequency of the tissue being imaged: brief RF pulses create an oscillating magnetic field perpendicular to the main field, and the resulting signal is captured, amplified, and digitized to extract the frequency and phase information from which brain images are constructed. The third component is the gradient coils, electromagnets that create controlled spatial variation in field strength — the key to encoding *where* in the brain each signal comes from, so that measurements can be reconstructed into 3-D images. The same hardware supports many kinds of measurement — anatomy, function, white-matter tracts, tissue elasticity, cerebrospinal fluid flow, and more; what differs is the software, the *pulse sequences* that orchestrate how RF energy is transmitted and measured.

fMRI is safe and noninvasive: millions of MRIs are performed each year in the U.S. alone without incident, and there are no known long-term effects of the magnetic field on biological tissue — which is why the same person can be scanned repeatedly across development, aging, learning, or treatment, with no regulatory limit on scan frequency (some studies have scanned one person on around 100 separate days). But safety depends on strict adherence to protocols, because the potential for serious harm is real. The single most important fact about the MR environment is that **the magnetic field is always on** — on weekends, on holidays, and when the console computers are powered down. The field is off only after a deliberate ramp-down (which takes days) or a "quench," the rapid loss of helium, itself a major adverse event. Accidents have happened precisely because workers assumed a powered-down console meant a powered-down magnet. The most tragic case occurred in 2001, when six-year-old Michael Colombini was killed by an oxygen tank that flew into the bore while he was being scanned — an event that prompted greatly increased scrutiny and standardization of safety procedures. Consensus guidelines from the American College of Radiology (ACR) and others now codify MR safety practice.

The most common risk is injury from **ferromagnetic** material — iron, nickel, cobalt, and alloys such as stainless steel — which the field can turn into projectiles or dislodge from inside the body. A subtle danger is that the pull on a magnetic object scales roughly with the inverse square of distance to the bore: a pen that tugs only faintly six feet away becomes impossible to hold at a foot and a half. The 5 Gauss line marked on the scanner-room floor shows where the field falls to a level that interacts only minimally with magnetic objects — but feeling nothing at that line says nothing about what happens closer in. Implants span the range from MR safe (e.g., titanium joints) through "MR conditional" (safe only under specific protocols) to strictly contraindicated (pacemakers and other electronic implants, or metal fragments in the eyes of metal workers). Careful screening is the main line of defense, with decisions about implants made by certified MR technologists or medical personnel. Even cosmetics matter: tattoos, eyeliner, and hair products containing ferromagnetic particles can distort images and, occasionally, cause burns.

Not all materials behave alike in the field. **Paramagnetic** materials become magnetized in an external field but do not stay magnetized when it is removed — deoxygenated hemoglobin is paramagnetic, which, as later chapters explain, is one of the reasons fMRI works at all. Moving paramagnetic material through the field induces currents and field distortions, one reason head motion is such a persistent problem for image quality. **Diamagnetic** materials — water, wood, plastic, and non-magnetic metals such as copper and gold — are weakly repelled by the field and generally MR safe. A second class of risk comes from RF energy itself: electromagnetic fields deposit heat in tissue, so scanner software estimates and limits the Specific Absorption Rate (SAR) given the pulse sequence and the participant's body mass. But metal conductors inside the RF coil can concentrate induced currents and cause burns even at otherwise harmless power levels — the same physics as sparks from metal in a microwave oven — which is why many centers ask participants to change into scrubs (underwire bras are a classic hazard). Rapidly switching gradients can occasionally cause peripheral nerve stimulation (mild twitching sensations), and they generate acoustic noise above 100 dB at 3T — comparable to a loud rock concert — making earplugs standard issue. For pregnancy, current evidence is reassuring: FDA and ACR guidance considers MRI within normal operating mode safe in pregnant patients for exams under 30 minutes at up to 3T, though 7T fetal imaging is less established and gadolinium contrast is not recommended.

Task fMRI adds another layer: delivering stimuli and recording behavior inside this hostile electromagnetic environment. Visual stimuli are typically projected onto a screen viewed through small mirrors mounted on the head coil (with a viewing angle of only about 15 degrees); audio requires MR-compatible headphones with noise attenuation; and specialized MR-safe equipment exists for thermal, pressure, olfactory, gustatory, and other stimulation. Responses are usually limited to button presses, eye movements, or joystick input, often accompanied by physiological recordings such as pulse and skin conductance. Every cable and device is a potential antenna: unshielded cables and RF leaks can produce large image artifacts, so connections are optical or shielded wherever possible, and it is wise to test image quality with and without any new device in the room. Because stimuli must be synchronized with the fMRI time series at sub-second precision — and operating-system timing interruptions can accumulate into substantial errors over a session — a robust practice is to log stimulus times, responses, physiology, and the start of every image acquisition (each TR) to a separate external recording system.

Finally, the participant is a human being lying in a narrow, noisy tube, asked to hold perfectly still for an hour. Head motion is the enemy of data quality: it changes the local magnetic field in ways that motion-correction algorithms only partially fix, so heads are stabilized with vacuum pillows, foam pads, straps, or bite bars, and tasks requiring overt speech or body movement demand special workarounds. The environment can provoke anxiety and cortisol release, and claustrophobia excludes some participants entirely; mock-scanner familiarization sessions — complete with recorded scanner noise — help acclimate children, patient populations, and anyone anxious about scanning, and also train stillness. The scanner differs from a typical testing room in other ways too: it is often chilly, participants lie supine (which changes cardiovascular function), and after an hour aches and stiffness set in. Yet the environment is not inherently stressful — many participants find it hypnotic, and a substantial fraction fall asleep during resting-state scans, which is why researchers keep sessions short (often under 1.5 hours), intersperse varied tasks, ensure participants are well rested, and commonly allow normal caffeine intake (which also enhances the BOLD response). Reassuringly, thousands of experiments comparing performance and emotional states inside and outside the scanner have found few if any systematic differences — but the best practice is still to pilot paradigms outside the magnet and verify that the phenomena you study survive the trip inside.

## Key ideas in pictures

:::{figure} images/ch11-scanner-components.png
:alt: Cutaway illustration of an MRI scanner showing the magnet, gradient coils, radiofrequency coil, patient table, and patient position
:width: 85%

**The components of an MR scanner.** The three main components are the magnet, which produces the strong static field; the radiofrequency coil, which transmits RF pulses and receives the signals used to construct images; and the gradient coils, which create controlled spatial variations in field strength that encode the location of each signal. The participant lies on the table with the head at the center of the bore. Reproduced under a CC BY license from Amanamba et al. (2020). *(Figure 11.1 from the book.)*
:::

:::{figure} images/ch11-cleaning-cart-bore.png
:alt: A floor cleaning machine pulled into the bore of an MRI scanner and stuck against the magnet opening
:width: 70%

**Cleaning equipment pulled into the bore of a 1.5T magnet.** The static magnetic field is *always on* — even when the console computers are powered down — and its pull grows roughly with the inverse square of distance to the bore, so an object that feels harmless at the door can become an unstoppable projectile near the magnet. Accidents like this one have happened when staff assumed a dark console meant the field was off. Reproduced under a CC BY-SA 4.0 license (Wikimedia Commons, user Xksev). *(Figure 11.2 from the book.)*
:::

**From arrival to scan: the safety workflow.** The diagram below summarizes how MR centers screen participants and manage the environment so that nothing ferromagnetic — on, in, or near a person — gets close to the bore.

```{mermaid}
flowchart TD
    A["Participant arrives"] --> B["Safety screening form + interview<br/>implants, surgeries, metal exposure,<br/>tattoos/cosmetics, pregnancy, claustrophobia"]
    B --> C{"Any implants or<br/>possible metal in body?"}
    C -- "Electronic implant (e.g., pacemaker)<br/>or metal fragments (e.g., in eyes)" --> D["Contraindicated:<br/>do NOT enter MR environment"]
    C -- "Known implant<br/>(e.g., joint, IUD, clips)" --> E["Review by certified MR<br/>technologist / medical personnel:<br/>MR safe vs. MR conditional<br/>(protocol-specific limits)"]
    C -- "No" --> F["Remove all metal: jewelry, watches,<br/>cards, pens; change into scrubs<br/>(no underwires, no magnetic cosmetics)"]
    E -- "Cleared" --> F
    F --> G["Earplugs in<br/>(gradient noise > 100 dB)"]
    G --> H["Enter scanner room past the 5 Gauss line<br/>only with screened, MR-safe items —<br/>remember: pull grows as 1/distance²<br/>and the field is ALWAYS on"]
    H --> I["Position on table: head restraints<br/>(vacuum pillow, pads, straps),<br/>mirror, response devices, squeeze ball"]
    I --> J["Scan: monitor comfort, motion,<br/>and timing synchronization"]
    style A fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    style D fill:#fee2e2,stroke:#ef4444,color:#7f1d1d
    style H fill:#fef3c7,stroke:#f59e0b,color:#713f12
    style J fill:#dcfce7,stroke:#22c55e,color:#14532d
```

## Thought questions

1. Several fatal or near-fatal MR accidents share a common root cause: a correct local observation ("the console is off," "I don't feel any pull here at the door") leading to a false global inference about the field. What features of the MR environment make it so prone to this class of error, and how would you design training, signage, room layout, and access procedures so that safety does not depend on every individual's physical intuition being correct?

2. Head motion corrupts fMRI data, and participants who move too much are often excluded from analysis — yet children, elderly adults, and many patient populations move more than healthy young adults, and anxious participants move more than comfortable ones. Trace the path from "exclusion for motion" to biased scientific conclusions in a concrete example (say, a developmental or psychiatric study), and propose a combination of before-scan (e.g., mock scanner), during-scan, and analysis-stage strategies that reduces the bias rather than just the motion.

3. The chapter argues that the scanner environment — supine posture, cold rooms, 100 dB noise, confinement, cortisol release, and even sleepiness — could in principle alter the psychological states being studied, yet comparisons of performance inside and outside the scanner find few differences. For an emotion or stress study you might design, which of these environmental factors would worry you most, what would "testing the paradigm outside the scanner" concretely look like, and what result would convince you that your phenomenon does *not* survive the trip into the magnet?

4. Precise timing synchronization between stimulus presentation and image acquisition seems like a mundane engineering detail, but the chapter recommends an entire independent recording computer for it. Explain how small, accumulating timing errors would propagate through an event-related fMRI analysis (think about what the model assumes about when each event occurred), why the damage might be invisible in the raw images, and why logging every TR pulse externally provides protection that trusting the stimulus software does not.

5. MRI has no regulatory limit on repeat scanning, no known long-term biological effects, and reassuring safety data in pregnancy up to 3T — yet many research centers still exclude potentially pregnant participants "to avoid unknown risks." Weigh this precautionary policy against its scientific and ethical costs (who is systematically excluded from neuroscience datasets, and what is not learned as a result). Where would you set the policy for a 3T research study, and would your answer change at 7T or with contrast agents?

## Quiz yourself

:::{dropdown} **Q1.** What are the three main hardware components of an MR scanner, and what does each do?
**Answer:** The magnet (a superconducting electromagnet) produces the strong static field; the radiofrequency (RF) coils transmit RF pulses and receive the signals used to construct images; and the gradient coils create controlled spatial variations in field strength that encode the spatial location of MR signals, allowing reconstruction into 3-D images.
:::

:::{dropdown} **Q2.** When is the scanner's main static magnetic field turned off?
**Answer:** Essentially never during normal operation — it stays on nights, weekends, and holidays, even when the console computers are powered down. It is off only after a deliberate ramp-down (which takes days and involves MR engineers) or a "quench," the rapid loss of helium, which is itself a major adverse event.
:::

:::{dropdown} **Q3.** How strong is the field of a typical scanner compared with the Earth's magnetic field?
**Answer:** A standard 1.5T clinical scanner is about 30,000 times stronger than the Earth's field (roughly 0.00005 T), and an ultra-high-field 7T scanner is about 140,000 times stronger — strong enough to pull ferromagnetic objects out of hands and into the bore.
:::

:::{dropdown} **Q4.** What is the difference between ferromagnetic, paramagnetic, and diamagnetic materials in the MR environment?
**Answer:** Ferromagnetic materials (iron, nickel, cobalt, many steels) interact strongly with magnets and can remain magnetized — these are the projectile and implant hazards. Paramagnetic materials (e.g., deoxyhemoglobin) become magnetized only while an external field is applied — a property that underlies the BOLD signal but also makes motion distort images. Diamagnetic materials (water, wood, plastic, copper, gold) are weakly repelled and are generally MR safe.
:::

:::{dropdown} **Q5.** What is the 5 Gauss line, and why can it give a false sense of security?
**Answer:** It is the boundary marked on the scanner-room floor outside which the field is at or below 5 Gauss (1 T = 10,000 Gauss) and interacts only minimally with magnetic objects. Because magnetic pull grows roughly with the inverse square of distance to the bore, an object that feels harmless at the line can become impossible to hold — or an implant can dislodge or malfunction — much closer to the magnet.
:::

:::{dropdown} **Q6.** What is SAR, and how are RF burns possible even when SAR limits are respected?
**Answer:** SAR (Specific Absorption Rate) is the rate at which RF energy is absorbed by tissue; scanner software estimates and limits it based on the pulse sequence and the participant's body mass. However, metal conductors inside the RF coil (implants, metallic cosmetic particles, underwires) can concentrate induced currents and heat to the point of causing burns even at RF power levels that are otherwise harmless — the same principle as metal sparking in a microwave oven.
:::

:::{dropdown} **Q7.** Why are earplugs standard issue for fMRI participants?
**Answer:** Rapid electrical pulses in the gradient coils produce loud tapping and buzzing noises exceeding 100 dB at 3T — comparable to a loud rock concert — which poses a risk of hearing damage over a scan session and adds a distracting element to the environment.
:::

:::{dropdown} **Q8.** What is a mock scanner, and why do researchers use one with children and patient populations?
**Answer:** A mock scanner is a non-magnetic replica of the MR environment, including an enclosed bore and simulated scanner noise. Familiarization sessions in it acclimate participants to the novelty, noise, and confinement — reducing anxiety and claustrophobic reactions — and provide training in lying still, which is especially valuable for children and patients, who tend to move more and are more prone to anxiety in the scanner.
:::
