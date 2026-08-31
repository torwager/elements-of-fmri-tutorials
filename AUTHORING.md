# Chapter authoring guide

Every chapter page follows this exact structure. Consistency across all 42 chapters matters more than local flourishes.

## File layout

- Chapter page: `partN/chNN-slug.md` (already stubbed; overwrite the stub)
- Python lab (code chapters only): `partN/labs/chNN-lab-python.ipynb` — **executed, outputs saved in the file** (build does not re-execute)
- MATLAB lab (code chapters only): `partN/labs/chNN_lab_matlab.m` — R2025a plain-text live-script format (`%%` sections, `%[text]` markdown comments); shown statically for now

Labs are added to `myst.yml` as `children` of their chapter page.

## Chapter page template (`.md`)

```markdown
---
title: "NN. Chapter Title"
subject: "Part N: Part Title"   # MUST be < 40 characters (MyST validation)
---

# Chapter Title

:::{admonition} What you will learn
:class: tip
3–5 bullet points: the skills/concepts this page teaches.
:::

## Overview

3–8 paragraphs distilling the chapter's key concepts. Written in the book's
voice: warm, precise, practitioner-focused. Use book figures where they help
(copied to `partN/images/`, referenced with a caption crediting the book).
Use real equations (`$$ ... $$`, MyST math) where the chapter has them.
NEVER paste long verbatim passages from the manuscript — distill and rephrase;
short quoted phrases are fine.

## Hands-on tutorial        ← code chapters only

Brief setup sentence, then synced tabs for every exercise block:

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab
(CANlab-tools code, kept runnable as shown; note required toolboxes)
:::
:::{tab-item} Python
:sync: python
(nilearn/numpy code mirroring the MATLAB version)
:::
::::

End the section with a card linking to the full labs:
"**Go deeper:** open the full Python lab notebook [→](./labs/chNN-lab-python.ipynb)
or download the MATLAB live script."
Below it, two badges on one line:
`[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/partN/labs/chNN-lab-python.ipynb)`
`[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=partN/labs/chNN_lab_matlab.m)`
(These resolve once the repo is public at launch.)

For conceptual chapters (2, 3, 4, 5, 9, 10, 11) replace this section with
`## Key ideas in pictures` — 2–4 book figures with rich captions, or a
mermaid/description diagram — no code.

## Thought questions

3–5 open-ended discussion questions as a numbered list. Each should require
integrating chapter concepts, not recall. No answers given.

## Quiz yourself

5–8 factual comprehension questions, each a dropdown:

:::{dropdown} **Q1.** Question text?
**Answer:** The answer, 1–3 sentences, self-contained.
:::

Order from easier to harder. Cover the chapter's main points.
```

## Python lab notebook rules

- First markdown cell: title + one-paragraph intro + "how to run this" note
  (browser/Colab/local + link back to the chapter page)
- If the notebook uses nilearn/nibabel, the FIRST code cell must be:
  `# Setup — needed when running in your browser or on Colab (safe to re-run anywhere)`
  `%pip install -q nilearn`
  (numpy/scipy/pandas/sklearn/statsmodels/networkx need no bootstrap — Pyodide auto-loads them.)
- Then imports only from: numpy, scipy, pandas, matplotlib,
  scikit-learn, statsmodels, nibabel, nilearn (the Pyodide-safe stack).
  `n_jobs=1` everywhere. No internet downloads unless from
  raw.githubusercontent/OSF (CORS-safe), data ≤ a few MB.
- Simulated or bundled small data only; each cell ≤ a few seconds runtime.
- Markdown cells narrate: state the concept, run it, interpret the output.
- Execute before committing:
  `resources/venv (private repo) → jupyter nbconvert --to notebook --execute --inplace <nb>`
- Notebook must run top-to-bottom with no errors.

## MATLAB lab rules

- Plain-text live script: `%%` section per exercise, `%[text]`-style rich
  comments are optional; standard `%%` + comments is fine.
- Use CANlab idioms (`fmri_data`, `load_image_set`, `plot`, `montage`,
  `onsets2fmridesign`, …); assume CanlabCore + SPM on the path.
- Keep runtime under a few minutes on small data (MATLAB Online free tier).

## Run-affordance elements (added site-wide 2026-08-31; include in any new chapter)

- Code chapter pages get, right before `## Overview`: a `🖥️ Ways to run this chapter's code`
  seealso-admonition (in-browser lab + power-icon instructions, Colab/MATLAB-Online badges,
  static-preview disclaimer) AND a `:::{div}` with `:class: run-quick` one-liner
  (floating bottom-right card via custom.css).
- A `{note}` goes immediately before the first tab-set: tabs are static previews; link to the lab.
- Lab notebooks: a markdown banner cell before the first code cell — interactive, ⏻ power icon,
  ~30–60 s first start, run top-to-bottom starting with the setup cell.

## Voice and attribution

- Figures from the book: caption ends with "*(Figure N.M from the book.)*"
- Code adapted from CANlab tutorials keeps a comment crediting the source repo.
- Do not reference "the manuscript", page numbers, or unpublished material.
