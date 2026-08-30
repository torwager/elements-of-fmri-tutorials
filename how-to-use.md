---
title: How to use this site
---

# How to use this site

Every chapter page follows the same structure, so you always know what to expect:

1. **Overview** — the chapter's key concepts, distilled, with figures from the book
2. **Hands-on tutorial** — runnable code demonstrating the main ideas (for chapters where code applies; some early chapters are conceptual and use figures and questions instead)
3. **Thought questions** — open-ended discussion questions, great for classes and journal clubs
4. **Quiz yourself** — factual comprehension questions; click any question to reveal the answer

## Running the code

You have three ways to run tutorial code, from zero-setup to full control:

:::{list-table}
:header-rows: 1

* - Option
  - What it is
  - When to use it
* - **▶ In-browser**
  - Python runs directly in your browser via WebAssembly — no account, no installation, nothing to maintain. The first cell takes ~30–90 s to set up, then everything is fast.
  - Simulations, GLM exercises, connectivity and machine-learning demos on small datasets.
* - **☁️ Open in Colab**
  - Opens the notebook on Google Colab (free Google account required).
  - Exercises using full-size datasets or heavier computation.
* - **💻 Run locally**
  - Download the notebook / script and run it with your own installation (see [Setup](./setup.md)).
  - Your own research workflow; MATLAB exercises with CANlab tools.
:::

Even if you never run anything, **every page shows the complete code with its actual output** — figures, statistics, and brain maps are rendered right on the page.

## MATLAB and Python tabs

Most exercises are provided in both languages using tabbed panels:

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% MATLAB version (CANlab tools)
obj = load_image_set('emotionreg');
plot(obj);
```
:::
:::{tab-item} Python
:sync: python

```python
# Python version (nilearn)
from nilearn import plotting
plotting.plot_stat_map(img)
```
:::
::::

The tabs are **synchronized**: choose MATLAB or Python once, and every tabbed example on the site switches to your choice.

MATLAB exercises come with an **Open in MATLAB Online** badge — one click opens the exercise in MATLAB Online with CANlab tools set up automatically. A free MathWorks account works (20 hours/month); academic users whose institution has a campus license get the full version.

## Quiz answers

Quiz questions hide their answers until you click:

:::{dropdown} **Q:** What does BOLD stand for?
**A:** Blood Oxygen Level-Dependent — the fMRI contrast mechanism sensitive to changes in blood oxygenation that accompany neural activity.
:::

Try it! Click the question above.
