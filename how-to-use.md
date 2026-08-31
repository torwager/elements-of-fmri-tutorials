---
title: How to use this site
---

# How to use this site

Every chapter page follows the same structure, so you always know what to expect:

1. **Overview** — the chapter's key concepts, distilled, with figures from the book
2. **Hands-on tutorial** — compact MATLAB/Python code examples for the main ideas (for chapters where code applies; some early chapters are conceptual and use figures and questions instead)
3. **Thought questions** — open-ended discussion questions, great for classes and journal clubs
4. **Quiz yourself** — factual comprehension questions; click any question to reveal the answer

## Chapter pages vs. labs — where code actually runs

Each code chapter has **two levels**:

:::{list-table}
:header-rows: 1

* - Page type
  - What it contains
  - Runnable?
* - **Chapter page** (e.g., *18. The General Linear Model*)
  - Concept overview + compact code **previews** in synced MATLAB/Python tabs, with copy buttons
  - No — the tabs are static previews of the key steps
* - **Lab page** (nested under the chapter in the sidebar, e.g., *Chapter 18 Lab — Python*)
  - The complete worked notebook, with **every cell's real output** (results and figures) shown
  - **Yes** — in your browser, on Colab, or locally
:::

Each code chapter page has a **"🖥️ Ways to run this chapter's code"** box near the top, linking to its lab and cloud options — and a quick-run shortcut card in the bottom-right corner on large screens.

## ▶ Running a lab in your browser (recommended — zero setup)

Python labs run **directly in your browser** via WebAssembly (JupyterLite/Pyodide). Nothing to install, no account needed:

1. Open a lab page and click the **⏻ power icon** at the **top right of the notebook content** (labeled *"Click for interactive code"*). On small windows the icon can be hidden — widen the browser window if you don't see it.
2. Wait for the kernel to start — the **first start takes ~30–60 seconds** while Python downloads into your browser (it's cached afterward, so later starts are fast).
3. **Run cells top to bottom, starting with the first code cell.** The first cell sets up imports (and for some labs installs `nilearn` via `%pip`). Later cells depend on earlier ones — if you jump ahead you'll get `NameError`s. If things get confused, use the restart button and run again from the top.

Once the kernel is live, each cell gets its own ▶ run button and a **run-all** control appears in the notebook toolbar.

**Cells are fully editable.** Once the kernel is running, click into any code cell, change parameters or code, and re-run it — this is the best way to build intuition (What happens with a smaller sample? A different threshold?). Your edits live only in your browser tab; reloading the page restores the original.

**If a cell hangs or you want to start over:** there is no per-cell stop button, but the **↺ restart control** in the notebook toolbar interrupts everything and gives you a fresh kernel — then run again from the top. Note that the *first* run of a lab can look "hung" for up to a minute while Python and packages download into your browser; the page may even be briefly unresponsive. That's normal — it's cached and fast afterward.

:::{tip}
Everything in the labs is simulated or tiny data — cells run in seconds, and nothing leaves your machine.
:::

## ☁️ Other ways to run

:::{list-table}
:header-rows: 1

* - Option
  - What it is
  - When to use it
* - **Open in Colab**
  - The badge on each chapter's run box opens the Python lab on Google Colab (free Google account required).
  - If you prefer a full Jupyter environment, or your browser struggles with WebAssembly.
* - **Open in MATLAB Online**
  - Opens the chapter's MATLAB live script in MATLAB Online. A free MathWorks account works (20 h/month); academics with a campus license get the full version.
  - For the MATLAB versions of the exercises with CANlab tools.
* - **💻 Run locally**
  - Download the notebook / live script (download icon in the page toolbar) and run with your own installation — see [Setup](./setup.md).
  - Your own research workflow, your own data.
:::

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

## Quiz answers

Quiz questions hide their answers until you click:

:::{dropdown} **Q:** What does BOLD stand for?
**A:** Blood Oxygen Level-Dependent — the fMRI contrast mechanism sensitive to changes in blood oxygenation that accompany neural activity.
:::

Try it! Click the question above.

:::{div}
:class: book-tile
![Cover of Elements of Functional Magnetic Resonance Imaging](cover-small.jpg)
**The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::
