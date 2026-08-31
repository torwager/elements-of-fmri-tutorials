---
title: Software setup
---

# Software setup

You can use most of this site with **no setup at all** — in-browser and Colab options cover the Python exercises. Set up software locally when you're ready to work with your own data.

## MATLAB + CANlab tools

The MATLAB exercises use the **CANlab Core Tools**, an object-oriented framework for fMRI analysis, along with SPM.

**Requirements:**

1. **MATLAB** (R2020a or newer recommended; R2025a+ for the best live-script experience) with the Statistics and Machine Learning Toolbox and Image Processing Toolbox
2. **SPM** — [SPM25](https://www.fil.ion.ucl.ac.uk/spm/) (or SPM12 with MATLAB ≤ R2023b)
3. **[CanlabCore](https://github.com/canlab/CanlabCore)** — the core object-oriented tools (`fmri_data`, `atlas`, `region`, `statistic_image`, …)
4. **[Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks)** *(optional, large)* — atlases, signature patterns, and masks used in the applied chapters
5. **[CANlab_help_examples](https://github.com/canlab/CANlab_help_examples)** *(optional)* — extended walkthroughs and sample data

**Install:**

```matlab
% Clone the repositories (or download ZIPs from GitHub), then:
addpath(genpath('/path/to/spm'));         % SPM
addpath(genpath('/path/to/CanlabCore'));  % CANlab Core
savepath;

% Verify your setup:
which fmri_data
obj = load_image_set('emotionreg');       % downloads a small sample dataset
plot(obj);
```

Full documentation, walkthroughs, and help are at **[canlab.github.io](https://canlab.github.io)**.

:::{tip}
Each MATLAB exercise page also has an **Open in MATLAB Online** badge that sets all of this up for you in the cloud — a good way to try the tools before installing anything.
:::

## Python + nilearn

The Python exercises primarily use **[nilearn](https://nilearn.github.io)** (statistical learning for neuroimaging), with [nibabel](https://nipy.org/nibabel/), numpy/scipy, pandas, matplotlib, and scikit-learn. Some tutorials also feature **[nltools](https://nltools.org)** from the [Cosan Lab](https://cosanlab.com), whose design is closely aligned with the CANlab MATLAB tools.

```bash
# With conda (recommended):
conda create -n elements-fmri python=3.12
conda activate elements-fmri
pip install nilearn nibabel matplotlib pandas scikit-learn statsmodels jupyterlab

# Optional:
pip install nltools
```

**Verify:**

```python
import nilearn
from nilearn import datasets, plotting
img = datasets.load_mni152_template()
plotting.plot_anat(img)
```

## Data used in the tutorials

Tutorial datasets are small, openly licensed, and downloaded automatically by the exercises themselves (from this site, OpenNeuro, Neurovault, or package data). No manual data setup is needed.

:::{div}
:class: book-tile
![Cover of Elements of Functional Magnetic Resonance Imaging](cover-small.jpg)
**The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::
