---
title: "32. Network Analysis"
subject: "Part 6: Brain Connectivity"
---

# Network Analysis

:::{admonition} What you will learn
:class: tip
- How to represent a brain as a graph — nodes and edges — and the difference between weighted connectivity matrices and binary adjacency matrices, and between undirected and directed networks
- How networks are constructed from functional connectomes by thresholding and binarizing, and what connection density measures
- How to characterize networks at local (degree, clustering, betweenness), mesoscale (communities, modularity), and global (characteristic path length, efficiency, assortativity) levels
- What makes regular, random, and small-world networks different, and how small-worldness is quantified against null networks
- Why threshold choice changes every graph metric, and how density-matched comparisons keep group analyses fair
:::

:::{admonition} 🖥️ Ways to run this chapter's code
:class: seealso
- **In your browser, no setup:** open the [interactive Python lab](./labs/ch32-lab-python.ipynb) and click the **⏻ power icon** at the top right of the notebook. Run cells top-to-bottom, starting with the first (setup/import) cell.
- **In the cloud:** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch32-lab-python.ipynb) · [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch32_lab_matlab.m)
- The code tabs on **this page** are static previews with copy buttons — the labs are where code runs.
:::

:::{div}
:class: run-quick
**Run this code:** [⚡ In-browser lab](./labs/ch32-lab-python.ipynb) · [Colab](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch32-lab-python.ipynb) · [MATLAB Online](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch32_lab_matlab.m)
:::

## Overview

Network analysis represents the brain as a set of interconnected **nodes** (brain regions or voxels) and **edges** (anatomical or functional connections), letting us study brain organization, network properties, and dynamics with the mathematical machinery of graph theory — a field whose roots reach back to Euler in the mid-1700s. In network neuroscience, functional connectomes are used to estimate network topology and higher-order graph-theoretic properties: global properties like efficiency and small-worldness, sub-network structure ("communities" or "modules"), and "hub" nodes that connect many regions or link different communities. All of these properties can vary over time and can be analyzed in relation to cognition, behavior, and disease — for example, comparing network topology between patient and control groups to reveal connectivity abnormalities.

A network is represented mathematically as a **graph** $G(V, E)$: a set of vertices $V$ and pairwise edges $E$ among them. Numerically, a graph is an $N \times N$ matrix, where $N$ is the number of nodes and each element represents an edge. In a **connectivity matrix** the elements are continuous ("weighted") values representing connection strength; in an **adjacency matrix** they are binarized to 1 (connected) or 0 (unconnected). Graphs may also be **undirected** (a symmetric matrix — the edge from $i$ to $j$ equals the edge from $j$ to $i$) or **directed** (asymmetric, with rows conventionally denoting sources and columns targets). Functional connectivity matrices are symmetric, implying undirected networks; effective connectivity methods such as Dynamic Causal Modeling or Granger causality (Chapters 35–36) yield directed ones.

:::{figure} images/ch32_fig2_graph_adjacency.png
:alt: A four-node example graph and its corresponding binary adjacency matrix
:width: 70%

A simple graph and its numerical representation. Four nodes (A–D) are linked by four undirected edges; the corresponding adjacency matrix contains a 1 for each connected pair and 0 otherwise, and is symmetric. From these ingredients we can compute each node's clustering coefficient (network average $C = 0.42$) and the shortest path length between every pair of nodes (characteristic path length $L = 1.33$). *(Figure 32.2 from the book.)*
:::

**Constructing** a brain network means filling in the entries of that matrix. In fMRI the matrix is typically a functional connectome (Chapter 30): correlations between the time series of all node pairs, using Pearson or Spearman correlations, partial correlations, coherence, or mutual information. Because many graph metrics require relatively sparse graphs, the next step is usually to **threshold** the matrix — removing weak connections — and often to binarize it into an adjacency matrix. Thresholds may be a fixed value, a fixed percentage of retained connections, or a statistical significance level (which raises a multiple comparisons problem across the $N(N-1)/2$ edge tests, analogous to mass-univariate GLM analysis). The **connection density** — surviving edges divided by possible edges — measures the resulting sparsity. Thresholding removes weak, potentially spurious connections, sharpens topological structure, and eases computation, but it has real disadvantages: the "right" threshold is rarely obvious, multiplicity grows rapidly with $N$, and — critically — subjects who differ in overall connectivity strength end up with different densities at a fixed threshold, which by itself changes every metric computed downstream. The standard fix is to vary the threshold per subject so that the number of edges (or the density) is held fixed.

Once built, a network can be **characterized** at several topological scales. *Local* measures describe individual nodes. A node's **degree** counts its neighbors,

$$
k_i = \sum_{j \neq i} a_{ij},
$$

and is the most common index of *centrality* — a node's importance for information transfer — and hence of "hub status" (with **strength** as the weighted analog). A node's **clustering coefficient** asks how many of its neighbors are also neighbors of each other: if nodes $i$, $j$, and $h$ are all interconnected they form a closed triangle, and with $t_i$ the number of triangles around node $i$,

$$
C_i = \frac{2\, t_i}{k_i (k_i - 1)}, \qquad C = \frac{1}{N} \sum_i C_i,
$$

where $C$ is the clustering coefficient of the whole network — a measure of *functional segregation*: specialized processing within densely interconnected groups. High clustering implies redundant connections, so losing one node matters less. **Betweenness centrality** counts the fraction of shortest paths between all node pairs that pass through a node, and the **participation coefficient** measures how evenly a node's edges are spread across communities — distinguishing *provincial hubs* (well connected within their own module) from *connector hubs* (linking different modules). *Global* measures summarize the whole network. The **shortest path length** $d_{ij}$ is the minimum number of edges (or minimum total weight) needed to travel between two nodes, and the **characteristic path length**

$$
L = \frac{1}{N} \sum_i L_i
$$

(where $L_i$ is the average distance from node $i$ to all others) indexes *functional integration*: short paths mean information can be combined rapidly across the network. Related global measures include global efficiency, assortativity (the tendency of high-degree nodes to connect to each other, forming densely interconnected "rich clubs"), and resilience — the ability to keep functioning as nodes and edges are removed. In between sit *mesoscale* properties: **communities** (modules) detected by algorithms such as Louvain, Girvan–Newman, or Clauset–Newman–Moore modularity optimization, hierarchical or spectral clustering, and stochastic block models; and, in multilayer networks spanning time points or data types, **recruitment** and **integration** of nodes across layers.

These metrics define characteristic **types of networks**. In a *regular* network (a lattice or ring) every node has the same degree: clustering $C$ is high but so is path length $L$ — redundant local communities, inefficient global transmission. In a *random* (Erdős–Rényi) network every node pair connects with equal probability $p$: both $C$ and $L$ are low — information travels easily, but the structure is vulnerable and unclustered. A **small-world** network resembles an ordered network with a few randomly rewired links: high $C$ *and* low $L$, making it simultaneously resilient and efficient. Since Watts and Strogatz's seminal observation that many social, biological, and technological networks share this property, the brain too has been argued to be small-world — maximizing segregation and integration while minimizing wiring cost. Small-worldness is quantified against a null network:

$$
\sigma = \frac{C / C_{rand}}{L / L_{rand}},
$$

with $\sigma > 1$ commonly taken as the cutoff. Other overall structures include core–periphery organization and disassortative networks, in which high-degree nodes preferentially connect to low-degree ones.

Finally, networks can be **compared** — patients versus controls, or network properties against behavior. So long as the graphs share the same set of nodes, any graph metric can serve as a dependent variable in a GLM regressed on group membership or cognitive performance, or as a feature in machine-learning models (Chapters 37–41). Group differences can be tested nonparametrically with permutation tests that shuffle group labels to build a null distribution for the difference in a metric. And to ask whether a network's organization is non-random — small-world, core–periphery — the benchmark is a population of simulated random graphs with the same connection density and number of connected nodes, yielding a null distribution for each topological metric of interest.

## Hands-on tutorial

In this tutorial you will build a graph from a simulated modular functional connectome, compute the core metrics from the chapter, and see first-hand why threshold choice matters. We simulate 30 nodes organized into three modules (plus a shared global signal that weakly links them), correlate their time series, threshold, and analyze the result.

**Step 1 — From connectivity matrix to graph and metrics.** We compute the correlation matrix, binarize it at $r > 0.3$, and compute density, degree, clustering, and characteristic path length.

:::{note}
The tabs below are **static previews** (with copy buttons) showing the key step in each language. To run and modify this code, use the [interactive in-browser lab](./labs/ch32-lab-python.ipynb) or the Colab / MATLAB Online links above.
:::

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
% Base MATLAB graph objects; CanlabCore on path for the fuller lab
rng(7);
n = 30; T = 200;
module  = repelem(1:3, 10)';                 % 3 modules of 10 nodes
g       = randn(1, T);                       % shared global signal
signals = randn(3, T);                       % one latent signal per module
ts = 0.8 * repmat(g, n, 1) + signals(module, :) + 1.2 * randn(n, T);

R = corr(ts');                               % 30 x 30 functional connectome
A = double(R > 0.3);  A(1:n+1:end) = 0;      % threshold + binarize
G = graph(A);

density = numedges(G) / nchoosek(n, 2)       % connection density
deg = degree(G);                             % node degree k_i

t_i = diag(A^3) / 2;  k = sum(A, 2);         % triangles around each node
Ci  = 2 * t_i ./ (k .* (k - 1));  Ci(k < 2) = 0;
C   = mean(Ci)                               % network clustering coefficient

D = distances(G);                            % shortest path lengths d_ij
comp = conncomp(G);  in_cc = comp == mode(comp);
Dcc  = D(in_cc, in_cc);                      % largest connected component
L    = mean(Dcc(Dcc > 0))                    % characteristic path length
```
:::
:::{tab-item} Python
:sync: python

```python
import numpy as np, networkx as nx

rng = np.random.default_rng(7)
n, T = 30, 200
module  = np.repeat([0, 1, 2], 10)            # 3 modules of 10 nodes
g       = rng.standard_normal(T)              # shared global signal
signals = rng.standard_normal((3, T))         # one latent signal per module
ts = 0.8 * g + signals[module] + 1.2 * rng.standard_normal((n, T))

R = np.corrcoef(ts)                           # 30 x 30 functional connectome
A = (R > 0.3).astype(int)                     # threshold + binarize
np.fill_diagonal(A, 0)
G = nx.from_numpy_array(A)

density = nx.density(G)                       # connection density
degree  = dict(G.degree())                    # node degree k_i
C       = nx.average_clustering(G)            # network clustering coefficient

Gcc = G.subgraph(max(nx.connected_components(G), key=len))
L   = nx.average_shortest_path_length(Gcc)    # characteristic path length
print(f"density={density:.2f}  C={C:.2f}  L={L:.2f}")
```
:::
::::

**Step 2 — The threshold changes everything.** Sweep the threshold and watch density, clustering, and path length co-vary — none of these metrics is interpretable without knowing the density it was computed at.

::::{tab-set}
:::{tab-item} MATLAB
:sync: matlab

```matlab
thr = 0.10:0.05:0.50;
[dens, Cs, Ls] = deal(zeros(size(thr)));
for i = 1:numel(thr)
    A = double(R > thr(i));  A(1:n+1:end) = 0;
    G = graph(A);
    dens(i) = numedges(G) / nchoosek(n, 2);
    t_i = diag(A^3)/2;  k = sum(A, 2);
    Ci = 2*t_i ./ (k .* (k-1));  Ci(k < 2) = 0;  Cs(i) = mean(Ci);
    D = distances(G);  comp = conncomp(G);
    Dcc = D(comp == mode(comp), comp == mode(comp));
    Ls(i) = mean(Dcc(Dcc > 0));
end
plot(thr, dens, '-o', thr, Cs, '-s', thr, Ls / max(Ls), '-^');
legend({'Density', 'Clustering C', 'Path length L (scaled)'});
xlabel('Correlation threshold');
```
:::
:::{tab-item} Python
:sync: python

```python
import matplotlib.pyplot as plt

thresholds = np.arange(0.10, 0.51, 0.05)
dens, Cs, Ls = [], [], []
for thr in thresholds:
    A = (R > thr).astype(int); np.fill_diagonal(A, 0)
    G = nx.from_numpy_array(A)
    dens.append(nx.density(G))
    Cs.append(nx.average_clustering(G))
    Gcc = G.subgraph(max(nx.connected_components(G), key=len))
    Ls.append(nx.average_shortest_path_length(Gcc))

plt.plot(thresholds, dens, "-o", label="Density")
plt.plot(thresholds, Cs, "-s", label="Clustering C")
plt.plot(thresholds, np.array(Ls) / max(Ls), "-^", label="Path length L (scaled)")
plt.xlabel("Correlation threshold"); plt.legend()
```
:::
::::

Every curve moves together: a "group difference in clustering" at a fixed threshold may be nothing more than a group difference in density. The full labs push further: they plant a **connector hub** in the network and find it with degree and betweenness centrality, detect communities with **greedy modularity optimization** and score them against the ground-truth modules, test **small-worldness** ($\sigma$) against density-matched random graphs, and demonstrate that comparing two "subjects" with different overall connectivity strength is confounded at a fixed threshold but fair after **density matching**.

:::{card} **Go deeper**
Open the full Python lab notebook [→](./labs/ch32-lab-python.ipynb) or download the [MATLAB live script](./labs/ch32_lab_matlab.m), which mirrors it using MATLAB graph objects and CANlab-style workflows.
:::

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/torwager/elements-of-fmri-tutorials/blob/main/part6/labs/ch32-lab-python.ipynb)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=torwager/elements-of-fmri-tutorials&file=part6/labs/ch32_lab_matlab.m)

## Thought questions

1. Node degree, betweenness centrality, and participation coefficient can all be used to call a region a "hub" — yet they can disagree sharply for the same node. Construct (or sketch) a network in which a node has modest degree but very high betweenness, and explain what functional role such a region might play in the brain. Why does the provincial/connector distinction matter for predicting the consequences of a lesion?
2. A study reports that patients have lower clustering coefficients than controls, using a fixed correlation threshold of $r > 0.25$ applied to every subject. The patients also show globally weaker functional connectivity (e.g., due to greater head motion or vascular differences). Walk through the chain of consequences from weaker correlations to the reported group difference. What analyses would convince you the topological difference is real?
3. Small-worldness $\sigma$ compares a network's clustering and path length to those of random graphs. What properties should the null networks be matched on (density? degree sequence? connectedness?), and how could a poor choice of null model manufacture — or hide — small-world structure?
4. Functional connectomes based on correlation are undirected, so all the metrics in this chapter ignore the direction of influence. For which network properties (hubs, communities, efficiency) do you think ignoring directionality is most misleading, and how might metrics computed on a directed effective-connectivity graph (Chapters 35–36) change the picture?
5. Community-detection algorithms will happily partition *any* network — including a random one — into "modules." Drawing on the network-comparison logic of the chapter, design a test for whether the modular structure detected in a resting-state connectome is statistically meaningful rather than an artifact of the algorithm.

## Quiz yourself

:::{dropdown} **Q1.** In a brain network, what do nodes and edges typically represent, and what distinguishes a connectivity matrix from an adjacency matrix?
**Answer:** Nodes are brain regions (or voxels) and edges are the functional or structural connections between them. A connectivity matrix holds continuous weighted values representing connection strength; an adjacency matrix is binarized, with 1 for connected and 0 for unconnected node pairs.
:::

:::{dropdown} **Q2.** Why are functional connectivity networks undirected, and which methods produce directed brain networks instead?
**Answer:** Functional connectivity is based on symmetric measures like correlation — the value linking $i$ to $j$ equals the one linking $j$ to $i$ — so the matrix is symmetric and the graph undirected. Effective connectivity techniques such as Dynamic Causal Modeling and Granger causality estimate directional influences and yield directed (asymmetric) networks.
:::

:::{dropdown} **Q3.** What is connection density, and why is thresholding usually applied before computing graph metrics?
**Answer:** Density is the number of edges surviving thresholding divided by the maximum possible number of edges — a measure of sparsity. Thresholding removes weak, potentially spurious connections, emphasizes topological structure, and reduces computational burden, and many graph metrics assume relatively sparse graphs.
:::

:::{dropdown} **Q4.** Define a node's degree and clustering coefficient, and say what aspect of brain function each is used to index.
**Answer:** Degree $k_i = \sum_{j \neq i} a_{ij}$ counts a node's neighbors and indexes centrality or hub status. The clustering coefficient $C_i = 2t_i / (k_i(k_i-1))$ measures the fraction of a node's neighbor pairs that are themselves connected (closed triangles); averaged over nodes it indexes functional segregation — specialized processing in densely interconnected groups.
:::

:::{dropdown} **Q5.** What is the characteristic path length, and what does a short value imply about a network?
**Answer:** It is the average shortest path length between all pairs of nodes, $L = \frac{1}{N}\sum_i L_i$. A short characteristic path length implies efficient functional integration — information can be combined rapidly across distributed regions because few edges must be traversed.
:::

:::{dropdown} **Q6.** Contrast regular, random, and small-world networks in terms of clustering coefficient $C$ and characteristic path length $L$.
**Answer:** Regular networks (lattices/rings) have high $C$ and high $L$ — strong local clustering but inefficient global transmission. Random (Erdős–Rényi) networks have low $C$ and low $L$. Small-world networks combine the best of both: high $C$ with low $L$, quantified by $\sigma = (C/C_{rand})/(L/L_{rand}) > 1$.
:::

:::{dropdown} **Q7.** What is the difference between a provincial hub and a connector hub?
**Answer:** A provincial hub's connections are mostly local, within its own community or module; a connector hub links nodes in different communities. The participation coefficient — how evenly a node's edges are distributed across communities — distinguishes them.
:::

:::{dropdown} **Q8.** When comparing graph metrics between patient and control groups, why is it problematic if the groups differ in connection density, and what are two remedies?
**Answer:** Essentially every graph metric depends on density, so a density difference alone can masquerade as a topological difference. Remedies: vary the threshold per subject to fix the number of edges or the connection density (density matching), and use nonparametric permutation tests — shuffling group labels to build a null distribution for the metric difference.
:::

:::{div}
:class: book-tile
📖 **The book:** [*Elements of Functional Magnetic Resonance Imaging*](https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/) — Wager & Lindquist, MIT Press
:::

---

[⌂ Back to home](https://torwager.github.io/elements-of-fmri-tutorials/) · [Table of contents](../contents.md) · [How to use this site](../how-to-use.md)
