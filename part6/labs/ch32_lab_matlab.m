%% Chapter 32 Lab - Graph Theory and Network Analysis (MATLAB)
% Elements of fMRI Analysis - Interactive Tutorials
%
% In this lab we build brain-style networks from a simulated modular
% functional connectome and characterize them with graph theory:
% degree, clustering, path length, betweenness, communities, modularity,
% small-worldness, and the effect of threshold choice on every metric.
%
% Uses base MATLAB graph objects (graph, degree, distances, conncomp,
% centrality). CANlab-style conventions; for real data, region-level time
% series can be extracted with CanlabCore's fmri_data/extract_roi_averages
% or the @brainpathway class, which computes region-by-region connectivity
% automatically. Adapted from CANlab tutorials (github.com/canlab).

%% 1. Simulate a modular functional connectome
% Thirty nodes in three modules of ten. Each node's time series mixes:
% a global signal shared by all nodes (weakly linking the modules),
% a module-specific latent signal (creating community structure), and
% independent noise. Node 1 is rewired into a CONNECTOR HUB: it carries
% all three module signals, so it correlates with every community.

rng(7);
n = 30; T = 200;
module  = repelem(1:3, 10)';           % ground-truth community labels
g       = randn(1, T);                 % shared global signal
signals = randn(3, T);                 % one latent signal per module

ts = 0.8 * repmat(g, n, 1) + signals(module, :) + 1.2 * randn(n, T);

% Plant a connector hub at node 1: equal parts of all three module signals
ts(1, :) = 0.8 * g + sum(signals, 1) / sqrt(3) + 0.6 * randn(1, T);

R = corr(ts');                         % 30 x 30 functional connectome

figure;
imagesc(R, [-1 1]); axis square; colorbar;
title('Simulated functional connectome (sorted by module)');
xlabel('Node'); ylabel('Node');

%% 2. Threshold, binarize, and build the graph
% Keep edges with r > 0.3 and binarize into an adjacency matrix.
% Connection density = surviving edges / possible edges.

thr0 = 0.3;
A = double(R > thr0);
A(1:n+1:end) = 0;                      % zero the diagonal (no self-loops)
G = graph(A);

density = numedges(G) / nchoosek(n, 2);
fprintf('Edges: %d of %d possible; density = %.2f\n', ...
    numedges(G), nchoosek(n, 2), density);

figure;
h = plot(G, 'Layout', 'force');        % spring-style layout
h.NodeCData = module;                  % color nodes by true module
colormap(lines(3));
title('Thresholded network (r > 0.3), colored by true module');

%% 3. Local metrics: degree, clustering, betweenness - find the hubs
% Degree k_i counts each node's neighbors (hub status). The clustering
% coefficient C_i = 2*t_i / (k_i*(k_i-1)) counts closed triangles.
% Betweenness centrality counts shortest paths through each node -
% connector hubs that bridge modules score high.

deg = degree(G);

t_i = diag(A^3) / 2;                   % triangles around each node
k   = sum(A, 2);
Ci  = 2 * t_i ./ (k .* (k - 1));
Ci(k < 2) = 0;

btw = centrality(G, 'betweenness');

[~, order] = sort(btw, 'descend');
fprintf('\nTop 5 nodes by betweenness centrality:\n');
fprintf('%6s %8s %8s %12s\n', 'node', 'degree', 'C_i', 'betweenness');
for i = 1:5
    j = order(i);
    fprintf('%6d %8d %8.2f %12.1f\n', j, deg(j), Ci(j), btw(j));
end
% Node 1 (the planted connector hub) should top the list: it bridges the
% three modules, so many shortest paths pass through it.

%% 4. Global metrics: clustering, path length, efficiency
% The network clustering coefficient C (segregation) is the mean of C_i.
% The characteristic path length L (integration) is the mean shortest
% path length; we compute it on the largest connected component.

C = mean(Ci);

D    = distances(G);                   % all shortest path lengths d_ij
comp = conncomp(G);
in_cc = comp == mode(comp);            % largest connected component
Dcc  = D(in_cc, in_cc);
L    = mean(Dcc(Dcc > 0));

% Global efficiency: mean of 1/d_ij over all pairs (0 for disconnected)
invD = 1 ./ D;  invD(1:n+1:end) = 0;  invD(isinf(D)) = 0;
Eglob = sum(invD(:)) / (n * (n - 1));

fprintf('\nClustering C = %.2f | Path length L = %.2f | Efficiency = %.2f\n', ...
    C, L, Eglob);

%% 5. Small-worldness: compare against density-matched random graphs
% A small-world network is more clustered than a random network with the
% same density, but has a similar path length:
%   sigma = (C/C_rand) / (L/L_rand) > 1.
% We build 20 Erdos-Renyi random graphs with the same number of nodes
% and edges and average their C and L.

n_rand = 20;
m = numedges(G);
[C_rand, L_rand] = deal(zeros(n_rand, 1));
for r = 1:n_rand
    % Random adjacency with exactly m edges among n nodes
    pairs = nchoosek(1:n, 2);
    sel = pairs(randperm(size(pairs, 1), m), :);
    Ar = zeros(n);
    Ar(sub2ind([n n], sel(:,1), sel(:,2))) = 1;
    Ar = Ar + Ar';
    Gr = graph(Ar);

    tr = diag(Ar^3) / 2;  kr = sum(Ar, 2);
    Cir = 2 * tr ./ (kr .* (kr - 1));  Cir(kr < 2) = 0;
    C_rand(r) = mean(Cir);

    Dr = distances(Gr);  compr = conncomp(Gr);
    Drc = Dr(compr == mode(compr), compr == mode(compr));
    L_rand(r) = mean(Drc(Drc > 0));
end

sigma = (C / mean(C_rand)) / (L / mean(L_rand));
fprintf('C/C_rand = %.2f, L/L_rand = %.2f, sigma = %.2f (>1 => small-world)\n', ...
    C / mean(C_rand), L / mean(L_rand), sigma);

%% 6. Community detection and modularity
% Base MATLAB has no built-in modularity optimizer, so we take a
% CANlab-style route: hierarchical clustering on the correlation-derived
% distance (1 - R), then score the partition with Newman's modularity Q.
% (In the Python lab we use greedy modularity optimization instead;
% dedicated tools include the Brain Connectivity Toolbox's Louvain.)

dist_vec = squareform(1 - R, 'tovector');
Z = linkage(dist_vec, 'average');
labels = cluster(Z, 'maxclust', 3);

Q_detected = modularity_q(A, labels);
Q_true     = modularity_q(A, module);
fprintf('\nModularity Q: detected partition = %.2f, true modules = %.2f\n', ...
    Q_detected, Q_true);

% Agreement with ground truth (up to label permutation)
fprintf('Cross-tab of detected communities x true modules:\n');
disp(crosstab(labels, module));

figure;
h = plot(G, 'Layout', 'force');
h.NodeCData = labels;
colormap(lines(3));
title(sprintf('Detected communities (Q = %.2f)', Q_detected));

%% 7. The threshold changes everything
% Sweep the threshold from 0.10 to 0.50: density, clustering, and path
% length all move together. No metric is interpretable without knowing
% the density at which it was computed.

thr = 0.10:0.05:0.50;
[dens, Cs, Ls] = deal(zeros(size(thr)));
for i = 1:numel(thr)
    Ai = double(R > thr(i));  Ai(1:n+1:end) = 0;
    Gi = graph(Ai);
    dens(i) = numedges(Gi) / nchoosek(n, 2);
    ti = diag(Ai^3)/2;  ki = sum(Ai, 2);
    Cii = 2*ti ./ (ki .* (ki-1));  Cii(ki < 2) = 0;
    Cs(i) = mean(Cii);
    Di = distances(Gi);  ci = conncomp(Gi);
    Dic = Di(ci == mode(ci), ci == mode(ci));
    Ls(i) = mean(Dic(Dic > 0));
end

figure;
plot(thr, dens, '-o', thr, Cs, '-s', thr, Ls / max(Ls), '-^', 'LineWidth', 1.5);
legend({'Density', 'Clustering C', 'Path length L (scaled)'}, 'Location', 'best');
xlabel('Correlation threshold'); ylabel('Metric value');
title('Every graph metric depends on the threshold');

%% 8. Density-matched comparison of two "subjects"
% Subject B has the SAME network architecture but weaker overall
% connectivity (more noise). At a fixed threshold, B's network is much
% sparser, and its metrics differ - a spurious "topological" difference.
% Matching density (keeping the same fraction of strongest edges in
% both subjects) makes the comparison fair.

tsB = 0.8 * repmat(g, n, 1) + signals(module, :) + 2.0 * randn(n, T);
tsB(1, :) = 0.8 * g + sum(signals, 1) / sqrt(3) + 1.2 * randn(1, T);
RB = corr(tsB');

% (a) Fixed threshold r > 0.3
AA = double(R  > thr0);  AA(1:n+1:end) = 0;
AB = double(RB > thr0);  AB(1:n+1:end) = 0;
fprintf('\nFixed threshold r > %.1f:\n', thr0);
fprintf('  Subject A: density %.2f, C = %.2f\n', ...
    sum(AA(:))/(n*(n-1)), local_C(AA));
fprintf('  Subject B: density %.2f, C = %.2f  <- sparser AND less clustered\n', ...
    sum(AB(:))/(n*(n-1)), local_C(AB));

% (b) Density-matched: keep the strongest 20% of edges in each subject
target_density = 0.20;
AAm = top_edges(R,  target_density);
ABm = top_edges(RB, target_density);
fprintf('Density matched at %.2f:\n', target_density);
fprintf('  Subject A: density %.2f, C = %.2f\n', ...
    sum(AAm(:))/(n*(n-1)), local_C(AAm));
fprintf('  Subject B: density %.2f, C = %.2f  <- now comparable\n', ...
    sum(ABm(:))/(n*(n-1)), local_C(ABm));

% With density matched, both subjects show similar topology - as they
% should, because their underlying architecture is identical.

%% Local functions

function Q = modularity_q(A, labels)
% Newman's modularity for a binary undirected adjacency matrix
m2 = sum(A(:));                        % 2m (each edge counted twice)
k  = sum(A, 2);
Q  = 0;
for c = unique(labels)'
    idx = labels == c;
    Q = Q + sum(sum(A(idx, idx))) / m2 - (sum(k(idx)) / m2)^2;
end
end

function C = local_C(A)
% Mean clustering coefficient of a binary undirected adjacency matrix
t = diag(A^3) / 2;  k = sum(A, 2);
Ci = 2 * t ./ (k .* (k - 1));  Ci(k < 2) = 0;
C = mean(Ci);
end

function A = top_edges(R, target_density)
% Binarize R keeping the top fraction of strongest off-diagonal edges
n = size(R, 1);
Ru = R;  Ru(logical(tril(ones(n)))) = -Inf;   % upper triangle only
vals = sort(Ru(isfinite(Ru)), 'descend');
n_keep = round(target_density * n * (n - 1) / 2);
thr = vals(n_keep);
A = double(R >= thr);  A(1:n+1:end) = 0;
end
