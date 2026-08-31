%% Chapter 42 lab — AI and Neuroscience (MATLAB)
%[text] In this lab you will build the two computational bridges between artificial
%[text] and biological networks described in Chapter 42, using base MATLAB only —
%[text] no deep learning toolboxes required. Part 1 simulates error-driven learning
%[text] with the Rescorla–Wagner and temporal difference (TD) models and shows that
%[text] the TD prediction error reproduces the classic phasic dopamine signal: a
%[text] burst to unexpected reward that transfers to the predictive cue with
%[text] learning, and a dip when an expected reward is omitted. Part 2 is a toy
%[text] representational similarity analysis (RSA): we train a tiny multilayer
%[text] network with hand-coded backpropagation and show that training for a task
%[text] is what makes a layer's representational geometry match a category-coding
%[text] "brain region" — the logic of the Yamins et al. experiment (Figure 42.3).
%[text]
%[text] Adapted from the CANlab Computational Foundations course (github.com/canlab).
%[text] Runs in base MATLAB (CanlabCore is not needed for this lab).
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part7/ch42-ai-and-neuroscience

%% Part 1.1 — Rescorla–Wagner: acquisition and extinction
%[text] The agent keeps a value estimate V for a cue and updates it from the
%[text] prediction error, the mismatch between reward received and reward expected:
%[text]
%[text] $$V_{t+1} = V_t + \alpha\,\delta_t, \qquad \delta_t = R_t - V_t$$
%[text]
%[text] where $V_t$ is the value (expected reward) on trial $t$, $R_t$ the reward
%[text] received (1 or 0), $\delta_t$ the prediction error, and $\alpha$ the
%[text] learning rate (0–1). A cue is rewarded on 80% of trials for 80 trials (acquisition), then
%[text] reward stops (extinction).

rng(42);                               % seed for reproducibility
alpha = 0.15;                          % alpha = learning rate (0-1)
n_trials = 120;                        % 80 acquisition + 40 extinction trials

r = double(rand(n_trials, 1) < 0.8);   % cue -> reward on 80% of trials
r(81:end) = 0;                         % extinction after trial 80

V = zeros(n_trials + 1, 1);            % value estimate
pe = zeros(n_trials, 1);               % prediction errors
for t = 1:n_trials
    pe(t) = r(t) - V(t);               % prediction error (dopamine-like)
    V(t + 1) = V(t) + alpha * pe(t);   % value update
end

figure('Color', 'w');
subplot(1, 2, 1);
plot(V, 'LineWidth', 2); hold on;
xline(80, '--', 'Color', [.5 .5 .5]);
title('Learned value V'); xlabel('Trial'); ylabel('V');
subplot(1, 2, 2);
plot(pe, 'r'); hold on; yline(0, 'k'); xline(80, '--', 'Color', [.5 .5 .5]);
title('Prediction error \delta'); xlabel('Trial'); ylabel('\delta');

%[text] Early in acquisition, rewards are unexpected: prediction errors are large
%[text] and positive and V climbs toward the true reward rate (0.8). As V
%[text] converges, the errors shrink. At extinction every trial delivers less than
%[text] expected, producing negative prediction errors that drive V back to zero.
%[text] If \delta were a neuron's firing relative to baseline, you would see bursts
%[text] early in learning and dips at extinction — already dopamine-like.

%% Part 1.2 — The learning rate shapes the value trajectory
%[text] A small alpha learns slowly but averages over many outcomes; a large alpha
%[text] chases each recent outcome (fast but volatile, and fast to extinguish).
%[text] Fitted to a participant's choices, this single parameter is a workhorse of
%[text] computational psychiatry.

figure('Color', 'w'); hold on;
for a = [0.05 0.2 0.6]
    Va = zeros(n_trials + 1, 1);
    for t = 1:n_trials
        Va(t + 1) = Va(t) + a * (r(t) - Va(t));
    end
    plot(Va, 'LineWidth', 2, 'DisplayName', sprintf('\\alpha = %.2f', a));
end
xline(80, '--', 'Color', [.5 .5 .5], 'HandleVisibility', 'off');
yline(0.8, ':', 'Color', [.5 .5 .5], 'HandleVisibility', 'off');
title('Value trajectories for different learning rates');
xlabel('Trial'); ylabel('V'); legend('show', 'Location', 'best');

%% Part 1.3 — TD learning: the dopamine signal moves to the cue
%[text] Real dopamine neurons respond at specific moments within a trial, and with
%[text] learning the reward burst disappears and a burst appears at the earliest
%[text] predictive cue. TD learning captures this by assigning a value V(s_t) to
%[text] each time step in the trial:
%[text]
%[text] $$\delta_t = r_t + \gamma V(s_{t+1}) - V(s_t)$$
%[text]
%[text] where $r_t$ is the reward at time step $t$, $V(s_t)$ the value of the
%[text] state occupied at step $t$, and $\gamma$ the temporal discount factor
%[text] (0–1) that down-weights value one step in the future.
%[text]
%[text] We simulate a trial as a chain of time steps: a brief inter-trial interval
%[text] (ITI), a cue at step 4, and a reward at step 14. Because cue onset is
%[text] unpredictable, ITI states keep value 0 (they are not updated); the jump
%[text] from 0 to V(cue) is itself a prediction error, producing the persistent
%[text] cue response.

T = 18; cue = 4; rew = 14;             % states within a trial (1-based indices)
gamma = 0.98;                          % gamma = discount factor (0-1)
alpha_td = 0.10; n_train = 600;        % TD learning rate; number of training trials

r_trial = zeros(T, 1); r_trial(rew) = 1;

V = zeros(T + 1, 1);                   % V(T+1) = 0 (end of trial)
deltas = zeros(n_train, T);
for trial = 1:n_train
    for t = 1:T
        d = r_trial(t) + gamma * V(t + 1) - V(t);
        deltas(trial, t) = d;
        if t >= cue                    % ITI states stay 0: cue onset unpredictable
            V(t) = V(t) + alpha_td * d;
        end
    end
end

% Probe trial with reward OMITTED (measure the error, no learning)
delta_omit = zeros(T, 1);
for t = 1:T
    delta_omit(t) = 0 + gamma * V(t + 1) - V(t);
end

figure('Color', 'w');
snapshots = {deltas(1, :)', 'Trial 1'; deltas(20, :)', 'Trial 20'; ...
             deltas(end, :)', 'Trial 600'; delta_omit, 'Omission probe'};
for k = 1:4
    subplot(1, 4, k);
    stem(1:T, snapshots{k, 1}, 'filled'); hold on;
    xline(cue, '--g'); xline(rew, '--', 'Color', [1 .6 0]);
    ylim([-1.1 1.1]); title(snapshots{k, 2}); xlabel('Time step in trial');
    if k == 1, ylabel('TD error \delta_t (dopamine-like)'); end
end

%[text] This is the Schultz result. Early in training the prediction error occurs
%[text] at the reward (unexpected reward -> burst). As value propagates backward
%[text] through the state chain, the error shrinks at the reward and marches
%[text] earlier in time. After training, the reward evokes almost no error — the
%[text] burst has transferred to the cue, whose unpredictable onset still carries
%[text] surprise (it appears at the ITI-to-cue transition, one step before the
%[text] dashed cue line: the moment the cue appears). When a fully predicted
%[text] reward is omitted, the error at the
%[text] expected reward time is negative — a dip below baseline, as dopamine
%[text] neurons pause when an expected reward fails to arrive.
%[text]
%[text] In model-based fMRI, trial-by-trial \delta values from a model fitted to
%[text] behavior serve as parametric regressors (Chapter 20), and ventral striatal
%[text] BOLD reliably tracks them. Remember the limits: BOLD is not dopamine, and
%[text] the mapping is an analogy at the level of computation.

%% Part 2.1 — Toy RSA: stimuli and simulated brain regions
%[text] RSA compares two systems at the level of representational geometry. For a
%[text] set of stimuli, each system's representational dissimilarity matrix (RDM)
%[text] holds the pairwise dissimilarities (here 1 - correlation) between the
%[text] response patterns the stimuli evoke. Systems with different units (nodes
%[text] vs. voxels) are compared by correlating their RDMs.
%[text]
%[text] We build 32 stimuli (4 categories x 8 exemplars, 20 features each) with a
%[text] deliberately weak category signal buried under exemplar variation, and two
%[text] simulated regions measured with noise: an early sensory region mixing the
%[text] raw features, and a category-selective region coding category identity.

n_cat = 4; n_ex = 8;                   % categories; exemplars per category
n_stim = n_cat * n_ex; n_feat = 20;    % 32 stimuli total, 20 features each
categ = repelem((1:n_cat)', n_ex);              % category labels

protos = randn(n_cat, n_feat);                  % category prototypes
X = 0.6 * protos(categ, :) + 1.6 * randn(n_stim, n_feat);
% weak category signal, buried under exemplar-specific variation

early_region = X * randn(n_feat, 60) + 2.0 * randn(n_stim, 60);
onehot = eye(n_cat); onehot = onehot(categ, :);
categ_region = 4.0 * onehot * randn(n_cat, 60) + 2.0 * randn(n_stim, 60);

rdm = @(A) 1 - corr(A');                        % correlation-distance RDM

figure('Color', 'w');
subplot(1, 2, 1); imagesc(rdm(early_region), [0 2]); axis square;
title('Early sensory region RDM'); xlabel('Stimulus'); ylabel('Stimulus');
subplot(1, 2, 2); imagesc(rdm(categ_region), [0 2]); axis square;
title('Category-selective region RDM'); xlabel('Stimulus'); colorbar;

%[text] The category-selective region shows the classic block-diagonal RDM: items
%[text] from the same category evoke similar patterns. The early region's RDM
%[text] reflects idiosyncratic exemplar similarity, because exemplar variation
%[text] dominates the raw features.

%% Part 2.2 — Train a tiny network with hand-coded backpropagation
%[text] One hidden layer of 12 ReLU units, softmax output over 4 categories,
%[text] full-batch gradient descent on the cross-entropy loss plus a little weight
%[text] decay (which encourages the network to discard input directions it does
%[text] not need) — the same algorithm that trains networks with billions of
%[text] parameters.

n_hidden = 12; n_epochs = 2000;        % hidden units; training epochs
lr = 0.5; wd = 3e-3;                   % lr = learning rate; wd = weight decay
Y = eye(n_cat); Y = Y(categ, :);                % one-hot targets

W1 = 0.5 * randn(n_feat, n_hidden); b1 = zeros(1, n_hidden);
W2 = 0.5 * randn(n_hidden, n_cat);  b2 = zeros(1, n_cat);

losses = zeros(n_epochs, 1);
for epoch = 1:n_epochs
    Z1 = X * W1 + b1;                           % forward pass
    A1 = max(Z1, 0);                            % ReLU hidden layer
    S = A1 * W2 + b2;
    P = exp(S - max(S, [], 2)); P = P ./ sum(P, 2);   % softmax
    losses(epoch) = -mean(sum(Y .* log(P + 1e-12), 2));

    dS = (P - Y) / n_stim;                      % backpropagation
    dW2 = A1' * dS + wd * W2; db2 = sum(dS, 1);
    dA1 = dS * W2';
    dZ1 = dA1 .* (Z1 > 0);
    dW1 = X' * dZ1 + wd * W1; db1 = sum(dZ1, 1);

    W2 = W2 - lr * dW2; b2 = b2 - lr * db2;     % gradient descent step
    W1 = W1 - lr * dW1; b1 = b1 - lr * db1;

    if epoch == 1                               % save the UNTRAINED hidden layer
        hidden_untrained = A1;
    end
end
Z1 = X * W1 + b1; hidden_trained = max(Z1, 0);  % layers AFTER training
S = hidden_trained * W2 + b2;
P = exp(S - max(S, [], 2)); P = P ./ sum(P, 2);
output_trained = P;

[~, pred] = max(P, [], 2);
fprintf('final training accuracy: %.0f%%   final loss: %.3f\n', ...
    100 * mean(pred == categ), losses(end));

figure('Color', 'w'); plot(losses, 'LineWidth', 2);
title('Training loss'); xlabel('Epoch'); ylabel('Cross-entropy');

%% Part 2.3 — Compare model RDMs to brain-region RDMs
%[text] Correlate the off-diagonal RDM entries (Spearman) for each model
%[text] representation against each brain region.

layers = {X, hidden_untrained, hidden_trained, output_trained};
layer_names = {'Input features', 'Hidden (untrained)', 'Hidden (trained)', 'Output (trained)'};
regions = {early_region, categ_region};
region_names = {'EarlySensory', 'CategorySelective'};

iu = find(triu(ones(n_stim), 1));               % off-diagonal entries only
rsa = zeros(4, 2);
for i = 1:4
    Ri = rdm(layers{i});
    for j = 1:2
        Rj = rdm(regions{j});
        rsa(i, j) = corr(Ri(iu), Rj(iu), 'type', 'Spearman');
    end
end

disp(array2table(rsa, 'VariableNames', region_names, 'RowNames', layer_names));

figure('Color', 'w');
bar(rsa);
set(gca, 'XTickLabel', layer_names);
ylabel('RDM correlation (Spearman)');
title('Which brain region does each representation resemble?');
legend(region_names, 'Location', 'northwest');

%[text] The raw input and the untrained hidden layer resemble the early sensory
%[text] region: a random ReLU projection largely preserves input geometry. After
%[text] training, the deeper into the network you look, the more categorical the
%[text] geometry becomes: the trained hidden layer sits in between, and the
%[text] trained output layer develops clear block structure, resembling the
%[text] category-selective region far more than the early one. Nothing about the
%[text] brain was used in training; brain-like geometry emerged from optimizing
%[text] task performance, growing more categorical with depth. That is the logic
%[text] of Yamins et al. (Figure 42.3): in CNNs optimized for object recognition,
%[text] intermediate layers match V4 while top layers match IT cortex.
%[text]
%[text] Caveats: a high RDM correlation shows shared geometry, not shared
%[text] mechanism; our "regions" were simulated with the very codes we tested for,
%[text] so real analyses must guard against circularity, use held-out stimuli, and
%[text] compare many candidate models; and backpropagation as used here has no
%[text] accepted biological implementation. See the chapter page for discussion
%[text] and thought questions.
