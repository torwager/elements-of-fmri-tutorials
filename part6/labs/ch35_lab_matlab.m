%% Chapter 35 Lab: Dynamic Causal Models (MATLAB)
% This lab accompanies Chapter 35, "Dynamic Causal Models". You will build
% the core logic of a DCM by hand: simulate a two-region bilinear neuronal
% system in which a modulatory input strengthens the z1 -> z2 connection,
% pass the latent neuronal states through a hemodynamic observation model,
% and compare a "modulation" model against a "no-modulation" model using
% BIC as a simple stand-in for log model evidence.
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part6/ch35-dynamic-causal-models
%
% Requirements: base MATLAB plus the Statistics and Machine Learning
% Toolbox (for gampdf). No SPM or DCM toolbox is needed — everything is
% simulated. Real DCM analyses use SPM (spm_dcm_*), which performs full
% Bayesian inversion of the neuronal + Balloon hemodynamic model.
%
% Runtime: under a minute. All data are simulated.

%% 1. A two-region bilinear neuronal model
% The bilinear DCM neuronal equation is
%
%   dz/dt = (A + sum_j u_j * B{j}) * z + C * u
%
% A  = intrinsic coupling among regions (here z1 -> z2 with weight 0.3,
%      plus self-decay terms on the diagonal)
% B2 = change in coupling induced by input u2 (u2 gates z1 -> z2)
% C  = direct driving influence of input u1 (drives region 1 only)
%
% Input u1 delivers brief driving pulses (like stimuli); input u2 is a
% sustained block "context" (like attention) that modulates the connection.

dt = 0.1;                 % integration step (s)
T  = 300;                 % total simulation time (s)
t  = (0:dt:T-dt)';        % time grid
n  = numel(t);

u1 = double(mod(t, 20) < 1);             % 1-s driving pulse every 20 s
u2 = double(mod(floor(t/60), 2) == 1);   % alternating 60-s modulatory blocks

A  = [-0.4  0; 0.3 -0.4];  % intrinsic coupling (z1 -> z2 = 0.3)
B2 = [0 0; 0.5 0];         % u2 adds 0.5 to the z1 -> z2 path
C  = [1; 0];               % u1 drives region 1 only

%% 2. Simulate the latent neuronal dynamics (Euler integration)
% We step the differential equation forward in time. During modulation-ON
% blocks the effective z1 -> z2 coupling is 0.3 + 0.5 = 0.8, so region 2
% responds much more strongly to the same driving pulses.

z = zeros(n, 2);
for i = 1:n-1
    dz = (A + u2(i) * B2) * z(i, :)' + C * u1(i);
    z(i+1, :) = z(i, :) + dt * dz';
end

figure;
subplot(3, 1, 1); plot(t, u1, 'k'); hold on; plot(t, u2, 'r', 'LineWidth', 2);
legend('u_1 (driving)', 'u_2 (modulatory)'); title('Experimental inputs');
subplot(3, 1, 2); plot(t, z(:, 1)); title('Latent neuronal state z_1');
subplot(3, 1, 3); plot(t, z(:, 2)); title('Latent neuronal state z_2');
xlabel('Time (s)');

%% 3. Hemodynamic observation model: from z to BOLD
% Real DCM uses an extended Balloon model (vasodilatory signal, inflow,
% volume, deoxyhemoglobin) with region-specific parameters. As a simple
% stand-in we convolve the neuronal states with a canonical double-gamma
% HRF, sample every TR = 2 s, and add measurement noise.

hrf_t = (0:dt:30)';                                  % 30-s HRF support
hrf   = gampdf(hrf_t, 6, 1) - gampdf(hrf_t, 16, 1) / 6;   % canonical double-gamma HRF (peak - undershoot/6)

step  = round(2 / dt);                               % fine samples per TR (TR = 2 s)
TR_t  = t(1:step:end);                               % acquisition times

y1_clean = conv(z(:, 1), hrf) * dt;  y1_clean = y1_clean(1:step:n);
y2_clean = conv(z(:, 2), hrf) * dt;  y2_clean = y2_clean(1:step:n);

rng(35);                                             % seed for reproducible noise
noise_sd = 0.05;                                     % measurement noise SD (a.u.)
y1 = y1_clean + noise_sd * randn(size(y1_clean));
y2 = y2_clean + noise_sd * randn(size(y2_clean));

figure;
plot(TR_t, y1, 'o-', TR_t, y2, 'o-');
legend('y_1 (observed BOLD)', 'y_2 (observed BOLD)');
xlabel('Time (s)'); title('Observed BOLD (TR = 2 s), modulation ON at 60-120 s and 180-240 s');

%% 4. Invert two candidate models
% We now play the analyst: given the noisy y2, which generative model
% explains it best?
%   Model 1 ("modulation"):     a21 and b21 free       (k = 2 parameters)
%   Model 0 ("no modulation"):  a21 free, b21 fixed 0  (k = 1 parameter)
% Because z1 does not depend on a21 or b21 (no z2 -> z1 connection), we can
% hold the simulated z1 fixed and refit only region 2's dynamics — a cheap
% least-squares stand-in for DCM's Bayesian (EM) inversion.

z1_true = z(:, 1);   % simulated z1, held fixed (no z2 -> z1 connection)
a22 = -0.4;          % region-2 self-decay, treated as known here

obj_mod   = @(th) dcm_rss(th(1), th(2), z1_true, u2, y2, hrf, dt, step, a22, n);
obj_nomod = @(th) dcm_rss(th(1), 0,     z1_true, u2, y2, hrf, dt, step, a22, n);

theta1 = fminsearch(obj_mod,   [0.2 0.2]);   rss1 = obj_mod(theta1);    k1 = 2;   % [a21 b21] start; k = free params
theta0 = fminsearch(obj_nomod, 0.2);         rss0 = obj_nomod(theta0);  k0 = 1;   % a21 start

fprintf('Modulation model:    a21 = %.3f, b21 = %.3f (truth: 0.30, 0.50), RSS = %.3f\n', ...
    theta1(1), theta1(2), rss1);
fprintf('No-modulation model: a21 = %.3f                                 RSS = %.3f\n', ...
    theta0(1), rss0);

%% 5. Compare models with BIC (a simple evidence proxy)
% BIC = n*log(RSS/n) + k*log(n) penalizes fit by complexity. The difference
% in BIC approximates twice the log Bayes factor; a large positive
% dBIC (no-mod minus mod) is strong evidence for the modulation model.
% (Real DCM uses a free-energy approximation to the log evidence instead.)

n_obs = numel(y2);
bic1  = n_obs * log(rss1 / n_obs) + k1 * log(n_obs);
bic0  = n_obs * log(rss0 / n_obs) + k0 * log(n_obs);

fprintf('BIC (modulation)    = %.1f\n', bic1);
fprintf('BIC (no modulation) = %.1f\n', bic0);
fprintf('dBIC = %.1f -> approximate Bayes factor exp(dBIC/2) = %.2g in favor of modulation\n', ...
    bic0 - bic1, exp((bic0 - bic1) / 2));

% Visualize why: the no-modulation model cannot produce different response
% amplitudes in ON vs OFF blocks.
pred1 = dcm_predict(theta1(1), theta1(2), z1_true, u2, hrf, dt, step, a22, n);
pred0 = dcm_predict(theta0(1), 0,         z1_true, u2, hrf, dt, step, a22, n);

figure;
plot(TR_t, y2, 'k.', TR_t, pred1, '-', TR_t, pred0, '-');
legend('observed y_2', 'modulation model fit', 'no-modulation model fit');
xlabel('Time (s)'); title('Model fits to region-2 BOLD');

%% 6. Control analysis: data generated WITHOUT modulation
% Model comparison should also protect us from over-claiming. Regenerate
% the data with b21 = 0 (no true modulation) and repeat: now the extra
% parameter buys almost no fit, and the BIC penalty makes the simpler
% model win.

z2_null = simulate_z2(0.3, 0, z1_true, u2, dt, a22, n);
y2_null_clean = conv(z2_null, hrf) * dt;
y2_null = y2_null_clean(1:step:n) + noise_sd * randn(n_obs, 1);

obj_mod_n   = @(th) dcm_rss(th(1), th(2), z1_true, u2, y2_null, hrf, dt, step, a22, n);
obj_nomod_n = @(th) dcm_rss(th(1), 0,     z1_true, u2, y2_null, hrf, dt, step, a22, n);

th1n = fminsearch(obj_mod_n,   [0.2 0.2]);  rss1n = obj_mod_n(th1n);   % [a21 b21] start
th0n = fminsearch(obj_nomod_n, 0.2);        rss0n = obj_nomod_n(th0n); % a21 start

bic1n = n_obs * log(rss1n / n_obs) + 2 * log(n_obs);
bic0n = n_obs * log(rss0n / n_obs) + 1 * log(n_obs);
fprintf('Null-truth data: dBIC (no-mod - mod) = %.1f (negative -> simpler model wins)\n', ...
    bic0n - bic1n);

%% 7. Where to go from here
% What real DCM adds beyond this toy version:
%  - A biophysical Balloon hemodynamic model per region, with parameters
%    estimated under shrinkage priors (not a fixed canonical HRF)
%  - Full Bayesian inversion (EM), yielding posterior distributions over
%    all coupling parameters, not just point estimates
%  - Model evidence via a free-energy approximation, and group inference
%    via Parametric Empirical Bayes (PEB)
% Caveats to carry with you: inferences are only as good as the specified
% model; omitting regions that influence the modeled ones can produce
% false conclusions about direction and strength of connections; and
% classical DCM is limited to small networks. Use DCM to identify systems
% and pathways rather than to make strong causal claims.
%
% In SPM, see spm_dcm_specify / spm_dcm_estimate / spm_dcm_peb to run the
% real thing on your own data.

%% Local functions

function z2 = simulate_z2(a21, b21, z1, u2, dt, a22, n)
% Simulate region 2's neuronal state given candidate coupling parameters.
% z1 is held fixed because no z2 -> z1 connection exists in this model.
z2 = zeros(n, 1);
drive = (a21 + b21 * u2) .* z1;
for i = 1:n-1
    z2(i+1) = z2(i) + dt * (drive(i) + a22 * z2(i));
end
end

function pred = dcm_predict(a21, b21, z1, u2, hrf, dt, step, a22, n)
% Predicted BOLD for region 2: simulate neuronal state, convolve with HRF,
% and sample at the TR.
z2 = simulate_z2(a21, b21, z1, u2, dt, a22, n);
pred = conv(z2, hrf) * dt;
pred = pred(1:step:n);
end

function r = dcm_rss(a21, b21, z1, u2, y2obs, hrf, dt, step, a22, n)
% Residual sum of squares of the candidate model's predicted BOLD.
pred = dcm_predict(a21, b21, z1, u2, hrf, dt, step, a22, n);
r = sum((y2obs - pred).^2);
end
