%% Chapter 14 Lab: BOLD Physiology and the Hemodynamic Response (MATLAB)
% This lab accompanies Chapter 14, "BOLD Physiology". You will build the
% canonical hemodynamic response function (HRF), simulate HRF variability,
% predict responses to brief events and sustained epochs by convolution
% under linear time-invariant (LTI) assumptions, and then break linearity:
% vascular saturation makes closely spaced events under-add and biases
% GLM amplitude estimates for densely presented conditions.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab). The saturation
% and bias demos are adapted from the CANlab nonlinear_saturation_bias_fmri
% simulation and hrf_saturation.m (CanlabCore/Model_building_tools).
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part3/ch14-bold-physiology
%
% Runtime: under a minute. All data are simulated.

%% 1. The canonical HRF
% SPM's canonical HRF is a difference of two gamma functions: one gamma
% creates the peak (~5 s after neural activity) and a second, later gamma
% is subtracted to create the post-stimulus undershoot. We sample it at
% high temporal resolution (0.1 s) and normalize so the peak = 1. From
% here on, "1 unit" is the peak response to a single brief event.
%
% Note: this canonical form omits the small "initial dip" (first ~0.5 s,
% roughly one fifth of the peak at 3 T), which requires fast sampling to
% detect.

dt = 0.1;                                % time resolution (s)
hrf = spm_hrf(dt);                       % canonical double-gamma HRF
hrf = hrf ./ max(hrf);                   % single-event peak = 1
t_hrf = (0:length(hrf) - 1)' .* dt;

[~, i_peak] = max(hrf);
[~, i_under] = min(hrf);
fprintf('Peak at %3.1f s; undershoot minimum at %3.1f s (%2.0f%% of peak)\n', ...
    t_hrf(i_peak), t_hrf(i_under), 100 * abs(min(hrf)));

create_figure('canonical hrf');
plot(t_hrf, hrf, 'LineWidth', 2); hold on;
plot(t_hrf, zeros(size(t_hrf)), 'Color', [.5 .5 .5]);
xlabel('Time (s)'); ylabel('Response (single-event peak = 1)');
title('Canonical double-gamma HRF');

%% 2. The HRF is not fixed: "young" vs. "elderly" responses
% The HRF varies across people, brain regions, and vascular health. With
% age (and with hypertension, diabetes, caffeine, etc.) it becomes lower
% in amplitude and more protracted. We mimic this by building the double
% gamma by hand with shifted parameters and reduced amplitude.

dgamma = @(t, a1, a2, c) gampdf(t, a1, 1) - c .* gampdf(t, a2, 1);

hrf_young = dgamma(t_hrf, 6, 16, 1/6);
hrf_young = hrf_young ./ max(hrf_young);           % canonical

hrf_old = dgamma(t_hrf, 8, 18, 1/6);
hrf_old = 0.7 .* hrf_old ./ max(hrf_old);          % delayed, blunted

create_figure('hrf variability');
plot(t_hrf, hrf_young, 'k-', 'LineWidth', 2); hold on;
plot(t_hrf, hrf_old, '--', 'Color', [.5 .5 .5], 'LineWidth', 2);
legend({'Typical young (canonical)' 'Typical elderly'});
xlabel('Time (s)'); ylabel('Response');
title('HRF amplitude and latency vary across people');

% If both groups are analyzed with the same canonical HRF, the mismatch
% for the elderly group reduces fit and estimated amplitudes -- a vascular
% difference that can masquerade as a neural one.

%% 3. The brain as an LTI system: brief events vs. sustained epochs
% Under LTI assumptions, the predicted BOLD response is the convolution of
% the neural stimulus function with the HRF: x(t) = (s * h)(t).
% onsets2fmridesign builds indicator functions, convolves them with the
% canonical HRF, and samples at the TR. A second column in each onset
% entry specifies event duration in seconds.

TR = 1;                                  % repetition time / sampling interval (s)
len = 60;                                % simulated run length (s)

ons = {};
ons{1} = [10 0.5];                        % brief event: onset 10 s, 0.5 s
ons{2} = [10 20];                         % epoch: onset 10 s, 20 s duration
X = onsets2fmridesign(ons, TR, len);      % intercept is the LAST column

create_figure('event vs epoch');
plot(0:TR:len - 1, X(:, 1), 'LineWidth', 2); hold on;
plot(0:TR:len - 1, X(:, 2), 'LineWidth', 2);
legend({'Brief event (0.5 s)' 'Epoch (20 s)'});
xlabel('Time (s)'); ylabel('Predicted BOLD');
title('Convolution predicts event and epoch responses');

% The brief event yields a transient HRF-shaped response; the epoch ramps
% up over ~10 s, plateaus while stimulation continues, then falls with an
% undershoot. This is why blocked designs yield large, sustained signals.

%% 4. Verify the LTI properties: scaling and superposition
% Scaling: doubling the input doubles the output. Superposition: the
% response to two events equals the sum of the single-event responses.
% We verify both with high-resolution convolution.

frame = (0:dt:len - dt)';
mkstim = @(onsets, dur) double(any( ...
    frame >= onsets(:)' & frame < onsets(:)' + dur, 2));
bold_linear = @(s) conv(s, hrf) .* dt;    % LTI prediction (then truncate)
truncate = @(x) x(1:length(frame));

s_a  = mkstim(10, 0.5);
s_b  = mkstim(22, 0.5);
s_ab = mkstim([10 22], 0.5);

scaling_ok = max(abs(truncate(bold_linear(2 .* s_a)) - ...
    2 .* truncate(bold_linear(s_a)))) < 1e-10;
superpos_ok = max(abs(truncate(bold_linear(s_ab)) - ...
    truncate(bold_linear(s_a)) - truncate(bold_linear(s_b)))) < 1e-10;

fprintf('Scaling holds: %d;  superposition holds: %d\n', scaling_ok, superpos_ok);

% Both hold exactly -- for the model. Real vasculature is another story.

%% 5. Breaking linearity: vascular saturation for closely spaced events
% Empirically, stimuli repeated within a few seconds evoke smaller
% responses than superposition predicts (neural and vascular refractory
% effects). Linearity is a good approximation for spacings of ~5 s or
% more, but nonlinearity is substantial below ~2 s.
%
% We model this with a compressive "squashing" function applied to the
% summed linear response: sat(x) = cap * tanh(x / cap). It is near-linear
% for small responses and flattens near the ceiling. CanlabCore implements
% the same idea as a piecewise-linear squash in hrf_saturation.m, applied
% via onsets2fmridesign(..., 'nonlinsaturation').

unit = max(truncate(bold_linear(mkstim(10, 0.5))));   % calibrate units
cap = 2;                                 % ceiling, in single-event peak units
sat = @(x) cap .* tanh(x ./ cap);
bold_sat = @(s) sat(truncate(bold_linear(s)) ./ unit);

s_pair = mkstim([20 21], 0.5);                        % two events 1 s apart
x_lin = truncate(bold_linear(s_pair)) ./ unit;

create_figure('saturation');
plot(frame, x_lin, '--', 'LineWidth', 2); hold on;
plot(frame, sat(x_lin), 'LineWidth', 2);
legend({'Linear (LTI) prediction' 'With vascular saturation'});
xlabel('Time (s)'); ylabel('Response (single-event units)');
title('Closely spaced events under-add');

fprintf('Peak, linear prediction: %3.2f; saturated: %3.2f\n', ...
    max(x_lin), max(sat(x_lin)));

%% 6. Refractory effects as a function of inter-stimulus interval
% How much response does a SECOND event add, as a function of ISI? We
% compute the marginal response (pair minus single) and express its
% integral as a fraction of the single-event response. Values near 1 mean
% linearity holds; values below 1 mean refractory loss.

frame = (0:dt:90 - dt)';                              % longer window
mkstim = @(onsets, dur) double(any( ...
    frame >= onsets(:)' & frame < onsets(:)' + dur, 2));
truncate = @(x) x(1:length(frame));
bold_sat = @(s) sat(truncate(bold_linear(s)) ./ unit);

isis = [1 2 3 4 5 6 8 10 12];         % inter-stimulus intervals to test (s)
ratios = zeros(size(isis));
r_one = bold_sat(mkstim(20, 0.5));

for i = 1:length(isis)
    r_two = bold_sat(mkstim([20 20+isis(i)], 0.5));
    marginal = r_two - r_one;                          % 2nd event's addition
    ratios(i) = sum(marginal) ./ sum(r_one);
end

create_figure('refractory curve');
plot(isis, ratios, 'o-', 'LineWidth', 2); hold on;
plot(isis, ones(size(isis)), '--', 'Color', [.5 .5 .5]);
xlabel('Inter-stimulus interval (s)');
ylabel('2nd-event / 1st-event response');
title('Refractory effects shrink with wider spacing');
set(gca, 'YLim', [0 1.1]);

disp(table(isis', ratios', 'VariableNames', {'ISI_s' 'RelativeResponse'}))

%% 7. Why this matters: saturation biases GLM estimates
% Two conditions evoke IDENTICAL per-event neural responses, but one is
% densely packed (trains of 5 events, 2 s apart) and one is sparse
% (isolated events). We generate data from the saturating system, then fit
% the standard LINEAR GLM. The dense regressor over-predicts during its
% trains, so OLS shrinks its beta: a systematic bias, not noise.

len2 = 280;                              % run length for the bias demo (s)
frame = (0:dt:len2 - dt)';
mkstim = @(onsets, dur) double(any( ...
    frame >= onsets(:)' & frame < onsets(:)' + dur, 2));
truncate = @(x) x(1:length(frame));
bold_sat = @(s) sat(truncate(bold_linear(s)) ./ unit);

sparse_ons = [20 80 140 200 250];                     % 5 isolated events
train_starts = [50 110 170 230];
dense_ons = [];
for t0 = train_starts
    dense_ons = [dense_ons, t0 + 2 .* (0:4)];         % 4 trains of 5 events
end

s_sparse = mkstim(sparse_ons, 0.5);
s_dense  = mkstim(dense_ons, 0.5);

% Linear regressors, as a standard GLM would build them
x_sparse = truncate(bold_linear(s_sparse)) ./ unit;
x_dense  = truncate(bold_linear(s_dense))  ./ unit;

% 'True' data: same per-event amplitude (1 unit), saturating system
y_true = bold_sat(s_sparse + s_dense);

% Sample at the TR and add measurement noise
rng(14);                                 % reproducible noise
idx = 1:round(TR / dt):length(frame);
Xd = [x_sparse(idx) x_dense(idx) ones(length(idx), 1)];
y = y_true(idx) + 0.1 .* randn(length(idx), 1);  % noise SD = 0.1 units

b = Xd \ y;                                           % OLS fit

fprintf('True per-event amplitude for BOTH conditions: 1.00\n');
fprintf('beta_sparse = %3.2f;  beta_dense = %3.2f\n', b(1), b(2));
fprintf('Spurious condition difference (sparse - dense): %3.2f\n', b(1) - b(2));

create_figure('bias demo');
plot(frame(idx), y, 'k.', 'MarkerSize', 6); hold on;
plot(frame, x_sparse + x_dense, 'r--', 'LineWidth', 1.5);
legend({'Observed (saturating system)' 'Linear prediction, true amplitude'});
xlabel('Time (s)'); ylabel('Signal (single-event units)');
title('Dense trains: linear model over-predicts, so OLS shrinks the beta');

% Remedies -- at the design stage: keep event density comparable across
% conditions you plan to compare, or keep ISIs >= ~5 s where linearity
% holds. At the modeling stage: compare the sparse condition to an equally
% sparse SUBSET of the dense events, or model the nonlinearity explicitly
% (e.g., onsets2fmridesign(..., 'nonlinsaturation')).

%% Explore on your own
% 1. Change the ceiling (cap = 4, or cap = 1.2) and re-run Sections 5-7.
%    How does the spurious condition difference change?
% 2. Re-run Section 7 with 6 s spacing within the dense trains. How much
%    of the bias disappears?
% 3. Build a regressor from only the FIRST event of each dense train
%    (matching the sparse condition's density) and compare betas.
% 4. Generate data with hrf_old (Section 2) but fit with regressors built
%    from the canonical HRF. How much amplitude is lost to shape mismatch,
%    and could it be mistaken for reduced neural activity in aging?
