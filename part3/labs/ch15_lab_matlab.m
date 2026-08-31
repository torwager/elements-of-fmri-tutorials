%% Lab 15: Sampling, Aliasing, and Smoothing
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part3/ch15-spatial-and-temporal-resolution
%[text] This lab builds hands-on intuition for the temporal and spatial resolution
%[text] limits discussed in Chapter 15 of *Elements of fMRI Analysis*. You will
%[text] construct sine waves and examine them with the fast Fourier transform (FFT),
%[text] watch a fast signal get *aliased* to a slow one when sampled below the
%[text] Nyquist limit, alias a simulated heartbeat by sampling it at typical fMRI
%[text] TRs, and measure what spatial smoothing does to small versus large activations.
%[text]
%[text] **Requirements:** base MATLAB. CANlab helper functions are used where noted
%[text] (CanlabCore on your path: github.com/canlab/CanlabCore), with base-MATLAB
%[text] equivalents shown so the script runs either way.
%[text]
%[text] Adapted from the CANlab tutorial "Signal basics, fft, and aliasing"
%[text] (canlab.github.io) and github.com/canlab/FMRI_simulations.

%% 1. Sine waves and the FFT
%[text] A sine wave is defined by its frequency F (Hz), amplitude, and phase theta:
%[text] $y(t) = \sin(2 \pi F t + \theta)$. We simulate one at a high sampling
%[text] frequency Fs — think of this as the "ground truth" signal before the scanner
%[text] samples it. The FFT re-expresses the signal as power at each frequency; the
%[text] highest representable frequency is the **Nyquist limit**, Fs/2.

Fs = 1000;            % Sampling frequency (samples/sec)
L = 1000;             % Length of signal (samples) -> 1 sec
t = (0:L-1) / Fs;     % Time vector (sec)

F = 10;               % Frequency of sine wave (Hz)
theta = 0;            % Phase (radians)
y = sin(2 * pi * F * t + theta);

yf = abs(fft(y)) / L;         % magnitude spectrum
f = Fs * (0:(L-1)) / L;       % frequency of each FFT coefficient

figure;
subplot(1, 2, 1); plot(t, y)
xlabel('Time (sec)'), title(sprintf('%d Hz sine wave (time domain)', F))
subplot(1, 2, 2); plot(f, yf)
xlim([0 50]), xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('Magnitude spectrum')

nyquist = Fs / 2;     % Any signal faster than this will be aliased
fprintf('Nyquist limit = Fs/2 = %d Hz\n', nyquist);

%[text] **Try it:** change F to 20, or set theta = pi/2 (turning sine into cosine).
%[text] The spectral peak moves with F; phase changes do not move the magnitude peak.
%[text] With CanlabCore on your path, try create_figure('sin wave') and
%[text] plot_vertical_line(nyquist) for CANlab-style figures.

%% 2. Aliasing: sampling a fast signal too slowly
%[text] The Nyquist theorem says we need more than 2 samples per cycle to represent
%[text] a periodic signal. Sampled more slowly, the signal is *reflected around the
%[text] Nyquist frequency* and appears as a spurious slow oscillation — aliasing.
%[text] Here a 10 Hz wave is re-sampled at 40, 15, and 12 Hz. Predicted apparent
%[text] frequency: |F - round(F/Fs_new)*Fs_new|.

Fs = 1000;              % Fs = dense "ground truth" sampling rate (Hz)
dur = 2;                % dur = signal duration (s)
t = 0:1/Fs:dur - 1/Fs;  % time vector (s)
F = 10;                 % F = true signal frequency (Hz)
y = sin(2 * pi * F * t);

new_rates = [40 15 12]; % slower sampling rates to test (Hz)
figure;
for i = 1:3
    fs_new = new_rates(i);
    step = round(Fs / fs_new);
    alias_f = abs(F - round(F / fs_new) * fs_new);

    subplot(3, 1, i);
    plot(t, y, 'Color', [.75 .75 .75]); hold on
    plot(t(1:step:end), y(1:step:end), 'k.-', 'LineWidth', 2);
    title(sprintf('Fs = %d Hz, Nyquist = %3.1f Hz -> apparent frequency ~ %3.1f Hz', ...
        fs_new, fs_new / 2, alias_f));
end
xlabel('Time (sec)')

%[text] At 40 Hz the samples trace the 10 Hz wave faithfully. At 15 Hz (Nyquist
%[text] 7.5 Hz) the signal folds to 5 Hz; at 12 Hz (Nyquist 6 Hz) it appears at
%[text] 2 Hz — a slow, entirely artifactual oscillation. Confirm in the frequency
%[text] domain:

fs_new = 12;                 % Nyquist = 6 Hz, below the 10 Hz signal
step = round(Fs / fs_new);   % keep every step-th sample of the dense signal
ysub = y(1:step:end);  n = numel(ysub);
fsub = (0:floor(n/2)) / (n / fs_new);
msub = abs(fft(ysub)) / n;  msub = msub(1:floor(n/2) + 1);

figure; plot(fsub, msub, 'k'); hold on
plot([fs_new/2 fs_new/2], ylim, 'k:')   % Nyquist limit (or: plot_vertical_line)
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('Spectrum of the 10 Hz signal sampled at 12 Hz')

[~, imax] = max(msub);
fprintf('Spectral peak at %3.2f Hz (true frequency: %d Hz)\n', fsub(imax), F);

%% 3. Aliasing a heartbeat at fMRI TRs
%[text] The cardiac cycle produces a periodic artifact at roughly 1 Hz (60 bpm),
%[text] far faster than task-related BOLD fluctuations (below ~0.1 Hz). Whether the
%[text] heartbeat stays separable from the task depends on the TR. We simulate a
%[text] spiky heartbeat with beat-to-beat variability, then sample it once per TR.

rng(2026);                          % seed for reproducible simulations
fs = 100;                           % fs = simulation rate (Hz) — plenty for a ~1 Hz signal
dur = 300;                          % dur = duration (s): 5 minutes, a short fMRI run
t = 0:1/fs:dur - 1/fs;

intervals = 1 + 0.06 * randn(round(dur * 1.2), 1);   % mean 1 s between beats (60 bpm), sd = 0.06 s jitter
intervals = min(max(intervals, 0.7), 1.4);
beat_times = cumsum(intervals);
beat_times = beat_times(beat_times < dur - 1);

heart = zeros(size(t));
for i = 1:numel(beat_times)
    heart = heart + exp(-(t - beat_times(i)).^2 / (2 * 0.05^2));  % brief pulse
end

figure; plot(t, heart, 'k'); xlim([0 20])
xlabel('Time (sec)'), title('Simulated heartbeat (~60 bpm, variable intervals)')

%[text] Sample at three TRs and inspect time and frequency domains. The task band
%[text] (~0.005–0.1 Hz) is where block and event-related designs live.

TRs = [0.5 1.0 2.0];    % TR = time between volumes (s)
figure;
for i = 1:3
    TR = TRs(i);
    step = round(TR * fs);
    samp = heart(1:step:end);
    ts = (0:numel(samp) - 1) * TR;
    nyq = 1 / (2 * TR);

    subplot(3, 2, 2*i - 1)
    plot(t, heart, 'Color', [.8 .8 .8]); hold on
    plot(ts, samp, '.-', 'Color', [.5 0 .5], 'LineWidth', 1.5)
    xlim([0 20]), xlabel('Time (sec)'), title(sprintf('TR = %3.1f s', TR))

    n = numel(samp);
    fsamp = (0:floor(n/2)) / (n * TR);
    msamp = abs(fft(samp - mean(samp))) / n;  msamp = msamp(1:floor(n/2) + 1);
    subplot(3, 2, 2*i)
    plot(fsamp, msamp, 'Color', [.5 0 .5]); hold on
    plot([nyq nyq], ylim, 'k:')
    xlim([0 1.1]), xlabel('Frequency (Hz)')
    title(sprintf('Nyquist = %3.2f Hz', nyq))
end

%[text] Quantify how much cardiac power lands in the task band at each TR:

fprintf('\nFraction of sampled cardiac power in the task band (0.005-0.1 Hz):\n');
for TR = [0.4 0.5 1.0 2.0 3.0]      % TRs from fast multiband to slow standard sampling
    samp = heart(1:round(TR * fs):end);
    n = numel(samp);
    fsamp = (0:floor(n/2)) / (n * TR);
    p = abs(fft(samp - mean(samp))).^2;  p = p(1:floor(n/2) + 1);
    in_band = sum(p(fsamp > 0.005 & fsamp < 0.1)) / sum(p(2:end));
    fprintf('  TR = %3.1f s (Nyquist %4.2f Hz): %5.1f %%\n', TR, 1/(2*TR), 100 * in_band);
end

%[text] Short TRs keep cardiac power near 1 Hz, out of the task band, where it can
%[text] be filtered or modeled. At TR = 2–3 s much of it aliases into task
%[text] frequencies and becomes inseparable from activation — one reason modern
%[text] multiband protocols push TR below ~0.5 s.

%% 4. Is smoothing lossy? Small vs. large activations
%[text] Spatial smoothing suppresses noise and absorbs some inter-subject
%[text] misalignment, but it acts as a matched filter: it attenuates activations
%[text] smaller than the kernel. We simulate a 1-D strip of cortex (1 voxel = 1 mm)
%[text] with two equal-amplitude activations — narrow (~3.5 mm FWHM) and broad
%[text] (~19 mm FWHM) — plus noise, then smooth with an 8 mm FWHM Gaussian.
%[text] Inspired by is_smoothing_lossy.m in github.com/canlab/FMRI_simulations.

rng(1);                                        % seed for reproducible noise
n_vox = 200;                                   % n_vox = number of voxels in the 1-D strip
x = (1:n_vox)';                                % position along the strip: 1 voxel = 1 mm
fwhm2sig = 1 / (2 * sqrt(2 * log(2)));         % FWHM -> Gaussian sigma

small = 2 * exp(-(x - 60).^2 / (2 * 1.5^2));   % ~3.5 mm FWHM activation
large = 2 * exp(-(x - 140).^2 / (2 * 8^2));    % ~19 mm FWHM activation
truth = small + large;
y = truth + randn(n_vox, 1);                   % noise sd = 1

kernel_fwhm = 8;                   % kernel_fwhm = smoothing kernel FWHM (mm), a typical choice
sig = kernel_fwhm * fwhm2sig;      % convert FWHM -> Gaussian sigma
kx = (-20:20)';                    % kernel support: +/- 20 mm
kern = exp(-kx.^2 / (2 * sig^2));  kern = kern / sum(kern);
ysmooth = conv(y, kern, 'same');

figure;
plot(x, y, 'Color', [.75 .75 .75]); hold on
plot(x, ysmooth, 'k', 'LineWidth', 2)
plot(x, truth, 'r--', 'LineWidth', 1.5)
legend({'Noisy data', 'Smoothed (8 mm FWHM)', 'True signal'})
xlabel('Position (mm)'), ylabel('Signal')
title('Equal-amplitude activations: narrow (left) vs. broad (right)')

fprintf('\nTrue peak amplitude, both blobs:      2.00\n');
fprintf('Smoothed peak, narrow blob (~3.5 mm): %3.2f\n', max(ysmooth(50:70)));
fprintf('Smoothed peak, broad blob (~19 mm):   %3.2f\n', max(ysmooth(130:150)));

%[text] The broad activation survives nearly intact; the narrow one — same true
%[text] amplitude — is flattened and could fall below a statistical threshold.
%[text] Sweep the activation width to see attenuation vs. detection sensitivity:

blob_fwhms = [2 3.5 5 8 12 19 30];   % activation FWHMs to sweep (mm)
n_sims = 200;                        % n_sims = noise realizations per width; more -> stabler estimates
center = n_vox / 2;                  % place each activation mid-strip

fprintf('\nBlob FWHM (mm) | smoothed peak amp | peak z (raw) | peak z (smoothed)\n');
for bf = blob_fwhms
    blob = 2 * exp(-(x - center).^2 / (2 * (bf * fwhm2sig)^2));
    amp_sm = conv(blob, kern, 'same');  amp_sm = amp_sm(center);

    z_raw = zeros(n_sims, 1);  z_sm = zeros(n_sims, 1);
    for s = 1:n_sims
        noise = randn(n_vox, 1);
        z_raw(s) = blob(center) + noise(center);
        smoothed = conv(blob + noise, kern, 'same');
        noise_sm = conv(noise, kern, 'same');
        z_sm(s) = smoothed(center) / std(noise_sm);
    end
    fprintf('%14.1f | %17.2f | %12.1f | %17.1f\n', ...
        bf, amp_sm, mean(z_raw), mean(z_sm));
end

%[text] Smoothing always attenuates peak amplitude for activations smaller than the
%[text] kernel — an analysis-imposed partial volume effect. But because it also
%[text] suppresses independent voxel noise, detection sensitivity (peak z) improves
%[text] for activations at or above the kernel scale. The kernel choice silently
%[text] decides which spatial scales your study can see — the "effective
%[text] resolution" argument of Chapter 15.

%% 5. Further explorations
%[text] * Change the heart rate to 72 bpm (1.2 Hz) and find the TRs at which it
%[text]   aliases closest to 0 Hz. Which common TRs are worst?
%[text] * Add a simulated respiratory signal (~0.3 Hz) and see which TRs alias it.
%[text] * With CanlabCore: use fmri_data and preprocess(dat, 'smooth', 8) to smooth
%[text]   real image volumes, and compare small subcortical vs. large cortical
%[text]   activation clusters.
