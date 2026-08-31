%% Chapter 23 Lab - Localizing and Interpreting Results (MATLAB)
% In this lab you will use CANlab atlas objects to localize and interpret
% fMRI results: load a modern combined atlas, select regions of interest
% (ROIs) a priori, extract per-subject ROI averages and test them, label a
% thresholded statistic map automatically, and render results on 3-D
% surfaces. A final section demonstrates why "circular" (post hoc) ROI
% selection biases effect estimates.
%
% Companion to: https://torwager.github.io/elements-of-fmri-tutorials/book/part4/ch23-localizing-and-interpreting-results
%
% Requirements: CanlabCore <https://github.com/canlab/CanlabCore> and
% SPM12 on your MATLAB path. All data are bundled with CanlabCore
% (no large downloads).
%
% Adapted from CANlab tutorials:
% github.com/canlab/CANlab_help_examples (canlab_help_3b_atlases_and_ROI_analysis,
% using_canlab_atlases, canlab_help_4b_3D_visualization)

%% About atlas objects
% An atlas-class object is a specialized subclass of fmri_data that stores
% a set of labeled parcels (and, when available, the probabilistic maps
% underlying the parcellation). Common uses include:
%
% * Labeling activated regions by best-matching atlas parcels (region.table)
% * ROI analysis on parcel averages
% * Defining nodes for connectivity and graph analyses
%
% The default in many CANlab functions is the "CANlab combined 2018"
% atlas, which merges the Glasser 2016 HCP cortical parcellation (180
% areas/hemisphere, volumetric MNI projection), Pauli 2016 basal ganglia,
% amygdala/hippocampal regions from the SPM Anatomy Toolbox, the Morel
% thalamus atlas, the Pauli CIT168 "reinforcement learning" subcortical
% atlas, the Diedrichsen cerebellar (SUIT) atlas, and named brainstem
% nuclei. References are stored in the object and printed by region.table.

help atlas

% List the named atlases you can load (Glasser, Schaefer, Yeo, Julich, ...):
help load_atlas

%% Load and visualize the CANlab combined 2018 atlas

atlas_obj = load_atlas('canlab2018_2mm');   % 'canlab2018_2mm' = combined 2018 atlas at 2 mm resolution

% How many parcels?
disp(atlas_obj)
fprintf('Number of labeled parcels: %d\n', num_regions(atlas_obj));

% Interactive orthogonal slices: click around to see labels
orthviews(atlas_obj);

% Slice montage of all parcels
o2 = montage(atlas_obj);
drawnow, snapnow

%% Select regions of interest -- a priori, by name
% Choosing ROIs from an atlas BEFORE looking at your test data is the key
% to a valid ROI analysis. Here we select thalamic parcels by name.

% All parcels with 'Thal' in the label:
thal = select_atlas_subset(atlas_obj, {'Thal'});
disp(thal.labels')

% A focused a priori set: sensory, association, and epithalamic regions
thal = select_atlas_subset(atlas_obj, ...
    {'Thal_Intra', 'Thal_VL', 'Thal_MD', 'Thal_LGN', 'Thal_MGN', 'Thal_Hb'});

% You can also collapse a set of parcels into one region with 'flatten':
whole_thal = select_atlas_subset(atlas_obj, {'Thal'}, 'flatten');

% Or select parcels near a coordinate (here vmPFC, within 20 mm):
vmpfc_set = select_regions_near_crosshairs(atlas_obj, ...
    'coords', [0 38 -11], 'thresh', 20);   % coords = vmPFC in MNI mm; thresh = 20 mm search radius
disp(vmpfc_set.labels')

%% Load a sample dataset
% Bundled with CanlabCore: first-level contrast images from 30 subjects
% performing an emotion regulation task ([reappraise negative - look
% negative] contrast). One image per subject, in MNI space.

image_obj = load_image_set('emotionreg');
descriptives(image_obj);

%% Extract ROI averages and run a group ROI analysis
% extract_roi_averages returns a region object array r; r(i).dat contains
% the average over voxels in region i, for each of the 30 subject images.

r = extract_roi_averages(image_obj, thal);

% Concatenate region averages into a subjects x regions matrix:
roi_avgs = cat(2, r.dat);

% Clean up labels for display:
thal_labels = format_strings_for_legend(thal.labels);

% barplot_columns plots each region (violin + points + SE bar) and runs a
% one-sample t-test per column -- one test per a priori region, instead of
% one test per voxel. This is the ROI analysis.
create_figure('Thalamus ROI analysis');
barplot_columns(roi_avgs, 'nofig', ...
    'colors', scn_standard_colors(size(roi_avgs, 2)), 'names', thal_labels);
xlabel('Thalamic region'); ylabel('[Reappraise - Look] contrast value');
drawnow, snapnow

% With few a priori regions, correcting across 6 tests (e.g., Bonferroni:
% p < .05/6) is far less punishing than correcting across ~100,000 voxels.

%% Label a thresholded map automatically with region.table
% The complementary direction: you have a whole-brain result and want to
% know WHERE it is. region() breaks a thresholded map into contiguous
% blobs, and table() labels each blob using atlas parcels, with references.

t = ttest(image_obj);                    % voxelwise one-sample t-test
t = threshold(t, .05, 'fdr', 'k', 10);   % q = .05 FDR threshold; k = 10-voxel minimum cluster extent

orthviews(t);
drawnow, snapnow

rois = region(t);                        % contiguous significant blobs
table(rois);                             % auto-labeled table with atlas refs

%% Render results on 3-D surfaces
% Surface projection aids interpretation and communication. The surface()
% method renders blobs on canonical cortical (and cutaway) surfaces.

surface(t);
drawnow, snapnow

% A registry of montages + surfaces in one canonical layout:
o2 = canlab_results_fmridisplay([], 'compact2', 'noverbose');
o2 = addblobs(o2, rois);
drawnow, snapnow

% Note: this projects volumetric (MNI) results onto a surface for display.
% Fully surface-based analysis (fsaverage / fs_LR alignment) is more
% accurate in cortex but requires surface preprocessing from the start.

%% A priori vs. circular ROI selection
% Selecting the ROI from the SAME data you then test ("the parcel that lit
% up most") is circular: selection by the maximum guarantees inflated
% effect sizes, even in pure noise. Demonstration on null data:

rng(23);                                       % seed for reproducibility
n_sub = 20; n_parcels = 60; n_sims = 2000;     % n_sub = subjects; n_parcels = candidate ROIs; n_sims = null studies
apriori_est  = zeros(n_sims, 1);
circular_est = zeros(n_sims, 1);

for i = 1:n_sims
    nulldata = randn(n_sub, n_parcels);        % no true effect anywhere
    m = mean(nulldata);
    tvals = m ./ (std(nulldata) ./ sqrt(n_sub));

    apriori_est(i)  = m(1);                    % parcel chosen in advance
    [~, best] = max(tvals);                    % parcel chosen post hoc
    circular_est(i) = m(best);
end

fprintf('Mean estimated effect, a priori ROI:  %.3f (truth = 0)\n', mean(apriori_est));
fprintf('Mean estimated effect, circular ROI:  %.3f (truth = 0)\n', mean(circular_est));

create_figure('Circularity bias');
histogram(apriori_est, 40); hold on;
histogram(circular_est, 40);
legend({'A priori ROI', 'Circular ROI (max of same data)'});
xlabel('Estimated effect (truth = 0)'); ylabel('Simulations');
drawnow, snapnow

% The circular estimates are strongly biased upward even though NO parcel
% has any true effect. Fixes: a priori atlas ROIs, independent data or
% contrasts for selection, or split-half / cross-validated selection.

%% Explore further
% * Try load_atlas('glasser') or load_atlas('schaefer400') and re-run the
%   ROI analysis with a different a priori region set.
% * Use select_atlas_subset with 'flatten' to test one average over a
%   whole structure instead of subregions.
% * canlab.github.io and github.com/canlab/Neuroimaging_Pattern_Masks host
%   many atlases, meta-analytic maps, and signature patterns.
