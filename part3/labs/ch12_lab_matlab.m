%% Chapter 12 Lab: fMRI Basics and Terminology (MATLAB)
% This lab accompanies Chapter 12, "fMRI Basics and Terminology". You will
% load a set of brain images into a CANlab fmri_data object, inspect its
% dimensions and voxel-to-world (affine) mapping, convert between voxel
% (matrix) coordinates and world (mm) coordinates, extract values from a
% single voxel, and display the images with orthviews and montage.
%
% Requirements: CanlabCore and SPM12 on your MATLAB path.
%   https://github.com/canlab/CanlabCore
% Code adapted from CANlab tutorials (github.com/canlab and
% CANlab_help_examples: canlab_help_2b_basic_image_visualization.m and
% Voxel_and_image_spaces_in_CANlab_fmri_data_objects.mlx).
%
% Runtime: 1-2 minutes. Uses a small sample dataset bundled with CanlabCore.

%% 1. Load a sample dataset into an fmri_data object
% load_image_set loads a pre-canned set of 30 contrast images (one per
% participant) from an emotion regulation study. When images are loaded, a
% default brain mask is applied -- values outside the brain are not stored.
%
% The object stores the image data in the .dat field, a matrix of
% [voxels x images]. Each COLUMN is one 3-D brain volume, unrolled into a
% vector of in-mask voxels. In a 4-D time-series dataset, the columns would
% be volumes acquired one per TR.

imgs = load_image_set('emotionreg', 'noverbose');

size(imgs.dat)          % [in-mask voxels x images]

%% 2. Visualize the dataset
% The plot() method produces an orthviews display plus diagnostic plots
% (means, histograms, outliers). This is a good first look at any dataset.

plot(imgs);
drawnow, snapnow

%% 3. Image dimensions and the affine matrix: .volInfo
% The .volInfo field stores the image space information:
%   .dim  -- the 3-D image dimensions in voxels, in x (left-right),
%            y (back-front), and z (bottom-top)
%   .mat  -- the 4x4 affine matrix mapping voxel (matrix) coordinates to
%            world (mm) coordinates, in the same format used by SPM

imgs.volInfo.dim        % e.g., [47 56 31] voxels
imgs.volInfo.mat        % affine: voxel -> mm

% Read the affine like this:
% - The absolute values of the diagonal are the voxel sizes in mm.
% - The last column is related to the origin: the world (mm) position of
%   the image, so that some voxel maps to [0 0 0] -- the anterior commissure
%   in standard space.
% - A NEGATIVE first diagonal element means the stored x-axis runs
%   right-to-left ("radiological" storage); the object flips images
%   left-right when loading so that they display in neurological format
%   (right side of image = right side of brain).
% - Zeros off the diagonal mean no rotations were applied.

voxel_size = abs(diag(imgs.volInfo.mat(1:3, 1:3)))'

%% 4. Converting between voxel space and world (mm) space
% Voxel coordinates [i j k] index into the 3-D data matrix. World
% coordinates [x y z] are in mm relative to the origin. The affine converts
% between them (using homogeneous coordinates, with a trailing 1):

vox = [25 30 16 1]';                 % a voxel coordinate, [i j k 1]'
mm  = imgs.volInfo.mat * vox         % -> [x y z 1]' in mm

% And the inverse mapping, mm -> voxel:
vox_back = imgs.volInfo.mat \ mm     % recovers [25 30 16 1]'

% CanlabCore also provides voxel2mm() and mm2voxel() convenience functions
% that wrap this same calculation.
%
% The .volInfo.xyzlist field lists the [i j k] voxel coordinates for every
% in-mask voxel (one row per voxel), in the same order as the rows of .dat:

imgs.volInfo.xyzlist(1:5, :)         % first five in-mask voxels

%% 5. The mean image, orthviews, and montage
% orthviews() shows one sagittal, one coronal, and one axial slice through
% a chosen point -- click around to explore in 3-D. montage() shows a
% canonical series of slices, the standard "whole picture" view for papers.

m = mean(imgs);                      % mean across the 30 images

orthviews(m);
drawnow, snapnow

create_figure('montage'); axis off;
montage(m);
drawnow, snapnow

%% 6. Extract values from a single voxel
% Each row of .dat is one voxel's value in every image. Here the 30 values
% are contrast estimates from 30 participants; in a 4-D time-series object
% the same row would be a time series, one value per TR (see the Python lab
% for a simulated example).

wh_voxel = 1000;                     % pick an arbitrary in-mask voxel

% Where is this voxel in the brain? Convert its coordinate to mm:
xyz_vox = imgs.volInfo.xyzlist(wh_voxel, :);
xyz_mm  = imgs.volInfo.mat * [xyz_vox 1]';
fprintf('Voxel [%d %d %d] is at [%3.0f %3.0f %3.0f] mm\n', ...
    xyz_vox, xyz_mm(1:3));

v = imgs.dat(wh_voxel, :)';          % this voxel's value in each image

create_figure('one voxel');
plot(v, 'o-');
xlabel('Image number (participant)');
ylabel('Contrast value (a.u.)');
title(sprintf('Voxel [%d %d %d]: one value per image', xyz_vox));
drawnow, snapnow

%% 7. Orientation check: which side is left?
% The brain is nearly symmetric, so flipping errors are easy to miss.
% Safeguards include using the NIfTI format consistently, checking the sign
% conventions in the affine, fiducial markers, and anatomical heuristics:
% the calcarine fissure usually deviates rightward because the left
% occipital lobe is larger ("left looms larger").
%
% CANlab display methods use neurological orientation: the right side of
% the image is the right side of the brain. You can verify with a region
% known to be lateralized, or by displaying an anatomically labeled atlas:

atl = load_atlas('canlab2018_2mm');
r = select_atlas_subset(atl, {'Ctx_LO1_L', 'Ctx_LO1_R'});  % L and R visual regions
montage(r, 'regioncenters');
drawnow, snapnow

%% Explore on your own
% 1. Pick a world coordinate from a paper (e.g., [-42 -58 -12] mm) and
%    convert it to a voxel coordinate in this dataset with the inverse
%    affine. Does it land inside the brain mask? (Check against
%    imgs.volInfo.xyzlist.)
% 2. Use histogram(imgs.dat(:)) to view the distribution of all values.
%    Why are absolute values arbitrary, and what would change if the
%    scanner scaled all intensities by 10?
% 3. Try orthviews on the standard deviation image (std(imgs)) instead of
%    the mean. Which brain areas are most variable across participants,
%    and what might drive that (anatomy, vasculature, registration)?
