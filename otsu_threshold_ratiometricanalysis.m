


close all, clear, clc
%% 1. Load Images
% Replace with your actual filenames
donor = double(imread('mtCa-KD68-6_c1_ORG.tif'));
fret = double(imread('mtCa-KD68-6_c2_ORG.tif'));

%% 2. Background Subtraction (Top-Hat Filtering)
% Top-hat is superior for mitochondria "dots" as it removes uneven 
% background while preserving small, bright structures.
se = strel('disk', 10); % '10' is the approximate radius of your largest mito
donor_sub = imtophat(donor, se);
fret_sub = imtophat(fret, se);

%% 3. Semi-Auto Thresholding (Masking)
% We create a binary mask from the Donor channel to identify mitochondria
thresh = graythresh(mat2gray(donor_sub)); % Otsu's method
mask = imbinarize(mat2gray(donor_sub), thresh);

% Clean up the mask: remove noise (objects smaller than 10 pixels)
mask = bwareaopen(mask, 10); 

%% 4. Pixel-by-Pixel Ratio Calculation
% Calculate ratio: FRET / Donor
% We use the corrected (subtracted) images
ratio_img = fret_sub ./ donor_sub;

%% 5. Precision Masking for Visualization
% By setting background pixels to NaN, MATLAB will not plot them in the heatmap,
% preventing background noise from distorting your color scale.
ratio_display = ratio_img;
ratio_display(~mask) = NaN; 

%% 6. ROI-Based Statistical Measurement
% Identify individual mitochondria as unique objects
CC = bwconncomp(mask); 
stats_donor = regionprops(CC, donor_sub, 'MeanIntensity');
stats_fret = regionprops(CC, fret_sub, 'MeanIntensity');

% Calculate the ratio per individual mitochondrion
mito_ratios = [stats_fret.MeanIntensity] ./ [stats_donor.MeanIntensity];
mn_ratios = mean(mito_ratios);
%% 7. Visualization (The Heatmap)
figure('Name', 'Ratiometric Analysis');

% Display the ratio image
h = imagesc(ratio_display); 
set(gcf, 'Position',  [0, 0, 1000, 1000])

% Apply styling
colormap("turbo"); % Or 'turbo' for better perceptual uniformity

% --- COLORBAR MODIFICATION ---
c = colorbar;           % Create colorbar handle
c.FontSize = 32;        % Set font size
c.Label.String = 'FRET Ratio'; % Optional: Add a label to the colorbar

title('Mitochondrial Calcium Heatmap', ...
    'FontSize',32,'FontName','ArialNarrow','FontWeight','normal');
axis image off; % Maintain aspect ratio and remove axis ticks

% CRITICAL: Set the color axis scale (e.g., 0.5 to 2.5) 
% This allows you to compare different cells/treatments on the same scale.
caxis([0.5, 2.5]); 

% Output summary to command window
fprintf('Detected %d mitochondria.\n', CC.NumObjects);
fprintf('Mean Ratio: %.3f\n', mean(mito_ratios));

%% save data

name= 'mtCa-KD68-6_ORG';
datFileName = ['FRET-DATA' name ];
save(datFileName, '-v7.3');