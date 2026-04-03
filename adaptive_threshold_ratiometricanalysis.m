% 1. Load and Clean
ch1 = double(imread('mtCa-KD68-4_c2_ORG.tif'));
ch2 = double(imread('mtCa-KD68-4_c3_ORG.tif'));

% 2. Background Subtraction (Morphological Top-hat is great for dots)
se = strel('disk', 10);
ch1_sub = imtophat(ch1, se);
ch2_sub = imtophat(ch2, se);

% 3. Semi-Auto Thresholding
mask = imbinarize(mat2gray(ch1_sub), 'adaptive'); % Adjusts to local intensity
mask = bwareaopen(mask, 10); % Remove noise smaller than 10 pixels

% 4. Connected Components (ROI identification)
CC = bwconncomp(mask); 

% 5. Measure Intensities from BOTH images using the SAME mask
stats1 = regionprops(CC, ch1_sub, 'MeanIntensity');
stats2 = regionprops(CC, ch2_sub, 'MeanIntensity');

% 6. Compute Ratios
ratios = [stats2.MeanIntensity] ./ [stats1.MeanIntensity];
mn_ratios = mean(ratios)
% 7. Visualize Overlay
imshow(label2rgb(labelmatrix(CC), 'jet', 'k'));
title('Automatically Detected Mitochondria ROIs');

%%
figure(2),clf
imagesc((ch2_sub ./ ch1_sub) .* mask); 
colormap("jet"); 
colorbar;
axis image; % Keeps mitochondria from looking stretched
title('Mitochondrial FRET Ratio (Acceptor/Donor)');

% Set the color axis (CRITICAL for comparing different cells)
caxis([0.5, 2.5]);