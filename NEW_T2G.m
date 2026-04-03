close all, clear, clc

%% 1. Select and Load Dual-Color TIFFs
% Open UI dialog to select the GREEN file
[file_name_green, path_name] = uigetfile({'*GREEN.tif;*GREEN.tiff', 'Green Image Files (*GREEN.tif)'}, ...
    'Select the Corrected GREEN File', 'MultiSelect', 'off');

if isequal(file_name_green, 0)
    disp('User selected Cancel. Exiting...');
    return;
end

% Automatically infer the RED filename
file_name_red = strrep(file_name_green, '_GREEN', '_RED');
current_file_green = fullfile(path_name, file_name_green);
current_file_red = fullfile(path_name, file_name_red);

if ~exist(current_file_red, 'file')
    error(['Could not find matching RED file: ', file_name_red, '. Make sure it is in the same folder.']);
end

disp(['Loading GREEN channel: ', file_name_green]);
M2_gcam = loadtiff(string(current_file_green));

disp(['Loading RED channel: ', file_name_red]);
M2_red = loadtiff(string(current_file_red));

[H, W, T] = size(M2_gcam);

%% 2. Extract Metadata (Timevec, Event Marker, and ROIs) from GREEN TIFF
disp('Extracting metadata from ImageDescription...');
tiff_info = imfinfo(current_file_green);

if isfield(tiff_info, 'ImageDescription') && ~isempty(tiff_info(1).ImageDescription)
    img_desc = tiff_info(1).ImageDescription;
    
    % Extract Time Vector
    t_tokens = regexp(img_desc, 'TimeVector_Seconds=\[(.*?)\]', 'tokens');
    timevec = str2num(t_tokens{1}{1}); 
    timevec = timevec(:) / 60; % Convert to minutes
    
    % Extract Event Marker (Interval Change)
    e_tokens = regexp(img_desc, 'EventMarker=Frame_(\d+)', 'tokens');
    if ~isempty(e_tokens)
        marks = str2double(e_tokens{1}{1});
    else
        marks = 35; 
    end
    
    % Extract ROI Geometries
    r_tokens = regexp(img_desc, 'ROI_\d+=\[(.*?),(.*?),(.*?)\]', 'tokens');
    roiCount = length(r_tokens);
    
    roi_data = zeros(roiCount, 3); % [X, Y, Radius]
    for r = 1:roiCount
        roi_data(r, 1) = str2double(r_tokens{r}{1}); 
        roi_data(r, 2) = str2double(r_tokens{r}{2}); 
        roi_data(r, 3) = str2double(r_tokens{r}{3}); 
    end
    disp(['Successfully extracted timevec and ', num2str(roiCount), ' ZEN ROIs.']);
else
    error('No ImageDescription metadata found in this TIFF file.');
end
%%
evnt = {'20 Thapsigargin'}; % Set your injection label here

%% 3. Create Spatial Masks from ROI Metadata
disp('Generating ROI masks...');
[Xgrid, Ygrid] = meshgrid(1:W, 1:H); 
islands_PixelIdxList = cell(roiCount, 1);

% Visualize the ROIs on the mean GREEN image
aveMap = mean(M2_gcam, 3);
figure(1); clf;
imagesc(aveMap); colormap bone; axis square; hold on;
title('Mean GCaMP Image with Extracted ZEN ROIs');

for c = 1:roiCount
    cx = roi_data(c, 1); cy = roi_data(c, 2); rad = roi_data(c, 3);
    mask = (Xgrid - cx).^2 + (Ygrid - cy).^2 <= rad^2;
    islands_PixelIdxList{c} = find(mask);
    
    viscircles([cx, cy], rad, 'Color', 'g', 'LineWidth', 1);
    text(cx, cy, num2str(c), 'Color', 'y', 'FontSize', 10, 'HorizontalAlignment', 'center');
end
hold off;

%% 4. Parallel Signal Extraction (Green and Red)
numCells = roiCount;
roiMean_gcam = zeros(numCells, T);
roiMean_red = zeros(numCells, T);

% Flatten PixelIdxList for lightning-fast accumarray
all_pixel_idx = []; all_cell_ids = [];
for c = 1:numCells
    curr_pixels = islands_PixelIdxList{c};
    all_pixel_idx = [all_pixel_idx; curr_pixels(:)]; 
    all_cell_ids = [all_cell_ids; repmat(c, length(curr_pixels), 1)];
end

disp('Starting parallel signal extraction across time points...');
tic;
parfor timei = 1:T
    % Read frames ONCE
    tmp_g = double(M2_gcam(:,:,timei)); 
    tmp_r = double(M2_red(:,:,timei));
    
    % Extract pixels
    vals_g = tmp_g(all_pixel_idx);
    vals_r = tmp_r(all_pixel_idx);
    
    % Calculate the mean for ALL cells instantly
    roiMean_gcam(:, timei) = accumarray(all_cell_ids, vals_g, [numCells, 1], @mean);
    roiMean_red(:, timei)  = accumarray(all_cell_ids, vals_r, [numCells, 1], @mean);
end
elapsedTime = toc; disp(['Extraction time: ' num2str(elapsedTime) ' seconds']);

%% 5. Calculate Ratio and dF/F0 (dR/R0)
disp('Calculating Ratio and standardizing to dF/F0...');

% Protect against division by near-zero in the Red channel
roiMean_red(roiMean_red < 0.01) = 0.01; 

% 1. Raw Ratio (Green / Red)
roiMean_ratio = roiMean_gcam ./ roiMean_red;

% CLEANUP: Clamp extreme artifacts in the Ratio to protect the Mean/SEM graph
roiMean_ratio(isnan(roiMean_ratio) | isinf(roiMean_ratio)) = 0;
roiMean_ratio(roiMean_ratio > 20) = 20; % Cap massive ratio noise spikes 
roiMean_ratio(roiMean_ratio < 0) = 0;   % Ratios cannot be negative

% 2. Calculate F0 (Baseline of the Ratio) before injection
baseline_start = max(1, marks(1)-11);
baseline_end = max(2, marks(1)-1);
F0_ratio = mean(roiMean_ratio(:, baseline_start:baseline_end), 2);

% Protect the baseline against division by near-zero 
F0_ratio(F0_ratio < 0.01) = 0.01; 

% 3. Calculate dF/F0 of the Ratio 
df_ratio = bsxfun(@rdivide, roiMean_ratio, F0_ratio) - 1;

% 4. CLEANUP: Clamp extreme artifacts in the dF/F0 
df_ratio(isnan(df_ratio) | isinf(df_ratio)) = 0; 
df_ratio(df_ratio > 10) = 10;   % Cap positive noise spikes 
df_ratio(df_ratio < -1.5) = -1.5; % Cap negative noise spikes
%% 6. Visualization: Raw GCaMP (Green)

disp('Plotting Raw GCaMP...');
mean_gcam = mean(roiMean_gcam, 1)';
sem_gcam  = std(roiMean_gcam, 0, 1)' / sqrt(numCells);

figure(2); clf; set(gcf, 'Position', [100, 100, 1200, 800]);

% --- Heatmap ---
subplot(2,2,1);
surf(timevec, 1:numCells, roiMean_gcam, 'LineStyle', 'none'); view(2);
axis tight; xlabel('Time (min.)'); ylabel('ROI number'); colorbar; 
title('Heatmap: Raw GCaMP');

% --- All Traces ---
subplot(2,2,2);
plot(timevec, roiMean_gcam); axis tight;
xlabel('Time (min.)'); ylabel('Intensity (a.u.)'); 
title('All Traces: Raw GCaMP');

% --- Mean Trace with SEM ---
subplot(2,2,[3 4]); hold on;
fill([timevec; flipud(timevec)], ...
     [mean_gcam+sem_gcam; flipud(mean_gcam-sem_gcam)], ...
     [0.8 0.95 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_gcam, 'Color', [0 0.6 0], 'LineWidth', 2);
ylabel('Mean Brightness (a.u.)'); xlabel('Time (min.)'); 
title('Mean Raw GCaMP Trace with SEM');

% Dynamic Y-Axis limits for Drug Bar
y_max = max(mean_gcam + sem_gcam);
y_range = y_max - min(mean_gcam - sem_gcam);
bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
set(gca, 'xlim', timevec([1 end]), 'ylim', [min(mean_gcam-sem_gcam)-0.1*y_range, text_y + 0.1*y_range], 'TickDir', 'none', 'fontsize', 12);

% Add drug treatment line and text
% Define the start and end time for the horizontal bar
t_start = timevec(marks(1));
t_end = timevec(min(T, T)); % Adjust '300' if your treatment duration differs


% Draw thick horizontal treatment bar
line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 3);

% Place text exactly in the middle and resting on top of the bar
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% 7. Visualization: Ratio (GCaMP / tdTomato)
disp('Plotting Ratio...');
mean_ratio = mean(roiMean_ratio, 1)';
sem_ratio  = std(roiMean_ratio, 0, 1)' / sqrt(numCells);

figure(3); clf; set(gcf, 'Position', [150, 150, 1200, 800]);

% --- Heatmap ---
subplot(2,2,1);
surf(timevec, 1:numCells, roiMean_ratio, 'LineStyle', 'none'); view(2);
axis tight; xlabel('Time (min.)'); ylabel('ROI number'); colorbar; 
title('Heatmap: Ratio (Green/Red)');

% --- All Traces ---
subplot(2,2,2);
plot(timevec, roiMean_ratio); axis tight;
xlabel('Time (min.)'); ylabel('Ratio (a.u.)'); 
title('All Traces: Ratio');

% --- Mean Trace with SEM ---
subplot(2,2,[3 4]); hold on;
fill([timevec; flipud(timevec)], ...
     [mean_ratio+sem_ratio; flipud(mean_ratio-sem_ratio)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_ratio, 'b', 'LineWidth', 2);
ylabel('Mean Ratio (a.u.)'); xlabel('Time (min.)'); 
title('Mean Ratio Trace with SEM');

% Dynamic Y-Axis limits for Drug Bar (Safe from NaNs and flatlines)
valid_idx = ~isnan(mean_ratio) & ~isnan(sem_ratio);

if any(valid_idx)
    y_max = max(mean_ratio(valid_idx) + sem_ratio(valid_idx));
    y_min = min(mean_ratio(valid_idx) - sem_ratio(valid_idx));
else
    y_max = 2; y_min = 0; % Fallback if completely corrupted
end

y_range = y_max - y_min;
if y_range == 0 % Fallback if data is a perfectly flat line
    y_range = 1; 
end

bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
plot_ymin = max(0, y_min - 0.1 * y_range); % Prevent negative bottom limits for ratio
plot_ymax = text_y + 0.1 * y_range;

set(gca, 'xlim', timevec([1 end]), 'ylim', [plot_ymin, plot_ymax], 'TickDir', 'none', 'fontsize', 12);

% Add drug treatment line and text
t_start = timevec(marks(1));
t_end = timevec(min(T, T)); % Adjust '300' if your treatment duration differs

line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 3);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 14, 'FontWeight', 'bold');
box off;
%% 8. Visualization: dF/F0 (Normalized Ratio) 
disp('Plotting dF/F0 of Ratio...');
mean_df = mean(df_ratio, 1)';
sem_df  = std(df_ratio, 0, 1)' / sqrt(numCells);

figure(4); clf; set(gcf, 'Position', [200, 200, 1200, 800]);

% --- Heatmap ---
subplot(2,2,1);
surf(timevec, 1:numCells, df_ratio, 'LineStyle', 'none'); view(2);
set(gca, 'ylim', [0.5 numCells+0.5], 'xlim', timevec([1 end]), ...
    'clim', [-0.5 10], 'YDir', 'reverse');
colorbar; title('Heatmap: \Delta(Ratio) / Ratio_0');
xlabel('Time (min.)'); ylabel('ROI number');

% --- All Traces ---
subplot(2,2,2);
plot(timevec, df_ratio); 
set(gca, 'xlim', timevec([1 end]), 'ylim', [-0.5 12]);
xlabel('Time (min.)'); ylabel('\DeltaF/F_0'); 
title('All Traces: \Delta(Ratio) / Ratio_0');

% --- Mean Trace with SEM ---
subplot(2,2,[3 4]); hold on;
fill([timevec; flipud(timevec)], ...
     [mean_df+sem_df; flipud(mean_df-sem_df)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_df, 'r', 'LineWidth', 2);
ylabel('Mean \DeltaF/F_0'); xlabel('Time (min.)'); 
title('Mean \DeltaF/F_0 Trace with SEM');

% Dynamic Y-Axis limits for Drug Bar (Safe from NaNs and flatlines)
valid_idx = ~isnan(mean_df) & ~isnan(sem_df);

if any(valid_idx)
    y_max = max(mean_df(valid_idx) + sem_df(valid_idx));
    y_min = min(mean_df(valid_idx) - sem_df(valid_idx));
else
    y_max = 1; y_min = -0.5; % Fallback if completely corrupted
end

y_range = y_max - y_min;
if y_range == 0 % Fallback if data is a perfectly flat line
    y_range = 1; 
end

bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
plot_ymin = y_min - 0.1 * y_range;
plot_ymax = text_y + 0.1 * y_range;

set(gca, 'xlim', timevec([1 end]), 'ylim', [plot_ymin, plot_ymax], 'TickDir', 'none', 'fontsize', 12);

% Add drug treatment line and text
t_start = timevec(marks(1));
t_end = timevec(min(T, T)); % Adjust '300' if your treatment duration differs


line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 3);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 14, 'FontWeight', 'bold');
box off;