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
    
    % --- Extract Time Vector ---
    time_tokens = regexp(img_desc, 'TimeVector_Seconds=\[(.*?)\]', 'tokens');
    if ~isempty(time_tokens)
        time_str = time_tokens{1}{1};
        timevec = str2num(time_str); % str2num automatically parses comma-separated strings
        timevec = timevec(:); % Ensure it is a column vector
        disp(['Successfully extracted timevec: ', num2str(length(timevec)), ' frames.']);
    else
        warning('TimeVector_Seconds format not found within the ImageDescription tag.');
        timevec = [];
    end
    
    % --- Extract ROI Geometries ---
    % Searches for ROI_XX=[X,Y,Radius]
    roi_tokens = regexp(img_desc, 'ROI_\d+=\[([\d\.]+),([\d\.]+),([\d\.]+)\]', 'tokens');
    numCells = length(roi_tokens);
    roiCount = numCells;
    disp(['Successfully extracted ', num2str(numCells), ' embedded ROIs.']);
else
    error('No ImageDescription metadata found in this TIFF file. Cannot extract ROIs.');
end

% Estimate sampling rate
marks = find(abs(diff(round(diff(timevec)))) >= 1.5) + 2;
timevec = timevec / 60;

%%
evnt = {'20 \muM Thap'}; % Set your injection label here
%evnt = {'Akh 10ng'}; % Set your injection label here
%% 3. Build Spatial Masks for ROIs from Metadata (Interactive & Save)
disp('Building interactive spatial masks for ROIs...');

% Define the path where adjusted ROIs would be saved
[~, base_name_noext, ~] = fileparts(file_name_green);
roi_save_path = fullfile(path_name, [base_name_noext, '_Adjusted_ROIs.mat']);

% To visualize the ROIs on an average projection
aveMap = squeeze(mean(M2_red, 3));
fig1 = figure(1); clf;

% Create axes explicitly to attach UI elements safely
ax = axes('Parent', fig1);
imagesc(ax, aveMap); colormap hot; axis square; 
title('Adjust ROIs: Drag to move, drag edge to resize, click & press "Delete" to remove.');
hold(ax, 'on');

% 1. Pre-draw imported ROIs (Check for saved file first, fallback to original metadata)
if exist(roi_save_path, 'file')
    disp(['Found previously adjusted ROIs. Loading from: ', roi_save_path]);
    
    % Load the previously saved coordinates and count
    load(roi_save_path, 'final_roi_coords', 'numCells');
    roiCount = numCells; % Keep counts synced
    
    for c = 1:numCells
        cx  = final_roi_coords(c, 1);
        cy  = final_roi_coords(c, 2);
        rad = final_roi_coords(c, 3);
        
        % Draw interactive ROI
        drawcircle(ax, 'Center', [cx, cy], 'Radius', rad, 'Label', num2str(c), ...
                   'Color', 'w', 'FaceAlpha', 0.1, 'LabelVisible', 'on');
    end
else
    disp('No previously adjusted ROIs found. Loading original metadata ROIs...');
    
    for c = 1:numCells
        % Parse coordinates from regex tokens
        cx  = str2double(roi_tokens{c}{1});
        cy  = str2double(roi_tokens{c}{2});
        rad = str2double(roi_tokens{c}{3});
        
        % Draw interactive ROI
        drawcircle(ax, 'Center', [cx, cy], 'Radius', rad, 'Label', num2str(c), ...
                   'Color', 'w', 'FaceAlpha', 0.1, 'LabelVisible', 'on');
    end
end

% 2. Initialize State for the While Loop
setappdata(fig1, 'action', 'wait');

% 3. Add UI Buttons for user interaction
% Buttons now securely pass states to the loop instead of executing heavy functions directly
btn_add = uicontrol('Parent', fig1, 'Style', 'pushbutton', 'String', 'Add New ROI', ...
                    'Position', [20 20 100 30], ...
                    'Callback', @(~,~) setappdata(fig1, 'action', 'add'));
                    
btn_done = uicontrol('Parent', fig1, 'Style', 'pushbutton', 'String', 'Done Adjusting', ...
                     'Position', [130 20 120 30], ...
                     'Callback', @(~,~) setappdata(fig1, 'action', 'done'));

% 4. The Interactive While Loop
disp('Waiting for user to adjust ROIs. Click "Add New ROI" to draw, or "Done Adjusting" when finished...');

while isvalid(fig1)
    action = getappdata(fig1, 'action');
    
    if strcmp(action, 'add')
        % Reset the state immediately so it doesn't infinite loop
        setappdata(fig1, 'action', 'wait');
        
        % Safely draw the new ROI in the main thread
        drawcircle(ax, 'Color', 'g', 'FaceAlpha', 0.1, 'Label', 'New');
        
    elseif strcmp(action, 'done')
        break; % Exit the loop when user clicks "Done Adjusting"
    end
    
    pause(0.1); % Small pause prevents MATLAB from freezing and polls for your clicks
end

% 5. --- Post-Interaction Processing ---
disp('Compiling and saving final ROI masks...');

% Find all circle objects currently on the axes
roi_objs = findobj(ax, 'Type', 'images.roi.circle');

% findobj returns newest objects first. Flip to roughly maintain original numbering order
roi_objs = flipud(roi_objs); 

% Update counts based on user additions/deletions
numCells = length(roi_objs);
roiCount = numCells;
disp(['Final ROI count: ', num2str(numCells)]);

all_pixel_idx = [];
all_cell_ids = [];
[X_grid, Y_grid] = meshgrid(1:W, 1:H); 

% Pre-allocate a matrix to save the final coordinates
final_roi_coords = zeros(numCells, 3); 

% 6. Build final masks from the updated objects and record coordinates
for c = 1:numCells
    % Clean up labels to reflect final sequential numbering
    roi_objs(c).Label = num2str(c);
    roi_objs(c).Color = 'y'; % Turn them yellow to show they are "locked in"
    
    % Get final adjusted coordinates
    cx = roi_objs(c).Center(1);
    cy = roi_objs(c).Center(2);
    rad = roi_objs(c).Radius;
    
    % Record to our export matrix
    final_roi_coords(c, :) = [cx, cy, rad];
    
    % Create spatial mask for this circle
    mask = (X_grid - cx).^2 + (Y_grid - cy).^2 <= rad.^2;
    curr_pixels = find(mask);
    
    % Flatten into 1D arrays for fast accumarray processing
    all_pixel_idx = [all_pixel_idx; curr_pixels(:)]; 
    all_cell_ids = [all_cell_ids; repmat(c, length(curr_pixels), 1)];
end

hold(ax, 'off');

% Clean up UI buttons so the figure looks clean for saving/review
delete(btn_add);
delete(btn_done);

% 7. Save the final coordinates to a .mat file
save(roi_save_path, 'final_roi_coords', 'numCells');
disp(['Saved adjusted ROI coordinates to: ', roi_save_path]);

%% 4. Parallel Signal Extraction (Green and Red)
disp('Starting parallel signal extraction across time points...');

roiMean_gcam = zeros(numCells, T);
roiMean_red = zeros(numCells, T);

% NOTE: all_pixel_idx and all_cell_ids are already built perfectly in Step 3!

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
%df_ratio(df_ratio > 10) = 10;     % Cap positive noise spikes 
%df_ratio(df_ratio < -1.5) = -1.5; % Cap negative noise spikes

% 5. --- NEW: GLOBAL SAVITZKY-GOLAY SMOOTHING ---
disp('Applying Savitzky-Golay filter (Order 3, Frame 11) to dF/F0...');
for i = 1:numCells
    % Smooth the individual cell's ratio trace
    roiMean_ratio(i, :) = sgolayfilt(roiMean_ratio(i, :), 3, 11);
    % Smooth the individual cell's dF/F0 trace
    df_ratio(i, :) = sgolayfilt(df_ratio(i, :), 3, 11);
end

%% 6. Visualization: Raw GCaMP (Green)
disp('Plotting Raw GCaMP...');
mean_gcam = mean(roiMean_gcam, 1)';
sem_gcam  = std(roiMean_gcam, 0, 1)' / sqrt(numCells);

% Define drug timing bounds once to use across all plots
t_start = timevec(marks(1));
t_end = timevec(end); 

figure(2); clf; set(gcf, 'Position', [100, 100, 1200, 800]);

% --- Heatmap ---
subplot(2,2,1);
surf(timevec, linspace(0.5,0.5+numCells,numCells), roiMean_gcam, 'LineStyle', 'none'); view(2);
axis tight; set(gca, 'ylim', [0.5 numCells+0.5], 'xlim', timevec([1 end]), 'YDir', 'reverse');
xlabel('Time (min.)'); ylabel('ROI number'); colorbar; 
title({'Heatmap: Raw GCaMP', ''});
box off;
% Add dashed line and text for heatmap
xline(t_start, '--w', 'LineWidth', 1); % White dashed line for visibility on heatmap
text(t_start, 0, evnt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');

% --- All Traces ---
subplot(2,2,2);
plot(timevec, roiMean_gcam); axis tight;
xlabel('Time (min.)'); ylabel('Intensity (a.u.)'); 
title({'All Traces: Raw GCaMP', ''});
box off;
% Add drug treatment line and text for All Traces
y_max_all = max(roiMean_gcam(:));
y_min_all = min(roiMean_gcam(:));
y_range_all = max(y_max_all - y_min_all, 1e-3);
bar_y_all = y_max_all + 0.15 * y_range_all;
text_y_all = y_max_all + 0.25 * y_range_all;
set(gca, 'ylim', [y_min_all - 0.1*y_range_all, text_y_all + 0.1*y_range_all]);
line([t_start t_end], [bar_y_all bar_y_all], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y_all + (0.02 * y_range_all), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

% --- Mean Trace with SEM ---
subplot(2,2,[3 4]); hold on;
fill([timevec; flipud(timevec)], ...
     [mean_gcam+sem_gcam; flipud(mean_gcam-sem_gcam)], ...
     [0.8 0.95 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_gcam, 'Color', [0 0.6 0], 'LineWidth', 2);
ylabel('Mean Brightness (a.u.)'); xlabel('Time (min.)'); 
title({'Mean Raw GCaMP Trace with SEM', ''});
box off;
% Dynamic Y-Axis limits for Drug Bar
y_max = max(mean_gcam + sem_gcam);
y_min = min(mean_gcam - sem_gcam);
y_range = max(y_max - y_min, 1e-3);
bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
set(gca, 'xlim', timevec([1 end]), 'ylim', [y_min - 0.1*y_range, text_y + 0.1*y_range], 'TickDir', 'none', 'FontSize', 12);
% Add drug treatment line and text
line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);

% Apply Arial and Normal Font Weight globally to Figure 2
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');

%% 7. Visualization: Ratio (GCaMP / tdTomato)
disp('Plotting Filtered Ratio...');
mean_ratio = mean(roiMean_ratio, 1)';
sem_ratio  = std(roiMean_ratio, 0, 1)' / sqrt(numCells);

figure(3); clf; set(gcf, 'Position', [150, 150, 1200, 800]);

% --- Heatmap ---
subplot(2,2,1);
surf(timevec, linspace(0.5,0.5+numCells,numCells), roiMean_ratio, 'LineStyle', 'none'); view(2);
axis tight; set(gca, 'ylim', [0.5 numCells+0.5], 'xlim', timevec([1 end]), 'YDir', 'reverse');
xlabel('Time (min.)'); ylabel('ROI number'); colorbar; 
title({'Heatmap: Filtered Ratio (Green/Red)', ''});
box off;
% Add dashed line and text for heatmap
xline(t_start, '--w', 'LineWidth', 2);
text(t_start, 0, evnt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');

% --- All Traces ---
subplot(2,2,2);
plot(timevec, roiMean_ratio); axis tight;
xlabel('Time (min.)'); ylabel('Ratio (a.u.)'); 
title({'All Traces: Filtered Ratio', ''});
box off;
% Add drug treatment line and text for All Traces
y_max_all = max(roiMean_ratio(:));
y_min_all = min(roiMean_ratio(:));
y_range_all = max(y_max_all - y_min_all, 1e-3);
bar_y_all = y_max_all + 0.15 * y_range_all;
text_y_all = y_max_all + 0.25 * y_range_all;
set(gca, 'ylim', [max(0, y_min_all - 0.1*y_range_all), text_y_all + 0.1*y_range_all]);
line([t_start t_end], [bar_y_all bar_y_all], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y_all + (0.02 * y_range_all), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

% --- Mean Trace with SEM ---
subplot(2,2,[3 4]); hold on;
fill([timevec; flipud(timevec)], ...
     [mean_ratio+sem_ratio; flipud(mean_ratio-sem_ratio)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_ratio, 'b', 'LineWidth', 1);
ylabel('Mean Ratio (a.u.)'); xlabel('Time (min.)'); 
title({'Mean Filtered Ratio Trace with SEM', ''});
box off;
% Dynamic Y-Axis limits for Drug Bar
valid_idx = ~isnan(mean_ratio) & ~isnan(sem_ratio);
if any(valid_idx)
    y_max = max(mean_ratio(valid_idx) + sem_ratio(valid_idx));
    y_min = min(mean_ratio(valid_idx) - sem_ratio(valid_idx));
else
    y_max = 2; y_min = 0; 
end
y_range = max(y_max - y_min, 1e-3); 
bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
plot_ymin = max(0, y_min - 0.1 * y_range); 
plot_ymax = text_y + 0.1 * y_range;
set(gca, 'xlim', timevec([1 end]), 'ylim', [plot_ymin, plot_ymax], 'TickDir', 'none', 'FontSize', 12);
% Add drug treatment line and text
line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);

% Apply Arial and Normal Font Weight globally to Figure 3
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');

%% 8. Visualization: dF/F0 (Normalized Ratio) 
disp('Plotting Filtered dF/F0 of Ratio...');
mean_df = mean(df_ratio, 1)';
sem_df  = std(df_ratio, 0, 1)' / sqrt(numCells);

figure(4); clf; set(gcf, 'Position', [200, 200, 1200, 800]);

% --- Heatmap ---
subplot(2,2,1);
surf(timevec, linspace(0.5,0.5+numCells,numCells), df_ratio, 'LineStyle', 'none'); view(2);
axis tight; set(gca, 'ylim', [0.5 numCells+0.5], 'xlim', timevec([1 end]), 'YDir', 'reverse');
colorbar; 
title({'Heatmap: Filtered \Delta(Ratio) / Ratio_0', ''});
xlabel('Time (min.)'); ylabel('ROI number');
box off;
% Add dashed line and text for heatmap
xline(t_start, '--w', 'LineWidth', 2);
text(t_start, 0.2, evnt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');

% --- All Traces ---
subplot(2,2,2);
plot(timevec, df_ratio); axis tight;
xlabel('Time (min.)'); ylabel('\DeltaF/F_0'); 
title({'All Traces: Filtered \Delta(Ratio) / Ratio_0', ''});
box off;
% Add drug treatment line and text for All Traces
y_max_all = max(df_ratio(:));
y_min_all = min(df_ratio(:));
y_range_all = max(y_max_all - y_min_all, 1e-3);
bar_y_all = y_max_all + 0.15 * y_range_all;
text_y_all = y_max_all + 0.25 * y_range_all;
set(gca, 'ylim', [y_min_all - 0.1*y_range_all, text_y_all + 0.1*y_range_all]);
line([t_start t_end], [bar_y_all bar_y_all], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y_all + (0.02 * y_range_all), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

% --- Mean Trace with SEM ---
subplot(2,2,[3 4]); hold on;
fill([timevec; flipud(timevec)], ...
     [mean_df+sem_df; flipud(mean_df-sem_df)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_df, 'r', 'LineWidth', 2);
ylabel('Mean \DeltaF/F_0'); xlabel('Time (min.)'); 
title({'Mean Filtered \DeltaF/F_0 Trace with SEM', ''});
box off;
% Dynamic Y-Axis limits for Drug Bar
valid_idx = ~isnan(mean_df) & ~isnan(sem_df);
if any(valid_idx)
    y_max = max(mean_df(valid_idx) + sem_df(valid_idx));
    y_min = min(mean_df(valid_idx) - sem_df(valid_idx));
else
    y_max = 1; y_min = -0.5;
end
y_range = max(y_max - y_min, 1e-3); 
bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
plot_ymin = y_min - 0.1 * y_range;
plot_ymax = text_y + 0.1 * y_range;
set(gca, 'xlim', timevec([1 end]), 'ylim', [plot_ymin, plot_ymax], 'TickDir', 'none', 'FontSize', 12);
% Add drug treatment line and text
line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);

% Apply Arial and Normal Font Weight globally to Figure 4
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');

%% 9. (Optional) All Cell Time Tracing
disp('Generating individual traces for all ROIs...');
for i=1:numCells
    figure(100+i); clf;
    set(gcf,'Units', 'inches', 'Position', [1 1 15 5]);
    plot(timevec, df_ratio(i,:), '-k', 'LineWidth', 1);
    
    ylabel('\DeltaF/F_0'); xlabel('Time (min.)');
    title({'', sprintf('Fluorescence of ROI %d (Filtered)', i), ''});
    
    set(gca, 'xlim', timevec([1 end]), 'ylim', [-0.8 10], 'TickDir', 'none');
    xline(timevec(marks(1)), '--b', evnt, 'FontWeight', 'bold', 'LineWidth', 1);
    box off;
    
    % Apply Arial formatting
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
    set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');
end

%% 10. (Optional) Single Cell Highlight Tracing
j = 4; % Define specific cell to highlight
if j <= numCells
    figure(200+j); clf;
    set(gcf,'Units', 'inches', 'Position', [1 1 15 5]);
    plot(timevec, df_ratio(j,:), '-k', 'LineWidth', 1);
    
    ylabel('\DeltaF/F_0'); xlabel('Time (min.)');
    title({'', sprintf('Fluorescence Highlight: ROI %d (Filtered)', j), ''});
    
    set(gca,'xlim', timevec([1 end]), 'ylim', [-1.5 12], 'YTick', (0:2:10), 'TickDir', 'none');
    box off;
    
    % Drug treatment line
    t_start = timevec(marks(1));
    t_end = timevec(end);
    line([t_start t_end], [10.5 10.5], 'color', 'k', 'Linewidth', 2);
    text(t_start + 0.5*(t_end - t_start), 11.5, evnt, 'HorizontalAlignment', 'center');
    
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
    set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');
end

%% 11. Autonomous Peak Analysis (Dynamic Tuning)
disp('Running Autonomous Peak Analysis on Filtered Data...');

% Initialize cell arrays
hdr_pks = cell(numCells,1);
hdr_locs = cell(numCells,1);

for i = 1:numCells
    % --- 1. SIGNAL SMOOTHING ---
    % Signal is already smoothed globally in Section 5. 
    smooth_signal = df_ratio(i, :);
    
    % Extract post-injection window for peak finding
    signal_post = smooth_signal(marks(1):end); 
    time_post = timevec(marks(1):end);
    
    % --- 2. DYNAMIC NOISE ESTIMATION ---
    % Evaluate pre-injection baseline stability (Exclude failed washings)
    baseline_check_idx = max(1, marks(1)-20);
    baseline_signal = smooth_signal(baseline_check_idx : marks(1)-1);
    
    if mean(baseline_signal) < 1 % Exclude failed washings
        
        % Calculate robust noise (sigma) using Median Absolute Deviation (MAD)
        baseline_median = median(baseline_signal);
        robust_sigma = median(abs(baseline_signal - baseline_median)) / 0.6745;
        robust_sigma = max(robust_sigma, 0.01); % Prevent sigma from being exactly 0
        
        % --- 3. DYNAMIC PARAMETER TUNING ---
        % Peak must stand out from its local surroundings by 2x the noise floor
        dyn_prom = max(0.2, 1.5 * robust_sigma); 
        
        % Peak must reach an absolute height of at least 1x the noise floor above 0
        dyn_height = max(0.2, 1 * robust_sigma); 
        
        % Minimum peak width (in minutes). 
        dt_min = mean(diff(timevec)); 
        dyn_width = dt_min * 13; % Must be at least ~x frames wide
        
        % --- 4. EXECUTE FINDPEAKS ---
        [pks_post, locs_post] = findpeaks(signal_post, time_post, ...
            'MinPeakProminence', dyn_prom, ...
            'MinPeakHeight', dyn_height, ...
            'MinPeakWidth', dyn_width);
        
        % --- 5. VISUALIZATION (Verification) ---
        figure(300+i); clf; set(gcf, 'Position', [150, 150, 1000, 600]);
        
        % Subplot 1: Filtered Ratio Check
        subplot(2,1,1); hold on;
        plot(time_post, roiMean_ratio(i, marks(1):end), 'k', 'LineWidth', 1.5); 
        set(gca, 'xlim', timevec([marks(1) end]), 'TickDir', 'none');
        title(sprintf('ROI %d: Filtered Ratio', i)); box off;
        
        % Subplot 2: dF/F0 Peak Detection
        subplot(2,1,2); hold on;
        plot(time_post, signal_post, 'r', 'LineWidth', 1.5); % Dark red smoothed dF/F0
        
        % Plot detected peaks
        if ~isempty(locs_post)
            plot(locs_post, pks_post + 0.35, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineStyle', 'none');
        end
        
        % Draw dynamic threshold line for visual reference
        yline(dyn_height, '--', 'Color', [0.5 0.5 0.5], 'Label', sprintf('Threshold: %.2f', dyn_height));
        
        set(gca, 'xlim', timevec([marks(1) end]), 'ylim', [-1 max(10, max(signal_post)+2)], 'TickDir', 'none');
        title(sprintf('ROI %d: Autonomous Peak Detection (Prominence: %.2f)', i, dyn_prom)); box off;
        
        set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
        set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');
        
        % Store entire peak array in one cell
        hdr_pks{i} = pks_post;
        hdr_locs{i} = locs_post;
    else
        % Excluded due to unstable baseline
        hdr_pks{i} = [];
        hdr_locs{i} = [];
    end
end

%% 12. Aggregate Peak Statistics
disp('Calculating Peak Statistics...');
peak_counts = cellfun(@numel, hdr_pks);
peak_counts(peak_counts == 0) = NaN; 
mean_oscil = mean(peak_counts, 1, 'omitnan');

peak_ave = cellfun(@mean, hdr_pks);
ave_amp = mean(peak_ave, 1,'omitmissing');

act_cell = sum(peak_counts > 0);
per_act_cell = (act_cell / numCells) * 100;

sem_amp  = std(peak_ave, 0, 1) / sqrt(act_cell);
sem_oscil  = std(peak_counts, 0, 1, 'omitnan') / sqrt(act_cell);

%% 13. Analyze 1st Peaks vs Remaining Peaks & QC Graphing
disp('Extracting 1st Peak Dynamics (Rise, Decay, Duration) and Generating QC Plots...');
st_hdr_locs = zeros(numCells,1);
st_hdr_pks = zeros(numCells,1);
st_int = zeros(numCells,1);
st_loss = zeros(numCells,1);

for k = 1:numCells
    if ~isempty(hdr_locs{k}) && ~isnan(hdr_locs{k}(1))
        
        st_hdr_locs(k) = hdr_locs{k}(1);
        st_hdr_pks(k)  = hdr_pks{k}(1);
        
        % Robustly find indices instead of using exact float == matches
        [~, peak_idx] = min(abs(timevec - st_hdr_locs(k)));
        
        % --- 1. RISE DYNAMICS ---
        pre_peak_signal = df_ratio(k, marks(1):peak_idx);
        st_int_bline = prctile(pre_peak_signal, 15);
        
        idx_rise_start = find(pre_peak_signal <= st_int_bline, 1, 'last');
        if isempty(idx_rise_start), idx_rise_start = 1; end
        st_int(k) = timevec(marks(1) + idx_rise_start - 1);
        
        % --- 2. DECAY DYNAMICS ---
        post_peak_signal = df_ratio(k, peak_idx:end);
        
        if isempty(find(post_peak_signal <= 0.3, 1))
            st_loss(k) = timevec(end); % Fallback to end of trace
            decay_baseline = 0.3;      % Fallback baseline
        else
            if length(hdr_pks{k}) == 1
                decay_baseline = prctile(post_peak_signal, 15);
                idx_decay_end = find(post_peak_signal <= decay_baseline, 1, 'first');
                if isempty(idx_decay_end), idx_decay_end = length(post_peak_signal); end
                st_loss(k) = timevec(peak_idx + idx_decay_end - 1);
            else
                [~, peak2_idx] = min(abs(timevec - hdr_locs{k}(2)));
                st_nd_dy = df_ratio(k, peak_idx:peak2_idx);
                decay_baseline = prctile(st_nd_dy, 18);
                idx_decay_end = find(st_nd_dy <= decay_baseline, 1, 'first');
                if isempty(idx_decay_end), idx_decay_end = length(st_nd_dy); end
                st_loss(k) = timevec(peak_idx + idx_decay_end - 1);
            end
        end
        
        % --- 3. QC GRAPH: 1st Peak Dynamics ---
        figure(400+k); clf;
        set(gcf, 'Position', [200, 200, 1000, 500]);
        hold on;
        
        % Plot the main trace (already filtered in Sec 5)
        plot(timevec, df_ratio(k,:), '-k', 'LineWidth', 1.5);
        
        % Plot all detected peaks for this cell
        plot(hdr_locs{k}, hdr_pks{k} + 0.3, 'rv', 'MarkerFaceColor', 'r', 'LineStyle', 'none');
        
        % Draw Rectangle for Rise Phase (Yellow)
        rise_w = max(0.01, st_hdr_locs(k) - st_int(k));
        rise_h = max(0.01, st_hdr_pks(k) + 0.35 - st_int_bline);
        rectangle('Position', [st_int(k), st_int_bline, rise_w, rise_h], ...
            'FaceColor', 'y', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            
        % Draw Rectangle for Decay Phase (Magenta)
        decay_w = max(0.01, st_loss(k) - st_hdr_locs(k));
        decay_h = max(0.01, st_hdr_pks(k) + 0.35 - decay_baseline);
        rectangle('Position', [st_hdr_locs(k), decay_baseline, decay_w, decay_h], ...
            'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        
        % Formatting
        ylabel('\DeltaF/F_0');
        xlabel('Time (min.)');
        title({'', sprintf('1st Peak Dynamics QC: ROI %d (Filtered)', k), ''});
        
        % Dynamic Y-lim to ensure large peaks don't get cut off
        y_max = max(10, max(df_ratio(k,:)) + 1);
        set(gca, 'xlim', timevec([1 end]), 'ylim', [-0.8 y_max], 'TickDir', 'none');
        box off;
        
        % Injection marker
        xline(timevec(marks(1)), '--b', evnt, 'FontWeight', 'bold', 'LineWidth', 1.5);
        
        % Apply Arial and Normal Font Weight globally
        set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
        set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');
        
        hold off;
    end
end

st_hdr_pks(st_hdr_pks == 0) = NaN;
st_hdr_locs(st_hdr_locs == 0) = NaN;
st_int(st_int == 0) = NaN;
st_loss(st_loss == 0) = NaN;

st_amp = mean(st_hdr_pks, 'omitnan');
st_loc = mean(st_hdr_locs, 'omitnan');

st_raise = st_hdr_locs - st_int;
st_decay = st_loss - st_hdr_locs;
st_dur = st_loss - st_int;

disp('Completed 1st peaks dynamics and generated QC plots.');

% --- Remaining Peaks ---
rest_locs = cell(numCells,1);
rest_pks = cell(numCells, 1);
for l = 1:numCells
    if length(hdr_locs{l}) > 1
        rest_locs{l} = hdr_locs{l}(2:end);
        rest_pks{l}  = hdr_pks{l}(2:end);
    end
end
rest_cnt_osc = cellfun(@numel, rest_pks);
rest_cnt_osc(rest_cnt_osc == 0) = NaN;
rest_cnt_amp = cellfun(@mean, rest_pks);

%% 14. Intracellular Calcium Dynamics (Thapsigargin, SOCE, VOCC)
disp('--- ER Depletion (Thapsigargin) AUC ---');

% Bring the Mean dF/F0 graph to the front for interaction
figure(4); 
disp('Click TWO points on the Mean dF/F0 graph to define the ER depletion window...');
[x_er, ~] = ginput(2); 
x_er = sort(x_er); % Ensure chronological order regardless of click order

% Map clicked time (min) to exact timevec indices
[~, ER_start] = min(abs(timevec - x_er(1)));
[~, ER_end]   = min(abs(timevec - x_er(2)));

df_thap = df_ratio(:, ER_start:ER_end);
timevec_thap = timevec(ER_start:ER_end);

figure(14); clf; set(gcf, 'Position', [250, 250, 800, 400]);
plot(timevec_thap, df_thap);
title({'ER Depletion Window (Thapsigargin)', ''});
xlabel('Time (min.)'); ylabel('\DeltaF/F_0');
set(gca, 'xlim', timevec([ER_start ER_end]), 'TickDir', 'none');
box off;

auc_er = nan(numCells, 1); % Preallocate AUC array with NaNs
for i = 1:numCells
    % findpeaks returns empty [], not NaN, if no peaks are found
    [~, tg_locs] = findpeaks(df_thap(i, :), 'MinPeakProminence', 1);
    
    if ~isempty(tg_locs)
        % Find where signal crosses 0 before the first peak
        tg_int = find(df_thap(i, 1:tg_locs(1)) <= 0, 1, 'last');
        if isempty(tg_int)
            tg_int = 1; % Failsafe if it never dropped below 0
        end
        
        % Calculate true Area Under Curve using exact time (minutes)
        auc_er(i) = trapz(timevec_thap(tg_int:end), df_thap(i, tg_int:end));
    end    
end

figure(15); clf;
bar(auc_er, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
title({'Thapsigargin AUC per ROI', ''});
xlabel('ROI Number'); ylabel('AUC (\DeltaF/F_0 \times min)');
box off;

auc_er_mean = mean(auc_er, 'omitnan');
disp(['-> Mean Thapsigargin AUC: ', num2str(auc_er_mean)]);

% Apply formatting
set(findall(figure(14), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
set(findall(figure(15), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');


%% --- SOCE Slope ---
disp('--- SOCE Slope ---');
figure(4);
disp('Click TWO points on the Mean dF/F0 graph to define the SOCE slope window...');
[x_soce, ~] = ginput(2);
x_soce = sort(x_soce);

[~, SOCE_start] = min(abs(timevec - x_soce(1)));
[~, SOCE_end]   = min(abs(timevec - x_soce(2)));

df_slpe = df_ratio(:, SOCE_end) - df_ratio(:, SOCE_start);
timevec_slpe = timevec(SOCE_end) - timevec(SOCE_start);

soce_slope = df_slpe ./ timevec_slpe; % Delta Y / Delta X
soce_slope_mean = mean(soce_slope, 'omitnan');
disp(['-> Mean SOCE Slope: ', num2str(soce_slope_mean)]);


%% --- VOCC + SOCE AUC ---
disp('--- VOCC + SOCE AUC ---');
figure(4);
disp('Click TWO points on the Mean dF/F0 graph to define the VOCC+SOCE window...');
[x_voso, ~] = ginput(2);
x_voso = sort(x_voso);

[~, bg_voso] = min(abs(timevec - x_voso(1)));
[~, end_voso] = min(abs(timevec - x_voso(2)));

fprintf('-> Mapped Frame Indices: Start = %d, End = %d\n', bg_voso, end_voso);

% Baseline correction based on the 20 frames right before the window
baseline_voso_idx = max(1, bg_voso-20);
baseline_voso = mean(df_ratio(:, baseline_voso_idx:bg_voso-1), 2);
df_voso = df_ratio(:, bg_voso:end_voso) - baseline_voso;
timevec_voso = timevec(bg_voso:end_voso);

figure(16); clf; set(gcf, 'Position', [300, 300, 800, 400]);
plot(timevec_voso, df_voso);
title({'VOCC+SOCE Window (Baseline Corrected)', ''});
xlabel('Time (min.)'); ylabel('\DeltaF/F_0');
set(gca, 'xlim', timevec([bg_voso end_voso]), 'TickDir', 'none');
box off;

% Calculate true Area Under Curve using exact time (minutes)
auc_voso = trapz(timevec_voso, df_voso, 2);

figure(17); clf;
bar(auc_voso, 'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none');
title({'VOCC+SOCE AUC per ROI', ''});
xlabel('ROI Number'); ylabel('AUC (\DeltaF/F_0 \times min)');
box off;

% Apply formatting
set(findall(figure(16), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
set(findall(figure(17), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');


%% 15. Memory Cleanup and Save
disp('Purging heavy image arrays and temporary variables to save memory...');

% 1. Delete the massive 3D image arrays and UI coordinate grids
clear M2_gcam M2_red aveMap X_grid Y_grid ax fig1 fig2 fig3 fig4 roi_objs mask curr_pixels;
clear tmp_g tmp_r vals_g vals_r all_pixel_idx all_cell_ids img_desc tiff_info;

% 2. Delete loop indices, UI handles, and raw metadata string tokens
clear c cx cy rad timei i j k l roi_tokens time_str time_tokens btn_add btn_done;

% 3. Delete plotting boundary calculations
clear t_start t_end y_max y_min y_range bar_y text_y plot_ymin plot_ymax;
clear y_max_all y_min_all y_range_all bar_y_all text_y_all valid_idx;

% 4. Delete intermediate baseline and peak finding thresholds
clear baseline_start baseline_end F0_ratio smooth_signal signal_post time_post;
clear baseline_check_idx baseline_signal baseline_median robust_sigma;
clear dyn_prom dyn_height dt_min dyn_width pks_post locs_post;

% 5. Delete intermediate peak dynamics and graphing variables
clear peak_idx pre_peak_signal st_int_bline idx_rise_start post_peak_signal;
clear decay_baseline idx_decay_end st_nd_dy peak2_idx rise_w rise_h decay_w decay_h;

% 6. Delete temporary AUC selection and mapping variables
clear x_er x_soce x_voso ER_start ER_end SOCE_start SOCE_end bg_voso end_voso;
clear baseline_voso_idx baseline_voso tg_locs tg_int df_thap timevec_thap;
clear df_slpe timevec_slpe df_voso timevec_voso;

disp('Saving lightweight analysis workspace...');
[~, base_name_green, ~] = fileparts(file_name_green);

% Strip out '_GREEN' if it exists to make a clean base name
save_name = strrep(base_name_green, '_GREEN', '');
datFileName = [save_name, '_AnalyzedData.mat'];
save(fullfile(path_name, datFileName), '-v7.3');

disp(['SUCCESS! Analysis saved as: ', datFileName]);
disp('Pipeline Complete.');


%% 
close all, clear, clc