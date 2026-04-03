close all; clear; clc;

%% Parallel Setup
if isempty(gcp('nocreate')) 
    parpool("Threads",14)
elseif contains(class(gcp('nocreate')), 'parallel.ProcessPool')
    % Clean up threads
    delete(gcp('nocreate'));
    parpool("Threads",14)
end
%% Read tiff file
% Open UI dialog to select one or multiple files
[file_list, path_name] = uigetfile({'*.tif;*.tiff', 'Image Files (*.tif, *.tiff)'}, ...
    'Select Corrected Files', 'MultiSelect', 'off');

% Check if the user clicked "Cancel"
if isequal(file_list, 0)
    disp('User selected Cancel. Exiting...');
    return;
end

% If only one file is selected, uigetfile returns a string. Convert to cell array for looping.
if ischar(file_list)
    file_list = {file_list};
end

% Extract the first file path (added {1} to properly index the cell array)
current_file = fullfile(path_name, file_list{1});

disp(['Loading: ', file_list{1}]);
M2 = loadtiff(string(current_file));
[img_h, img_w, T] = size(M2);

%% Extract Time Vector and ROI Metadata from TIFF
disp('Extracting metadata from ImageDescription...');

% Read the TIFF header information from the first frame
tiff_info = imfinfo(current_file);

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

%% define events
evnt= {'Akh 10 ng'};

%% Build Spatial Masks for ROIs from Metadata (Interactive & Save)
disp('Building interactive spatial masks from embedded ROI coordinates...');

% To visualize the ROIs on an average projection
aveMap = squeeze(mean(M2, 3));
fig3 = figure(3); clf;

% Create axes explicitly to attach UI elements safely
ax = axes('Parent', fig3);
imagesc(ax, aveMap); colormap hot; axis square; 
title('Adjust ROIs: Drag to move, drag edge to resize, click & press "Delete" to remove.');
hold(ax, 'on');

% 1. Pre-draw imported ROIs using interactive drawcircle objects
for c = 1:numCells
    % Parse coordinates from regex tokens
    cx  = str2double(roi_tokens{c}{1});
    cy  = str2double(roi_tokens{c}{2});
    rad = str2double(roi_tokens{c}{3});
    
    % Draw interactive ROI
    drawcircle(ax, 'Center', [cx, cy], 'Radius', rad, 'Label', num2str(c), ...
               'Color', 'w', 'FaceAlpha', 0.1, 'LabelVisible', 'on');
end

% 2. Add UI Buttons for user interaction
btn_add = uicontrol('Parent', fig3, 'Style', 'pushbutton', 'String', 'Add New ROI', ...
                    'Position', [20 20 100 30], ...
                    'Callback', @(src, event) drawcircle(ax, 'Color', 'g', 'FaceAlpha', 0.1, 'Label', 'New'));
                    
btn_done = uicontrol('Parent', fig3, 'Style', 'pushbutton', 'String', 'Done Adjusting', ...
                     'Position', [130 20 120 30], ...
                     'Callback', @(src, event) uiresume(fig3));

% 3. Pause script execution until the user clicks "Done Adjusting"
disp('Waiting for user to adjust ROIs. Click "Done Adjusting" on the figure when finished...');
uiwait(fig3); 

% 4. --- Post-Interaction Processing ---
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
[X_grid, Y_grid] = meshgrid(1:img_w, 1:img_h);

% Pre-allocate a matrix to save the final coordinates
final_roi_coords = zeros(numCells, 3); 

% 5. Build final masks from the updated objects and record coordinates
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

% 6. Save the final coordinates to a .mat file
[~, base_name_noext, ~] = fileparts(file_list{1});
roi_save_path = fullfile(path_name, [base_name_noext, '_Adjusted_ROIs.mat']);
save(roi_save_path, 'final_roi_coords', 'numCells');
disp(['Saved adjusted ROI coordinates to: ', roi_save_path]);

%% preallocate 2D array for computing traces
roiMeanTraces = zeros(numCells, T);

%% 2. Optimized Parallel Time Loop
disp('Starting parallel processing across time points...');
tic;

% Process timepoints in parallel. 
parfor timei = 1:T
    % --- I/O Phase ---
    M2_tmp = double(M2(:,:,timei)); 
    
    % --- Computation Phase ---
    % Extract only the pixels that belong to metadata cells
    M2_vals = M2_tmp(all_pixel_idx);
    
    % Calculate the mean for ALL cells instantly using accumarray
    roiMeanTraces(:, timei) = accumarray(all_cell_ids, M2_vals, [numCells, 1], @mean);
end

elapsedTime = toc;
disp(['Load & Process time: ' num2str(elapsedTime) ' seconds']);

%% Visualization

% Show heatmap
figure(5), clf
set(gcf,'Unit','pixel','position',[100 100 1400 400])
subplot(131)
surf(timevec,linspace(0.5,roiCount+0.5,roiCount),roiMeanTraces,'LineStyle','none')
set(gca,'ylim',[0.5 roiCount+0.5], ...
    'xlim',timevec([1 end]))
view(2)
xlabel('Time (min.)')
ylabel('ROI number')
colorbar
axis square

subplot(132)
% Plot all traces

plot(timevec, roiMeanTraces)
ylabel('Brightness (a.u.)')
title('Fluorescence of all ROIs','fontname','Arial','FontWeight','normal')
xlabel('Time (min.)')
set(gca, 'xlim', timevec([1 end]))
axis square

%Mean Raw Trace with SEM ---
subplot(133)
mean_raw = mean(roiMeanTraces, 1, 'omitnan')';
sem_raw = std(roiMeanTraces, 0, 1, 'omitnan')' / sqrt(roiCount);

hold on;
fill([timevec; flipud(timevec)], ...
     [mean_raw+sem_raw; flipud(mean_raw-sem_raw)], ...
     [0.8 0.95 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light green shading
plot(timevec, mean_raw, 'Color', [0 0.6 0], 'LineWidth', 2);
ylabel('Mean Brightness (a.u.)')
xlabel('Time (min.)')
title('Mean Raw Trace with SEM','fontname','Arial','FontWeight','normal')
set(gca, 'xlim', timevec([1 end]), 'TickDir', 'none')
axis square
box off

%% Convert to dF/F
disp('Standardizing to dF');
% Define baseline frames (before injection)
baseline_start = max(1, marks(1)-11);
baseline_end = max(2, marks(1)-1);
F0 = mean(roiMeanTraces(:, baseline_start:baseline_end), 2);

% Protect the baseline against division by near-zero 
F0(F0 < 0.01) = 0.01; 
df_ROIs = bsxfun(@rdivide, roiMeanTraces, F0) - 1;

% Clean extreme artifacts to protect the Mean/SEM graph
df_ROIs(isnan(df_ROIs) | isinf(df_ROIs)) = 0; 
yaxis_max = 1.05* max(df_ROIs,[],'all');
%% dF/F visualization
figure(6), clf
set(gcf,'Unit','pixel','position',[100 100 1400 400])
subplot(131)
surf(timevec,linspace(0.5,roiCount+0.5,roiCount),df_ROIs,'LineStyle','none')
set(gca,'ylim',[0.5 roiCount+0.5], ...
    'ytick',[1 roiCount],...
        'xlim', timevec([1 end]), ...
    'clim',[-0.75 0.75], ...
    'YDir','reverse', ...
    'fontsize',11, ...
    'fontname','Arial Narrow')
xline(timevec(marks),'LineWidth',1.5,'Color',[1 1 1],'LineStyle','--')
text(timevec(marks(1)),-0.2,evnt,"HorizontalAlignment","left",'FontSize',11)
colorbar
axis square
box off
grid off
view(2)

subplot(132)
plot(timevec, df_ROIs);
ylabel('\DeltaF/F_0')
xlabel('Time (min.)')
title('dF of all ROIs', ...
    'FontWeight','normal', ...
    'fontsize',11, ...
    'FontName','Arial Narrow')
set(gca,'xlim', timevec([1 end]), ...
    'ylim', [-1 1.1*yaxis_max], ...
    'YTick',(-1:0.25:yaxis_max), ...
    'TickDir','none', ...
    'Fontsize',11, ...
    'Fontname','Arial Narrow')
axis square
box off

%drawing 1 lines to indicate drugs treatment on top of graph
t_end = timevec(min(300, T));
line([timevec(marks(1)) t_end],[1.02*yaxis_max 1.01*yaxis_max],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(t_end - timevec(marks(1))),1.1*yaxis_max,evnt, ...
    'HorizontalAlignment','center','FontSize',11)

% Mean dF/F0 Trace with SEM ---
subplot(133)
mean_df = mean(df_ROIs, 1, 'omitnan')';
sem_df = std(df_ROIs, 0, 1, 'omitnan')' / sqrt(roiCount);

hold on;
fill([timevec; flipud(timevec)], ...
     [mean_df+sem_df; flipud(mean_df-sem_df)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light red shading
plot(timevec, mean_df, 'r', 'LineWidth', 2);

ylabel('Mean \DeltaF/F_0')
xlabel('Time (min.)')
title('Mean \DeltaF/F_0 with SEM', ...
    'FontWeight','normal', ...
    'fontsize',11, ...
    'FontName','Arial Narrow')

% Dynamic Y-Axis limits for Drug Bar (Safe from NaNs and flatlines)
valid_idx = ~isnan(mean_df) & ~isnan(sem_df);
if any(valid_idx)
    y_max = max(mean_df(valid_idx) + sem_df(valid_idx));
    y_min = min(mean_df(valid_idx) - sem_df(valid_idx));
else
    y_max = 1; y_min = -0.5;
end
y_range = max(y_max - y_min, 1);

bar_y = y_max + 0.15 * y_range;
text_y = y_max + 0.25 * y_range;
plot_ymin = y_min - 0.1 * y_range;
plot_ymax = text_y + 0.1 * y_range;

set(gca,'xlim', timevec([1 end]), ...
    'ylim', [plot_ymin, plot_ymax], ...
    'TickDir','none', ...
    'Fontsize',11, ...
    'Fontname','Arial Narrow')
axis square
box off

% Draw vertical injection marker
t_start = timevec(marks(1));

% Draw thick horizontal treatment bar and centered text
line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 3);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11);

%% Clear heavy image data
clear('M1','M2','Y','img_data','img_plane_image','img_data_planes','shifts1','shifts2', 'M2_tmp', 'X_grid', 'Y_grid');

%% 6.save data
% Split the string at the first dot and take the first part
parts = split(file_list{1}, '_');
file_name = parts{1};
save(file_name,'-v7.3');

%% Clean up threads
delete(gcp('nocreate'));