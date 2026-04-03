close all; clear; clc;


%% Parallel Setup
if isempty(gcp('nocreate'))
    parpool('Threads', 14);
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
T = size(M2,ndims(M2));

%% Extract Time Vector from TIFF Metadata
disp('Extracting time vector from ImageDescription metadata...');

% Read the TIFF header information from the first frame
tiff_info = imfinfo(current_file);

% Check if the ImageDescription tag exists in the file header
if isfield(tiff_info, 'ImageDescription') && ~isempty(tiff_info(1).ImageDescription)
    img_desc = tiff_info(1).ImageDescription;
    
    % Use regular expression to find the string inside TimeVector_Seconds=[...]
    tokens = regexp(img_desc, 'TimeVector_Seconds=\[(.*?)\]', 'tokens');
    
    if ~isempty(tokens)
        % Extract the comma-separated string and convert to a numeric array
        time_str = tokens{1}{1};
        timevec = str2num(time_str); % str2num automatically parses comma-separated strings
        
        % Ensure it is a column vector
        timevec = timevec(:); 
        
        disp(['Successfully extracted timevec: ', num2str(length(timevec)), ' frames spanning ', num2str(max(timevec)), ' seconds.']);
    else
        warning('TimeVector_Seconds format not found within the ImageDescription tag.');
        timevec = [];
    end
else
    warning('No ImageDescription metadata found in this TIFF file.');
    timevec = [];
end

% Estimate sampling rate

marks= find(abs(diff(round(diff(timevec))))>=1.5) +2;

timevec = timevec /60;
%% define events
evnt= {'20 Thapsigargin'};
%% Identify and Remove "Background Noise"
disp('Processing background subtraction and de-noising...');
% compute the average map from the non-rigid corrected stack
aveMap = squeeze(mean(M2, 3));

% Normalize the initial average map to a [0 1] range
aveMap = aveMap - min(aveMap(:));
aveMap = aveMap ./ max(aveMap, [], 'all');

% STEP 1: Estimate the background using morphological opening
% NOTE: The disk radius (30) should be slightly larger than your largest cell. 
se = strel('disk', 1); 
backgroundMap = imopen(aveMap, se);

% STEP 2: Subtract the background to isolate the foreground
foregroundMap = aveMap - backgroundMap;

% STEP 3: De-noise the resulting image
% Using a 3x3 median filter to remove speckle noise while preserving ROI edges
%foregroundMap = medfilt2(foregroundMap, [3 3]);

% Normalize the final foreground map to [0 1]
foregroundMap = foregroundMap - min(foregroundMap(:));
foregroundMap = foregroundMap ./ max(foregroundMap, [], 'all');
% STEP 4: threshold the foreground map
threshval = .05;
threshimg = foregroundMap > threshval;

%% pause to see what we've done so far

figure(3), clf
colormap hot

% the average map (many code statements in one line!)
subplot(221), imagesc(aveMap), axis square, title('Mean')

% the foreground image
subplot(222)
imagesc(foregroundMap)
axis square
title('Isolated')
set(gca,'clim',[0 .3])
colorbar

% and visualize that
subplot(223)
imagesc(threshimg)
axis square
title('binarized')

%% clustering
% get cluster information
islands = bwconncomp(threshimg);

% identify the cluster sizes
cellsizes = cellfun(@length,islands.PixelIdxList);

% find small and large cells
cells2cut = cellsizes<150| cellsizes>2000;

% remove those cells
islands.PixelIdxList(cells2cut) = [];

% update the number of remaining clusters ("neurons")
islands.NumObjects = numel(islands.PixelIdxList);


% finally, recreate the threshold image without rejected clusters 
threshimgFilt = false(size(aveMap));
for i=1:islands.NumObjects
    threshimgFilt(islands.PixelIdxList{i}) = true;
end
roiCount = islands.NumObjects;

% and visualize that
subplot(224)
imagesc(threshimgFilt)
axis square
title('binarized')
%% preallocate 2D and 1D array for computing
numCells = islands.NumObjects;
 roiMeanTraces = zeros(numCells, T);

% Flatten the PixelIdxList into 1D arrays for accumarray
% This entirely removes the need for an inner "cell" loop
all_pixel_idx = [];
all_cell_ids = [];

for c = 1:numCells
    curr_pixels = islands.PixelIdxList{c};
    % Ensure column vectors
    all_pixel_idx = [all_pixel_idx; curr_pixels(:)]; 
    all_cell_ids = [all_cell_ids; repmat(c, length(curr_pixels), 1)];
end

%% 2. Optimized Parallel Time Loop
disp('Starting parallel processing across time points...');
tic;

% Process timepoints in parallel. 
% Each worker gets a time point, reads the frame ONCE, and processes all cells.
parfor timei = 1:T
    
    % --- I/O Phase ---
    % Squeeze is usually unnecessary if indexing 2D explicitly, but kept for safety
    M2_tmp = double(M2(:,:,timei)); 
    
    % --- Computation Phase ---
    % 1. Extract only the pixels that belong to cells (ignores background)
    M2_vals = M2_tmp(all_pixel_idx);
    
    % 2. Calculate the mean for ALL cells instantly using accumarray
    roiMeanTraces(:, timei) = accumarray(all_cell_ids, M2_vals, [numCells, 1], @mean);
    
end

elapsedTime = toc;
disp(['Load & Process time: ' num2str(elapsedTime) ' seconds']);

%% Visualization

% Show heatmap
figure(5), clf
surf(timevec,linspace(0.5,roiCount+0.5,roiCount),roiMeanTraces,'LineStyle','none')
set(gca,'ylim',[0.5 roiCount+0.5], ...
    'xlim',timevec([1 end]))
view(2)
xlabel('Time (min.)')
ylabel('ROI number')
colorbar

% Plot all traces
figure(6), clf
plot(timevec, roiMeanTraces)
ylabel('Brightness (a.u.)')
title('Fluorescence of all ROIs')
xlabel('Time (min.)')
set(gca, 'xlim', timevec([1 end]), ...
    'YLim',[0 256])

%% Convert to dF/F
disp('Standardizing to dF');
df_ROIs = bsxfun(@rdivide, roiMeanTraces, mean(roiMeanTraces(:,marks(1)-11:marks(1)-1), 2))-1;

%% dF/F visualization
figure(11), clf
surf(timevec,linspace(0.5,roiCount+0.5,roiCount),df_ROIs,'LineStyle','none')
set(gca,'ylim',[0.5 roiCount+0.5], ...
    'ytick',[1 roiCount],...
        'xlim', timevec([1 end]), ...
    'clim',[-0.75 0.75], ...
    'YDir','reverse', ...
    'fontsize',20, ...
    'fontname','Arial Narrow')
xline(timevec(marks),'LineWidth',1.5,'Color',[1 1 1],'LineStyle','--')
text(timevec(marks(1)),-0.2,evnt,"HorizontalAlignment","left",'FontSize',20)
colorbar
axis square
box off
grid off
view(2)
%%
figure(12), clf
set(gcf,'Units', 'inches')
set(gcf,'Position', [2 2 10 7])
plot(timevec, df_ROIs);
ylabel('\DeltaF/F_0')
xlabel('Time (min.)')
title('Fluorescence of all ROIs (normalized)', ...
    'FontWeight','normal', ...
    'fontsize',22, ...
    'FontName','Arial Narrow')
set(gca,'xlim', timevec([1 end]), ...
    'ylim', [-0.75 1], ...
    'YTick',(-1:0.25:0.5), ...
    'TickDir','none', ...
    'Fontsize',20, ...
    'Fontname','Arial Narrow')

box off
%drawing 1 lines to indicate drugs treatment on top of graph
line([timevec(marks(1)) timevec(300)],[0.75 0.75],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(timevec(300)- timevec(marks(1))),0.85,evnt, ...
    'HorizontalAlignment','center','FontSize',20)

%% Clear heavy image data
clear('M1','M2','Y','img_data','img_plane_image','img_data_planes','shifts1','shifts2');

%% 6.save data
% Split the string at the first dot and take the first part
parts = split(file_list{1}, '_');
file_name = parts{1};
save(file_name,'-v7.3');

%% Clean up threads
delete(gcp('nocreate'));