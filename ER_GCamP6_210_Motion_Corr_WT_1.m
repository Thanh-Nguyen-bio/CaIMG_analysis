close all; clear; clc;


%% Parallel Setup
if isempty(gcp('nocreate'))
    parpool('Threads', 14);
end


%% Load Excel time series data
filename = ['ER-GCamP6210-KD68-3.csv'];

% Read the table
dataTbl = readtable(filename);


% Identify the time column
timeColIdx = find(contains(dataTbl.Properties.VariableNames, 'Time', 'IgnoreCase', true), 1);

if isempty(timeColIdx)
    error('No time column found!');
end

% Extract time vector
timevec = dataTbl{:, timeColIdx};
timevec = rmmissing(timevec);
timevec = timevec * 60;
%npnts = length(timevec);

% Estimate sampling rate

marks= find(abs(diff(round(diff(timevec))))>=1.5) +2;

timevec = dataTbl{:, timeColIdx};
timevec = rmmissing(timevec);
%% define events
evnt= {'Akh 10 ng'};

%% Read tiff file
% Open UI dialog to select one or multiple files
[file_list, path_name] = uigetfile({'*.tif;*.tiff', 'Image Files (, *.tif, *.tiff)'}, ...
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

current_file = fullfile(path_name, file_list);

M2 = loadtiff(string(current_file),1,510);
T = size(M2,ndims(M2));
%% Identify and Remove "Background Noise"
disp('Processing background subtraction and de-noising...');
% compute the average map from the non-rigid corrected stack
aveMap = squeeze(mean(M2, 3));

% Normalize the initial average map to a [0 1] range
aveMap = aveMap - min(aveMap(:));
aveMap = aveMap ./ max(aveMap, [], 'all');

% STEP 1: Estimate the background using morphological opening
% NOTE: The disk radius (30) should be slightly larger than your largest cell. 
se = strel('disk', 30); 
backgroundMap = imopen(aveMap, se);

% STEP 2: Subtract the background to isolate the foreground
foregroundMap = aveMap - backgroundMap;

% STEP 3: De-noise the resulting image
% Using a 3x3 median filter to remove speckle noise while preserving ROI edges
%foregroundMap = medfilt2(foregroundMap, [3 3]);

% Normalize the final foreground map to [0 1]
foregroundMap = foregroundMap - min(foregroundMap(:));
foregroundMap = foregroundMap ./ max(foregroundMap, [], 'all');

%% Draw Multiple Round ROIs and Calculate Means Through Stack
figure(1); clf;
imagesc(foregroundMap,[-0.1 0.5]); colormap hot; axis image;

title('Draw Circular ROIs');
hold on;

keepDrawing = true;
roiCount = 0;
roiMasks = {}; % Cell array to store masks

while keepDrawing
    roiCount = roiCount + 1;
    
    % Prompt user to draw a circle
    title(sprintf('Draw ROI #%d (Click and drag to draw)', roiCount));
    hROI = drawcircle('Color', 'r');
    
    % Pause to allow manual adjustments
    input('Adjust the ROI, then press ENTER in the command window to confirm...');
    
    % Add the mask to our cell array
    roiMasks{roiCount} = createMask(hROI);
    
    % Add a text label to the figure so you know which ROI is which
    text(hROI.Center(1), hROI.Center(2), num2str(roiCount), ...
        'Color', 'y', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % Ask if the user wants to draw another one
    reply = questdlg('Would you like to draw another ROI?', ...
        'Continue Drawing?', 'Yes', 'No', 'Yes');
    
    if strcmp(reply, 'No')
        keepDrawing = false;
    end
end

title(sprintf('Finished drawing %d ROIs.', roiCount));

% Pre-allocate matrix for the mean traces (Frames x ROIs)
roiMeanTraces = zeros(T, roiCount);

disp('Calculating mean intensities for all ROIs...');
% Calculate the mean of each ROI for each frame
parfor i = 1:roiCount
    currentMask = roiMasks{i};
    for t = 1:T
        currentFrame = M2(:, :, t);
        roiMeanTraces(t, i) = mean(currentFrame(currentMask));
    end
end

roiMeanTraces =roiMeanTraces';

%% Clean up threads
delete(gcp('nocreate'));
%% Visualization

% Show heatmap
figure(5), clf
surf(timevec,linspace(0.5,roiCount+0.5,roiCount),roiMeanTraces,'LineStyle','none')
set(gca,'ylim',[0.5 roiCount+0.5])
set(gca, 'xlim', timevec([1 end]))
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
set(gca, 'xlim', timevec([1 end]))

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
parts = split(filename, '.');
file_name = parts{1};
save(file_name,'-v7.3');