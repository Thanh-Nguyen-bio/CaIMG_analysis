
close all, clear, clc
%% Load and Prepare Data
%Fluo4-time-lapse
img_data = bfopen('mtGCamP5GDsRed-WT-1-Akh1ng.czi');
img_data_planes = img_data{1, 1};
% Bio-Formats returns a cell array where column 1 is the pixels
img_plane_image = img_data_planes(:,1); 
total_image = size(img_plane_image,1);
% Concatenate into 3D matrix and convert to double immediately
slide_odd= img_plane_image(1:2:end); 
slide_even = img_plane_image(2:2:end);
gcamp_img = double(cat(3, slide_odd{:}));

tdtom_img = double(cat(3, slide_even{:}));

% Verify dimensions
npnts = 0.5*total_image;


%% Load Excel time series data
filename = ['mtGCamP5GDsRed-WT-1-Akh1ng.csv'];

% Read the table
dataTbl = readtable(filename);

% Display column names
disp('Available column names:');
disp(dataTbl.Properties.VariableNames')

% Identify the time column
timeColIdx = find(contains(dataTbl.Properties.VariableNames, 'Time', 'IgnoreCase', true), 1);

if isempty(timeColIdx)
    error('No time column found!');
end

% Extract time vector
timevec = dataTbl{:, timeColIdx};
timevec = rmmissing(timevec);
timevec = timevec * 60;
npnts = length(timevec);

% Estimate sampling rate

marks= find(abs(diff(round(diff(timevec))))>=1.5) +2;

timevec = dataTbl{:, timeColIdx};
timevec = rmmissing(timevec);

%% define events
evnt= {'Akh 1 ng'};

%% Visualization Loop
depth_green = 4096;
green_clmp = [zeros(depth_green,1), linspace(0,1,depth_green)',zeros(depth_green,1) ];
fig1= figure(1); clf;
% Initialize image with the first frame
tl1 = tiledlayout(1,2)
tith = title(tl1,'Initializing...')
ax1= nexttile;
imgh1 = imagesc(gcamp_img(:,:,1)); 
title('GCamP5G');
axis square
colormap(ax1,green_clmp)
set(gca, 'clim', [0 1024]) % Standardizing contrast
colorbar

ax2= nexttile;
imgh2 = imagesc(tdtom_img(:,:,1)); 
title('TdTomato');
axis square
colormap(ax2,'parula')
set(gca,'clim',[0 4096]) % Standardizing contrast
colorbar

%%
nFrames = min(120, npnts); % Ensure we don't exceed available frames

for framei = 60:nFrames
     %Check if the figure and image handle still exist before updating
    if ~ishandle(imgh1)
        warning('Figure was closed. Stopping loop.');
        break; 
    end
    
    % Update the image data
    
    set(imgh1, 'CData', gcamp_img(:,:,framei));
    set(imgh2, 'CData', tdtom_img(:,:,framei));
    % Update the title
    set(tith, 'String', sprintf('Frame: %d / %d', framei, nFrames));
    
    % Force MATLAB to draw the update before pausing
    drawnow; 
    pause(0.002);
end


%%
rati_img = gcamp_img ./ tdtom_img;

%% Visualization Loop
fig2 = figure(2); clf;
% Initialize image with the first frame
imgh_rat = imagesc(rati_img(:,:,1)); 
tith2 = title('Initializing...');
axis square
set(gca,'xtick', [], 'YTick', [])
colormap jet
set(gca, 'clim', [0 1]) % Standardizing contrast
colorbar
%%
nFrames = min(1500, npnts); % Ensure we don't exceed available frames

for framei = 1:nFrames
    % Check if the figure and image handle still exist before updating
    %if ~ishandle(imgh)
    %    warning('Figure was closed. Stopping loop.');
    %    break; 
    %end
    
    % Update the image data
    set(imgh_rat, 'CData', rati_img(:,:,framei));
    
    % Update the title
    set(tith2, 'String', sprintf('Frame: %d / %d', framei, nFrames));
    
    % Force MATLAB to draw the update before pausing
    drawnow; 
    pause(0.002);
end

%% identify and remove "background noise"

% STEP 1a: compute the average map
avemap = squeeze(mean(tdtom_img,3));

% normalize to a range of [0 1]
avemap = avemap-min(avemap(:));
avemap = avemap./max(avemap,[],'all');


% STEP 1b: apply MATLAB's local adaptive histogram equalization
avemap = adapthisteq(avemap);


% STEP 2: estimate the background as a fuzzy version of the image
background = imgaussfilt(avemap,5);


% STEP 3: create the boosted-SNR map by subtracting the 'background'
foreground = avemap-background;

%% pause to see what we've done so far

figure(3), clf
colormap hot

% the average map (many code statements in one line!)
subplot(221), imagesc(avemap), axis square, title('Mean')


% the background image
subplot(222)
imagesc(background)
axis square
title('Background')

% the foreground image
subplot(223)
imagesc(foreground)
axis square
title('Isolated')
set(gca,'clim',[0 .3])
colorbar

%% continuing...

% STEP 4: threshold the foreground map
threshval = .05;
threshimg = foreground > threshval;


% and visualize that
subplot(224)
imagesc(threshimg)
axis square
title('binarized')

%% clustering
% get cluster information
islands = bwconncomp(threshimg);

% identify the cluster sizes
cellsizes = cellfun(@length,islands.PixelIdxList);

% find small and large cells
cells2cut = cellsizes<10 | cellsizes>Inf;

% remove those cells
islands.PixelIdxList(cells2cut) = [];

% update the number of remaining clusters ("neurons")
islands.NumObjects = numel(islands.PixelIdxList);


% finally, recreate the threshold image without rejected clusters 
threshimgFilt = false(size(avemap));
for i=1:islands.NumObjects
    threshimgFilt(islands.PixelIdxList{i}) = true;
end
nROIs = islands.NumObjects;

%% visualize

% same as in previous video, redrawn for the before/after show
figure(4), clf
subplot(121)
imagesc(threshimg)
axis square
title('binarized (original)')
colormap gray


% show again for comparison
subplot(122)
imagesc(threshimgFilt)
axis square
title('binarized (filtered)')

%% get time courses from all "cells"

% initialize time series matrix
time_trace_gcamp = zeros(islands.NumObjects,npnts);
time_trace_tdtom = zeros(islands.NumObjects,npnts);
% extract data from each cell over time
for celli=1:islands.NumObjects
    
    % done per time point because cells are 2D
    for timei=1:npnts
        
        % get the entire map from this time point
        tmp_gcamp = squeeze(gcamp_img(:,:,timei));
        tmp_tdtom = squeeze(tdtom_img(:,:,timei));
        
        % compute the average of all pixels in this time point
        time_trace_gcamp(celli,timei) = mean( tmp_gcamp(islands.PixelIdxList{celli}) );
        time_trace_tdtom(celli,timei) = mean( tmp_tdtom(islands.PixelIdxList{celli}) );
    end
end

%% Visualization

% Plot an example trace
figure(4), clf
%rprest = randi(nROIs,1);
j= 5 ;
plot(timevec, time_trace_gcamp(j,:))
ylabel('Brightness (a.u.)')
set(gca, 'xlim', timevec([1 end]), ...
    'ylim', [0 1000])
title(sprintf('Fluorescence of ROIs %d (normalized)', j))
box off

% Show heatmap
figure(5), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),time_trace_gcamp,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca, 'xlim', timevec([1 end]))
view(2)
xlabel('Time (min.)')
ylabel('ROI number')
colorbar

% Plot all traces
figure(6), clf
plot(timevec, time_trace_gcamp)
ylabel('Brightness (a.u.)')
title('Fluorescence of all ROIs')
xlabel('Time (min.)')
set(gca, 'xlim', timevec([1 end]), ...
    'YLim',[0 1500])


%% Visualization Tomato


% Show heatmap
figure(7), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),time_trace_tdtom,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca, 'xlim', timevec([1 end]))
view(2)
xlabel('Time (min.)')
ylabel('ROI number')
colorbar

% Plot all traces
figure(8), clf
plot(timevec, time_trace_tdtom)
ylabel('Brightness (a.u.)')
title('Tomato fluorescence of all ROIs')
xlabel('Time (min.)')
set(gca, 'xlim', timevec([1 end]))

%%
rati_ROIs = time_trace_gcamp ./time_trace_tdtom;

%%
figure(9), clf
title('Normalized to Tomato')
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),rati_ROIs,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca,'xlim', timevec([1 end]))
xlabel('Time (min.)')
ylabel('ROI number')
colorbar
view(2)
figure(10), clf
plot(timevec, rati_ROIs)
ylabel('Normalized to tdTomato')
xlabel('Time (min.)')
title('Fluorescence of all ROIs (normalized)')
set(gca, 'xlim',timevec([1 end]))

%% Convert to dF/F

df_ROIs = bsxfun(@rdivide, rati_ROIs, mean(rati_ROIs(:,marks(1)-20:marks(1)-1), 2))-1;

%% dF/F visualization
figure(11), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),df_ROIs,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5], ...
    'ytick',[1 nROIs],...
        'xlim', timevec([1 end]), ...
    'clim',[-1 10], ...
    'YDir','reverse', ...
    'fontsize',20, ...
    'fontname','Arial Narrow')
xline(timevec(marks),'LineWidth',1.5,'Color',[1 1 1],'LineStyle','--')
text(timevec(marks(1)),-2.4,evnt,"HorizontalAlignment","left",'FontSize',20)
colorbar
axis square
box off
grid off
view(2)

figure(12), clf
set(gcf,'Units', 'inches')
set(gcf,'Position', [1 1 15 5])
plot(timevec, df_ROIs);
ylabel('\DeltaF/F_0')
xlabel('Time (min.)')
title('Fluorescence of all ROIs (normalized)', ...
    'FontWeight','normal', ...
    'fontsize',22, ...
    'FontName','Arial Narrow')
set(gca,'xlim', timevec([1 end]), ...
    'ylim', [-1.5 12], ...
    'YTick',(0:2:10), ...
    'TickDir','none', ...
    'Fontsize',20, ...
    'Fontname','Arial Narrow')

box off
%drawing 1 lines to indicate drugs treatment on top of graph
line([timevec(marks(1)) timevec(end)],[10.5 10.5],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),11.5,evnt, ...
    'HorizontalAlignment','center','FontSize',20)