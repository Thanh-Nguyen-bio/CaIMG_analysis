
close all, clear,clc
%%
% 1. Specify the folder path
myFolder = '/media/sjmoons/Lab_RNA_seq_data/Thanh/20260209-basal-WT-DEL/20260229-basal-WT-DELta';

% 2. Get a list of all TIFF files in that folder
filePattern_340 = fullfile(myFolder, '0911-11.*'); % Handles .tif and .tiff
tifFiles_340 = dir(filePattern_340);

filePattern_380 = fullfile(myFolder, '0911-12.*'); % Handles .tif and .tiff
tifFiles_380 = dir(filePattern_380);
frame_number = length(tifFiles_340);

% 3. Loop through each file
for k = 1:frame_number
    baseFileName_340 = tifFiles_340(k).name;
    fullFileName_340 = fullfile(myFolder, baseFileName_340);
    
    baseFileName_380 = tifFiles_380(k).name;
    fullFileName_380 = fullfile(myFolder, baseFileName_380);
    
    fprintf('Now reading %s\n', baseFileName_340);
    fprintf('Now reading %s\n', baseFileName_380);
    % 4. Read the image
    image_340(:,:,k) = imread(fullFileName_340);
    image_380(:,:,k) = imread(fullFileName_380);
  
end


%%
blank_image = uint16(ones(size(image_340,1),size(image_340,2)));
bckgrnd_340 = uint16(ceil(max(image_340(500:end,1:12,:),[],'all')));
bckgrnd_340_img = uint16(bckgrnd_340 .* blank_image);
bckgrnd_380 = uint16(ceil(max(image_380(500:end,1:12,:),[],'all')));
bckgrnd_380_img = uint16(bckgrnd_380 .* blank_image);

%%
for i=1:frame_number
    image_340(:,:,i) = imsubtract(image_340(:,:,i),bckgrnd_340_img);
    image_380(:,:,i) = imsubtract(image_380(:,:,i),bckgrnd_380_img);
end
%%
image_340 = double(image_340);
image_380 = double(image_380);
ratio_340_380 = image_340 ./ image_380;
%%
fig_1 = figure(1),clf;

fig_1_til1 = tiledlayout(1,3);

fig_1_title= title(fig_1_til1,'Initializing...');
ax1= nexttile;
fig_1_340 = imagesc(image_340(:,:,1));
title("340")
set(ax1, 'clim', [0 1]*5000,'xtick',[],'ytick',[])
colormap(ax1,"gray")

axis square
colorbar

ax2=nexttile;
fig_1_380 = imagesc(image_380(:,:,1));
title("380")
colormap(ax2,"gray")
set(ax2, 'clim', [0 1]*5000,'xtick',[],'ytick',[])
axis square
colorbar

ax3=nexttile;
fig_1_ratio = imagesc(ratio_340_380(:,:,1));
title("340/380")
axis square
colormap(ax3,"jet")
set(ax3, 'clim', [0.2 1]*2,'xtick',[],'ytick',[])
colorbar

%%
vid_1=VideoWriter('WT-1b-Akh10.avi','Motion JPEG AVI');
vid_1.FrameRate = 24;
open(vid_1)

%%
nFrames = min(720, frame_number); % Ensure we don't exceed available frames

for framei = 300:nFrames
     %Check if the figure and image handle still exist before updating
    if ~ishandle(fig_1)
        warning('Figure was closed. Stopping loop.');
        break; 
    end
     % Update the title
    set(fig_1_title, 'String', sprintf('Frame: %d / %d', framei, nFrames));
    % Update the image data
    
    set(fig_1_340,'CData',image_340(:,:,framei));
    set(fig_1_380,'CData',image_380(:,:,framei));
    set(fig_1_ratio,'CData', ratio_340_380(:,:,framei));
   
    
    vid_1_frame = getframe(gcf);
    writeVideo(vid_1,vid_1_frame);
    % Force MATLAB to draw the update before pausing
    drawnow; 
    pause(0.02);
end


close(vid_1)

%% identify and remove "background noise"

% STEP 1a: compute the average map
avemap = squeeze(mean(image_380,3));

% normalize to a range of [0 1]
%avemap = avemap-min(avemap(:));
avemap = avemap./max(avemap,[],'all');


% STEP 1b: apply MATLAB's local adaptive histogram equalization
avemap = adapthisteq(avemap);


% STEP 2: estimate the background as a fuzzy version of the image
background = imgaussfilt(avemap,80);


% STEP 3: create the boosted-SNR map by subtracting the 'background'
foreground = avemap-background;

%% pause to see what we've done so far

figure(3), clf
colormap hot

% the average map (many code statements in one line!)
subplot(221), 
imagesc(avemap)
axis square 
title('Mean')



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
set(gca,'clim',[0 0.3])
colorbar


%% continuing...

% STEP 4: threshold the foreground map
threshval = 0.089;
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
cells2cut = cellsizes<9 | cellsizes>Inf;

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
time_trace_340 = zeros(islands.NumObjects,frame_number);
time_trace_380 = zeros(islands.NumObjects,frame_number);

% extract data from each cell over time
for celli=1:islands.NumObjects
    
    % done per time point because cells are 2D
    for timei=1:frame_number
        
        % get the entire map from this time point
        tmp_340 = squeeze(image_340(:,:,timei));
        tmp_380 = squeeze(image_380(:,:,timei));
        
        % compute the average of all pixels in this time point
        time_trace_340(celli,timei) = mean( tmp_340(islands.PixelIdxList{celli}) );
        time_trace_380(celli,timei) = mean( tmp_380(islands.PixelIdxList{celli}) );
    end
end


%%
timevec = (1:1:frame_number);
% Show heatmap
figure(5), clf

tiledlayout(1,3)
nexttile
plot( timevec,time_trace_340)
set(gca,'XLim', timevec([1 end]))
axis square
title('340')

nexttile
plot( timevec,time_trace_380)
set(gca,'XLim', timevec([1 end]))
axis square
title('380')

nexttile
plot( timevec,time_trace_340 ./ time_trace_380)
set(gca,'XLim', timevec([1 end]),'ylim',[0.5 1.2])
axis square
title('Ratio')
xlabel('Time (min.)')
ylabel('F340/F380')
colorbar


%%
figure(6), clf
% Show heatmap

tiledlayout(1,3)
nexttile
imagesc( timevec,[],time_trace_340)
axis square
nexttile
imagesc ( timevec,[],time_trace_380)
axis square
nexttile
imagesc( timevec,[],time_trace_340 ./ time_trace_380)
axis square

xlabel('Time (min.)')
ylabel('ROI number')
colorbar


%% show image with ROIs
figure(7),clf
img_anatomy = cat(3,avemap,avemap,avemap);
bx_1=imagesc(img_anatomy)

hold on


fig_7_roi = imagesc(threshimgFilt)
set(fig_7_roi,'AlphaData',threshimgFilt*0.1)
set(gca,'clim',[-1 1],'xtick',[],'ytick',[])
colormap parula
axis square


%% save