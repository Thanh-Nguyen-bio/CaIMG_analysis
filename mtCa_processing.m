%% clear

close all, clear, clc



%% Load and Prepare Data
%Fluo4-time-lapse
img_data = bfopen('WT-3-5mmCAF-2Ion.czi');
img_data_planes = img_data{1, 1};
% Bio-Formats returns a cell array where column 1 is the pixels
img_plane_image = img_data_planes(:,1); 
total_image = size(img_plane_image,1);
% Concatenate into 3D matrix and convert to double immediately
slide_odd= img_plane_image(1:2:end); 
slide_even = img_plane_image(2:2:end);
cfp_matObj = double(cat(3, slide_odd{:}));

fret_matObj = double(cat(3, slide_even{:}));

% Verify dimensions
npnts = 0.5*total_image;

%% Load Excel time series data
filename = ['KD68-1-Akh10ng.csv'];

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

%% make up timevec
timevec = 1:npnts;
%% reagent treatment and time

marks = [60 300];
evnt = {'CAF 5 mM','Ionomycin'};
%% Visualization Loop
fig1 = figure(1); clf;
% Initialize image with the first frame
imgh = imagesc(fret_matObj(:,:,1)); 
tith = title('Initializing...');
axis square
colormap parula
set(gca, 'clim', [0 128]) % Standardizing contrast
colorbar

nFrames = min(200, npnts); % Ensure we don't exceed available frames

for framei = 1:nFrames
    % Check if the figure and image handle still exist before updating
    %if ~ishandle(imgh)
    %    warning('Figure was closed. Stopping loop.');
    %    break; 
    %end
    
    % Update the image data
    set(imgh, 'CData', fret_matObj(:,:,framei));
    
    % Update the title
    set(tith, 'String', sprintf('Frame: %d / %d', framei, nFrames));
    
    % Force MATLAB to draw the update before pausing
    drawnow; 
    pause(0.001);
end

%% identify and remove "background noise"

% STEP 1a: compute the average map
avemap = mean (cfp_matObj,3);

%%
% 2. Background Subtraction (Top-Hat Filtering)
% Top-hat is superior for mitochondria "dots" as it removes uneven 
% background while preserving small, bright structures.
se = strel('disk', 5); % '10' is the approximate radius of your largest mito
fore_grnd = imtophat(avemap, se);

% We create a binary mask from the Donor channel to identify mitochondria
thresh = graythresh(mat2gray(fore_grnd)); % Otsu's method
mask = imbinarize(mat2gray(fore_grnd), thresh);

% Clean up the mask: remove noise (objects smaller than 10 pixels)
mask = bwareaopen(mask, 5); 

figure(3),clf
tiledlayout(1,3)
ax1 = nexttile
imagesc(ax1, avemap)
axis square
colormap(ax1,"jet")
colorbar

ax2 = nexttile
imagesc(ax2,fore_grnd)
colormap(ax2,"jet")
axis square
colorbar

ax3 = nexttile
imagesc(ax3,mask)
colormap (ax3,'gray')
axis square
colorbar

%%

% get cluster information
islands = bwconncomp(mask);

% identify the cluster sizes
mt_sizes = cellfun(@length,islands.PixelIdxList);

% find small and large cells
mts_2cut = mt_sizes<5 | mt_sizes> Inf;

% and remove those cells
islands.PixelIdxList(mts_2cut) = [];

% update the number of remaining clusters ("neurons")
islands.NumObjects = numel(islands.PixelIdxList);

% finally, recreate the threshold image without rejected clusters 
threshimgFilt = false(size(avemap));
for i=1:islands.NumObjects
    threshimgFilt(islands.PixelIdxList{i}) = true; % what value shall we assign these pixels?
end

figure(4), clf
subplot(121)
imagesc(mask)
axis square
title('binarized (original)')
colormap gray



% show again for comparison
subplot(122)
imagesc(threshimgFilt)
axis square
title('binarized (filtered)')

%% 1. Pre-allocation and Vectorization Setup (Runs Once)
% Set up parallel pool for the Intel Xeon (14 threads)
if isempty(gcp('nocreate'))
    parpool('Threads', 14);
end
%% preallocate 2D and 1D array for computing
numCells = islands.NumObjects;
fret_mitos = zeros(numCells, npnts);
cfp_mitos  = zeros(numCells, npnts);

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
parfor timei = 1:npnts
    
    % --- I/O Phase ---
    % Squeeze is usually unnecessary if indexing 2D explicitly, but kept for safety
    fret_tmp = double(fret_matObj(:,:,timei)); 
    cfp_tmp  = double(cfp_matObj(:,:,timei));
    
    % --- Computation Phase ---
    % 1. Extract only the pixels that belong to cells (ignores background)
    fret_vals = fret_tmp(all_pixel_idx);
    cfp_vals  = cfp_tmp(all_pixel_idx);
    
    % 2. Calculate the mean for ALL cells instantly using accumarray
    fret_mitos(:, timei) = accumarray(all_cell_ids, fret_vals, [numCells, 1], @mean);
    cfp_mitos(:, timei)  = accumarray(all_cell_ids, cfp_vals, [numCells, 1], @mean);
    
end

elapsedTime = toc;
disp(['Load & Process time: ' num2str(elapsedTime) ' seconds']);

%% Clean up
delete(gcp('nocreate'));

%% visualize some time courses
figure(5), clf

tiledlayout(1,2)
% show all neurons at a time
bx1 = nexttile
imagesc(timevec,[],fret_mitos)
xlabel('Time (sec.)')
ylabel('Mito number')
title('FRET channel')
axis square
colorbar

bx2 = nexttile
imagesc(timevec,[],cfp_mitos)
xlabel('Time (sec.)')
ylabel('Mito number')
title('CFP channel')
axis square
colorbar


% show all neurons at the same time
figure(6), clf
tiledlayout(2,1)
cx1= nexttile
plot(timevec,fret_mitos)
ylabel('Brightness (a.u.)')
title('FRET of all mito cluster')
set(gca,'xlim',timevec([1 end]))
xlabel('Time (sec.)')

cx2 = nexttile
plot(timevec,cfp_mitos)
ylabel('Brightness (a.u.)')
title('CFP of all mito cluster')
set(gca,'xlim',timevec([1 end]))
xlabel('Time (sec.)')


%% ratiometric time tracing
ratio_mitos = fret_mitos./ cfp_mitos;

%% Savitzky-Golay filtering
odr = 5; % define order
f_lgth = 11; %define filter length 
flp_ratios = ratio_mitos';
fil_ratios = sgolayfilt(flp_ratios,odr,f_lgth);
fil_ratios = fil_ratios';
%% visualize ratiometric tracing
figure(7), clf

imagesc(timevec,[],ratio_mitos)
xlabel('Time (sec.)')
ylabel('Mito number')
title('Ratiometric FRET/CFP')
axis square
set(gca,'clim', [-0.5 5])
colorbar

figure(71), clf
imagesc(timevec,[],fil_ratios)
xlabel('Time (sec.)')
ylabel('Mito number')
title('Ratiometric FRET/CFP (Filtered)')
axis square
set(gca,'clim', [-0.5 5])
colorbar

figure(8), clf
plot(timevec,fil_ratios)
ylabel('FRET/CFP')
title('Ratiometric of all mito cluster')
set(gca,'xlim',timevec([1 end]),'YLim',[-0.5 5])
xlabel('Time (sec.)')
xline(timevec([marks]),'--r',evnt,'fontsize',11,'FontWeight','bold')


%%
ra= randi(islands.NumObjects,1);
    figure(2000+ra), clf
    set(gcf,'Units', 'centimeters');
    set(gcf,'Position', [0 0 16 16])
    tiledlayout(3,1);

    fx1=nexttile;
    plot(timevec, fret_mitos(ra,:),'-g','LineWidth',1)
    ylabel('F acceptor')
    xlabel('Time (min.)')
    title(sprintf('Calcium content of mitochondria %d', ra))
    set(gca, 'xlim', timevec([1 end]))
    %set(gca, 'ylim', [0.15 0.7])
    xline(timevec([marks]),'--r',evnt,'fontsize',11,'FontWeight','bold')

    fx2= nexttile;
    plot(timevec, cfp_mitos(ra,:),'-c','LineWidth',1)
    ylabel('F donor')
    xlabel('Time (min.)')
    title(sprintf('Calcium content of mitochondria %d', ra))
    set(gca, 'xlim', timevec([1 end]))
    %set(gca, 'ylim', [0.15 0.7])
    xline(timevec([marks]),'--r',evnt,'fontsize',11,'FontWeight','bold')

    fx3= nexttile;
    plot(timevec, ratio_mitos(ra,:),':k','LineWidth',1.5)
    hold on
    plot(timevec, fil_ratios(ra,:),'-b','LineWidth',2)
    ylabel('F acceptor/donor')
    xlabel('Time (min.)')
    title(sprintf('Calcium content of mitochondria %d', ra))
    set(gca, 'xlim', timevec([1 end]))
    %set(gca, 'ylim', [0.15 0.7])
    xline(timevec([marks]),'--r',evnt,'fontsize',11,'FontWeight','bold')
%% (Optional) all cell time tracing
for i=1:size(ratio_mitos,1)
    figure(1000+i), clf
    plot(timevec, ratio_mitos(i,:),':k','LineWidth',1)
    hold on
    plot(timevec, fil_ratios(i,:),'-b','LineWidth',2)
    ylabel('FRET/YFP')
    xlabel('Time (min.)')
    title(sprintf('Calcium content of mitochondria %d', i))
    set(gca, 'xlim', timevec([1 end]))
    set(gca, 'ylim', [0.15 0.7])
    xline(timevec([marks]),'--r',evnt,'fontsize',11,'FontWeight','bold')

end

%% clear image data
clear img_data img_data_planes img_plane_image slide_even slide_odd fret_matObj cfp_matObj
whos
%% (Option) Save data
[~, name, ~] = fileparts(filename);
datFileName = ['DATA' name ];
save(datFileName, '-v7.3');