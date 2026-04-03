close all; clear; clc;

%% Parallel Setup
if isempty(gcp('nocreate'))
    parpool('local', 14);
end


%% Select Multiple Files
% Open UI dialog to select one or multiple files
[file_list, path_name] = uigetfile({'*.czi;*.tif;*.tiff', 'Image Files (*.czi, *.tif, *.tiff)'}, ...
    'Select Image Files', 'MultiSelect', 'on');

% Check if the user clicked "Cancel"
if isequal(file_list, 0)
    disp('User selected Cancel. Exiting...');
    return;
end

% If only one file is selected, uigetfile returns a string. Convert to cell array for looping.
if ischar(file_list)
    file_list = {file_list};
end

%% Loop Through Each Selected File
for f_idx = 1:length(file_list)
    % Get current file name and base name for saving outputs
    current_file = fullfile(path_name, file_list{f_idx});
    [~, base_name, ~] = fileparts(file_list{f_idx});
    
    fprintf('\n=== Processing File %d of %d: %s ===\n', f_idx, length(file_list), base_name);
    
    %% Load and Prepare Data
    tic;
    disp('Loading image data...');
    img_data = bfopen(current_file);
    img_data_planes = img_data{1, 1};
    img_data_signal = img_data_planes(:,1);

    load_toc = toc; disp(['Loading image data time is:' num2str(load_toc) 'secs']);
    % Extract base channels (Do not ratio them yet)
    Y_red = single(cat(3,img_data_signal{1:2:end})); 
    Y_gcam = single(cat(3,img_data_signal{2:2:end})); 
    T = size(Y_red, ndims(Y_red));
    
    % Create a raw ratio image purely for the downstream metrics plots
    Y_red_safe = Y_red; Y_red_safe(Y_red_safe == 0) = eps;
    Y_ratio = Y_gcam ./ Y_red_safe;
    Y_ratio(isnan(Y_ratio) | isinf(Y_ratio)) = 0;
    Y_ratio(Y_ratio > quantile(Y_ratio(:), 0.999)) = quantile(Y_ratio(:), 0.999);
    
    %% set parameters (first try out rigid motion correction)
    disp('Running Rigid Motion Correction on RED channel...');
    options_rigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
        'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
    
    %% perform motion correction ON RED CHANNEL
    tic; [M1_red, shifts1, template1, options_rigid] = normcorre(Y_red, options_rigid); 
    rigid_toc = toc; disp(['Rigid correction time is:' num2str(rigid_toc)]);
    
    % Apply rigid shifts to Green channel and create Rigid Ratio (for metrics)
    M1_gcam = apply_shifts(Y_gcam, shifts1, options_rigid);
    M1_red_safe = M1_red; M1_red_safe(M1_red_safe == 0) = eps;
    M1_ratio = M1_gcam ./ M1_red_safe;
    M1_ratio(isnan(M1_ratio) | isinf(M1_ratio)) = 0;
    M1_ratio(M1_ratio > quantile(M1_ratio(:), 0.999)) = quantile(M1_ratio(:), 0.999);
    
    %% now try non-rigid motion correction
    disp('Running Non-Rigid Motion Correction on RED channel...');
    options_nonrigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
        'grid_size',[32,32], 'mot_uf',4, 'max_shift',15, 'max_dev',3, 'us_fac',50, ...
        'bin_width', 50, 'init_batch', 100, 'mem_batch_size', 200);
        
    % Calculate non-rigid shifts on Red
    tic; [M2_rough_red, shifts2, template2, options_nonrigid] = normcorre_batch(Y_red, options_nonrigid); 
    non_toc = toc; disp(['Non-rigid correction time is:',num2str(non_toc)]);
    
    % Apply identical non-rigid shifts to Green
    disp('Applying calculated shifts to GCAMP channel...');
    M2_rough_gcam = apply_shifts(Y_gcam, shifts2, options_nonrigid);
    
    %% Denoise M2_rough arrays
    disp('Denoising the non-rigid corrected stacks...');
    tic;
    M2_red = zeros(size(M2_rough_red), 'single');
    M2_gcam = zeros(size(M2_rough_gcam), 'single');
    
    parfor t = 1:T
        M2_red(:,:,t) = medfilt2(M2_rough_red(:,:,t), [3 3]);
        M2_gcam(:,:,t) = medfilt2(M2_rough_gcam(:,:,t), [3 3]);
    end
    dnoise_toc = toc; disp(['Denoising time is:',num2str(dnoise_toc)]);
    
    %% NOW calculate the clean Final Ratio Image (M2)
    disp('Calculating Final Ratio Image...');
    M2_red_safe = M2_red; M2_red_safe(M2_red_safe == 0) = eps;
    M2 = M2_gcam ./ M2_red_safe;
    
    % Clean output for plotting and TIFF saving
    M2(isnan(M2) | isinf(M2)) = 0;
    ratio_cutoff = quantile(M2(:), 0.999);
    M2(M2 > ratio_cutoff) = ratio_cutoff;
    
    %% compute metrics (Mapping back to your original variable names)
    disp('Computing Metrics...');
    Y = Y_ratio; 
    M1 = M1_ratio;
    
    nnY = quantile(Y(:),0.005);
    mmY = quantile(Y(:),0.995);
    [cY,mY,vY] = motion_metrics(Y,10);
    [cM1,mM1,vM1] = motion_metrics(M1,10);
    [cM2,mM2,vM2] = motion_metrics(M2,10);
    T = length(cY);
    %% plot metrics (Invisible Figure)
    fig_metrics = figure('Name', 'Metrics', 'Position', [100, 100, 1200, 800], 'Visible', 'off');
        ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
        ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
        ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
        subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
        subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
        subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2,ax3],'xy')
    
    % SAVE METRICS FIGURE
    disp('Exporting Metrics Figure...');
    exportgraphics(fig_metrics, fullfile(path_name, [base_name, '_metrics.png']), 'Resolution', 300);
        
    %% plot shifts (Invisible Figure)
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
    shifts_x = squeeze(shifts_nr(:,1,:))';
    shifts_y = squeeze(shifts_nr(:,2,:))';
    patch_id = 1:size(shifts_x,2);
    str = strtrim(cellstr(int2str(patch_id.')));
    str = cellfun(@(x) ['patch # ',x],str,'un',0);
    
    fig_shifts = figure('Name', 'Shifts', 'Position', [150, 150, 800, 800], 'Visible', 'off');
        ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
                xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax1,ax2,ax3],'x')
        
    % SAVE SHIFTS FIGURE
    disp('Exporting Shifts Figure...');
    exportgraphics(fig_shifts, fullfile(path_name, [base_name, '_shifts.png']), 'Resolution', 300);
        
    %% EXTREMELY FAST MOVIE RENDERING (Direct Matrix to AVI)
    disp('Rendering and Saving Movie directly from matrices...');
    
    % Initialize VideoWriter (Switched to Motion JPEG AVI for Linux compatibility)
    v = VideoWriter(fullfile(path_name, [base_name, '_movie.avi']), 'Motion JPEG AVI');
    v.FrameRate = 30; 
    open(v);
    
    % Pre-load the 'bone' colormap
    cmap = bone(256);
    tic;
    for t = 1:T
        % 1. Extract current frames
        frameY = Y(:,:,t);
        frameM2 = M2(:,:,t);
        
        % 2. Normalize both frames to [0, 1] using your existing nnY and mmY limits
        imgY  = max(0, min(1, (frameY - nnY) / (mmY - nnY)));
        imgM2 = max(0, min(1, (frameM2 - nnY) / (mmY - nnY)));
        
        % 3. Concatenate them side-by-side (Left: Raw, Right: Corrected & Denoised)
        combined_img = [imgY, imgM2];
        
        % 4. Convert the normalized 2D matrix into a 3D RGB image
        rgb_frame = ind2rgb(round(combined_img * 255) + 1, cmap);
        
        % 5. Write the raw matrix directly to the video file
        writeVideo(v, rgb_frame);
    end
    close(v);
    v_toc = toc; disp(['Rendering and Saving Movie in:',v_toc]);
    
  %% Generate Time Vector from Extracted Metadata Parameters
    disp('Extracting timing metadata from Bio-Formats...');
    
    % Bio-Formats stores the metadata hashtable in the second column
    meta_hash = img_data{1, 2}; 
    
    % 1. Extract the acquisition intervals dynamically
    interval_1 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #1'));
    interval_2 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #2'));
    
    % 2. Extract the UTC Timestamps for Start and Marker events
    startTime_str = char(meta_hash.get('Global Information|Image|T|StartTime'));
    markerTime_str = char(meta_hash.get('Global Information|TimelineTrack|TimelineElement|Time #1'));
    
    % 3. Convert ISO 8601 strings to MATLAB datetime to calculate elapsed time
    t_start = datetime(startTime_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS''Z''', 'TimeZone', 'UTC');
    t_marker = datetime(markerTime_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS''Z''', 'TimeZone', 'UTC');
    
    % Calculate exact seconds between acquisition start and the interval change
    timeOfChange = seconds(t_marker - t_start); 
    
    % 4. Calculate frames per phase automatically
    totalFrames = T; 
    frames_phase1 = round(timeOfChange / interval_1);
    frames_phase2 = totalFrames - frames_phase1;
    
    % 5. Construct the final time vector
    timevec_phase1 = (0 : frames_phase1-1) * interval_1;
    timevec_phase2 = timevec_phase1(end) + (1 : frames_phase2) * interval_2;
    timevec = [timevec_phase1, timevec_phase2]; 
    
    % 6. Format the metadata block so Fiji/ImageJ parses it easily
    time_str = sprintf('%.2f,', timevec);
    time_str = time_str(1:end-1); % Remove the trailing comma
    metadata_string = sprintf('TimeVector_Seconds=[%s]\nTotalFrames=%d\nEventMarker=Frame_%d', time_str, totalFrames, frames_phase1);
    
    disp(['Time vector generated: ', num2str(totalFrames), ' frames spanning ', num2str(timevec(end)), ' seconds.']);

    %% Save Variable M2 as multi-page TIFF using saveastiff
    disp('Saving Non-Rigid Corrected & Denoised Stack (M2) as BigTIFF with embedded time metadata...');
    tiff_filename = fullfile(path_name, [base_name, '_M2_corrected.tif']);
    
    % Convert to uint8 (preserving 8-bit PMT depth)
    M2_uint8 = uint8(M2); 
    
    % Configure options for fast, lossless, BigTIFF saving
    options_tiff.color = false;
    options_tiff.compress = 'lzw'; 
    options_tiff.message = true;
    options_tiff.append = false;
    options_tiff.overwrite = true;
    options_tiff.big = true; 
    
    % Inject the custom time vector metadata here
    options_tiff.ImageDescription = metadata_string; 
    
    % Execute the save function 
    saveastiff(M2_uint8, tiff_filename, options_tiff);
    
    disp(['Finished processing: ', base_name]);
    close all;
end

disp('All selected files have been successfully processed and saved!');

%% Clean up threads
delete(gcp('nocreate'));