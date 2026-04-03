close all; clear; clc;
%% Parallel Setup
if isempty(gcp('nocreate'))
    parpool('local', 14);
end

%% Select Multiple Files
[file_list, path_name] = uigetfile({'*.czi;*.tif;*.tiff', 'Image Files (*.czi, *.tif, *.tiff)'}, ...
    'Select Image Files', 'MultiSelect', 'on');

if isequal(file_list, 0)
    disp('User selected Cancel. Exiting...');
    return;
end

if ischar(file_list)
    file_list = {file_list};
end

%% Loop Through Each Selected File
for f_idx = 1:length(file_list)
    current_file = fullfile(path_name, file_list{f_idx});
    [~, base_name, ~] = fileparts(file_list{f_idx});
    
    fprintf('\n=== Processing File %d of %d: %s ===\n', f_idx, length(file_list), base_name);
    
    %% Load and Prepare Data
    tic;
    disp('Loading image data...');
    img_data = bfopen(current_file);
    img_data_planes = img_data{1, 1};
    load_toc = toc; disp(['Loading image data time is: ' num2str(load_toc) ' secs']);
    img_plane_image = img_data_planes(:,1);
    % Extract base channels independently
    Y_red = single(cat(3,img_plane_image{1:2:end})); 
    Y_gcam = single(cat(3,img_plane_image{2:2:end})); 
    T = size(Y_red, ndims(Y_red));
    
    %% Set parameters (Rigid motion correction)
    disp('Running Rigid Motion Correction on RED channel...');
    options_rigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
        'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
    
    %% Perform Rigid motion correction
    tic; [M1_red, shifts1, template1, options_rigid] = normcorre(Y_red, options_rigid); 
    rigid_toc = toc; disp(['Rigid correction time is: ' num2str(rigid_toc)]);
    
    % Apply rigid shifts to Green channel
    M1_gcam = apply_shifts(Y_gcam, shifts1, options_rigid);
    
    %% Set parameters (Non-Rigid motion correction)
    disp('Running Non-Rigid Motion Correction on RED channel...');
    options_nonrigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
        'grid_size',[32,32], 'mot_uf',4, 'max_shift',15, 'max_dev',3, 'us_fac',50, ...
        'bin_width', 50, 'init_batch', 100, 'mem_batch_size', 200);
        
    % Calculate non-rigid shifts on Red
    tic; [M2_rough_red, shifts2, template2, options_nonrigid] = normcorre_batch(Y_red, options_nonrigid); 
    non_toc = toc; disp(['Non-rigid correction time is: ',num2str(non_toc)]);
    
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
    dnoise_toc = toc; disp(['Denoising time is: ',num2str(dnoise_toc)]);
    
    %% Compute Metrics (Tracking the Green Channel for Functional Signal Quality)
    disp('Computing Metrics (GCaMP Channel)...');
    nnY_gcam = quantile(Y_gcam(:),0.005);
    mmY_gcam = quantile(Y_gcam(:),0.995);
    [cY,mY,vY] = motion_metrics(Y_gcam,10);
    [cM1,mM1,vM1] = motion_metrics(M1_gcam,10);
    [cM2,mM2,vM2] = motion_metrics(M2_gcam,10);
    
    %% Plot metrics (Invisible Figure)
    fig_metrics = figure('Name', 'Metrics', 'Position', [100, 100, 1200, 800], 'Visible', 'off');
        ax1 = subplot(2,3,1); imagesc(mY,[nnY_gcam,mmY_gcam]);  axis equal; axis tight; axis off; title('mean raw GCAMP','fontsize',14,'fontweight','bold')
        ax2 = subplot(2,3,2); imagesc(mM1,[nnY_gcam,mmY_gcam]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
        ax3 = subplot(2,3,3); imagesc(mM2,[nnY_gcam,mmY_gcam]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
        subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw GCAMP','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
        subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
        subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2,ax3],'xy')
    
    disp('Exporting Metrics Figure...');
    exportgraphics(fig_metrics, fullfile(path_name, [base_name, '_metrics_GCaMP.png']), 'Resolution', 300);
        
    %% Plot shifts (Invisible Figure)
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(Y_red)-1,T);
    shifts_x = squeeze(shifts_nr(:,1,:))';
    shifts_y = squeeze(shifts_nr(:,2,:))';
    
    fig_shifts = figure('Name', 'Shifts', 'Position', [150, 150, 800, 800], 'Visible', 'off');
        ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw GCAMP','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
                xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax1,ax2,ax3],'x')
        
    disp('Exporting Shifts Figure...');
    exportgraphics(fig_shifts, fullfile(path_name, [base_name, '_shifts.png']), 'Resolution', 300);
        
    %% FAST MOVIE RENDERING (Direct Dual-Color RGB Matrix to AVI)
    disp('Rendering Dual-Color RGB Movie directly from matrices...');
    v = VideoWriter(fullfile(path_name, [base_name, '_movie_DualColor.avi']), 'Motion JPEG AVI');
    v.FrameRate = 30; 
    open(v);
    
    % Establish normalization bounds for movie rendering
    nnY_red = quantile(Y_red(:),0.005); mmY_red = quantile(Y_red(:),0.995);
    
    tic;
    for t = 1:T
        % Normalize raw and corrected frames to [0, 1] bounds
        img_raw_R = max(0, min(1, (Y_red(:,:,t) - nnY_red) / (mmY_red - nnY_red)));
        img_raw_G = max(0, min(1, (Y_gcam(:,:,t) - nnY_gcam) / (mmY_gcam - nnY_gcam)));
        
        img_cor_R = max(0, min(1, (M2_red(:,:,t) - nnY_red) / (mmY_red - nnY_red)));
        img_cor_G = max(0, min(1, (M2_gcam(:,:,t) - nnY_gcam) / (mmY_gcam - nnY_gcam)));
        
        % Build TrueColor RGB images (Red, Green, Blue=zeros)
        raw_rgb = cat(3, img_raw_R, img_raw_G, zeros(size(img_raw_R)));
        cor_rgb = cat(3, img_cor_R, img_cor_G, zeros(size(img_cor_R)));
        
        % Concatenate Side-by-Side (Left: Raw, Right: Corrected)
        combined_img = [raw_rgb, cor_rgb];
        
        writeVideo(v, combined_img);
    end
    close(v);
    v_toc = toc; disp(['Rendering and Saving Movie in: ', num2str(v_toc)]);
    
 
    %% Generate Time Vector and ROI Data from Bio-Formats Metadata
    disp('Extracting timing and ROI metadata from Bio-Formats...');
    meta_hash = img_data{1, 2}; 
    
    % --- 1. EXTRACT TIME VECTOR ---
    interval_1 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #1'));
    interval_2 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #2'));
    
    startTime_str = char(meta_hash.get('Global Information|Image|T|StartTime'));
    markerTime_str = char(meta_hash.get('Global Information|TimelineTrack|TimelineElement|Time #1'));
    
    t_start = datetime(startTime_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS''Z''', 'TimeZone', 'UTC');
    t_marker = datetime(markerTime_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS''Z''', 'TimeZone', 'UTC');
    
    timeOfChange = seconds(t_marker - t_start); 
    totalFrames = T; 
    frames_phase1 = round(timeOfChange / interval_1);
    frames_phase2 = totalFrames - frames_phase1;
    
    timevec_phase1 = (0 : frames_phase1-1) * interval_1;
    timevec_phase2 = timevec_phase1(end) + (1 : frames_phase2) * interval_2;
    timevec = [timevec_phase1, timevec_phase2]; 
    
    time_str = sprintf('%.2f,', timevec);
    time_str = time_str(1:end-1); 
    
    % --- 2. EXTRACT ROI GEOMETRIES ---
    num_ROIs = 20;
    roi_metadata_str = ''; % Blank string to append ROI text
    
    for r = 1:num_ROIs
        idx_str = sprintf('%02d', r); % Formats 1 as '01'
        
        str_X = ['Global Layer|Circle|Geometry|CenterX #', idx_str];
        str_Y = ['Global Layer|Circle|Geometry|CenterY #', idx_str];
        str_R = ['Global Layer|Circle|Geometry|Radius #', idx_str];
        
        % Extract values 
        cx = str2double(meta_hash.get(str_X));
        cy = str2double(meta_hash.get(str_Y));
        rad = str2double(meta_hash.get(str_R));
        
        % Build a clean string format: ROI_01=[X,Y,Radius]
        roi_metadata_str = sprintf('%sROI_%02d=[%.2f,%.2f,%.2f]\n', roi_metadata_str, r, cx, cy, rad);
    end
    
    % --- 3. BUILD FINAL METADATA HEADER ---
    % Combine Time parameters and ROI parameters into one massive string
    metadata_string = sprintf('TimeVector_Seconds=[%s]\nTotalFrames=%d\nEventMarker=Frame_%d\n%s', ...
        time_str, totalFrames, frames_phase1, roi_metadata_str);
    
    disp(['Time vector generated: ', num2str(totalFrames), ' frames spanning ', num2str(timevec(end)), ' seconds.']);
    disp(['Successfully embedded ', num2str(num_ROIs), ' ROIs into TIFF metadata.']);
    
    %% Save Red and Green Channels as Separate BigTIFFs
    disp('Saving separate Red and Green Channels as BigTIFFs...');
    
    % Configure options for fast, lossless, BigTIFF saving
    options_tiff.color = false;
    options_tiff.compress = 'lzw'; 
    options_tiff.message = true;
    options_tiff.append = false;
    options_tiff.overwrite = true;
    options_tiff.big = true; 
    
    % Inject the massive custom metadata string here
    options_tiff.ImageDescription = metadata_string; 
    
    % 1. Export RED Channel 
    tiff_red = fullfile(path_name, [base_name, '_M2_RED.tif']);
    saveastiff(uint16(M2_red), tiff_red, options_tiff);
    
    % 2. Export GREEN Channel 
    tiff_gcam = fullfile(path_name, [base_name, '_M2_GREEN.tif']);
    saveastiff(uint16(M2_gcam), tiff_gcam, options_tiff);
    
    disp(['Finished processing: ', base_name]);
    close all;
end
disp('All selected files have been successfully processed and saved!');

%% Clean up threads
delete(gcp('nocreate'));