
close all; clear; clc;
%% Parallel Setup

if isempty(gcp('nocreate')) 
    parpool("Processes",14)
elseif contains(class(gcp('nocreate')), 'parallel.ThreadPool')
    % Clean up threads
    delete(gcp('nocreate'));
    parpool("Processes",14)
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
    meta_hash = img_data{1, 2}; 
   
%% --- DYNAMIC CHANNEL DETECTION (METADATA FIRST, LUT FALLBACK) ---
        
    % 1. Extract Channel Names
    chan1_name = lower(char(meta_hash.get('Global Information|Image|Channel|Name #1')));
    chan2_name = lower(char(meta_hash.get('Global Information|Image|Channel|Name #2')));
    
    disp(['Metadata Names - Ch1: ', chan1_name, ' | Ch2: ', chan2_name]);
    
    % Define priority keywords for Red channel
    red_keys = {'red', 'tomato', 'tdtom', 'mcherry', 'dsred', '555', '568', 'a555', 'a568'};
    
    % --- STEP 1: KEYWORD CHECK (High Priority) ---
    is_ch1_red = any(contains(chan1_name, red_keys));
    is_ch2_red = any(contains(chan2_name, red_keys));
    
    if is_ch1_red && ~is_ch2_red
        red_idx = 1; green_idx = 2;
        disp('Red channel identified by Metadata Name: Channel 1');
    elseif is_ch2_red && ~is_ch1_red
        red_idx = 2; green_idx = 1;
        disp('Red channel identified by Metadata Name: Channel 2');
    else
        % --- STEP 2: LUT MATRIX CHECK (Fallback) ---
        img_col_cel = img_data{1,3};
        disp('Metadata names ambiguous. Falling back to LUT Matrix check...');
        
        lut_ch2 = img_col_cel{2};
        
        % Compare Red column (1) vs Green column (2) at the max intensity index

        r2 = lut_ch2(end, 1); g2 = lut_ch2(end, 2);
        
        if g2 > r2
            red_idx = 1; green_idx = 2;
            disp('Red channel identified by LUT Color: Channel 1');
        elseif r2 > g2
            red_idx = 2; green_idx = 1;
            disp('Red channel identified by LUT Color: Channel 2');
        else
            % --- STEP 3: ULTIMATE FALLBACK ---
            warning('All detection methods failed. Defaulting to Ch2=Red, Ch1=Green.');
            red_idx = 2; green_idx = 1;
        end
    end
    
    fprintf('FINAL ASSIGNMENT: Red = Channel %d, Green = Channel %d\n', red_idx, green_idx);

    %% --- EXTRACT CHANNELS BASED ON DETECTED INDICES ---
    % Use the indices to step through the interleaved plane array correctly
    Y_red  = single(cat(3, img_plane_image{red_idx:2:end})); 
    Y_gcam = single(cat(3, img_plane_image{green_idx:2:end})); 
    
    T = size(Y_red, ndims(Y_red));
    
    % Detect Bit-Depth for later export
    orig_class = class(img_plane_image{1});
    disp(['Detected original image format: ', orig_class]);
    
%% --- DYNAMIC PARAMETER TUNING (Heuristic approach) ---
    disp('Dynamically analyzing image to tune NoRMCorre parameters...');
    
    % 1. Dynamically calculate max_shift based on actual sample drift
    % Take the first frame and a frame from the middle of the recording
    frame_first = Y_red(:,:,1);
    frame_mid = Y_red(:,:, round(T/2));
    
    % Perform a fast 2D cross-correlation to find the global pixel drift
    c = normxcorr2(frame_first, frame_mid);
    [~, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c), imax);
    
    % Calculate exact pixel drift and add a 5-pixel safety buffer
    drift_y = abs(ypeak - size(frame_first, 1));
    drift_x = abs(xpeak - size(frame_first, 2));
    dyn_max_shift = ceil(max(drift_y, drift_x)) + 5; 
    
    % Clamp the shift to reasonable limits (e.g., between 5 and 30 pixels)
    dyn_max_shift = max(5, min(dyn_max_shift, 30));
    disp(['-> Dynamic max_shift set to: ', num2str(dyn_max_shift), ' pixels.']);
    
    
    % 2. Dynamically calculate grid_size based on Metadata ROIs
    % We assume roi_data (from your metadata extraction) exists. 
    % If average cell radius is 5 pixels, diameter is 10. Optimal grid is ~2.5x diameter.
    if exist('roi_data', 'var') && ~isempty(roi_data)
        avg_radius = mean(roi_data(:, 3));
        opt_grid_dim = round((avg_radius * 2) * 2.5);
    else
        opt_grid_dim = 32; % Fallback
    end
    
    % Snap the grid dimension to the nearest multiple of 16 for FFT efficiency
    opt_grid_dim = max(16, round(opt_grid_dim / 16) * 16);
    dyn_grid_size = [opt_grid_dim, opt_grid_dim];
    disp(['-> Dynamic non-rigid grid_size set to: [', num2str(dyn_grid_size(1)), ',', num2str(dyn_grid_size(2)), '].']);
    
     
    %% set parameters (Using Dynamic Variables)
    disp('Running Rigid Motion Correction on RED channel...');
    options_rigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
        'bin_width', 200, ...
        'max_shift', dyn_max_shift, ...   % <--- INJECTED DYNAMIC SHIFT
        'us_fac', 50, 'init_batch', 200);

   
    %% Perform Rigid motion correction
    tic; [M1_red, shifts1, template1, options_rigid] = normcorre(Y_red, options_rigid); 
    rigid_toc = toc; disp(['Rigid correction time is: ' num2str(rigid_toc)]);
    
    % Apply rigid shifts to Green channel
    M1_gcam = apply_shifts(Y_gcam, shifts1, options_rigid);
    
    %% now try non-rigid motion correction
    disp('Running Non-Rigid Motion Correction on RED channel...');
    options_nonrigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
        'grid_size', dyn_grid_size, ...   % <--- INJECTED DYNAMIC GRID SIZE
        'mot_uf', 4, ...
        'max_shift', dyn_max_shift, ...   % <--- INJECTED DYNAMIC SHIFT
        'max_dev', 3, 'us_fac', 50, ...
        'bin_width', 50, 'init_batch', 100, 'mem_batch_size', 200);

    % Calculate non-rigid shifts on Red
    tic; [M2_rough_red, shifts2, template2, options_nonrigid] = normcorre_batch(Y_red, options_nonrigid); 
    non_toc = toc; disp(['Non-rigid correction time is: ',num2str(non_toc)]);
    
    % Apply identical non-rigid shifts to Green
    disp('Applying calculated shifts to GCAMP channel...');
    M2_rough_gcam = apply_shifts(Y_gcam, shifts2, options_nonrigid);
    
    %% --- DYNAMIC DENOISING TUNING ---
    disp('Analyzing SNR to tune adaptive denoising...');
    
    % Use the final non-rigid template to estimate the sample's SNR
    % because the template has the best signal representation.
    temp_vals = template2(:);
    signal_level = quantile(temp_vals, 0.98); % Intensity of cells
    noise_level = quantile(temp_vals, 0.20);  % Intensity of background
    
    % Calculate a simple SNR proxy
    current_snr = signal_level / max(noise_level, eps);
    disp(['-> Estimated SNR: ', num2str(current_snr)]);
    
    % Determine Median Filter size based on SNR thresholds
    if current_snr < 2
        % Very Noisy: Use a larger 5x5 filter
        dyn_filt_size = [5 5];
        disp('-> Action: High noise detected. Applying aggressive 5x5 median filter.');
    elseif current_snr < 5
        % Moderate Noise: Standard 3x3 filter
        dyn_filt_size = [3 3];
        disp('-> Action: Normal SNR. Applying standard 3x3 median filter.');
    else
        % Very Clean: Use a small 2x2 filter to preserve fine cell edges
        dyn_filt_size = [2 2];
        disp('-> Action: High SNR. Applying conservative 2x2 median filter.');
    end

    %% Denoise using Dynamic Filter Size
    disp('Denoising the non-rigid corrected stacks...');
    tic;
    M2_red = zeros(size(M2_rough_red), 'single');
    M2_gcam = zeros(size(M2_rough_gcam), 'single');
    
    parfor t = 1:T
        % Use the dynamically selected filter size
        M2_red(:,:,t) = medfilt2(M2_rough_red(:,:,t), dyn_filt_size);
        M2_gcam(:,:,t) = medfilt2(M2_rough_gcam(:,:,t), dyn_filt_size);
    end
    dnoise_toc = toc; disp(['Denoising time is: ', num2str(dnoise_toc)]);
    
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
    disp(['Saving separate Red and Green Channels as ', orig_class, ' BigTIFFs...']);
    
    % Configure options for fast, lossless, BigTIFF saving
    options_tiff.color = false;
    options_tiff.compress = 'lzw'; 
    options_tiff.message = true;
    options_tiff.append = false;
    options_tiff.overwrite = true;
    options_tiff.big = true; 
    options_tiff.ImageDescription = metadata_string; 
    
    % dynamically cast the matrices back to their original bit-depth
    M2_red_export = cast(M2_red, orig_class);
    M2_gcam_export = cast(M2_gcam, orig_class);
    
    % 1. Export RED Channel 
    tiff_red = fullfile(path_name, [base_name, '_M2_RED.tif']);
    saveastiff(M2_red_export, tiff_red, options_tiff);
    
    % 2. Export GREEN Channel 
    tiff_gcam = fullfile(path_name, [base_name, '_M2_GREEN.tif']);
    saveastiff(M2_gcam_export, tiff_gcam, options_tiff);
    
    disp(['Finished processing: ', base_name]);
    close all;
end
disp('All selected files have been successfully processed and saved!');

%% Clean up threads
delete(gcp('nocreate'));