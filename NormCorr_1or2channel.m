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
    meta_hash = img_data{1, 2}; 
    omeMeta = img_data{1, 4};
   
    load_toc = toc; disp(['Loading image data time is: ' num2str(load_toc) ' secs']);
    
    % Detect how many channels are in the file using OME metadata
    try
        numChannels = omeMeta.getPixelsSizeC(0).getValue();
    catch
        numChannels = omeMeta.getChannelCount(0);
    end
    fprintf('Detected %d channel(s) in file.\n', numChannels);
    
    img_plane_image = img_data_planes(:,1);
    orig_class = class(img_plane_image{1});
    disp(['Detected original image format: ', orig_class]);
    
    if numChannels == 1
        %% --- PROCESS SINGLE COLOR ---
        disp('Processing as Single Color recording (Functional Channel only).');
        Y_gcam = single(cat(3, img_plane_image{1:1:end}));
        T = size(Y_gcam, 3);
        
    else
        %% --- PROCESS DUAL COLOR ---
        disp('Processing as Dual Color recording.');
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
        
        % Extract interleaved channels
        Y_red  = single(cat(3, img_plane_image{red_idx:2:end})); 
        Y_gcam = single(cat(3, img_plane_image{green_idx:2:end})); 
        T = size(Y_red, 3);
    end
    
    
%% --- DYNAMIC PARAMETER TUNING (Heuristic approach) ---
    disp('Dynamically analyzing image to tune NoRMCorre parameters...');
    
    % Use GCaMP for single color drift tracking, Red for dual color
    if numChannels == 1
        frame_first = Y_gcam(:,:,1);
        frame_mid = Y_gcam(:,:, round(T/2));
    else
        frame_first = Y_red(:,:,1);
        frame_mid = Y_red(:,:, round(T/2));
    end
    
    % Perform a fast 2D cross-correlation to find the global pixel drift
    c = normxcorr2(frame_first, frame_mid);
    [~, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c), imax);
    
    drift_y = abs(ypeak - size(frame_first, 1));
    drift_x = abs(xpeak - size(frame_first, 2));
    dyn_max_shift = ceil(max(drift_y, drift_x)) + 5; 
    dyn_max_shift = max(5, min(dyn_max_shift, 30));
    disp(['-> Dynamic max_shift set to: ', num2str(dyn_max_shift), ' pixels.']);
    
    if exist('roi_data', 'var') && ~isempty(roi_data)
        avg_radius = mean(roi_data(:, 3));
        opt_grid_dim = round((avg_radius * 2) * 2.5);
    else
        opt_grid_dim = 32; % Fallback
    end
    
    opt_grid_dim = max(16, round(opt_grid_dim / 16) * 16);
    dyn_grid_size = [opt_grid_dim, opt_grid_dim];
    disp(['-> Dynamic non-rigid grid_size set to: [', num2str(dyn_grid_size(1)), ',', num2str(dyn_grid_size(2)), '].']);
    
     
%% --- MOTION CORRECTION ---
    if numChannels == 1
        %% Single Color MC
        disp('Running Rigid Motion Correction on GCaMP channel...');
        options_rigid = NoRMCorreSetParms('d1',size(Y_gcam,1),'d2',size(Y_gcam,2), ...
            'bin_width', 200, 'max_shift', dyn_max_shift, 'us_fac', 50, 'init_batch', 200);
        tic; [M1_gcam, shifts1, template1, options_rigid] = normcorre(Y_gcam, options_rigid); 
        rigid_toc = toc; disp(['Rigid correction time: ' num2str(rigid_toc)]);
        
        disp('Running Non-Rigid Motion Correction on GCaMP channel...');
        options_nonrigid = NoRMCorreSetParms('d1',size(Y_gcam,1),'d2',size(Y_gcam,2), ...
            'grid_size', dyn_grid_size, 'mot_uf', 4, 'max_shift', dyn_max_shift, ...
            'max_dev', 3, 'us_fac', 50, 'bin_width', 50, 'init_batch', 100, 'mem_batch_size', 200);
        tic; [M2_rough_gcam, shifts2, template2, options_nonrigid] = normcorre_batch(Y_gcam, options_nonrigid); 
        non_toc = toc; disp(['Non-rigid correction time: ',num2str(non_toc)]);
        
    else
        %% Dual Color MC
        disp('Running Rigid Motion Correction on RED channel...');
        options_rigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
            'bin_width', 200, 'max_shift', dyn_max_shift, 'us_fac', 50, 'init_batch', 200);
        tic; [M1_red, shifts1, template1, options_rigid] = normcorre(Y_red, options_rigid); 
        rigid_toc = toc; disp(['Rigid correction time: ' num2str(rigid_toc)]);
        
        M1_gcam = apply_shifts(Y_gcam, shifts1, options_rigid);
        
        disp('Running Non-Rigid Motion Correction on RED channel...');
        options_nonrigid = NoRMCorreSetParms('d1',size(Y_red,1),'d2',size(Y_red,2), ...
            'grid_size', dyn_grid_size, 'mot_uf', 4, 'max_shift', dyn_max_shift, ...
            'max_dev', 3, 'us_fac', 50, 'bin_width', 50, 'init_batch', 100, 'mem_batch_size', 200);
        tic; [M2_rough_red, shifts2, template2, options_nonrigid] = normcorre_batch(Y_red, options_nonrigid); 
        non_toc = toc; disp(['Non-rigid correction time: ',num2str(non_toc)]);
        
        disp('Applying calculated shifts to GCAMP channel...');
        M2_rough_gcam = apply_shifts(Y_gcam, shifts2, options_nonrigid);
    end
    
%% --- DYNAMIC DENOISING TUNING ---
    disp('Analyzing SNR to tune adaptive denoising and deconvolution...');
    temp_vals = template2(:);
    signal_level = quantile(temp_vals, 0.98); 
    noise_level = quantile(temp_vals, 0.20);  
    current_snr = signal_level / max(noise_level, eps);
    disp(['-> Estimated SNR: ', num2str(current_snr)]);
    
    if current_snr < 2
        dyn_filt_size = [5 5]; 
        dyn_num_it = 3; % Minimal iterations to avoid noise amplification
        disp('-> Action: High noise. 5x5 filter, Deconv iterations limited to 3.');
    elseif current_snr < 5
        dyn_filt_size = [3 3]; 
        dyn_num_it = 5; % Standard iterations
        disp('-> Action: Normal SNR. 3x3 filter, Deconv iterations set to 5.');
    else
        dyn_filt_size = [2 2]; 
        dyn_num_it = 10; % Safe to push deconvolution further for maximum sharpness
        disp('-> Action: High SNR. 2x2 filter, Deconv iterations pushed to 10.');
    end
    
    disp('Denoising the non-rigid corrected stacks...');
    tic;
    M2_gcam = zeros(size(M2_rough_gcam), 'single');
    
    if numChannels == 1
        parfor t = 1:T
            M2_gcam(:,:,t) = medfilt2(M2_rough_gcam(:,:,t), dyn_filt_size);
        end
    else
        M2_red = zeros(size(M2_rough_red), 'single');
        parfor t = 1:T
            M2_red(:,:,t) = medfilt2(M2_rough_red(:,:,t), dyn_filt_size);
            M2_gcam(:,:,t) = medfilt2(M2_rough_gcam(:,:,t), dyn_filt_size);
        end
    end
    dnoise_toc = toc; disp(['Denoising time is: ', num2str(dnoise_toc)]);

%% --- DYNAMIC DECONVOLUTION (METADATA-DRIVEN) ---
    disp('Extracting optical parameters for dynamic deconvolution...');
    
    % Extract NA (Objective Numerical Aperture)
    NA = str2double(char(meta_hash.get('Global Information|Image|Objective|LensNA')));
    if isnan(NA), NA = 0.8; end % Fallback for common objectives
    
    % Extract Pixel Size in microns
    pixel_size = str2double(char(meta_hash.get('Global Information|Image|Scaling|Scaling|X'))) * 1e6; 
    if isnan(pixel_size), pixel_size = 0.5; end % Fallback
    
    if numChannels == 1
        % Single Color Wavelength Extraction
        wave_green = str2double(char(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|MultiTrackSetup|Track|Channel|EmissionWavelength #1')));
        if isnan(wave_green), wave_green = 520; end
        
        fprintf('Optical Profile: NA=%.2f, Green_λ=%dnm, Pixel=%.3fµm\n', NA, round(wave_green), pixel_size);
        
        % Calculate Theoretical PSF
        sigma_green = (0.61 * (wave_green/1000) / NA) / pixel_size;
        psf_green = fspecial('gaussian', round(sigma_green*4), sigma_green);
        
        disp(['Performing Richardson-Lucy deconvolution with ', num2str(dyn_num_it), ' iterations...']);
        tic;
        parfor t = 1:T
            M2_gcam(:,:,t) = deconvlucy(M2_gcam(:,:,t), psf_green, dyn_num_it);
        end
        deconv_toc = toc; disp(['Deconvolution time is: ', num2str(deconv_toc)]);
        
    else
        % Dual Color Wavelength Extraction
        wave_red = str2double(char(meta_hash.get(['Global Experiment|AcquisitionBlock|TimeSeriesSetup|MultiTrackSetup|Track|Channel|EmissionWavelength #' num2str(red_idx)])));
        wave_green = str2double(char(meta_hash.get(['Global Experiment|AcquisitionBlock|TimeSeriesSetup|MultiTrackSetup|Track|Channel|EmissionWavelength #' num2str(green_idx)])));
        
        if isnan(wave_red), wave_red = 600; end
        if isnan(wave_green), wave_green = 520; end
        
        fprintf('Optical Profile: NA=%.2f, Red_λ=%dnm, Green_λ=%dnm, Pixel=%.3fµm\n', NA, round(wave_red), round(wave_green), pixel_size);

        % Calculate Theoretical PSFs
        sigma_red = (0.61 * (wave_red/1000) / NA) / pixel_size;
        sigma_green = (0.61 * (wave_green/1000) / NA) / pixel_size;
        
        psf_red = fspecial('gaussian', round(sigma_red*4), sigma_red);
        psf_green = fspecial('gaussian', round(sigma_green*4), sigma_green);
        
        disp(['Performing Richardson-Lucy deconvolution with ', num2str(dyn_num_it), ' iterations...']);
        tic;
        parfor t = 1:T
            M2_red(:,:,t) = deconvlucy(M2_red(:,:,t), psf_red, dyn_num_it);
            M2_gcam(:,:,t) = deconvlucy(M2_gcam(:,:,t), psf_green, dyn_num_it);
        end
        deconv_toc = toc; disp(['Deconvolution time is: ', num2str(deconv_toc)]);
    end


    %% --- CREATE PREPROCESSING FOLDER ---
disp('Setting up preprocessing directory...');
% Define the new subfolder path inside the original directory
prep_folder = fullfile(path_name, 'preprocessing');

% Check if the folder already exists; if not, create it
if ~exist(prep_folder, 'dir')
    mkdir(prep_folder);
    disp(['Created new directory: ', prep_folder]);
else
    disp(['Directory already exists: ', prep_folder]);
end
%% Compute Metrics (Tracking the Green Channel)
    disp('Computing Metrics (GCaMP Channel)...');
    nnY_gcam = quantile(Y_gcam(:),0.005);
    mmY_gcam = quantile(Y_gcam(:),0.995);
    [cY,mY,vY] = motion_metrics(Y_gcam,10);
    [cM1,mM1,vM1] = motion_metrics(M1_gcam,10);
    [cM2,mM2,vM2] = motion_metrics(M2_gcam,10);
    
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
    exportgraphics(fig_metrics, fullfile(prep_folder, [base_name, '_metrics_GCaMP.png']), 'Resolution', 300);
        
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(Y_gcam)-1,T);
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
    exportgraphics(fig_shifts, fullfile(prep_folder, [base_name, '_shifts.png']), 'Resolution', 300);
        
%% FAST MOVIE RENDERING
    disp('Rendering RGB Movie directly from matrices...');
    tic;
    if numChannels == 1
        v = VideoWriter(fullfile(prep_folder, [base_name, '_movie_SingleColor.avi']), 'Motion JPEG AVI');
        v.FrameRate = 30; open(v);
        
        for t = 1:T
            img_raw_G = max(0, min(1, (Y_gcam(:,:,t) - nnY_gcam) / (mmY_gcam - nnY_gcam)));
            img_cor_G = max(0, min(1, (M2_gcam(:,:,t) - nnY_gcam) / (mmY_gcam - nnY_gcam)));
            
            raw_rgb = cat(3, zeros(size(img_raw_G)), img_raw_G, zeros(size(img_raw_G)));
            cor_rgb = cat(3, zeros(size(img_cor_G)), img_cor_G, zeros(size(img_cor_G)));
            
            writeVideo(v, [raw_rgb, cor_rgb]);
        end
    else
        v = VideoWriter(fullfile(prep_folder, [base_name, '_movie_DualColor.avi']), 'Motion JPEG AVI');
        v.FrameRate = 30; open(v);
        nnY_red = quantile(Y_red(:),0.005); mmY_red = quantile(Y_red(:),0.995);
        
        for t = 1:T
            img_raw_R = max(0, min(1, (Y_red(:,:,t) - nnY_red) / (mmY_red - nnY_red)));
            img_raw_G = max(0, min(1, (Y_gcam(:,:,t) - nnY_gcam) / (mmY_gcam - nnY_gcam)));
            img_cor_R = max(0, min(1, (M2_red(:,:,t) - nnY_red) / (mmY_red - nnY_red)));
            img_cor_G = max(0, min(1, (M2_gcam(:,:,t) - nnY_gcam) / (mmY_gcam - nnY_gcam)));
            
            raw_rgb = cat(3, img_raw_R, img_raw_G, zeros(size(img_raw_R)));
            cor_rgb = cat(3, img_cor_R, img_cor_G, zeros(size(img_cor_R)));
            
            writeVideo(v, [raw_rgb, cor_rgb]);
        end
    end
    close(v);
    v_toc = toc; disp(['Rendering and Saving Movie in: ', num2str(v_toc)]);
    
%% Generate Time Vector and ROI Data from Bio-Formats Metadata
    disp('Extracting timing and ROI metadata from Bio-Formats...');
    
    totalFrames = T; 
    timeOfChange = NaN;
    frames_phase1 = 0;
    
    % --- 1. EXTRACT EVENT MARKER TIME ---
    % Find the time of drug injection / stimulation
    startTime_str = char(meta_hash.get('Global Information|Image|T|StartTime'));
    markerTime_str = char(meta_hash.get('Global Information|TimelineTrack|TimelineElement|Time #1'));
    
    if ~isempty(startTime_str) && ~isempty(markerTime_str)
        try
            % Truncate to 19 characters to avoid MATLAB fractional-second parsing errors
            t_start = datetime(startTime_str(1:19), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
            t_marker = datetime(markerTime_str(1:19), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
            timeOfChange = seconds(t_marker - t_start); 
        catch
            try
                % Fallback to automatic parsing
                t_start = datetime(startTime_str);
                t_marker = datetime(markerTime_str);
                timeOfChange = seconds(t_marker - t_start);
            catch
                disp('Warning: Could not parse Event Marker timestamps.');
            end
        end
    end
    
    % --- 2. EXTRACT TIME VECTOR ---
    timevec = zeros(1, T);
    valid_ome_time = true;
    
    % ATTEMPT A: Frame-by-Frame DeltaT from OME-XML (Handles "Original" compressed format natively)
    try
        for i = 1:T
            % OME-XML Plane index is 0-based. Get time for the first channel of each timepoint.
            plane_index = (i - 1) * numChannels; 
            dt_obj = omeMeta.getPlaneDeltaT(0, plane_index);
            if ~isempty(dt_obj)
                timevec(i) = double(dt_obj.value());
            else
                valid_ome_time = false; break;
            end
        end
    catch
        valid_ome_time = false;
    end
    
    % ATTEMPT B: Legacy Hashtable parsing (Your original method for "Uncompressed")
    if ~valid_ome_time || any(isnan(timevec))
        disp('OME-XML DeltaT unavailable. Falling back to Hashtable Interval parsing...');
        interval_1 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #1'));
        interval_2 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #2'));
        
        if ~isnan(interval_1) && ~isnan(interval_2) && ~isnan(timeOfChange)
            frames_phase1 = round(timeOfChange / interval_1);
            frames_phase2 = totalFrames - frames_phase1;
            
            timevec_phase1 = (0 : frames_phase1-1) * interval_1;
            timevec_phase2 = timevec_phase1(end) + (1 : frames_phase2) * interval_2;
            timevec = [timevec_phase1, timevec_phase2];
        else
            % ATTEMPT C: Universal Fallback using standard Time Increment
            try
                dt = double(omeMeta.getPixelsTimeIncrement(0).value());
            catch
                dt = 1.0; % Default to 1 sec
                warning('SAFE FAIL: All timing metadata missing. Defaulting to 1s/frame.');
            end
            timevec = (0:T-1) * dt;
        end
    end
    
    % Sync the marker frame index if it wasn't explicitly calculated in Attempt B
    if frames_phase1 == 0 
        if ~isnan(timeOfChange)
            [~, frames_phase1] = min(abs(timevec - timeOfChange));
        else
            frames_phase1 = round(T/3); % Fallback injection frame
            disp('No injection marker found. Defaulting event marker to 1/3 of the recording.');
        end
    end
    
    % Convert to comma-separated string for TIFF header
    time_str = sprintf('%.2f,', timevec);
    time_str = time_str(1:end-1); 
    
    % --- 3. EXTRACT ROI GEOMETRIES DYNAMICALLY ---
    roi_metadata_str = ''; 
    num_ROIs = 0;
    
    while true
        % Check for both zero-padded (e.g., #01) and normal (e.g., #1) formats
        idx_pad = sprintf('%02d', num_ROIs + 1); 
        idx_nopad = sprintf('%d', num_ROIs + 1); 
        
        str_X_pad = ['Global Layer|Circle|Geometry|CenterX #', idx_pad];
        str_X_nopad = ['Global Layer|Circle|Geometry|CenterX #', idx_nopad];
        
        % Try to fetch the X coordinate
        val_X = meta_hash.get(str_X_pad);
        active_idx = idx_pad;
        
        if isempty(val_X)
            val_X = meta_hash.get(str_X_nopad);
            active_idx = idx_nopad;
        end
        
        % If both formats return empty, we have reached the end of the ROIs
        if isempty(val_X)
            break;
        end
        
        % If an ROI is found, increment the count and fetch Y and Radius
        num_ROIs = num_ROIs + 1;
        
        str_Y = ['Global Layer|Circle|Geometry|CenterY #', active_idx];
        str_R = ['Global Layer|Circle|Geometry|Radius #', active_idx];
        
        cx = str2double(val_X);
        cy = str2double(meta_hash.get(str_Y));
        rad = str2double(meta_hash.get(str_R));
        
        % Append to the master string
        roi_metadata_str = sprintf('%sROI_%02d=[%.2f,%.2f,%.2f]\n', roi_metadata_str, num_ROIs, cx, cy, rad);
    end
    
    metadata_string = sprintf('TimeVector_Seconds=[%s]\nTotalFrames=%d\nEventMarker=Frame_%d\n%s', ...
        time_str, totalFrames, frames_phase1, roi_metadata_str);
        
    disp(['Time vector generated: ', num2str(totalFrames), ' frames spanning ', num2str(max(timevec)), ' seconds.']);
    disp(['Successfully embedded ', num2str(num_ROIs), ' ROIs into TIFF metadata.']);
    
    
%% Save Channels as Separate BigTIFFs

disp('Exporting corrected images to preprocessing folder...');

% Ensure we have the base name of the current file
[~, base_name, ~] = fileparts(file_list{f_idx}); % Assuming 'f' is your file loop index

% Define the NEW save paths using 'prep_folder' instead of 'path_name'
save_name_green = fullfile(prep_folder,[base_name, '_M2_GREEN.tif']);

    options_tiff.color = false;
    options_tiff.compress = 'lzw'; 
    options_tiff.message = true;
    options_tiff.append = false;
    options_tiff.overwrite = true;
    options_tiff.big = true; 
    options_tiff.ImageDescription = metadata_string; 
    % Save the files

    % Always export Green Channel
    M2_gcam_export = cast(M2_gcam, orig_class);
    saveastiff(M2_gcam_export, save_name_green, options_tiff);
    
    % If Dual Channel, also export Red
    if numChannels > 1
        save_name_red = fullfile(prep_folder, [base_name, '_M2_RED.tif']);
        M2_red_export = cast(M2_red, orig_class);
        saveastiff(M2_red_export, save_name_red, options_tiff);
    end
    
    disp(['Finished processing: ', base_name]);
    close all;
end
disp('All selected files have been successfully processed and saved!');

%% Clean up threads
delete(gcp('nocreate'));