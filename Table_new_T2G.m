close all; clear; clc;

%% 1. Load Excel Time Series Data with Multiple ROIs
disp('Select the CSV file containing the extracted fluorescence data...');
[filename, filepath] = uigetfile('*.csv', 'Select Fluorescence CSV');
if isequal(filename, 0)
    disp('User canceled. Exiting...'); return;
end

full_filepath = fullfile(filepath, filename);
disp(['Loading: ', filename]);

% Read the table
dataTbl = readtable(full_filepath);

% Identify the time column
timeColIdx = find(contains(dataTbl.Properties.VariableNames, 'Time', 'IgnoreCase', true), 1);
if isempty(timeColIdx)
    error('No time column found!');
end

% Extract and process time vector
timevec_raw = dataTbl{:, timeColIdx};
timevec_raw = rmmissing(timevec_raw);
timevec_sec = timevec_raw * 60; 

% Estimate sampling rate and injection mark
marks = find(abs(diff(round(diff(timevec_sec)))) >= 1.5) + 2;
if size(marks, 1) > 2
    marks = marks(2);
end

timevec = timevec_raw(:); % Ensure column vector in minutes

% --- Extract GCaMP (Green) Channel ---
fluorColsIdx = find(contains(dataTbl.Properties.VariableNames, 'IntensityMeanThrs__GCMP6_T2', 'IgnoreCase', true));
fluorColsIdx(fluorColsIdx == timeColIdx) = [];

if isempty(fluorColsIdx)
    error('No GCaMP fluorescence data columns found!');
end

neuronts = dataTbl{:, fluorColsIdx}';
neuronts = rmmissing(neuronts, 2);

% --- Extract tdTomato (Red) Reference Channel ---
reffluorColsIdx = find(contains(dataTbl.Properties.VariableNames, 'IntensityMeanThrs__tdTom_T1', 'IgnoreCase', true));
reffluorColsIdx(reffluorColsIdx == timeColIdx) = [];

if isempty(reffluorColsIdx)
    disp('No tdTomato reference columns found. Proceeding with single-channel GCaMP only.');
    refneuronts = [];
else
    refneuronts = dataTbl{:, reffluorColsIdx}';
    refneuronts = rmmissing(refneuronts, 2);
end

numCells = size(neuronts, 1);
nROIs = numCells;

%% Define Event
%evnt = {'1mM CuSO_{4}'};
evnt = {'20 \muM Thap'};
%evnt = {'Akh 10 ng'};
disp(['Successfully loaded ', num2str(numCells), ' ROIs over ', num2str(length(timevec)), ' timepoints.']);

%% 2. Calculate Ratio, dF/F0, and Apply Savitzky-Golay Filter
disp('Calculating Ratio, dF/F0, and applying Savitzky-Golay filter...');

% --- Calculate Ratio (if dual-channel) or Use Raw (if single-channel) ---
if ~isempty(refneuronts)
    disp('Dual-channel data detected. Calculating GCaMP / tdTomato ratio...');
    % Protect against division by near-zero in the Red channel
    refneuronts(refneuronts < 0.01) = 0.01; 
    
    roiMean_raw = neuronts ./ refneuronts;
    
    % CLEANUP: Clamp extreme artifacts in the Ratio
    roiMean_raw(isnan(roiMean_raw) | isinf(roiMean_raw)) = 0;
    roiMean_raw(roiMean_raw > 20) = 20; % Cap massive ratio noise spikes 
    roiMean_raw(roiMean_raw < 0) = 0;   % Ratios cannot be negative
else
    disp('Single-channel data detected. Proceeding without ratio normalization...');
    roiMean_raw = neuronts;
end

% --- Calculate dF/F0 ---
% Calculate F0 (Baseline of the Raw signal/Ratio) before injection
baseline_start = max(1, marks(1)-11);
baseline_end = max(2, marks(1)-1);
F0_raw = mean(roiMean_raw(:, baseline_start:baseline_end), 2);

% Protect the baseline against division by near-zero 
F0_raw(F0_raw < 0.01) = 0.01; 

% Calculate dF/F0 
df_ratio = bsxfun(@rdivide, roiMean_raw, F0_raw) - 1;

% CLEANUP: Clamp extreme artifacts in dF/F0
df_ratio(isnan(df_ratio) | isinf(df_ratio)) = 0; 

% --- GLOBAL SAVITZKY-GOLAY SMOOTHING ---
smo_frm = 15;
fprintf('Applying Savitzky-Golay filter (Order 3, Frame %d ) to data...\n',smo_frm);
for i = 1:numCells
    roiMean_raw(i, :) = sgolayfilt(roiMean_raw(i, :), 3, smo_frm);
    df_ratio(i, :) = sgolayfilt(df_ratio(i, :), 3, smo_frm);
end

%% 3. Visualization: Mean trace 
disp('Plotting Filtered GCaMP/TdTomato...');
mean_raw = mean(roiMean_raw, 1)';
sem_raw  = std(roiMean_raw, 0, 1)' / sqrt(numCells);
mean_df = mean(df_ratio, 1)';
sem_df  = std(df_ratio, 0, 1)' / sqrt(numCells);

% Define drug timing bounds once to use across all plots
t_start = timevec(marks(1));
t_end = timevec(end); 
figure(1),clf; set(gcf, 'Position', [100, 100, 1200, 800]);
% --- Mean Trace with SEM ---
subplot(2,1,1); hold on;
fill([timevec; flipud(timevec)], [mean_raw+sem_raw; flipud(mean_raw-sem_raw)], ...
     [0.8 0.95 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_raw, 'Color', [0 0.6 0], 'LineWidth', 2);
ylabel('Mean GCaMP/TdTomato'); xlabel('Time (min.)'); 
title({'Mean Filtered GCaMP/TdTomato', ''});
box off;
y_max = max(mean_raw + sem_raw); y_min = min(mean_raw - sem_raw);
y_range = max(y_max - y_min, 1e-3);
bar_y = y_max + 0.15 * y_range; text_y = y_max + 0.25 * y_range;
set(gca, 'xlim', timevec([1 end]), 'ylim', [y_min - 0.1*y_range, text_y + 0.1*y_range], 'TickDir', 'none', 'FontSize', 12);
line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
% --- Mean Trace with SEM ---
subplot(2,1,2); hold on;
fill([timevec; flipud(timevec)], [mean_df+sem_df; flipud(mean_df-sem_df)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
plot(timevec, mean_df, 'r', 'LineWidth', 2);
ylabel('Mean \DeltaF/F_0'); xlabel('Time (min.)'); 
title({'Mean Filtered \DeltaF/F_0 Trace with SEM', ''});
box off;

valid_idx = ~isnan(mean_df) & ~isnan(sem_df);
if any(valid_idx)
    y_max = max(mean_df(valid_idx) + sem_df(valid_idx)); y_min = min(mean_df(valid_idx) - sem_df(valid_idx));
else
    y_max = 1; y_min = -0.5;
end
y_range = max(y_max - y_min, 1e-3); 
bar_y = y_max + 0.15 * y_range; text_y = y_max + 0.25 * y_range;
plot_ymin = y_min - 0.1 * y_range; plot_ymax = text_y + 0.1 * y_range;
set(gca, 'xlim', timevec([1 end]), 'ylim', [plot_ymin, plot_ymax], 'TickDir', 'none', 'FontSize', 12);

line([t_start t_end], [bar_y bar_y], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y + (0.02 * y_range), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');

%% 4. Visualization: Normalized GCaMP/TdTomato  
figure(2); clf; set(gcf, 'Position', [100, 100, 1200, 800]);
% --- Heatmap ---
subplot(2,1,1);
surf(timevec, linspace(0.5,0.5+numCells,numCells), roiMean_raw, 'LineStyle', 'none'); view(2);
axis tight; set(gca, 'ylim', [0.5 numCells+0.5], 'xlim', timevec([1 end]), 'YDir', 'reverse');
xlabel('Time (min.)'); ylabel('ROI number'); colorbar; 
title({'Heatmap: Filtered GCaMP/TdTomato', ''});
box off;
xline(t_start, '--w', 'LineWidth', 1); 
text(t_start, 0, evnt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');

% --- All Traces ---
subplot(2,1,2);
plot(timevec, roiMean_raw); axis tight;
xlabel('Time (min.)'); ylabel('GCaMP/TdTomato'); 
title({'All Traces: Filtered GCaMP/TdTomato', ''});
box off;
y_max_all = max(roiMean_raw(:)); y_min_all = min(roiMean_raw(:));
y_range_all = max(y_max_all - y_min_all, 1e-3);
bar_y_all = y_max_all + 0.15 * y_range_all; text_y_all = y_max_all + 0.25 * y_range_all;
set(gca, 'ylim', [y_min_all - 0.1*y_range_all, text_y_all + 0.1*y_range_all]);
line([t_start t_end], [bar_y_all bar_y_all], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y_all + (0.02 * y_range_all), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');

%% 5. Visualization: dF/F0
disp('Plotting Filtered dF/F0...');


figure(3); clf; set(gcf, 'Position', [200, 200, 1200, 800]);

% --- Heatmap ---
subplot(2,1,1);
surf(timevec, linspace(0.5,0.5+numCells,numCells), df_ratio, 'LineStyle', 'none'); view(2);
axis tight; set(gca, 'ylim', [0.5 numCells+0.5], 'xlim', timevec([1 end]), 'YDir', 'reverse');
colorbar; 
title({'Heatmap: Filtered \DeltaF/F_0', ''});
xlabel('Time (min.)'); ylabel('ROI number');
box off;
xline(t_start, '--w', 'LineWidth', 2);
text(t_start, 0.2, evnt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');

% --- All Traces ---
subplot(2,1,2);
plot(timevec, df_ratio); axis tight;
xlabel('Time (min.)'); ylabel('\DeltaF/F_0'); 
title({'All Traces: Filtered \DeltaF/F_0', ''});
box off;
y_max_all = max(df_ratio(:)); y_min_all = min(df_ratio(:));
y_range_all = max(y_max_all - y_min_all, 1e-3);
bar_y_all = y_max_all + 0.15 * y_range_all; text_y_all = y_max_all + 0.25 * y_range_all;
set(gca, 'ylim', [y_min_all - 0.1*y_range_all, text_y_all + 0.1*y_range_all]);
line([t_start t_end], [bar_y_all bar_y_all], 'Color', 'k', 'LineWidth', 1);
text(t_start + (t_end - t_start)/2, bar_y_all + (0.02 * y_range_all), evnt, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);


set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');

%% 6. (Optional) All Cell Time Tracing
disp('Generating individual traces for all ROIs...');
figure(100); clf; set(gcf, 'Position', [0, 50, 960, 945]);
t = tiledlayout(ceil(numCells/2),2,'TileSpacing','compact','Padding','compact');
title(t,'Normalization of all ROIs');
for i=1:numCells
    nexttile
    plot(timevec, roiMean_raw(i,:), '-k', 'LineWidth', 1);
    
    
    if  i >= numCells-1
    xlabel('Time (min.)');
    end
    ylabel(sprintf('ROI %d \nGCamP/TdTomato',i)); 
    set(gca, 'xlim', timevec([1 end]), 'ylim', [0 max(roiMean_raw,[],'all')], 'TickDir', 'none');
    xline(timevec(marks(1)), '--b', evnt, 'FontWeight', 'bold', 'LineWidth', 1, 'LabelHorizontalAlignment', 'center');
    box off;
end


set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');


%% 7. (Optional) All Cell Time Tracing
disp('Generating individual traces for all ROIs...');
figure(200); clf; set(gcf, 'Position', [960, 50, 960, 945]);
t = tiledlayout(ceil(numCells/2),2,'TileSpacing','compact','Padding','compact')
title(t,'Standarlized fluorescence of all ROIs');
for i=1:numCells
    nexttile
    
    plot(timevec, df_ratio(i,:), '-k', 'LineWidth', 1);
    
    
    if  i >= numCells-1
    xlabel('Time (min.)');
    end
    ylabel(sprintf('ROI %d \n\\DeltaF/F_0',i)); 
    set(gca, 'xlim', timevec([1 end]), 'ylim', [-0.8 10], 'TickDir', 'none');
    xline(timevec(marks(1)), '--b', evnt, 'FontWeight', 'bold', 'LineWidth', 1, 'LabelHorizontalAlignment', 'center');
box off;
end


set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');

%% 8. (Optional) Single Cell Highlight Tracing
j = 4; % Define specific cell to highlight
if j <= numCells
    figure(200+j); clf;
    set(gcf,'Units', 'inches', 'Position', [1 1 15 5]);
    plot(timevec, df_ratio(j,:), '-k', 'LineWidth', 1);
    
    ylabel('\DeltaF/F_0'); xlabel('Time (min.)');
    title({'', sprintf('Fluorescence Highlight: ROI %d (Filtered)', j), ''});
    
    set(gca,'xlim', timevec([1 end]), 'ylim', [-1.5 12], 'YTick', (0:2:10), 'TickDir', 'none');
    box off;
    
    line([t_start t_end], [10.5 10.5], 'color', 'k', 'Linewidth', 2);
    text(t_start + 0.5*(t_end - t_start), 11.5, evnt, 'HorizontalAlignment', 'center');
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
end

%% 9. Autonomous Peak Analysis (Dynamic Tuning)
disp('Running Autonomous Peak Analysis on Filtered Data...');

hdr_pks = cell(numCells,1);
hdr_locs = cell(numCells,1);

for i = 1:numCells
    % Signal is already smoothed globally in Section 2. 
    smooth_signal = df_ratio(i, :);
    signal_post = smooth_signal(marks(1):end); 
    time_post = timevec(marks(1):end);
    
    % Evaluate pre-injection baseline stability
    baseline_check_idx = max(1, marks(1)-20);
    baseline_signal = smooth_signal(baseline_check_idx : marks(1)-1);
    
    if mean(baseline_signal) < 1 
        % Median Absolute Deviation (MAD)
        baseline_median = median(baseline_signal);
        robust_sigma = median(abs(baseline_signal - baseline_median)) / 0.6745;
        robust_sigma = max(robust_sigma, 0.01); 
        
        dyn_prom = max(0.8, 1.5 * robust_sigma); 
        dyn_height = max(0.3, 1 * robust_sigma); 
        dt_min = mean(diff(timevec)); 
        dyn_width = dt_min * 30; 
        dtance = 1;
        [pks_post, locs_post] = findpeaks(signal_post, time_post, ...
            'MinPeakProminence', dyn_prom, 'MinPeakHeight', dyn_height, ...
            'MinPeakWidth', dyn_width,'MinPeakDistance',dtance);
        
        figure(300+i); clf; set(gcf, 'Position', [150, 150, 1000, 600]);
        
        % Subplot 1: Filtered Raw Check
        subplot(2,1,1); hold on;
        plot(time_post, roiMean_raw(i, marks(1):end), 'k', 'LineWidth', 1.5); 
        set(gca, 'xlim', timevec([marks(1) end]), 'TickDir', 'none');
        title(sprintf('ROI %d: Filtered Raw Signal', i)); box off;
        
        % Subplot 2: dF/F0 Peak Detection
        subplot(2,1,2); hold on;
        plot(time_post, signal_post, 'r', 'LineWidth', 1.5); 
        
        if ~isempty(locs_post)
            plot(locs_post, pks_post + 0.35, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineStyle', 'none');
        end
        
        yline(dyn_height, '--', 'Color', [0.5 0.5 0.5], 'Label', sprintf('Threshold: %.2f', dyn_height));
        set(gca, 'xlim', timevec([marks(1) end]), 'ylim', [-1 max(10, max(signal_post)+2)], 'TickDir', 'none');
        title(sprintf('ROI %d: Autonomous Peak Detection (Prominence: %.2f)', i, dyn_prom)); box off;
        
        set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
        
        hdr_pks{i} = pks_post; hdr_locs{i} = locs_post;
    else
        hdr_pks{i} = []; hdr_locs{i} = [];
    end
end

%% Aggregate Peak Statistics
disp('Calculating Peak Statistics...');
peak_counts = cellfun(@numel, hdr_pks);
peak_counts(peak_counts == 0) = NaN; 
mean_oscil = mean(peak_counts, 1, 'omitnan');

peak_ave = cellfun(@mean, hdr_pks);
ave_amp = mean(peak_ave, 1, 'omitnan');

act_cell = sum(peak_counts > 0, 'omitnan');
per_act_cell = (act_cell / numCells) * 100;

sem_amp  = std(peak_ave, 0, 1, 'omitnan') / sqrt(act_cell);
sem_oscil  = std(peak_counts, 0, 1, 'omitnan') / sqrt(act_cell);


%% 10. Analyze 1st Peaks vs Remaining Peaks & Comprehensive AUC (INTERACTIVE)
disp('Extracting Peak Dynamics. Adjust Rectangles interactively for AUC accuracy...');
st_hdr_locs = zeros(numCells,1); st_hdr_pks = zeros(numCells,1);
st_int = zeros(numCells,1); st_loss = zeros(numCells,1);

% --- Initialize New AUC Arrays ---
auc_1st_peak   = nan(numCells, 1);
auc_5min_post  = nan(numCells, 1);
auc_all_peaks  = nan(numCells, 1);

for k = 1:numCells
    if ~isempty(hdr_locs{k}) && ~isnan(hdr_locs{k}(1))
        
        st_hdr_locs(k) = hdr_locs{k}(1);
        st_hdr_pks(k)  = hdr_pks{k}(1);
        [~, peak_idx] = min(abs(timevec - st_hdr_locs(k)));
        
        % AUTODETECT RISE DYNAMICS
        pre_peak_signal = df_ratio(k, marks(1):peak_idx);
        st_int_bline = prctile(pre_peak_signal, 15);
        idx_rise_start = find(pre_peak_signal <= st_int_bline, 1, 'last');
        if isempty(idx_rise_start), idx_rise_start = 1; end
        st_int(k) = timevec(marks(1) + idx_rise_start - 1);
        
        % AUTODETECT DECAY DYNAMICS
        post_peak_signal = df_ratio(k, peak_idx:end);
        if isempty(find(post_peak_signal <= 0.3, 1))
            st_loss(k) = timevec(end); 
            decay_baseline = 0.3;      
        else
            if length(hdr_pks{k}) == 1
                decay_baseline = prctile(post_peak_signal, 15);
                idx_decay_end = find(post_peak_signal <= decay_baseline, 1, 'first');
                if isempty(idx_decay_end), idx_decay_end = length(post_peak_signal); end
                st_loss(k) = timevec(peak_idx + idx_decay_end - 1);
            else
                [~, peak2_idx] = min(abs(timevec - hdr_locs{k}(2)));
                st_nd_dy = df_ratio(k, peak_idx:peak2_idx);
                decay_baseline = prctile(st_nd_dy, 18);
                idx_decay_end = find(st_nd_dy <= decay_baseline, 1, 'first');
                if isempty(idx_decay_end), idx_decay_end = length(st_nd_dy); end
                st_loss(k) = timevec(peak_idx + idx_decay_end - 1);
            end
        end
        
        % ==========================================
        % QC GRAPH & INTERACTIVE ADJUSTMENT
        % ==========================================
        figQC = figure(400+k); clf; set(figQC, 'Position', [200, 200, 1000, 500]); hold on;
        plot(timevec, df_ratio(k,:), '-k', 'LineWidth', 1.5);
        plot(hdr_locs{k}, hdr_pks{k} + 0.3, 'rv', 'MarkerFaceColor', 'r', 'LineStyle', 'none');
        
        % Interactive DrawRectangles
        rise_w = max(0.01, st_hdr_locs(k) - st_int(k));
        rise_h = max(0.01, st_hdr_pks(k) + 0.35 - st_int_bline);
        rect_rise = drawrectangle('Position', [st_int(k), st_int_bline, rise_w, rise_h], ...
                                  'Color', 'y', 'FaceAlpha', 0.5);
            
        decay_w = max(0.01, st_loss(k) - st_hdr_locs(k));
        decay_h = max(0.01, st_hdr_pks(k) + 0.35 - decay_baseline);
        rect_decay = drawrectangle('Position', [st_hdr_locs(k), decay_baseline, decay_w, decay_h], ...
                                   'Color', 'm', 'FaceAlpha', 0.5);
        
        ylabel('\DeltaF/F_0'); xlabel('Time (min.)');
        title({'', sprintf('ROI %d: Adjust Left edge of Yellow for Start, Right edge of Magenta for End', k), ''});
        y_max = max(10, max(df_ratio(k,:)) + 1);
        set(gca, 'xlim', timevec([1 end]), 'ylim', [-0.8 y_max], 'TickDir', 'none'); box off;
        xline(timevec(marks(1)), '--b', evnt, 'FontWeight', 'bold', 'LineWidth', 1.5);
        set(findall(figQC, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
        
        % Add Confirm Button to Pause execution
        btn_confirm = uicontrol('Parent', figQC, 'Style', 'pushbutton', 'String', 'Confirm & Next', ...
                                'Position', [20 20 120 40], 'FontSize', 10, 'FontWeight', 'bold', ...
                                'Callback', 'uiresume(gcbf)');
        
        disp(['Waiting for user to adjust ROI ', num2str(k), ' boundaries...']);
        uiwait(figQC);
        
        % Extract New User-Adjusted Boundaries
        if isvalid(figQC)
            if isvalid(rect_rise)
                pos_rise = rect_rise.Position;
                st_int(k) = pos_rise(1); % New Start Time (Left Edge)
            end
            if isvalid(rect_decay)
                pos_decay = rect_decay.Position;
                st_loss(k) = pos_decay(1) + pos_decay(3); % New End Time (Left Edge + Width)
            end
            delete(btn_confirm); % Cleanup button so graph is clean
            hold off;
        end

        % ==========================================
        % --- CALCULATE AUC METRICS AFTER CONFIRM ---
        % ==========================================
        [~, idx_start_1st] = min(abs(timevec - st_int(k)));
        [~, idx_end_1st]   = min(abs(timevec - st_loss(k)));
        
        % 1. AUC of 1st Peak (Start to Loss)
        time_1st = timevec(idx_start_1st:idx_end_1st);
        sig_1st  = df_ratio(k, idx_start_1st:idx_end_1st);
        auc_1st_peak(k) = trapz(time_1st, max(0, sig_1st - st_int_bline));
        
        % 2. AUC within 5 minutes after induced (detected by adjusted st_int)
        [~, idx_5min] = min(abs(timevec - (st_int(k) + 5)));
        idx_5min = min(idx_5min, length(timevec)); 
        time_5min = timevec(idx_start_1st:idx_5min);
        sig_5min  = df_ratio(k, idx_start_1st:idx_5min);
        auc_5min_post(k) = trapz(time_5min, max(0, sig_5min - st_int_bline));
    end
    
    % 3. AUC of ALL peaks at all time
    time_post = timevec(marks(1):end);
    sig_post  = df_ratio(k, marks(1):end);
    auc_all_peaks(k) = trapz(time_post, max(0, sig_post));
end

st_hdr_pks(st_hdr_pks == 0) = NaN; st_hdr_locs(st_hdr_locs == 0) = NaN;
st_int(st_int == 0) = NaN; st_loss(st_loss == 0) = NaN;

st_amp = mean(st_hdr_pks, 'omitnan'); st_loc = mean(st_hdr_locs, 'omitnan');
st_raise = st_hdr_locs - st_int;
st_decay = st_loss - st_hdr_locs;
st_dur = st_loss - st_int;

% --- Aggregate Mean AUCs ---
mean_auc_1st  = mean(auc_1st_peak, 'omitnan');
mean_auc_5min = mean(auc_5min_post, 'omitnan');
mean_auc_all  = mean(auc_all_peaks, 'omitnan');

disp('--- Comprehensive AUC Metrics ---');
disp(['-> Mean 1st Peak AUC: ', num2str(mean_auc_1st)]);
disp(['-> Mean AUC (5 mins post-induction): ', num2str(mean_auc_5min)]);
disp(['-> Mean Total AUC (All peaks): ', num2str(mean_auc_all)]);

% Remaining Peaks
rest_locs = cell(numCells,1); rest_pks = cell(numCells, 1);
for l = 1:numCells
    if length(hdr_locs{l}) > 1
        rest_locs{l} = hdr_locs{l}(2:end); rest_pks{l}  = hdr_pks{l}(2:end);
    end
end
rest_cnt_osc = cellfun(@numel, rest_pks); rest_cnt_osc(rest_cnt_osc == 0) = NaN;
rest_cnt_amp = cellfun(@mean, rest_pks);

%% 11. Intracellular Calcium Dynamics (Thapsigargin, SOCE, VOCC)
disp('--- ER Depletion (Thapsigargin) AUC ---');
figure(3); 
disp('Click TWO points on the Mean dF/F0 graph to define the ER depletion window...');
[x_er, ~] = ginput(2); x_er = sort(x_er); 
[~, ER_start] = min(abs(timevec - x_er(1)));
[~, ER_end]   = min(abs(timevec - x_er(2)));

df_thap = df_ratio(:, ER_start:ER_end);
timevec_thap = timevec(ER_start:ER_end);

figure(14); clf; set(gcf, 'Position', [250, 250, 800, 400]);
plot(timevec_thap, df_thap);
title({'ER Depletion Window (Thapsigargin)', ''}); xlabel('Time (min.)'); ylabel('\DeltaF/F_0');
set(gca, 'xlim', timevec([ER_start ER_end]), 'TickDir', 'none'); box off;

auc_er = nan(numCells, 1); 
for i = 1:numCells
    [~, tg_locs] = findpeaks(df_thap(i, :), 'MinPeakProminence', 1);
    if ~isempty(tg_locs)
        tg_int = find(df_thap(i, 1:tg_locs(1)) <= 0, 1, 'last');
        if isempty(tg_int), tg_int = 1; end
        auc_er(i) = trapz(timevec_thap(tg_int:end), df_thap(i, tg_int:end));
    end    
end

figure(15); clf; bar(auc_er, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
title({'Thapsigargin AUC per ROI', ''}); xlabel('ROI Number'); ylabel('AUC (\DeltaF/F_0 \times min)'); box off;
auc_er_mean = mean(auc_er, 'omitnan');
disp(['-> Mean Thapsigargin AUC: ', num2str(auc_er_mean)]);
set(findall(figure(14), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
set(findall(figure(15), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');

%%  12.--- SOCE Slope (Autonomous Max Slope Detection) ---
disp('--- SOCE Slope ---');
figure(3); 
disp('Click ONE point on the Mean dF/F0 graph to define the START of the SOCE induction...');
[x_soce, ~] = ginput(1);

% Find the closest frame index for the clicked start time
[~, SOCE_start] = min(abs(timevec - x_soce));

% Define a search window (e.g., 5 minutes) to look for the maximum slope
% (Adjust this value if your rise phase takes longer to reach the plateau)
search_duration = 5; 
[~, SOCE_end] = min(abs(timevec - (x_soce + search_duration)));

timevec_search = timevec(SOCE_start:SOCE_end);
dt = mean(diff(timevec_search)); % Average time step in minutes

soce_slope = zeros(numCells, 1);
max_slope_idx = zeros(numCells, 1);

% Prepare QC Figure
figure(18); clf; set(gcf, 'Position', [350, 350, 1000, 600]); hold on;
title({'Autonomous SOCE Slope Detection', 'Black lines indicate the point of maximum steepness for each ROI'}); 
xlabel('Time (min.)'); ylabel('\DeltaF/F_0');
colors = lines(numCells);

for i = 1:numCells
    signal_search = df_ratio(i, SOCE_start:SOCE_end);
    
    % Calculate the derivative (rate of change: dF/dt)
    % Using 'gradient' provides a slightly smoother derivative than 'diff'
    dF_dt = gradient(signal_search) ./ dt;
    
    % Find the maximum slope (steepest part of the rise)
    [max_val, max_idx] = max(dF_dt);
    
    soce_slope(i) = max_val;
    max_slope_idx(i) = max_idx;
    
    % --- Plotting for QC ---
    plot(timevec_search, signal_search, 'Color', [colors(i,:) 0.6], 'LineWidth', 1.5);
    
    % Draw a tangent line at the exact point of maximum slope for visual confirmation
    t_max = timevec_search(max_idx);
    y_max = signal_search(max_idx);
    
    % Define a small line segment to represent the tangent visually
    line_length = 0.1; % minutes
    t_tangent = [t_max - line_length/2, t_max + line_length/2];
    y_tangent = y_max + max_val * (t_tangent - t_max);
    
    plot(t_tangent, y_tangent, 'Color', 'k', 'LineWidth', 1.75);
    plot(t_max, y_max, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4); 
end

set(gca, 'xlim', timevec([SOCE_start SOCE_end]), 'TickDir', 'none'); box off;
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
hold off;

soce_slope_mean = mean(soce_slope, 'omitnan');
disp(['-> Mean Autonomous SOCE Slope: ', num2str(soce_slope_mean)]);

% Generate Bar Chart
figure(19); clf; bar(soce_slope, 'FaceColor', [0.4 0.8 0.4], 'EdgeColor', 'none');
title({'Maximum SOCE Slope per ROI', ''}); 
xlabel('ROI Number'); ylabel('Max Slope (\DeltaF/F_0 / min)'); box off;
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
%% --- VOCC + SOCE AUC (Autonomous 3-Metrics) ---
disp('--- VOCC + SOCE AUC ---');
disp('Autonomously detecting true onset and saturation points using the previous click...');

% Preallocate arrays for the three different AUC calculations
auc_selected_3min = nan(numCells, 1);
auc_onset_3min    = nan(numCells, 1);
auc_onset_sat     = nan(numCells, 1);

% Arrays to store the autonomously found indices for QC graphing
all_onset_idx = ones(numCells, 1);
all_sat_idx   = ones(numCells, 1);

% Metric 1 End Index: 3 minutes after the globally selected point
[~, selected_end_idx] = min(abs(timevec - (x_soce + 3)));

for i = 1:numCells
    % 1. DEFINE BASELINE
    % Use 20 frames right before the globally clicked start time
    if SOCE_start > 1
        baseline_idx_start = max(1, SOCE_start - 20);
        baseline_val = mean(df_ratio(i, baseline_idx_start : SOCE_start - 1));
    else
        baseline_val = df_ratio(i, 1);
    end
    
    % --- METRIC 1: AUC in 3 mins from selected point ---
    time_1 = timevec(SOCE_start : selected_end_idx);
    sig_1  = df_ratio(i, SOCE_start : selected_end_idx) - baseline_val;
    auc_selected_3min(i) = trapz(time_1, max(0, sig_1)); % max(0) prevents negative area from noise
    
    % 2. AUTODETECT ACTUAL RISING TIMEPOINT (ONSET)
    % Convert relative max_slope_idx from the previous section to an absolute index
    abs_max_idx = SOCE_start + max_slope_idx(i) - 1;
    
    % Step backward from the steepest slope until the signal hits the local baseline
    idx_onset = abs_max_idx;
    while idx_onset > SOCE_start && df_ratio(i, idx_onset) > baseline_val+0.05
        idx_onset = idx_onset - 1;
    end
    all_onset_idx(i) = idx_onset;
    
    % --- METRIC 2: AUC in 3 mins from ACTUAL rising point ---
    [~, idx_onset_3min] = min(abs(timevec - (timevec(idx_onset) + 3)));
    idx_onset_3min = min(idx_onset_3min, length(timevec)); % Failsafe
    
    time_2 = timevec(idx_onset : idx_onset_3min);
    sig_2  = df_ratio(i, idx_onset : idx_onset_3min) - baseline_val;
    auc_onset_3min(i) = trapz(time_2, max(0, sig_2));
    
    % 3. AUTODETECT FIRST SATURATE TIMEPOINT (PLATEAU)
    % Look forward from the steepest slope until the derivative drops to <= 0 (stops rising)
    % Cap the search window at 2 minutes after the max slope to prevent runaway
    search_sat_end = min(length(timevec), abs_max_idx + round(2 / mean(diff(timevec)))); 
    dF_dt_forward = gradient(df_ratio(i, abs_max_idx : search_sat_end));
    
    stop_rising = find(dF_dt_forward <= 0, 1, 'first');
    if isempty(stop_rising)
        idx_sat = search_sat_end; % Fallback if it never perfectly plateaus
    else
        idx_sat = abs_max_idx + stop_rising - 1;
    end
    all_sat_idx(i) = idx_sat;
    
    % --- METRIC 3: AUC between Actual Rising and Saturation ---
    time_3 = timevec(idx_onset : idx_sat);
    sig_3  = df_ratio(i, idx_onset : idx_sat) - baseline_val;
    auc_onset_sat(i) = trapz(time_3, max(0, sig_3));
end

%% --- Visualizations for VOCC + SOCE Metrics ---

% QC FIGURE: Traces with Detected Onset & Saturation Markers
figure(16); clf; set(gcf, 'Position', [200, 300, 1000, 500]); hold on;
title({'VOCC+SOCE Dynamics (Baseline Corrected)', 'Triangles = Detected Onset (Rise) | Squares = Detected Saturation (Plateau)'}); 
xlabel('Time (min.)'); ylabel('\DeltaF/F_0 (Baseline Subtracted)');

colors = lines(numCells);
plot_end = min(length(timevec), max(all_sat_idx) + round(2 / mean(diff(timevec)))); % 2 mins past last saturation

for i = 1:numCells
    baseline_val = mean(df_ratio(i, max(1, SOCE_start - 20) : max(1, SOCE_start - 1)));
    
    % Plot normalized trace
    sig_plot = df_ratio(i, SOCE_start:plot_end) - baseline_val;
    plot(timevec(SOCE_start:plot_end), sig_plot, 'Color', [colors(i,:) 0.5], 'LineWidth', 1);
    
    % Plot Onset Marker
    t_on = timevec(all_onset_idx(i));
    y_on = df_ratio(i, all_onset_idx(i)) - baseline_val;
    plot(t_on, y_on, '^', 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1);
    
    % Plot Saturation Marker
    t_sat = timevec(all_sat_idx(i));
    y_sat = df_ratio(i, all_sat_idx(i)) - baseline_val;
    plot(t_sat, y_sat, 's', 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1);
end
set(gca, 'xlim', timevec([SOCE_start, plot_end]), 'TickDir', 'none'); box off;

% GROUPED BAR CHART: Comparing the 3 AUC Metrics
figure(17); clf; set(gcf, 'Position', [400, 300, 900, 500]);
bar_data = [auc_selected_3min, auc_onset_3min, auc_onset_sat];
b = bar(bar_data, 'grouped', 'EdgeColor', 'none');

% Style the grouped bars
b(1).FaceColor = [0.7 0.7 0.7]; % Gray: Global 3 Min
b(2).FaceColor = [0.2 0.6 0.8]; % Blue: Actual Rise 3 Min
b(3).FaceColor = [0.8 0.3 0.3]; % Red: Rise to Saturation

legend({'3 Min from Selected Point', '3 Min from Actual Rise', 'Between Rise & Saturation'}, 'Location', 'northwest');
title({'VOCC+SOCE AUC Comparisons per ROI', ''}); 
xlabel('ROI Number'); ylabel('AUC (\DeltaF/F_0 \times min)'); box off;

% Apply Formatting
set(findall(figure(16), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
set(findall(figure(17), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');

% Console Summary
disp(['-> Mean AUC (3 mins from Global Click): ', num2str(mean(auc_selected_3min, 'omitnan'))]);
disp(['-> Mean AUC (3 mins from Actual Rise):  ', num2str(mean(auc_onset_3min, 'omitnan'))]);
disp(['-> Mean AUC (Actual Rise to Saturation): ', num2str(mean(auc_onset_sat, 'omitnan'))]);


%% 13. Find the trace closest to the MEDIAN
% Calculate the median value across all ROIs at each time point
% Assuming your matrix is named 'traces' 
% (Rows = ROIs, Columns = Time points)
traces = df_ratio;
theoreticalMedianTrace = median(traces, 1);

% Calculate the sum of squared distances from each ROI to the theoretical median
distancesToMedian = sum((traces - theoreticalMedianTrace).^2, 2);

% Find the row index with the minimum distance
[~, closestMedianIdx] = min(distancesToMedian);

% Extract the actual ROI trace closest to the median
medianTraceROI = traces(closestMedianIdx, :);

%% 14. Find the trace closest to the MODE
% Calculate the mode across all ROIs at each time point
% Note: If your data is highly continuous, you may need to round it slightly 
% first (e.g., mode(round(traces, 2), 1)) for 'mode' to find meaningful overlaps.
theoreticalModeTrace = mode(traces, 1);

% Calculate the sum of squared distances from each ROI to the theoretical mode
distancesToMode = sum((traces - theoreticalModeTrace).^2, 2);

% Find the row index with the minimum distance
[~, closestModeIdx] = min(distancesToMode);

% Extract the actual ROI trace closest to the mode
modeTraceROI = traces(closestModeIdx, :);
fprintf('ROI closest to the Median is Row: %d\n', closestMedianIdx);
fprintf('ROI closest to the Mode is Row: %d\n', closestModeIdx);
% Display the results to verify
% Optional: Plot them to visualize
figure;
tiledlayout(2,1)
nexttile
plot(medianTraceROI, 'LineWidth', 1.5, 'DisplayName', 'Closest to Median');
legend('Location', 'best');
xlabel('Time points');
ylabel('\DeltaF/F_0');
title('Representative Traces');
nexttile
plot(modeTraceROI, 'LineWidth', 1.5, 'DisplayName', 'Closest to Mode');
legend('Location', 'best');
xlabel('Time points');
ylabel('\DeltaF/F_0');
title('Representative Traces');
%% 15. Memory Cleanup and Save
disp('Purging temporary buffered variables to save memory...');

% 1. Clean up Initial Loading & Base Metric Loop Indices
clear dataTbl timevec_raw timevec_sec timeColIdx fluorColsIdx;
clear c cx cy rad timei i j k l;
%clear t_start t_end y_max y_min y_range bar_y text_y plot_ymin plot_ymax;
%clear y_max_all y_min_all y_range_all bar_y_all text_y_all valid_idx;

% 2. Clean up Peak Detection & Savitzky-Golay Buffers
clear baseline_start baseline_end F0_raw smooth_signal signal_post time_post;
clear baseline_check_idx baseline_signal baseline_median robust_sigma;
clear dyn_prom dyn_height dt_min dyn_width pks_post locs_post;
clear peak_idx pre_peak_signal st_int_bline idx_rise_start post_peak_signal;
clear decay_baseline idx_decay_end st_nd_dy peak2_idx rise_w rise_h decay_w decay_h;
clear rect_rise rect_decay pos_rise pos_decay btn_confirm figQC;

% 3. Clean up General AUC & Thapsigargin Variables
clear x_er ER_start ER_end df_thap timevec_thap tg_locs tg_int;
clear time_1st sig_1st time_5min sig_5min sig_post idx_start_1st idx_end_1st idx_5min;

% 4. Clean up SOCE Slope Detection Buffers
clear x_soce SOCE_start SOCE_end search_duration dt signal_search dF_dt;
clear max_slope_idx max_val max_idx t_max y_max t_tangent y_tangent line_length colors;

% 5. Clean up VOCC+SOCE Autonomous 3-Metric Buffers
clear selected_end_idx baseline_idx_start baseline_val abs_max_idx;
clear time_1 sig_1 idx_onset idx_onset_3min time_2 sig_2;
clear search_sat_end dF_dt_forward stop_rising idx_sat time_3 sig_3;
clear plot_end sig_plot t_on y_on t_sat y_sat bar_data b;
clear all_onset_idx all_sat_idx;

disp('Saving lightweight analysis workspace...');
[~, base_name, ~] = fileparts(filename);
datFileName = [base_name, '_AnalyzedData.mat'];
save(fullfile(filepath, datFileName), '-v7.3');

disp(['SUCCESS! Analysis saved as: ', datFileName]);
disp('Pipeline Complete.');


%%
close all, clear, clc
