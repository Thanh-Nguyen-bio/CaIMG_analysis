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
%evnt = {'20 \muM Thap'};
evnt = {'Akh 10 ng'};
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

%% 4. Visualization: Normalized GCaMP/TdTomato  
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

%% 10. Aggregate Peak Statistics
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


%% 11. Analyze 1st Peaks vs Remaining Peaks & Comprehensive AUC (INTERACTIVE)
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

%% 12. Intracellular Calcium Dynamics (Thapsigargin, SOCE, VOCC)
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
%% 13.SOCE Slope
disp('--- SOCE Slope ---');
figure(3); disp('Click TWO points on the Mean dF/F0 graph to define the SOCE slope window...');
[x_soce, ~] = ginput(2); x_soce = sort(x_soce);
[~, SOCE_start] = min(abs(timevec - x_soce(1)));
[~, SOCE_end]   = min(abs(timevec - x_soce(2)));
df_slpe = df_ratio(:, SOCE_end) - df_ratio(:, SOCE_start);
timevec_slpe = timevec(SOCE_end) - timevec(SOCE_start);
soce_slope = df_slpe ./ timevec_slpe; 
soce_slope_mean = mean(soce_slope, 'omitnan');
disp(['-> Mean SOCE Slope: ', num2str(soce_slope_mean)]);

disp('--- VOCC + SOCE AUC (3-Minute Window) ---');
figure(3); 
disp('Click ONE point on the Mean dF/F0 graph to define the START of the VOCC+SOCE window...');
[x_voso, ~] = ginput(1); % Only requires 1 click now

% Find the closest frame index for the clicked start time
[~, bg_voso] = min(abs(timevec - x_voso));

% Calculate the target end time (Start + 3 minutes)
target_end_time = timevec(bg_voso) + 3;

% Find the closest frame index for the 3-minute end mark 
% (This automatically snaps to the last frame if 3 mins exceeds the recording length)
[~, end_voso] = min(abs(timevec - target_end_time));

fprintf('-> Mapped Frame Indices: Start = %d, End = %d (approx 3 mins)\n', bg_voso, end_voso);

% Baseline correction based on the 20 frames right before the window
baseline_voso_idx = max(1, bg_voso-20);
baseline_voso = mean(df_ratio(:, baseline_voso_idx:bg_voso-1), 2);
df_voso = df_ratio(:, bg_voso:end_voso) - baseline_voso;
timevec_voso = timevec(bg_voso:end_voso);

figure(16); clf; set(gcf, 'Position', [300, 300, 800, 400]);
plot(timevec_voso, df_voso);
title({'VOCC+SOCE Window (Baseline Corrected, 3 Min)', ''}); xlabel('Time (min.)'); ylabel('\DeltaF/F_0');
set(gca, 'xlim', timevec([bg_voso end_voso]), 'TickDir', 'none'); box off;

% Calculate true Area Under Curve using exact time (minutes)
auc_voso = trapz(timevec_voso, df_voso, 2);

figure(17); clf; bar(auc_voso, 'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none');
title({'VOCC+SOCE AUC per ROI (3 Min Window)', ''}); xlabel('ROI Number'); ylabel('AUC (\DeltaF/F_0 \times min)'); box off;

set(findall(figure(16), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');
set(findall(figure(17), '-property', 'FontName'), 'FontName', 'Arial', 'FontWeight', 'normal');


%% 14. Find the trace closest to the MEDIAN
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

%% 2. Find the trace closest to the MODE
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
%% 14. Memory Cleanup and Save
disp('Purging temporary variables to save memory...');
% Clean up tables, loop indices, string tokens, and graphing boundaries
clear dataTbl timevec_raw timevec_sec timeColIdx fluorColsIdx;
clear c cx cy rad timei i j k l;
%clear t_start t_end y_max y_min y_range bar_y text_y plot_ymin plot_ymax;
%clear y_max_all y_min_all y_range_all bar_y_all text_y_all valid_idx;
%clear baseline_start baseline_end F0_raw smooth_signal signal_post time_post;
%clear baseline_check_idx baseline_signal baseline_median robust_sigma;
clear dyn_prom dyn_height dt_min dyn_width pks_post locs_post;
clear peak_idx pre_peak_signal st_int_bline idx_rise_start post_peak_signal;
clear decay_baseline idx_decay_end st_nd_dy peak2_idx rise_w rise_h decay_w decay_h;
clear x_er x_soce x_voso ER_start ER_end SOCE_start SOCE_end bg_voso end_voso;
clear baseline_voso_idx baseline_voso tg_locs tg_int df_thap timevec_thap;
clear df_slpe timevec_slpe df_voso timevec_voso;

% Purge temporary AUC loop variables
clear time_1st sig_1st time_5min sig_5min sig_post idx_start_1st idx_end_1st idx_5min;
clear rect_rise rect_decay pos_rise pos_decay btn_confirm figQC;

disp('Saving lightweight analysis workspace...');
[~, base_name, ~] = fileparts(filename);
datFileName = [base_name, '_AnalyzedData.mat'];
save(fullfile(filepath, datFileName), '-v7.3');

disp(['SUCCESS! Analysis saved as: ', datFileName]);
disp('Pipeline Complete.');

%%
close all, clear, clc