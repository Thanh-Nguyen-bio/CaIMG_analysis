%% Comprehensive Calcium Imaging Analysis Pipeline (v4.1 - Final Cleaned Version)
% Author: AI-generated (Time unit correction applied)
% Date: October 3, 2025
%
% Key Update: Code modified to remove Z-score calculation, findpeaks, and all peak-related analysis/plotting.
% 1. Dynamic Frame Rate (Fs) and frame step (dt_min) calculation.
% 2. Motion Correction using Ch2 (TdTomato/Reference) fluorescence.
% 3. All time metrics and plots are labeled in minutes.
%% 1. Interactive File Selection and Initialization
clear, clc, close all;
% Use uigetfile to allow the user to select multiple CSV files
disp('Please select one or more CSV files to analyze...');
[filenames, pathname] = uigetfile('*.csv', 'Select Calcium Data Files', 'MultiSelect', 'on');
% Handle case where user cancels
if isequal(filenames, 0)
    disp('User canceled file selection.');
    return;
end
if ~iscell(filenames)
    filenames = {filenames};
end
% Pre-allocate a cell array to store results for each file
all_results = cell(1, length(filenames));
%% 2. Loop Through Each File for Analysis
for i_file = 1:length(filenames)
    filename = fullfile(pathname, filenames{i_file});
    disp(['Processing file: ' filename]);
    
    dataTbl = readtable(filename, 'VariableNamingRule', 'preserve');
    
    % --- Identify Time and Fluorescence Columns ---
    timeColIdx = find(contains(dataTbl.Properties.VariableNames, 'Time', 'IgnoreCase', true), 1);
    ch1ColIdx = find(contains(dataTbl.Properties.VariableNames, 'Ch1', 'IgnoreCase', true));
    ch2ColIdx = find(contains(dataTbl.Properties.VariableNames, 'Ch2', 'IgnoreCase', true));
    if isempty(timeColIdx) || isempty(ch1ColIdx) || isempty(ch2ColIdx) || length(ch1ColIdx) ~= length(ch2ColIdx)
        warning(['Skipping file ' filenames{i_file} ': Missing Time, Ch1, or Ch2 data, or unequal number of ROIs.']);
        continue;
    end
    
    % Extract data and transpose to [numROIs x numTimePoints]
    timevec = dataTbl{:, timeColIdx}; % Units: Minutes
    neuronts_ch1 = dataTbl{:, ch1ColIdx}'; % Signal (GCaMP)
    neuronts_ch2 = dataTbl{:, ch2ColIdx}'; % Reference (TdTomato/Motion)
    
    % --- NaN Handling ---
    nan_cols = any(isnan(neuronts_ch1), 1) | any(isnan(neuronts_ch2), 1);
    valid_cols = ~nan_cols;
    
    neuronts_ch1 = neuronts_ch1(:, valid_cols);
    neuronts_ch2 = neuronts_ch2(:, valid_cols);
    timevec = timevec(valid_cols);
    
    nROIs = size(neuronts_ch1, 1);
    nFrames = size(neuronts_ch1, 2);
    
    % --- Dynamic Frame Rate Calculation (Updated for Minutes) ---
    time_steps = diff(timevec);
    dt_min = mean(time_steps); % Time step in minutes
    
    % Calculate Fs in Hz (Frames per second) for informational display
    frame_rate_Hz = 1 / (dt_min * 60); 
    disp(['Calculated Frame Rate (Fs): ' num2str(frame_rate_Hz, '%.2f') ' Hz']);
    disp(['Time Step (dt): ' num2str(dt_min, '%.4f') ' minutes']);
    
    % --- Normalization and Motion Correction ---
    baseline_frames = min(20, nFrames);
    baseline_ch1 = mean(neuronts_ch1(:, 1:baseline_frames), 2);
    baseline_ch2 = mean(neuronts_ch2(:, 1:baseline_frames), 2);
    
    dFF_ch1 = bsxfun(@rdivide, neuronts_ch1, baseline_ch1) - 1;
    dFF_ch2 = bsxfun(@rdivide, neuronts_ch2, baseline_ch2) - 1;
    
    % Motion Correction: dFF_corrected = dFF_signal - dFF_reference
    dFF_corr = dFF_ch1 - dFF_ch2;
    
    % Pre-allocate results table for this file (Will be empty since peak analysis is removed)
    analyzed_data = table(); 
    
    %% 3 & 4. Peak Detection, Filtering, and Property Analysis (REMOVED)
    
    all_results{i_file} = analyzed_data; 
    
    %% 5A. Visualization - Summary Figure (Plots dFF_corr only)
    hFigSummary = figure('Name', ['Summary: ' filenames{i_file}], 'NumberTitle', 'off');
    set(hFigSummary, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    
    plot(timevec, dFF_corr);
    
    xlabel('Time (min)'); % Corrected label
    ylabel('Corrected \DeltaF/F_0');
    title(['All ROIs (Ch2 Corrected) | Fs: ' num2str(frame_rate_Hz, '%.2f') ' Hz']);
    
    set(gca,'xlim', timevec([1 end]));
    
    % Save the Summary figure
    [~, name] = fileparts(filenames{i_file});
    saveas(gcf, [pathname name '_summary_corrected.png']);
    disp('Summary Figure saved successfully.');
    
    %% 5B. Visualization - Random 5 ROI Traces (dF/F0 only)
    hFigRandom = figure('Name', ['Random Traces Comparison: ' filenames{i_file}], 'NumberTitle', 'off');
    set(hFigRandom, 'units', 'normalized', 'outerposition', [0.1 0.1 0.8 0.8]);
    num_random_plots = min(5, nROIs);
    rand_indices = randperm(nROIs, num_random_plots);
    
    % Tiled layout changed to 1 column 
    tiledlayout(num_random_plots, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 
    
    for idx = 1:num_random_plots
        roi_index = rand_indices(idx);
        
        % Plot 1: Motion-Corrected dF/F0
        nexttile;
        plot(timevec, dFF_corr(roi_index,:), 'b', 'LineWidth', 1.5);
        
        set(gca,'xlim', timevec([1 end]));
        
        title(['ROI ' num2str(roi_index) ' - Corrected \DeltaF/F_0'], 'FontSize', 8);
        if idx == num_random_plots
            xlabel('Time (min)'); % Corrected label
        end
        ylabel('\DeltaF/F_0');
        grid on;
    end
    
    sgtitle(['Randomly Selected ' num2str(num_random_plots) ' ROI Corrected \DeltaF/F_0 Traces'], 'FontWeight', 'bold');
    
    % Save the Random Traces figure
    saveas(hFigRandom, [pathname name '_random_traces.png']);
    disp('Random Traces Comparison Figure saved successfully.');

    %% 5C. Visualization - Heatmap of Corrected dF/F0 (NEW FEATURE)
    figure(11), clf; % Use figure(11) as requested
    
    % Plot the surface using dFF_corr (Corrected Signal)
    surf(timevec, linspace(0.5, nROIs+0.5, nROIs), dFF_corr, 'LineStyle', 'none');
    
    % Apply formatting settings
    set(gca, 'ylim', [-1.5 nROIs+2.5], 'FontSize', 16);
    set(gca, 'ytick', (1:5:nROIs));
    set(gca, 'xlim', timevec([1 end]), 'FontSize', 16);
    
    colorbar;
    set(gca, 'clim', [-1 10], 'fontsize', 16); % Set color limits (adjust as needed for dFF)
    set(gca, 'YDir', 'reverse');
    
    grid off;
    view(2); % Switch to 2D view (top-down) for a standard heatmap look

    xlabel('Time (min)');
    ylabel('ROI Number');
    title(['Heatmap of Corrected \DeltaF/F_0 | Fs: ' num2str(frame_rate_Hz, '%.2f') ' Hz']);

    % Save the Heatmap figure
    saveas(gcf, [pathname name '_heatmap.png']);
    disp('Heatmap Figure saved successfully.');

end
%% 6. Final Data Export
% Combine all results into a single table
final_results_table = vertcat(all_results{:});
% Export the final table to a single CSV file
if ~isempty(final_results_table)
    output_filename = fullfile(pathname, 'All_Analyzed_Peak_Data.csv');
    writetable(final_results_table, output_filename);
    disp(['All analyzed data saved to: ' output_filename]);
else
    disp('Peak analysis was disabled. No data was exported to CSV.');
end