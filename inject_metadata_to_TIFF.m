close all; clear; clc;

%% 1. Select the RAW source file (.czi)
disp('Select the RAW (.czi) file containing the correct metadata...');
[czi_file, czi_path] = uigetfile({'*.czi', 'Zeiss CZI Files (*.czi)'}, 'Select Original RAW File');
if isequal(czi_file, 0)
    disp('Canceled. Exiting...'); return;
end
raw_filepath = fullfile(czi_path, czi_file);

%% 2. Select the TARGET TIFF files to update
disp('Select the exported TIFF file(s) you want to inject the metadata into...');
[tif_list, tif_path] = uigetfile({'*.tif;*.tiff', 'TIFF Files (*.tif, *.tiff)'}, ...
    'Select Target TIFF(s)', 'MultiSelect', 'on');
if isequal(tif_list, 0)
    disp('Canceled. Exiting...'); return;
end
if ischar(tif_list)
    tif_list = {tif_list};
end

%% 3. Extract Metadata from RAW File
disp(['Loading metadata from: ', czi_file]);
tic;
img_data = bfopen(raw_filepath);
meta_hash = img_data{1, 2}; 
omeMeta = img_data{1, 4};
% We only need the total frame count from the OME data
T = omeMeta.getPixelsSizeT(0).getValue();
numChannels = omeMeta.getChannelCount(0);
disp(['Extracted basic info: ', num2str(T), ' frames, ', num2str(numChannels), ' channels.']);
toc;

disp('Extracting robust timing and ROI metadata...');
totalFrames = T; 
timeOfChange = NaN;
frames_phase1 = 0;

% --- EXTRACT EVENT MARKER TIME ---
startTime_str = char(meta_hash.get('Global Information|Image|T|StartTime'));
markerTime_str = char(meta_hash.get('Global Information|TimelineTrack|TimelineElement|Time #1'));

if ~isempty(startTime_str) && ~isempty(markerTime_str)
    try
        t_start = datetime(startTime_str(1:19), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
        t_marker = datetime(markerTime_str(1:19), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
        timeOfChange = seconds(t_marker - t_start); 
    catch
        try
            t_start = datetime(startTime_str); t_marker = datetime(markerTime_str);
            timeOfChange = seconds(t_marker - t_start);
        catch
            disp('Warning: Could not parse Event Marker timestamps.');
        end
    end
end

% --- EXTRACT TIME VECTOR ---
timevec = zeros(1, T);
valid_ome_time = true;

try
    for i = 1:T
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

if ~valid_ome_time || any(isnan(timevec))
    disp('Falling back to Hashtable Interval parsing...');
    interval_1 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #1'));
    interval_2 = str2double(meta_hash.get('Global Experiment|AcquisitionBlock|TimeSeriesSetup|TimeSeriesTriggerIntervalModeCollection|TimeSeriesTriggerIntervalMode|OnTrigger|Interval #2'));
    
    if ~isnan(interval_1) && ~isnan(interval_2) && ~isnan(timeOfChange)
        frames_phase1 = round(timeOfChange / interval_1);
        frames_phase2 = totalFrames - frames_phase1;
        timevec_phase1 = (0 : frames_phase1-1) * interval_1;
        timevec_phase2 = timevec_phase1(end) + (1 : frames_phase2) * interval_2;
        timevec = [timevec_phase1, timevec_phase2];
    else
        try
            dt = double(omeMeta.getPixelsTimeIncrement(0).value());
        catch
            dt = 1.0; 
        end
        timevec = (0:T-1) * dt;
    end
end

if frames_phase1 == 0 
    if ~isnan(timeOfChange)
        [~, frames_phase1] = min(abs(timevec - timeOfChange));
    else
        frames_phase1 = round(T/3); 
    end
end

time_str = sprintf('%.2f,', timevec);
time_str = time_str(1:end-1); 

% --- EXTRACT ROI GEOMETRIES DYNAMICALLY ---
roi_metadata_str = ''; 
num_ROIs = 0;

while true
    idx_pad = sprintf('%02d', num_ROIs + 1); 
    idx_nopad = sprintf('%d', num_ROIs + 1); 
    
    val_X = meta_hash.get(['Global Layer|Circle|Geometry|CenterX #', idx_pad]);
    active_idx = idx_pad;
    
    if isempty(val_X)
        val_X = meta_hash.get(['Global Layer|Circle|Geometry|CenterX #', idx_nopad]);
        active_idx = idx_nopad;
    end
    
    if isempty(val_X), break; end
    
    num_ROIs = num_ROIs + 1;
    cx = str2double(val_X);
    cy = str2double(meta_hash.get(['Global Layer|Circle|Geometry|CenterY #', active_idx]));
    rad = str2double(meta_hash.get(['Global Layer|Circle|Geometry|Radius #', active_idx]));
    
    roi_metadata_str = sprintf('%sROI_%02d=[%.2f,%.2f,%.2f]\n', roi_metadata_str, num_ROIs, cx, cy, rad);
end

% Construct final payload
metadata_string = sprintf('TimeVector_Seconds=[%s]\nTotalFrames=%d\nEventMarker=Frame_%d\n%s', ...
    time_str, totalFrames, frames_phase1, roi_metadata_str);

disp('--- METADATA PAYLOAD READY ---');
disp(['Time Vector: ', num2str(T), ' frames spanning ', num2str(max(timevec)), 's']);
disp(['ROIs Detected: ', num2str(num_ROIs)]);

%% 4. Inject Metadata into TIFFs
disp('Starting TIFF injection process...');

% Setup BigTIFF options
options_tiff.color = false;
options_tiff.compress = 'lzw'; 
options_tiff.message = true;
options_tiff.append = false;
options_tiff.overwrite = true;
options_tiff.big = true; 
options_tiff.ImageDescription = metadata_string; 

for f_idx = 1:length(tif_list)
    current_tif = fullfile(tif_path, tif_list{f_idx});
    [~, base_name, ext] = fileparts(current_tif);
    backup_tif = fullfile(tif_path, [base_name, '_backup', ext]);
    
    fprintf('\nUpdating file %d of %d: %s\n', f_idx, length(tif_list), tif_list{f_idx});
    
    % Load the existing TIFF data
    disp('  -> Loading image matrix...');
    img_stack = loadtiff(string(current_tif));
    
    % Create a backup of the original TIFF
    disp(['  -> Creating backup: ', base_name, '_backup', ext]);
    movefile(current_tif, backup_tif);
    
    % Re-save it with the new header under the original name
    disp('  -> Saving with updated ImageDescription header...');
    saveastiff(img_stack, current_tif, options_tiff);
    
    disp('  -> Done.');
end

disp('=======================================');
disp('All TIFF files successfully updated and backed up!');