close all; clear; clc;


%% Select Multiple Files
% Open UI dialog to select one or multiple files
[file_list, path_name] = uigetfile({'*.czi;*.tif;*.tiff', 'Image Files (*.czi, *.tif, *.tiff)'}, ...
    'Select Image Files', 'MultiSelect', 'off');

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
    