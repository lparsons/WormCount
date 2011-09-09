function summary_results = count_worms_directory(varargin)
% COUNT_WORMS_DIRECTORY analyzes sets of images and estimates worm count
%   
%   SUMMARY_RESULTS = count_worms_directory() allows gui
%       selection of the directory to analyze.  The directory must contain
%       sets of whole plate images for each strain, time point, and
%       replicate. File names must conform for the following specification:
%
%           [strain]_T[time]_[replicate].png
%
%       SUMMARY_RESULTS is a Nx9 cell array with the counts for each trial
%           This is saved in a file named 'worm_counts.csv' in the
%           directory specified.
%
%   SUMMARY_RESULTS = count_worms_directory(directory)
%
%   SUMMARY_RESULTS = count_worms_directory(directory, param_name, param_value, ...)
%
%   Parameters:
%       minsize - Regions smaller than min_size will be discarded
%           default = 10
%       maxsize - Regions smaller than max_size will be used to determine 
%            the size of a single worm
%           default = 40
%       origin_diameter - Diameter of the area considered the origin (in pixels)
%           default = 150
%       treatment_diameter - Diameter of the area considered for each treatment (in pixels)
%           default = 300
%       split_total - If true, split total image into four smaller images
%           default = 1
%       use_previous = If true, attempt to use previously generated results
%           default = 0
%       debug - if true then output debug info


p = inputParser;
p.FunctionName = 'count_worms_directory';
p.addOptional('inputDir', '', @isdir);
p.addOptional('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p.addOptional('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p.addParamValue('origin_diameter',150,@isnumeric); % Diameter of area on plate for origin (in pixels)
p.addParamValue('treatment_diameter',300,@isnumeric); % Diameter of area on plate for each treatment (in pixels)
p.addParamValue('split_total',1,@isnumeric); % If true, split total image into four smaller images
p.addParamValue('use_previous',0,@isnumeric); % If true, use previous results if found
p.addParamValue('debug',0,@isnumeric);
p.parse(varargin{:});

% If directory not specified, allow user to choose the directory
if ( (isfield(p.Results,'inputDir')) && ~strcmp(p.Results.inputDir,''))
    input_dir = p.Results.inputDir;
else
    input_dir = uigetdir([],'Select Directory');
end

% Get list of files
plate_images = dir([input_dir filesep '*.png']);
num_images = size(plate_images,1);
col_headings = {'Strain', 'Time', 'Replicate', 'Est. Worm Size', 'Eth', 'But', 'Ori', 'Tot', 'CI'};
summary_results{num_images+1,size(col_headings,2)} = [];
summary_results(1,:) = col_headings;

for i=1:num_images
    t = regexpi(plate_images(i).name, '([A-Za-z0-9]+)_T([A-Za-z0-9-]+)_([A-Za-z0-9]+).png', 'tokens');
    [PATHSTR,NAME,EXT,VERSN] = fileparts(plate_images(i).name);  %#ok<ASGLU>
    data_path = [input_dir filesep PATHSTR filesep NAME '_results'];
    [s,mess,messid] = mkdir(data_path);
    data_filename = fullfile(data_path, [NAME '_data.mat' VERSN]); 

    %% Get Worm Counts
    if exist(data_filename, 'file') && p.Results.use_previous
        load(data_filename);
        display(sprintf('Using previously generated counts for %s', NAME))
    else
        plate_results = count_worms_plate(...
            [input_dir filesep plate_images(i).name], ...
            'minsize', p.Results.minsize, ...
            'maxsize', p.Results.maxsize, ...
            'origin_diameter', p.Results.origin_diameter, ...
            'treatment_diameter', p.Results.treatment_diameter, ...
            'split_total', p.Results.split_total, ...
            'debug', p.Results.debug);
        save(data_filename, 'plate_results');
    end
    if isempty(t)
        exp_info = [{i} {''} {''}];
    else
        exp_info = t{1};
    end   
    summary_results(i+1,:) = [exp_info plate_results.worm_size plate_results.eth plate_results.but plate_results.ori plate_results.tot plate_results.ci];
end

cellwrite([input_dir filesep 'worm_counts.csv'],summary_results,',','wt');
end