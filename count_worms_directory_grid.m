%   Copyright 2011 Lance R. Parsons, Lewis-Sigler Institute for Integrative Genomics, Princeton University
%
%   Licensed under the BSD 2-Clause License 
%   http://www.opensource.org/licenses/BSD-2-Clause

function summary_results = count_worms_directory_grid(varargin)
% COUNT_WORMS_DIRECTORY_GRID analyzes sets of images and estimates worm count
%   
%   SUMMARY_RESULTS = count_worms_directory() allows gui
%       selection of the directory to analyze.  The directory must contain
%       sets of whole plate images for each strain, time point, and
%       replicate. File names must conform for the following specification:
%
%           [strain]_T[time]_[replicate].png
%
%       SUMMARY_RESULTS is a NxM cell array with the counts for each trial
%           This is saved in a file named 'worm_counts.csv' in the
%           directory specified.
%
%   SUMMARY_RESULTS = count_worms_directory(directory)
%
%   SUMMARY_RESULTS = count_worms_directory(directory, param_name, param_value, ...)
%
%   Parameters:
%       minsize - Regions smaller than min_size will be discarded
%           default = 5
%       maxsize - Regions smaller than max_size will be used to determine 
%            the size of a single worm
%           default = 30
%       grid_cols - Number of columns in the grid
%           default = 6
%       grid_rows - Number of rows in the grid
%           default = 6
%       category_grid - Struct with fields for each category
%           Each field contains N x 2 matrix with x and y coords for
%               the grid positions belonging to that category
%           default = struct(...
%                 'left', ...
%                 [1,1; 2,1; 3,1; 4,1; 5,1; 6,1 ...
%                 ;1,2; 2,2; 3,2; 4,2; 5,2; 6,2], ...
%                 'center_top_bottom', ...
%                 [1,3; 2,3; 5,3; 6,3 ...
%                 ;1,4; 2,4; 5,4; 6,4], ...
%                 'center_center', ...
%                 [3,3; 4,3 ...
%                 ;3,4; 4,4], ...
%                 'right', ...
%                 [1,5; 2,5; 3,5; 4,5; 5,5; 6,5 ...
%                 ;1,6; 2,6; 3,6; 4,6; 5,6; 6,6] ...
%                 );
%       use_previous = If true, attempt to use previously generated results
%           default = 0
%       debug - if true then output debug info

[pathstr, name, ext] = fileparts(mfilename('fullpath')) ; %#ok<NASGU,ASGLU>
addpath(genpath([pathstr filesep 'lib']));

p = inputParser;
p.FunctionName = 'count_worms_directory';
p.addOptional('inputDir', '', @isdir);
p.addOptional('minsize',5,@isnumeric); % Regions smaller than this will be discarded
p.addOptional('maxsize',30,@isnumeric); % Regions smaller than this will determine single worm size
p.addParamValue('grid_cols',6,@isnumeric); % Number of columns in the grid
p.addParamValue('grid_rows',6,@isnumeric); % Number of rows in the grid
p.addParamValue('category_grid',default_category_grid(),@isstruct); % Struct with fields for each category
p.addParamValue('split_total',1,@isnumeric); % If true, split total image into four smaller images
p.addParamValue('use_previous',0,@isnumeric); % If true, use previous results if found
p.addParamValue('debug',false,@islogical);
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
col_headings = [{'Strain', 'Time', 'Replicate', 'Est. Worm Size'} fieldnames(p.Results.category_grid)'];
summary_results{num_images+1,size(col_headings,2)} = [];
summary_results(1,:) = col_headings;

for i=1:num_images
    t = regexpi(plate_images(i).name, '([A-Za-z0-9]+)_T([A-Za-z0-9-]+)_([A-Za-z0-9]+).png', 'tokens');
    [PATHSTR,NAME,EXT] = fileparts(plate_images(i).name);   %#ok<NASGU>
    data_path = [input_dir filesep PATHSTR filesep NAME '_results'];
    [s,mess,messid] = mkdir(data_path); %#ok<ASGLU,NASGU>
    data_filename = fullfile(data_path, [NAME '_data.mat']); 

    %% Get Worm Counts
    if exist(data_filename, 'file') && p.Results.use_previous
        load(data_filename);
        display(sprintf('Using previously generated counts for %s', NAME))
    else
        [worm_count, worm_size] = count_worms_grid([input_dir filesep plate_images(i).name], ...
            p.Results.grid_cols, p.Results.grid_rows, ...
            'minsize', p.Results.minsize, ...
            'maxsize', p.Results.maxsize, ...
            'debug', p.Results.debug);
        plate_results = sum_grid(worm_count, p.Results.category_grid);
        plate_results.worm_size = worm_size;
        %plate_results.ci = (plate_results.left - plate_results.right) / ...
        %    (plate_results.left + plate_results.right + plate_results.center_top_bottom);
        save(data_filename, 'plate_results');
    end
    if isempty(t)
        exp_info = [{i} {''} {''}];
    else
        exp_info = t{1};
    end
    results = [exp_info plate_results.worm_size];
    categories = fieldnames(p.Results.category_grid);
    for c = 1:size(categories,1)
        category = categories{c};
        results = [results plate_results.(category)];
    end
    summary_results(i+1,:) = results;
end

cellwrite([input_dir filesep 'worm_counts.csv'],summary_results,',','wt');
end