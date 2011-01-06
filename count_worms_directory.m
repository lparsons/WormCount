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
%           default = 5
%       maxsize - Regions smaller than max_size will be used to determine 
%            the size of a single worm
%           default = 30
%       area_width - Width of the area around each treatment (in pixels)
%           default = 300
%       split_total - If true, split total image into four smaller images
%       debug - if true then output debug info


p = inputParser;
p.FunctionName = 'count_worms_directory';
p.addOptional('inputDir', '', @isdir);
p.addOptional('minsize',5,@isnumeric); % Regions smaller than this will be discarded
p.addOptional('maxsize',30,@isnumeric); % Regions smaller than this will determine single worm size
p.addParamValue('area_width',300,@isnumeric); % Width of area on plate for each treatment (in pixels)
p.addParamValue('split_total',1,@isnumeric); % If true, split total image into four smaller images
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
    plate_results = count_worms_plate(...
        [input_dir filesep plate_images(i).name], ...
        'minsize', p.Results.minsize, ...
        'maxsize', p.Results.maxsize, ...
        'area_width', p.Results.area_width, ...
        'split_total', p.Results.split_total, ...
        'debug', p.Results.debug);
    ci = (plate_results{3} - plate_results{2}) / (plate_results{5} - plate_results{4});
    summary_results(i+1,:) = [t{1} plate_results ci];
end

cellwrite([input_dir filesep 'worm_counts.csv'],summary_results,',','wt');
end


%% OLD CODE

% but_images = dir([input_dir filesep 'but*.png']);
% trials = [];
% for i=1:size(but_images,1)
%     t = regexpi(but_images(i).name, 'but([0-9]+).png', 'tokens');
%     n = str2double(t{1});
%     trials = [trials, n];
% end
% 
% disp(['Found ' num2str(size(trials,2)) ' trials in ''' input_dir '''']);
% 
% all_results = {};
% summary_results = {'Trial', 'Eth', 'But', 'Ori', 'Tot'};
% types = {'eth', 'but', 'ori', 'tot'};
% for i=1:size(trials,2)
%     trial_results = {trials(i)};
%     for t=1:size(types,2)
%         disp([types{t} num2str(trials(i))]);
%         [num_worms, worm_size] = count_worms_image([input_dir filesep types{t} num2str(trials(i)) '.png'], 'minsize', min_worm_size, 'maxsize', max_worm_size);
%         all_results = vertcat(all_results, {types{t}, i, worm_size, num_worms});
%         trial_results = horzcat(trial_results, num_worms);
%     end
%     summary_results = vertcat(summary_results, trial_results);
% end
% 
% cellwrite([input_dir filesep 'worm_counts_stats.csv'],all_results,',','wt');
% cellwrite([input_dir filesep 'worm_counts_summary.csv'],summary_results,',','wt');
% end