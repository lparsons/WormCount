function [summary_results, all_results] = count_worms_directory(varargin)
% COUNT_WORMS_DIRECTORY analyzes sets of images and estimates worm count
%   
%   [SUMMARY_RESULTS, ALL_RESULTS] = count_worms_directory() allows gui
%       selection of the directory to analyze.  The directory must contain
%       sets of images numbered 1 to N where N is number of trials.
%
%           e.g. eth1.png, but1.png, ori1.png, tot1.png
%
%       SUMMARY_RESULTS is a Nx5 cell array with the counts for each trial
%           This is saved in a file named 'worm_counts_summary.csv' in the
%           directory specified.
%
%       ALL_RESULTS outputs one row per image, indicating the treatment,
%           trial number, the size of a worm, and the worm count.  This is
%           saved in a file named 'worm_counts_stats.csv'.
%
%   [SUMMARY_RESULTS, ALL_RESULTS] = count_worms_directory(directory)
%
%   [SUMMARY_RESULTS, ALL_RESULTS] = count_worms_directory(directory, minsize, maxsize)
%       minsize - Regions smaller than min_size will be discarded
%           default = 10
%       maxsize - Regions smaller than max_size will be used to determine 
%            the size of a single worm
%           default = 80

p = inputParser;
p.FunctionName = 'count_worms_directory';
p.addOptional('inputDir', '', @isdir);
p.addOptional('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p.addOptional('maxsize',80,@isnumeric); % Regions smaller than this will determine single worm size
p.parse(varargin{:});
min_worm_size = p.Results.minsize; 
max_worm_size = p.Results.maxsize; 

if ( (isfield(p.Results,'inputDir')) && ~strcmp(p.Results.inputDir,''))
    input_dir = p.Results.inputDir;
else
    input_dir = uigetdir([],'Select Directory');
end

but_images = dir([input_dir filesep 'but*.png']);
trials = [];
for i=1:size(but_images,1)
    t = regexpi(but_images(i).name, 'but([0-9]+).png', 'tokens');
    n = str2double(t{1});
    trials = [trials, n];
end

disp(['Found ' num2str(size(trials,2)) ' trials in ''' input_dir '''']);

all_results = {};
summary_results = {'Trial', 'Eth', 'But', 'Ori', 'Tot'};
types = {'eth', 'but', 'ori', 'tot'};
for i=1:size(trials,2)
    trial_results = {trials(i)};
    for t=1:size(types,2)
        disp([types{t} num2str(trials(i))]);
        [worm_size, num_worms] = count_worms_image([input_dir filesep types{t} num2str(trials(i)) '.png'], 'minsize', min_worm_size, 'maxsize', max_worm_size);
        all_results = vertcat(all_results, {types{t}, i, worm_size, num_worms});
        trial_results = horzcat(trial_results, num_worms);
    end
    summary_results = vertcat(summary_results, trial_results);
end

cellwrite([input_dir filesep 'worm_counts_stats.csv'],all_results,',','wt');
cellwrite([input_dir filesep 'worm_counts_summary.csv'],summary_results,',','wt');
end