function [summary_results, all_results] = count_worms_directory(varargin)

p = inputParser;
p.FunctionName = 'count_worms_directory';
p.addOptional('inputDir', '', @isdir);
p.parse(varargin{:});

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
        [worm_size, num_worms] = count_worms_image([input_dir filesep types{t} num2str(trials(i)) '.png']);
        all_results = vertcat(all_results, {types{t}, i, worm_size, num_worms});
        trial_results = horzcat(trial_results, num_worms);
    end
    summary_results = vertcat(summary_results, trial_results);
end

cellwrite([input_dir filesep 'worm_counts_stats.csv'],all_results,',','wt');
cellwrite([input_dir filesep 'worm_counts_summary.csv'],summary_results,',','wt');
end