function [results] = count_worms_directory(varargin)

i_p = inputParser;
i_p.FunctionName = 'count_worms_directory';
i_p.addOptional('input_dir',@isdir);
i_p.parse(varargin{:});

if (isfield(i_p,'input_dir'))
    input_dir = i_p.Results.input_dir;
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

disp(['Found ' num2str(size(trials,2)) ' trials.']);

results = {};
types = {'but', 'eth', 'ori', 'tot'};
for i=1:size(trials,2)
    for t=1:size(types,2)
        disp([types{t} num2str(trials(i))]);
        [worm_size, num_worms] = count_worms_image([input_dir filesep types{t} num2str(trials(i)) '.png']);
        results = vertcat(results, {types{t}, i, worm_size, num_worms});
    end
end

cellwrite([input_dir filesep 'worm_counts.csv'],results);

end