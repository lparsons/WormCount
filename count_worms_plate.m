function [plate_results, all_results] = count_worms_plate(varargin)
% COUNT_WORMS_PLATE analyzes an image of a full plate estimates worm count
% for each of four locations on the plate
%   
%   [PLATE_RESULTS, ALL_RESULTS] = count_worms_plate() allows gui
%       selection of the pate image file to analyze.
%
%       PLATE_RESULTS is a 1x4 cell array with the counts for each
%       treatment
%
%       ALL_RESULTS outputs one row per treatment, indicating the treatment,
%           the size of a worm, and the worm count.
%
%   [SUMMARY_RESULTS, ALL_RESULTS] = count_worms_plate(filename)
%
%   [SUMMARY_RESULTS, ALL_RESULTS] = count_worms_directory(filename, minsize, maxsize, area_width)
%       minsize - Regions smaller than min_size will be discarded
%           default = 10
%       maxsize - Regions smaller than max_size will be used to determine 
%            the size of a single worm
%           default = 80
%       area_width - Width of the area around each treatment (in pixels)
%           default = 300

% Image specified
p1 = inputParser;
p1.FunctionName = 'count_worms_plate';
p1.addOptional('filename', '', @ischar);
p1.addParamValue('minsize',15,@isnumeric); % Regions smaller than this will be discarded
p1.addParamValue('maxsize',100,@isnumeric); % Regions smaller than this will determine single worm size
p1.addParamValue('area_width',300,@isnumeric); % Width of area on plate for each treatment (in pixels)
p1.addParamValue('debug',0,@isnumeric);

% No image specified, select using GUI
p2 = inputParser;
p2.addParamValue('minsize',15,@isnumeric); % Regions smaller than this will be discarded
p2.addParamValue('maxsize',100,@isnumeric); % Regions smaller than this will determine single worm size
p2.addParamValue('area_width',300,@isnumeric); % Width of area on plate for each treatment (in pixels)
p2.addParamValue('debug',0,@isnumeric);

try
    p1.parse(varargin{:})
    p = p1;
catch e1
    try
        p2.parse(varargin{:})
        p = p2;
    catch e2
        exception = MException(...
            'count_worms_plate:arglist',...
            'Error in input argument list');
        exception = addCause(exception, e1);
        exception = addCause(exception, e2);
        throw(exception)
    end
end

p.parse(varargin{:});
min_worm_size = p.Results.minsize; 
max_worm_size = p.Results.maxsize; 

%% Select File
if ( (isfield(p.Results,'filename')) && ~strcmp(p.Results.filename,''))
    fullfilename = p.Results.filename;
else
    [FileName,PathName,FilterIndex] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
        '*.*','All Files' },'Select Image File');
    fullfilename = [PathName, filesep, FileName];
end

%% Read in image
image.info = imfinfo( fullfilename );
image.data = imread(fullfilename);

%% Mask dish
plate_fig = figure();
image_handle = imshow(image.data);
title('Select region to analyze and double click when finished')

% Select region
h = imellipse(gca, [10 10 1000 1000]);
setFixedAspectRatioMode(h, true)
position = wait(h);
mask = createMask(h,image_handle);

% Complement of image
I_comp = imcomplement(image.data);
% BG Substract
I_bsub = imtophat(I_comp,strel('disk',15));
bg = imcomplement(I_comp - I_bsub);

% Set total image
cropped_images.tot = image.data;
cropped_images.tot(~mask) = median(cropped_images.tot(mask));
imshow(cropped_images.tot)


%% Choose area for each treatment from plate image
image_handle = imshow(cropped_images.tot);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure

types = {'eth', 'but', 'ori', 'tot'};
selected_points = zeros(size(types,2),2);
roi_width = p.Results.area_width;
for t=1:size(types,2)-1
    disp(types{t})
    title(sprintf('Select center of %s area', types{t}))
    [x,y] = ginput(1);
    selected_points(t,:) = [x,y];
    box = [x-(roi_width/2) y-(roi_width/2) roi_width roi_width];
    hold on
    rectangle('Position', box, 'EdgeColor', 'r')
    text(box(1), box(2), types{t}, 'BackgroundColor', [.8 .7 .7], 'VerticalAlignment', 'bottom')
    hold off
    cropped_images.(types{t}) = imcrop(cropped_images.tot, box);
end
close(plate_fig)


%% For each treatment, analyze area around each selected center
plate_results = {};
all_results = {};
for t=1:size(types,2)
    img = cropped_images.(types{t});
    [num_worms, worm_size] = count_worms_image(img, 'minsize', min_worm_size, 'maxsize', max_worm_size, 'debug', p.Results.debug);
    all_results = vertcat(all_results, {types{t}, worm_size, num_worms});
    plate_results = horzcat(plate_results, num_worms);
end

plate_results = vertcat({'Eth', 'But', 'Ori', 'Tot'}, plate_results);