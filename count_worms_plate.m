function plate_results = count_worms_plate(varargin)
% COUNT_WORMS_PLATE analyzes an image of a full plate estimates worm count
% for each of four locations on the plate
%
%   PLATE_RESULTS = count_worms_plate() allows gui
%       selection of the pate image file to analyze.
%
%       PLATE_RESULTS is a 1x5 cell array with the estimated worm size and
%       the worm counts for each treatment
%
%   PLATE_RESULTS = count_worms_plate(filename)
%
%   PLATE_RESULTS = count_worms_directory(filename, param_name, param_value, ...)
%
%   Parameters:
%       minsize - Regions smaller than min_size will be discarded
%           default = 10
%       maxsize - Regions smaller than max_size will be used to determine
%            the size of a single worm
%           default = 40
%       area_width - Width of the area around each treatment (in pixels)
%           default = 300
%       split_total - If true, split total image into four smaller images
%       debug - if true then output debug info

% Image specified
p1 = inputParser;
p1.FunctionName = 'count_worms_plate';
p1.addOptional('filename', '', @ischar);
p1.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p1.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p1.addParamValue('area_width',300,@isnumeric); % Width of area on plate for each treatment (in pixels)
p1.addParamValue('split_total',1,@isnumeric); % If true, split total image into four smaller images
p1.addParamValue('debug',0,@isnumeric);

% No image specified, select using GUI
p2 = inputParser;
p2.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p2.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p2.addParamValue('area_width',300,@isnumeric); % Width of area on plate for each treatment (in pixels)
p2.addParamValue('split_total',1,@isnumeric); % If true, split total image into four smaller images
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
[PATHSTR,NAME,EXT,VERSN] = fileparts(fullfilename);
data_path = [PATHSTR filesep NAME '_results'];
[s,mess,messid] = mkdir(data_path);
image_overlay_filename = fullfile(data_path, [NAME '_overlay.png' VERSN]);
data_filename = fullfile(data_path, [NAME '_data.mat' VERSN]);

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

% Grey out unmasked area
masked_total = image.data;
masked_total(~mask) = median(image.data(mask));

% Crop image
mask_props = regionprops(bwlabel(mask), 'BoundingBox');
masked_total = imcrop(masked_total, mask_props.BoundingBox);
image_handle = imshow(masked_total);
set(plate_fig, 'Position', get(0,'Screensize')); % Maximize figure

%% Setup total image(s)
cropped_images.tot = {};
if p.Results.split_total
    w = floor(size(masked_total,1)/2);
    h = floor(size(masked_total,2)/2);
    cropped_images.tot{1} = masked_total(1:h,1:w);
    cropped_images.tot{2} = masked_total(1:h,w+1:end);
    cropped_images.tot{3} = masked_total(h+1:end,1:w);
    cropped_images.tot{4} = masked_total(h+1:end,w+1:end);
else
    cropped_images.tot{1} = masked_total;
end


%% Choose area for each treatment from plate image
types = {'eth', 'but', 'ori'};
selected_points = zeros(size(types,2),2);
roi_width = p.Results.area_width;
for t=1:size(types,2)
    disp(types{t})
    title(sprintf('Select center of %s area', types{t}))
    [x,y] = ginput(1);
    selected_points(t,:) = [x,y];
    box = [x-(roi_width/2) y-(roi_width/2) roi_width roi_width];
    hold on
    rectangle('Position', box, 'Curvature', [1,1], 'EdgeColor', 'r')
    text(x, y, types{t}, 'BackgroundColor', [.8 .7 .7], 'VerticalAlignment', 'bottom')
    hold off
    cropped_images.(types{t}) = imcrop(masked_total, box);
end
close(plate_fig)


%% Mask areas on total image(s)
masks = {};
for i=1:size(cropped_images.tot,2)
    img = cropped_images.tot{i};
    masks{i} = find_worms_image(img, 'minsize', min_worm_size, 'maxsize', max_worm_size, 'debug', p.Results.debug);
end
full_mask = false(size(masked_total));
if p.Results.split_total
    w = floor(size(masked_total,1)/2);
    h = floor(size(masked_total,2)/2);
    full_mask(1:h,1:w) = masks{1};
    full_mask(1:h,w+1:end) = masks{2};
    full_mask(h+1:end,1:w) = masks{3};
    full_mask(h+1:end,w+1:end) = masks{4};
else
    full_mask = masks{1};
end

%% Count worms on total image
[total_num_worms, total_worm_size] = count_worms_image('worm_mask', full_mask, 'minsize', min_worm_size, 'maxsize', max_worm_size, 'debug', p.Results.debug);
plate_results.worm_size = total_worm_size;
plate_results.tot = total_num_worms;

%% For each treatment, analyze area around each selected center, generate overlay summary
overlay_figure = figure('Visible', 'off');
imshow(imoverlay(masked_total, bwperim(full_mask), [0 .5 0]));
for t=1:size(types,2)
    xy = selected_points(t,:);
    box = [xy(1)-(roi_width/2) xy(2)-(roi_width/2) roi_width roi_width];
    region_handle = imellipse(gca, box);
    region_mask = createMask(region_handle);
    delete(region_handle)
    region_worm_mask = region_mask & full_mask;
    [num_worms, worm_size] = count_worms_image('worm_mask', region_worm_mask, 'avg_worm_size', total_worm_size, 'minsize', min_worm_size, 'maxsize', max_worm_size);
    plate_results.(types{t}) = num_worms;
    rectangle('Position', box, 'Curvature', [1,1], 'EdgeColor', 'r')
    text(xy(1), xy(2), [types{t}, ' = ', num2str(num_worms)], 'BackgroundColor', 'none', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center')
end

plate_results.ci = chemotaxis_index(plate_results);

title_text = sprintf('Total = %s, CI = %s, Worm size = %s', ...
    num2str(total_num_worms), ...
    num2str(plate_results.ci), ...
    num2str(plate_results.worm_size));

text(floor(size(masked_total,1)/2), ...
    5, ...
    title_text, ...
    'BackgroundColor', [.8 .7 .7], ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center' ...
    )

% Print Overlay Summary Figure
screen_DPI = get(0, 'ScreenPixelsPerInch');
set(overlay_figure, 'Units', 'pixels', 'Position', [32, 32, size(masked_total, 2), size(masked_total, 1)]);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gca, 'Units', 'points');
set(overlay_figure, 'Units', 'points', 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
print(overlay_figure, '-dpng', sprintf('-r%d', screen_DPI), image_overlay_filename);
close(overlay_figure);
end











