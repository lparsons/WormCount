function [num_worms, worm_size] = count_worms_image(varargin)
% COUNT_WORMS_IMAGE function analyzes an image and estimates worm count
%
%   [WORM_SIZE, NUM_WORMS] = count_worms_images() allows gui selection of
%       image to analyze and allows the user to deselect erroneously
%       identified regions.  When finished selecting regions, press escape
%       to continue.
%
%       WORM_SIZE is the estimated size of a single worm (in pixels)
%
%       NUM_WORMS is the estimated number of worms in the image based on
%           the worm_size.
%
%   [WORM_SIZE, NUM_WORMS] = count_worms_images([filename OR image_data])
%
%   [WORM_SIZE, NUM_WORMS] = count_worms_images([filename OR image_data], minsize, maxsize)
%       minsize - Regions smaller than min_size will be discarded
%           default = 10
%       maxsize - Regions smaller than max_size will be used to determine
%            the size of a single worm
%           default = 40
%
%   [WORM_SIZE, NUM_WORMS] = count_worms_images([filename OR image_data], minsize, maxsize, debug)
%       debug [0/1] flag outputs various image overlays
%            default = 0 (off)

%% Parse arguments
p1 = inputParser;
p1.FunctionName = 'count_worms_image';
p1.addRequired('image_data',@isnumeric);
p1.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p1.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p1.addParamValue('avg_worm_size',0,@isnumeric); % Use predetermined worm size
p1.addParamValue('debug',0,@isnumeric);

p2 = inputParser;
p2.FunctionName = 'count_worms_image';
p2.addRequired('filename',@ischar);
p2.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p2.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p2.addParamValue('avg_worm_size',0,@isnumeric); % Use predetermined worm size
p2.addParamValue('debug',0,@isnumeric);

p3 = inputParser;
p3.FunctionName = 'count_worms_image';
p3.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p3.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p3.addParamValue('debug',0,@isnumeric);
p3.addParamValue('avg_worm_size',0,@isnumeric); % Use predetermined worm size
p3.addParamValue('worm_mask',0,@islogical);

try
    p1.parse(varargin{:})
    i_p = p1;
catch e1
    try
        p2.parse(varargin{:})
        i_p = p2;
    catch e2
        try
            p3.parse(varargin{:})
            i_p = p3;
        catch e3
            
            exception = MException(...
                'count_worms_plate:arglist',...
                'Error in input argument list');
            exception = addCause(exception, e1);
            exception = addCause(exception, e2);
            exception = addCause(exception, e3);
            throw(exception)
        end
    end
end

%% Get image data and find worms
debug = i_p.Results.debug;
min_worm_size = i_p.Results.minsize;
max_worm_size = i_p.Results.maxsize;

if isfield(i_p.Results,'worm_mask') && islogical(i_p.Results.worm_mask)

    worm_mask = i_p.Results.worm_mask;
else
    if ( isfield(i_p.Results,'image_data') && ~strcmp(i_p.Results.image_data,0))
        image.data = i_p.Results.image_data;
    else
        if ( (isfield(i_p.Results,'filename')) && ~strcmp(i_p.Results.filename,''))
            fullfilename = i_p.Results.filename;
        else
            [FileName,PathName,FilterIndex] = uigetfile(...
                {'*.jpg;*.tif;*.png;*.gif',...
                'All Image Files';...
                '*.*',...
                'All Files' }, ...
                'Select Image File');
            fullfilename = [PathName, filesep, FileName];
        end
        image.info = imfinfo( fullfilename );
        image.data = imread(fullfilename);
    end
    worm_mask = find_worms_image(image.data, ...
        'minsize', min_worm_size, ...
        'maxsize', max_worm_size, ...
        'debug', i_p.Results.debug);
end


%% Estimate single worm size
wormdata = regionprops(bwlabel(worm_mask, 8), 'Area', 'PixelIdxList');
worm_areas = [wormdata.Area];
worm_size = i_p.Results.avg_worm_size;
if worm_size == 0
    worm_size = median(worm_areas(worm_areas<=max_worm_size));
end

%% Estimate number of worms in image
num_worms = round(sum(worm_mask(:))/worm_size);

%% Debug output
if (debug)
    
    % Debug image of all pixels considered worms
    if exist('image', 'var')
        overlay2 = imoverlay(image.data, worm_mask, [.3 1 .3]);
    else
        overlay2 = worm_mask;
    end
    
    % Single worm image for debugging
    single_worms = false(size(worm_mask));
    single_worms(vertcat(wormdata(worm_areas<max_worm_size).PixelIdxList)) = true;
    RGB_label = label2rgb(bwlabel(single_worms,8), @lines, 'k', 'shuffle');
    
    % Display images and debug info
    fprintf('Estimated size of one worm: %.2f\n', worm_size);
    fprintf('Estimated number of worms: %.0f\n', num_worms);
    figure, imshow(RGB_label);
    figure, imshow(overlay2);
    if exist('image', 'var')
        figure, imshow(image.data);
    end
end

end

