function worm_mask = find_worms_image(varargin)
% FIND_WORMS_IMAGE function analyzes an image and returns a mask indicating
% worm positions
%
%   WORM_MASK = find_worms_image() allows gui selection of
%       image to analyze and allows the user to deselect erroneously
%       identified regions.  When finished selecting regions, press escape
%       to continue.
%
%   WORM_MASK = count_worms_images([filename OR image_data])
%
%   WORM_MASK = count_worms_images([filename OR image_data], minsize, maxsize)
%       minsize - Regions smaller than min_size will be discarded
%           default = 5
%       maxsize - Regions smaller than max_size will be used to determine
%            the size of a single worm
%           default = 40
%
%   WORM_MASK = count_worms_images([filename OR image_data], minsize, maxsize, debug)
%       debug [0/1] flag outputs various image overlays
%            default = 0 (off)

%% Parse arguments
p1 = inputParser;
p1.FunctionName = 'count_worms_image';
p1.addOptional('image_data',0,@isnumeric);
p1.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p1.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p1.addParamValue('debug',0,@isnumeric);

p2 = inputParser;
p2.FunctionName = 'count_worms_image';
p2.addOptional('filename','',@ischar);
p2.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p2.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p2.addParamValue('debug',0,@isnumeric);

p3 = inputParser;
p3.FunctionName = 'count_worms_image';
p3.addParamValue('minsize',10,@isnumeric); % Regions smaller than this will be discarded
p3.addParamValue('maxsize',40,@isnumeric); % Regions smaller than this will determine single worm size
p3.addParamValue('debug',0,@isnumeric);

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

if ( isfield(i_p.Results,'image_data') && ~strcmp(i_p.Results.image_data,0))
    image.data = i_p.Results.image_data;
else
    if ( (isfield(i_p.Results,'filename')) && ~strcmp(i_p.Results.filename,''))
        fullfilename = i_p.Results.filename;
    else
        [FileName,PathName,FilterIndex] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
            '*.*','All Files' },'Select Image File');
        fullfilename = [PathName, filesep, FileName];
    end
    image.info = imfinfo( fullfilename );
    image.data = imread(fullfilename);
end

debug = i_p.Results.debug;
min_worm_size = i_p.Results.minsize;
max_worm_size = i_p.Results.maxsize;

% Set Highlight Color
highlightColor = [0 .5 0];

%% Convert to grayscale
if ndims(image.data) > 2
    I_gray = rgb2gray(image.data);
else
    I_gray = mat2gray(image.data);
end

%% Setup figure for manual review
%names = regexp(image.info.Filename,'(?<path>.*)/(?<filename>.*)','names');
%review_fig = figure('Name', names.filename,'MenuBar','none','ToolBar','none');
review_fig = figure('Visible', 'off');
review_ax = gca;
%plot_fig(image.data, 1);
%title('Original Image');

% Find worms
[worm_mask, bg] = find_worms(I_gray, min_worm_size);

% Show review image
reviewimg = imoverlay(image.data, bwperim(worm_mask), highlightColor);
%subplot(1,2,2);
%[of, h_im] = plot_fig(reviewimg, 2);
%[of, h_im] = plot_fig_2(reviewimg);
h_im = imshow(reviewimg, 'Parent', review_ax);
set(gcf, 'Position', get(0,'Screensize'));  % Maximize view
%truesize(review_fig); % 100% size view
%movegui(review_fig,'center')
title('Select regions to ignore, press <ESC> when done');
set(review_fig, 'Visible', 'on')
e = imrect(gca);

%% Looping manual review
while ~isempty(e)
    
    % Replace selected regions with estimation of background
    mask = createMask(e,h_im);
    I_gray(mask) = bg(mask);
    
    % Reanalyze - esp useful when very dark regions are removed from the
    % image, allowing the algorithm to find the lighter worms
    [worm_mask, bg, I_smooth] = find_worms(I_gray, min_worm_size);
    
    % Setup and review results
    % The function IMOVERLAY creates a mask-based image overlay. It takes input
    % image and a binary mask, and it produces an output image whose masked
    % pixels have been replaced by a specified color.
    % MatLab Central -
    % http://www.mathworks.com/matlabcentral/fileexchange/10502
    reviewimg = imoverlay(image.data, bwperim(worm_mask), highlightColor);
    %[of, h_im] = plot_fig(reviewimg, 2);
    h_im = imshow(reviewimg, 'Parent', review_ax);
    title('Select regions to ignore, press <ESC> when done');
    e = imrect(gca);
end
close gcf;


% %% Debug output
% if (debug)
%     
%     % Debug image of all pixels considered worms
%     overlay2 = imoverlay(image.data, worm_mask, [.3 1 .3]);
%     
%     % Single worm image for debugging
%     single_worms = false(size(worm_mask));
%     single_worms(vertcat(wormdata(worm_areas<max_worm_size).PixelIdxList)) = true;
%     RGB_label = label2rgb(bwlabel(single_worms,4), @lines, 'k', 'shuffle');
%     
%     % Display images and debug info
%     fprintf('Estimated size of one worm: %.2f\n', worm_size);
%     fprintf('Estimated number of worms: %.0f\n', num_worms);
%     figure, imshow(RGB_label);
%     figure, imshow(overlay2);
%     figure, imshow(image.data);
% end

end


% find_worms function identifies dark worms in image
%
%   [MASK, BG] = plot_fig(image, min_worm_size) 
%       Finds worms larger than min_worm_size in image
%       MASK = Logical matrix indicating pixels corresponding to dark areas
%           (worms)
%       BG = Esitmated background of image (from tophat transform)
function [mask, bg, I_smooth] = find_worms(image, min_worm_size)

%% Complement of image
I_comp = imcomplement(image);

%% BG Substract
I_bsub = imtophat(I_comp,strel('disk',5));
bg = imcomplement(I_comp - I_bsub);

%% Noise removal
I_smooth = wiener2(I_bsub,[5 5]);

%% Thresholding
% Grayscale image conversion
I_gray = mat2gray(I_smooth);
threshold = median(I_gray(:))*3;
bw = im2bw(I_gray, threshold);

%% Morphologically open binary image (remove small objects) < min_worm_size
%   Determine connected components (8 pixel neighborhood)
%   Compute area of each component
%   Remove those below specified value
mask = bwareaopen(bw, min_worm_size, 8);

end
