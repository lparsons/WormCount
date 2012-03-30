%   Copyright 2011 Lance R. Parsons, Lewis-Sigler Institute for Integrative Genomics, Princeton University
%
%   Licensed under the BSD 2-Clause License 
%   http://www.opensource.org/licenses/BSD-2-Clause

function worm_mask = find_worms_image(varargin)
% FIND_WORMS_IMAGE function analyzes an image and returns a mask indicating
% worm positions
%
%   WORM_MASK = find_worms_image(minsize, maxsize) allows gui selection of
%       image to analyze and allows the user to deselect erroneously
%       identified regions.  When finished selecting regions, press escape
%       to continue.
%
%       minsize - Regions smaller than min_size will be discarded
%       maxsize - Regions smaller than max_size will be used to determine
%            the size of a single worm
%
%   WORM_MASK = count_worms_images([filename OR image_data], minsize, maxsize)
%       Specify either filename OR image_data
%
%       filename - specifies the name of an image file to process
%       image_data - specifies numeric image data to process
%
%   WORM_MASK = count_worms_images([filename OR image_data], minsize, maxsize, debug)
%       debug [0/1] flag outputs various image overlays
%            default = 0 (off)

%%
[pathstr, name, ext] = fileparts(mfilename('fullpath')) ; %#ok<NASGU,ASGLU>
addpath(genpath([pathstr filesep 'lib']));

%% Parse arguments
p1 = inputParser;
p1.FunctionName = 'find_worms_image';
p1.addRequired('image_data',@isnumeric);
p1.addRequired('minsize',@isnumeric); % Regions smaller than this will be discarded
p1.addRequired('maxsize',@isnumeric); % Regions smaller than this will determine single worm size
p1.addParamValue('debug',0,@isnumeric);

p2 = inputParser;
p2.FunctionName = 'find_worms_image';
p2.addRequired('filename',@ischar);
p2.addRequired('minsize',@isnumeric); % Regions smaller than this will be discarded
p2.addRequired('maxsize',@isnumeric); % Regions smaller than this will determine single worm size
p2.addParamValue('debug',0,@isnumeric);

p3 = inputParser;
p3.FunctionName = 'find_worms_image';
p3.addRequired('minsize',@isnumeric); % Regions smaller than this will be discarded
p3.addRequired('maxsize',@isnumeric); % Regions smaller than this will determine single worm size
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
highlightColor = [1 0 0];

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
plot_fig(image.data, 1);
title('Original Image');

% Enhance image
I_enhanced = enhance_image(I_gray);

% Find worms
worm_mask = find_worms(I_enhanced, min_worm_size);

% Show review image
reviewimg = imoverlay(image.data, bwperim(worm_mask), highlightColor);
[of, h_im] = plot_fig(reviewimg, 2);
title('Select regions to ignore, press <ESC> when done');
truesize(review_fig); % 100% size view
movegui(review_fig,'center')

%h_im = imshow(reviewimg, 'Parent', review_ax);
%set(review_fig, 'Position', get(0,'Screensize'));  % Maximize view
%title('Select regions to ignore, press <ESC> when done');
%set(review_fig, 'Visible', 'on')

e = imrect(gca);


%% Looping manual review
while ~isempty(e)
    
    % Replace selected regions with estimation of background
    mask = createMask(e,h_im);
    %I_gray(mask) = bg(mask);
    I_enhanced(mask) = median(I_enhanced(~mask));
    
    % Reanalyze - esp useful when very dark regions are removed from the
    % image, allowing the algorithm to find the lighter worms
    %[worm_mask, bg, I_smooth] = find_worms(I_gray, min_worm_size);
    worm_mask= find_worms(I_enhanced, min_worm_size);
    
    % Remove roi region
    worm_mask = worm_mask & ~mask;
    
    % Setup and review results
    % The function IMOVERLAY creates a mask-based image overlay. It takes input
    % image and a binary mask, and it produces an output image whose masked
    % pixels have been replaced by a specified color.
    % MatLab Central -
    % http://www.mathworks.com/matlabcentral/fileexchange/10502
    reviewimg = imoverlay(image.data, bwperim(worm_mask), highlightColor);
    
    [of, h_im] = plot_fig(reviewimg, 2);
    
    %h_im = imshow(reviewimg, 'Parent', review_ax);
    %title('Select regions to ignore, press <ESC> when done');
    
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

% enchance_image function enhances images for worm detection
%
%   image = enhance_image(image) 
function enhanced_image = enhance_image(image)

%% Enchance
%I_enhance = imadjust(image, [], [], 4);

%% Complement of image
I_comp = imcomplement(image);

%% BG Substract
I_bsub = imtophat(I_comp,strel('disk',10));
%bg = imcomplement(I_comp - I_bsub);


%% Noise removal
%I_smooth = wiener2(I_enhance,[5 5]);

enhanced_image = I_bsub;

end


% find_worms function identifies dark worms in image
%
%   [MASK, BG] = plot_fig(image, min_worm_size) 
%       Finds worms larger than min_worm_size in image
%       MASK = Logical matrix indicating pixels corresponding to dark areas
%           (worms)
%       BG = Esitmated background of image (from tophat transform)
function mask = find_worms(image, min_worm_size)


%% Thresholding
% Grayscale image conversion
I_gray = mat2gray(image);
threshold = median(I_gray(:))*3;
%threshold = median(I_gray(I_gray>0)) + std(I_gray(I_gray>0)) * 2.5;
bw = im2bw(I_gray, threshold);

%% Morphologically open binary image (remove small objects) < min_worm_size
%   Determine connected components (8 pixel neighborhood)
%   Compute area of each component
%   Remove those below specified value
mask = bwareaopen(bw, min_worm_size, 8);

end


% plot_fig function plots a subimage for manual review
%
%   [SUBPLOT_HANDLE, IMAGE_HANDLE] = plot_fig(image, loc) 
%       Plots 'image' in a subplot at location loc
function [of, im] = plot_fig(image, loc)
of = subplot(1,2,loc); im = subimage(image);
set(of,'xtick',[],'ytick',[]);
p = get(of, 'pos');
if loc==1
    p(1) = p(1) - 0.05;
    p(3) = p(3) + 0.1;
    
else
    p(1) = p(1) - 0.05;
    p(3) = p(3) + 0.1;
    
end
set(of, 'pos', p);
end