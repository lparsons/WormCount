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
%   [WORM_SIZE, NUM_WORMS] = count_worms_images(filename)
%
%   [WORM_SIZE, NUM_WORMS] = count_worms_images(filename, minsize, maxsize)
%       minsize - Regions smaller than min_size will be discarded
%           default = 10
%       maxsize - Regions smaller than max_size will be used to determine
%            the size of a single worm
%           default = 80
%
%   [WORM_SIZE, NUM_WORMS] = count_worms_images(filename, minsize, maxsize, debug)
%       debug [0/1] flag outputs various image overlays
%            default = 0 (off)

%% Parse arguments
i_p = inputParser;
i_p.FunctionName = 'count_worms_image';
i_p.addOptional('filename','',@ischar);
i_p.addOptional('minsize',15,@isnumeric); % Regions smaller than this will be discarded
i_p.addOptional('maxsize',100,@isnumeric); % Regions smaller than this will determine single worm size
i_p.addOptional('debug',0,@isnumeric);
i_p.parse(varargin{:});

if ( (isfield(i_p.Results,'filename')) && ~strcmp(i_p.Results.filename,''))
    fullfilename = i_p.Results.filename;
else
    [FileName,PathName,FilterIndex] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
        '*.*','All Files' },'Select Image File');
    fullfilename = [PathName, filesep, FileName];
end

debug = i_p.Results.debug;
min_worm_size = i_p.Results.minsize;
max_worm_size = i_p.Results.maxsize;

% Set Highlight Color
highlightColor = [0 .5 0];

%% Read in image, convert to grayscale
image.info = imfinfo( fullfilename );
image.data = imread(fullfilename);
if ndims(image.data) > 2
    I_gray = rgb2gray(image.data);
else
    I_gray = mat2gray(image.data);
end

%% Setup figure for manual review
names = regexp(image.info.Filename,'(?<path>.*)/(?<filename>.*)','names');
%review_fig = figure('Name', names.filename,'MenuBar','none','ToolBar','none');
review_fig = figure();
plot_fig(image.data, 1);
title('Original Image');

% Find worms
[worm_mask, bg] = find_worms(I_gray, min_worm_size);

% Show review image
reviewimg = imoverlay(image.data, bwperim(worm_mask), highlightColor);
subplot(1,2,2);
[of, h_im] = plot_fig(reviewimg, 2);
truesize(review_fig);
movegui(review_fig,'center')
title('Select regions to ignore, press <ESC> when done');
e = imrect(gca);

%% Looping manual review
while ~isempty(e)
    
    % Replace selected regions with estimation of background
    mask = createMask(e,h_im);
    I_gray(mask) = bg(mask);
    
    % Reanalyze - esp useful when very dark regions are removed from the
    % image, allowing the algorithm to find the lighter worms
    [worm_mask, bg] = find_worms(I_gray, min_worm_size);
    
    % Setup and review results
    % The function IMOVERLAY creates a mask-based image overlay. It takes input
    % image and a binary mask, and it produces an output image whose masked
    % pixels have been replaced by a specified color.
    % MatLab Central -
    % http://www.mathworks.com/matlabcentral/fileexchange/10502
    reviewimg = imoverlay(image.data, bwperim(worm_mask), highlightColor);
    [of, h_im] = plot_fig(reviewimg, 2);
    title('Select regions to ignore, press <ESC> when done');
    e = imrect(gca);
end
close gcf;



%% Estimate single worm size
wormdata = regionprops(bwlabel(worm_mask, 4), 'Area', 'PixelIdxList');
worm_areas = [wormdata.Area];
worm_size = median(worm_areas(worm_areas<=max_worm_size));

%% Estimate number of worms in image
num_worms = round(sum(worm_mask(:))/worm_size);

%% Debug output
if (debug)
    
    % Debug image of all pixels considered worms
    overlay2 = imoverlay(image.data, worm_mask, [.3 1 .3]);
    
    % Single worm image for debugging
    single_worms = false(size(worm_mask));
    single_worms(vertcat(wormdata(worm_areas<max_worm_size).PixelIdxList)) = true;
    RGB_label = label2rgb(bwlabel(single_worms,4), @lines, 'k', 'shuffle');
    
    % Display images and debug info
    fprintf('Estimated size of one worm: %.2f\n', worm_size);
    fprintf('Estimated number of worms: %.0f\n', num_worms);
    figure, imshow(RGB_label);
    figure, imshow(overlay2);
    figure, imshow(image.data);
end

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

% find_worms function identifies dark worms in image
%
%   [MASK, BG] = plot_fig(image, min_worm_size) 
%       Finds worms larger than min_worm_size in image
%       MASK = Logical matrix indicating pixels corresponding to dark areas
%           (worms)
%       BG = Esitmated background of image (from tophat transform)
function [mask, bg] = find_worms(image, min_worm_size)

%% Complement of image
I_comp = imcomplement(image);

%% BG Substract
I_bsub = imtophat(I_comp,strel('disk',15));
bg = imcomplement(I_comp - I_bsub);

%% Enhance

% Noise removal
I_adj = wiener2(I_bsub,[5 5]);

% Gamma correction
I_adj = imadjust(I_adj, [], [], 1.2);

%% Global image threshold using Otsu's method
threshold = graythresh(I_adj);
threshold = max(threshold,.04); % Sanity check on threshold
bw = im2bw(I_adj, threshold);

%% Morphologically open binary image (remove small objects) < min_worm_size
%   Determine connected components (4 pixel neighborhood)
%   Compute area of each component
%   Remove those below specified value
mask = bwareaopen(bw, min_worm_size, 4);


end
