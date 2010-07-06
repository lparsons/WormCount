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

% Read in image
image.info = imfinfo( fullfilename );
image.data = imread(fullfilename);

if ndims(image.data) > 2
    I_gray = rgb2gray(image.data);
else
    I_gray = mat2gray(image.data);
end

%I_sc = imadjust(Ig,stretchlim(Ig),[0,1],.5);

%% Setup figure
names = regexp(image.info.Filename,'(?<path>.*)/(?<filename>.*)','names');
review_fig = figure('Name', names.filename,'MenuBar','none','ToolBar','none');
plot_fig(image.data, 1);
title('Original Image');


%% Manual review
[worm_mask, bg] = find_worms(I_gray, min_worm_size);
reviewimg = imoverlay(image.data, bwperim(bwmorph(worm_mask,'thicken',1)), [.3 .8 .3]);

subplot(1,2,2);
[of, h_im] = plot_fig(reviewimg, 2);
truesize(review_fig);
movegui(review_fig,'center')


title('Select regions to ignore, press <ESC> when done');
e = imrect(gca);

while ~isempty(e)
    mask = createMask(e,h_im);
    
    I_gray(mask)=bg(mask);
    
    [worm_mask, bg] = find_worms(I_gray, min_worm_size);
    
    %reviewimg = imoverlay(I_gray, bwperim(bwmorph(worm_mask,'thicken',1)), [.3 .8 .3]);
    reviewimg = imoverlay(image.data, bwperim(bwmorph(worm_mask,'thicken',1)), [.3 .8 .3]);
    
    [of, h_im] = plot_fig(reviewimg, 2);
    title('Select regions to ignore, press <ESC> when done');
    e = imrect(gca);
end
close gcf;

%% Morphological closing (dilation followed by erosion).
% bw5 = bwmorph(bw4, 'close');

%% With n = Inf, thickens objects by adding pixels to the exterior of
% objects until doing so would result in previously unconnected objects
% being 8-connected. This option preserves the Euler number.
% bw6 = bwmorph(bw5, 'thicken', cellBorderThicken);

% Get a binary image containing only the perimeter pixels of objects in
% the input image BW1. A pixel is part of the perimeter if it is nonzero
% and it is connected to at least one zero-valued pixel. The default
% connectivity is 4.
%bw5_perim = bwperim(bw5);

% The function IMOVERLAY creates a mask-based image overlay. It takes input
% image and a binary mask, and it produces an output image whose masked
% pixels have been replaced by a specified color.
% MatLab Central -
% http://www.mathworks.com/matlabcentral/fileexchange/10502
%overlay1 = imoverlay(I_bsub, bw, [.3 1 .3]);
overlay2 = imoverlay(image.data, worm_mask, [.3 1 .3]);


%% Estimate worm size
%cc = bwconncomp(worm_mask, 4);
wormdata = regionprops(bwlabel(worm_mask, 4), 'Area', 'PixelIdxList');
worm_areas = [wormdata.Area];
worm_size = median(worm_areas(worm_areas<=max_worm_size));

single_worms = false(size(worm_mask));
single_worms(vertcat(wormdata(worm_areas<max_worm_size).PixelIdxList)) = true;
RGB_label = label2rgb(bwlabel(single_worms,4), @lines, 'k', 'shuffle');

%figure, imshow(single_worms);

num_worms = round(sum(worm_mask(:))/worm_size);

% Debug output
if (debug)
    fprintf('Estimated size of one worm: %.2f\n', worm_size);
    fprintf('Estimated number of worms: %.0f\n', num_worms);
    
    figure, imshow(RGB_label);
    figure, imshow(overlay2);
    %figure, imshow(overlay1);
    %figure, imshow(I_bsub);
    figure, imshow(image.data);
end

%figure, imshow(overlay2);
end

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

function [mask, bg] = find_worms(image, min_worm_size)

I_comp = imcomplement(image);

%se = strel('disk',15);
%I_comp = imsubtract(imadd(I_comp,imtophat(I_comp,se)), imbothat(I_comp,se));

%% BG Substract
% background = imopen(I_comp,strel('disk',15));
% I_bsub = mat2gray(I_comp - background);
I_bsub = imtophat(I_comp,strel('disk',15));
bg =imcomplement(I_comp - I_bsub);

%% Contrast Enchance

% Stretch Limits
%I_adj = imadjust(I_bsub,stretchlim(I_bsub),[min(I_bsub(:)),max(I_bsub(:))],2);
%I_sm = medfilt2(I_adj);
%I_adj = imadjust(I_bsub,[],[0,1],1);

% Tophat/Bottomhat
%se = strel(ones(size(I_bsub)));
%I_adj = imsubtract(imadd(I_bsub,imtophat(I_bsub,se)), imbothat(I_bsub,se));

% Adaptive contrast enhancement
%I_adj = adapthisteq(I_bsub);


%I_adj = imadjust(I_bsub);
I_adj = I_bsub;

%% Global image threshold using Otsu's method
threshold = graythresh(I_adj);
%threshold = max(threshold,.2); % Sanity check on threshold
threshold = max(threshold,.04); % Sanity check on threshold

bw = im2bw(I_adj, threshold);
%
% [bw, threshold] = edge(I_adj,'canny');
%
% bw2 = imfill(bw,'holes');
%bw2 = imclose(bw, strel('disk', 10));

%bw = imextendedmax(I_adj,.99,8);
%bw = im2bw(I_adj, .2);

%bw2 = bwmorph(bw,'close');

%% Cleanup thresholded image
%   Fill in holes - pixes that cannot be reached by filling in the
%   background from the edge of the image (using 4 pixel neighborhood)
% bw2 = imfill(bw,'holes');

%% Morphological opening using a 5x5 block.  The morphological open
% operation is an erosion followed by a dilation, using the same
% structuring element for both operations.
% morphOpenStruct = ones(2,2);
% bw3 = imopen(bw, morphOpenStruct);

%% Morphologically open binary image (remove small objects) < min_worm_size
%   Determine connected components (4 pixel neighborhood)
%   Compute area of each component
%   Remove those below specified value
mask = bwareaopen(bw, min_worm_size, 4);


end