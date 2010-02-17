function [worm_size, num_worms] = count_worms_image(varargin)

i_p = inputParser;
i_p.FunctionName = 'count_worms_image';
i_p.addRequired('filename',@ischar);
i_p.addOptional('debug',0,@isnumeric);
i_p.parse(varargin{:});

filename = i_p.Results.filename;
debug = i_p.Results.debug;


min_worm_size = 10; % Regions smaller than this will be discarded
max_worm_size = 100; % Regions smaller than this will determine single worm size

image.info = imfinfo( filename );
image.data = imread(filename);

I = image.data;
I_sc = mat2gray(I);
I_comp = imcomplement(I_sc);

%% BG Substract
background = imopen(I_comp,strel('disk',15));
I_bsub = I_comp - background;


%% Global image threshold using Otsu's method
threshold = graythresh(I_bsub);
threshold = max(threshold,.2); % Sanity check on threshold
bw = im2bw(I_bsub, threshold); 


%% Cleanup thresholded image
%   Fill in holes - pixes that cannot be reached by filling in the
%   background from the edge of the image (using 4 pixel neighborhood)
% bw2 = imfill(bw,'holes');

%% Morphological opening using a 5x5 block.  The morphological open
% operation is an erosion followed by a dilation, using the same
% structuring element for both operations.
% morphOpenStruct = ones(2,2);
% bw3 = imopen(bw, morphOpenStruct);

%% Morphologically open binary image (remove small objects) < 40
%   Determine connected components (4 pixel neighborhood)
%   Compute area of each component
%   Remove those below specified value
bw4 = bwareaopen(bw, min_worm_size, 4);

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
worm_mask = bw4;
overlay1 = imoverlay(I_bsub, bw, [.3 1 .3]);
overlay2 = imoverlay(I, worm_mask, [.3 1 .3]);


%% Estimate worm size
%cc = bwconncomp(worm_mask, 4);
wormdata = regionprops(bwlabel(worm_mask, 4), 'Area', 'PixelIdxList');
worm_areas = [wormdata.Area];
worm_size = median(worm_areas(worm_areas<max_worm_size));

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
    figure, imshow(overlay1);
    figure, imshow(I_bsub);
    figure, imshow(image.data);
end

