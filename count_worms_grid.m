%   Copyright 2011 Lance R. Parsons, Lewis-Sigler Institute for Integrative Genomics, Princeton University
%
%   Licensed under the BSD 2-Clause License 
%   http://www.opensource.org/licenses/BSD-2-Clause

function [worm_counts, worm_size] = count_worms_grid(image, numrows, numcols, varargin)
% TODO count_worms_grid documentation
%
%

[pathstr, name, ext] = fileparts(mfilename('fullpath')) ; %#ok<NASGU,ASGLU>
addpath(genpath([pathstr filesep 'lib']));

ip = inputParser;
ip.FunctionName = 'count_worms_grid';
ip.addRequired('image', @(x)exist(x,'file') );
ip.addRequired('numcols', @isnumeric);
ip.addRequired('numrows', @isnumeric);
ip.addOptional('minsize',5,@isnumeric); % Regions smaller than this will be discarded
ip.addOptional('maxsize',30,@isnumeric); % Regions smaller than this will determine single worm size
ip.addParamValue('debug',false,@islogical);
ip.parse(image, numrows, numcols, varargin{:});


%% Even background and find dark regions
I  = imread(ip.Results.image);
Ibh = imbothat(I, strel('disk',5,0));
BW = im2bw(Ibh,graythresh(Ibh));
if ip.Results.debug
    figure, imshow(imoverlay(I,BW,[255 0 0]), 'Border','tight')
end

%% Find Top Left and Bottom Right of Grid and crop image
cc = bwconncomp(BW);
stats = regionprops(cc);
idx = find([stats.Area] == max([stats.Area]));
grid = ismember(labelmatrix(cc), idx);
if ip.Results.debug
    figure, imshow(imoverlay(I,grid,[0 255 0]), 'Border','tight')
end

grid_stats = regionprops(grid, 'BoundingBox', 'Extrema');

grid_crop_image = imcrop(I,grid_stats.BoundingBox);

grid_crop_bw = imcrop(grid,grid_stats.BoundingBox);
%grid_crop_bh = imbothat(grid_crop_image, strel('disk',5,0));
%grid_crop_bw = im2bw(grid_crop_bh,graythresh(grid_crop_bh));

if ip.Results.debug
    figure, imshow(imoverlay(grid_crop_image,grid_crop_bw,[255 0 0]), 'Border','tight')
end

%% Hough Transform
[H,T,R] = hough(grid_crop_bw);
if ip.Results.debug
    figure, imshow(H,[],'XData',T,'YData',R,...
        'InitialMagnification','fit', 'Border','tight');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
end
P  = houghpeaks(H,50,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
if ip.Results.debug
    plot(x,y,'s','color','white');
end

%% Find lines and plot them
lines = houghlines(grid_crop_bw,T,R,P,'FillGap',50,'MinLength',50);
if ip.Results.debug
    figure, imshow(grid_crop_image, 'Border','tight'), hold on
    max_len = 0;
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
        if abs(xy(1,1)-xy(2,1)) > abs(xy(1,2)-xy(2,2))
            c = 'green';
        else
            c = 'red';
        end
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',c);
        
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
        
        % Determine the endpoints of the longest line segment
        len = norm(lines(k).point1 - lines(k).point2);
        if ( len > max_len)
            max_len = len;
            xy_long = xy;
        end
    end
    
    % highlight the longest line segment
    plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
end

%% Find horizonal and vertical lines, and midpoints, order by midpoints (top->bottom)
approximate_grid_height = size(grid_crop_image,1)/numrows;
approximate_grid_width = size(grid_crop_image,2)/numcols;

vertlines = struct('point1',{},'point2',{},'theta',{},'rho',{},'sortposition',{},'length',{});
horzlines = struct('point1',{},'point2',{},'theta',{},'rho',{},'sortposition',{},'length',{});
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    if abs(xy(1,1)-xy(2,1)) < abs(xy(1,2)-xy(2,2)) % Vertical line
        lines(k).sortposition = (xy(1,1)+xy(2,1))/2;
        lines(k).length = norm(lines(k).point1 - lines(k).point2);
        vertlines(end+1) = lines(k);
        %vertlines(k).position = (xy(1,2)+xy(2,2))/2;
    else                                           % Horizontal line
        lines(k).sortposition = (xy(1,2)+xy(2,2))/2;
        lines(k).length = norm(lines(k).point1 - lines(k).point2);
        horzlines(end+1) = lines(k);
        %horzlines(k).position = (xy(1,1)+xy(2,1))/2;
    end
end
horzlines = nestedSortStruct(horzlines, 'sortposition');
vertlines = nestedSortStruct(vertlines, 'sortposition');

%% Split lines into vertical and horizonal, removing duplicates
lastpos = -approximate_grid_height;
grid_horzlines = struct('point1',{},'point2',{});
for k = 1:length(horzlines)
    if (horzlines(k).sortposition - lastpos) < approximate_grid_height * 0.75
        % Not sufficient distance, likely not grid line
        continue
    else
        grid_horzlines(end+1).point1 = horzlines(k).point1;
        grid_horzlines(end).point2 = horzlines(k).point2;
    end
    lastpos = horzlines(k).sortposition;
end

lastpos = -approximate_grid_height;
grid_vertlines = struct('point1',{},'point2',{});
for k = 1:length(vertlines)
    if (vertlines(k).sortposition - lastpos) < approximate_grid_width * 0.75
        % Not sufficient distance, likely duplicate grid line
        continue
    else
        grid_vertlines(end+1).point1 = vertlines(k).point1;
        grid_vertlines(end).point2 = vertlines(k).point2;
    end
    lastpos = vertlines(k).sortposition;
end

if ip.Results.debug
    figure, imshow(grid_crop_image, 'Border','tight'), hold on
    plot_lines(grid_horzlines, 'green');
    plot_lines(grid_vertlines, 'blue');
    hold off
end

%% Extend lines (horizontal)
extension_factor = 0.1; % Amount to extend lines beyone image to ensure intersection
width_extension = size(grid_crop_image,2) * extension_factor;
max_width = size(grid_crop_image,2) + width_extension;
min_width = 1 - width_extension;
fullgrid_horzlines = struct('point1',{},'point2',{});

for k = 1:length(grid_horzlines)
    this_line = grid_horzlines(k);
    p1 = this_line.point1;
    p2 = this_line.point2;
    deltax = p2(1) - p1(1);
    if deltax < 0
        p1 = this_line.point2;
        p2 = this_line.point1;
        deltax = p2(1) - p1(1);
    end
    deltay = p2(2) - p1(2);
    % Extend left
    left_end_point_deltax = p1(1) - min_width;
    left_end_point_deltay = deltay * (left_end_point_deltax/deltax);
    left_end_point = [min_width, p1(2) - left_end_point_deltay];
    % Extend right
    right_end_point_deltax = max_width - p1(1);
    right_end_point_deltay = deltay * (right_end_point_deltax/deltax);
    right_end_point = [max_width, p1(2) + right_end_point_deltay];
    new_line.point1 = left_end_point;
    new_line.point2 = right_end_point;
    
    fullgrid_horzlines(end+1) = new_line;
end

%% Extend Vertical
height_extension = size(grid_crop_image,1) * extension_factor;
max_height = size(grid_crop_image,1) + height_extension;
min_height = 1 - height_extension;
fullgrid_vertlines = struct('point1',{},'point2',{});

for k = 1:length(grid_vertlines)
    this_line = grid_vertlines(k);
    p1 = this_line.point1;
    p2 = this_line.point2;
    deltay = p2(2) - p1(2);
    if deltay < 0
        p1 = this_line.point2;
        p2 = this_line.point1;
        deltay = p2(2) - p1(2);
    end
    deltax = p2(1) - p1(1);
    % Extend up
    top_end_point_deltay = p1(2) - min_height;
    top_end_point_deltax = deltax * (top_end_point_deltay/deltay);
    top_end_point = [p1(1) - top_end_point_deltax, min_height];
    % Extend right
    bottom_end_point_deltay = max_height - p1(2);
    bottom_end_point_deltax = deltax * (bottom_end_point_deltay/deltay);
    bottom_end_point = [p1(1) + bottom_end_point_deltax, max_height];
    % new point
    new_line.point1 = top_end_point;
    new_line.point2 = bottom_end_point;
    fullgrid_vertlines(end+1) = new_line;
end

%% Show extended lines
if ip.Results.debug
    figure, imshow(grid_crop_image, 'Border','tight'), hold on
    plot_lines(fullgrid_horzlines, 'green');
    plot_lines(fullgrid_vertlines, 'blue');
    hold off
end

%% Find grid intersections
grid_corners = cell((numrows+1),(numcols+1));
for j = 1:numrows+1
    for k = 1:numcols+1
        [y,x] = intersections([fullgrid_horzlines(j).point1(1), fullgrid_horzlines(j).point2(1), ...
            fullgrid_vertlines(k).point1(1), fullgrid_vertlines(k).point2(1)],...
            [fullgrid_horzlines(j).point1(2), fullgrid_horzlines(j).point2(2), ...
            fullgrid_vertlines(k).point1(2), fullgrid_vertlines(k).point2(2)]);
        grid_corners{j,k} = [y,x];
    end
end

%% Make mask from grid corners
f = figure('Visible', 'off');
imshow(grid_crop_image);
masks = cell((numrows),(numcols));
for j = 1:numrows
    for k = 1:numcols
        topleft = grid_corners{j,k};
        topright = grid_corners{j,k+1};
        bottomleft = grid_corners{j+1,k};
        bottomright = grid_corners{j+1,k+1};
        roi = impoly(gca,[topleft;bottomleft;bottomright;topright]);
        mask = createMask(roi);
        masks{j,k}=mask;
    end
end
close(f);

%% Adjust mask to be inside grid lines (hopefully not touching)
inside_grid_masks = cell((numrows),(numcols));
inside_all_grid_mask = false(size(grid_crop_image));
for j = 1:numrows
    for k = 1:numcols
        %inside_grid_mask = masks{j,k} & (~grid_crop_bw);
        inside_grid_mask = bwmorph(masks{j,k},'shrink', 2);
        gm1 = inside_grid_mask & (~grid_crop_bw);
        gm2 = bwareaopen(gm1,round((approximate_grid_height*approximate_grid_width)*.2),4);
        %gm3 = bwmorph(gm2,'fill');
        gm3 = bwmorph(gm2,'spur');
        %ch = bwconvhull(gm2);
        %gm5 = bwmorph(gm3,'open', Inf);
        inside_grid_masks{j,k} = gm3 & inside_grid_mask;
        inside_all_grid_mask(inside_grid_masks{j,k}) = true;
    end
end
if ip.Results.debug
    figure; imshow(imoverlay(grid_crop_image,inside_all_grid_mask,[0,255,0]),'Border','tight');
end

%% Count worms
worm_counts = zeros((numrows),(numcols));
threshold_strength = 0.8;
%worm_mask = find_worms_image(grid_crop_image, 40, 100);
grid_crop_bh = imbothat(grid_crop_image, strel('disk',5,0));

% thresh_image = grid_crop_bh;
% thresh_image(grid_crop_bw)=median(grid_crop_bh(:));
% thresh = graythresh(thresh_image);
% worm_mask = im2bw(grid_crop_bh,thresh);
% worm_mask = worm_mask & inside_all_grid_mask;
% %worm_mask = bwareaopen(worm_mask, 3);

worm_mask = im2bw(grid_crop_bh,graythresh(grid_crop_bh) * threshold_strength);
worm_mask_disp = bwareaopen(inside_all_grid_mask & worm_mask, ip.Results.minsize);
if ip.Results.debug
    figure; imshow(imoverlay(grid_crop_image,worm_mask_disp,[0,255,0]),'Border','tight');
end
[num_worms, worm_size] = count_worms_image('worm_mask', worm_mask, 'minsize', ip.Results.minsize, 'maxsize', ip.Results.maxsize); %#ok<ASGLU>
for j = 1:numrows
    for k = 1:numcols
        this_worm_mask = worm_mask & inside_grid_masks{j,k};
        %figure, imshow(this_worm_mask); pause;
        num_worms = count_worms_image('worm_mask', this_worm_mask, 'avg_worm_size', worm_size);
        worm_counts(j,k) = num_worms;
    end
end


%% Save overlay image

[PATHSTR,NAME,EXT] = fileparts(ip.Results.image); %#ok<NASGU>
data_path = [PATHSTR filesep NAME '_results'];
[s,mess,messid] = mkdir(data_path); %#ok<NASGU,ASGLU>
image_overlay_filename = fullfile(data_path, [NAME '_overlay.png']);

overlay_figure = figure('Visible', 'off');
imshow(imoverlay(grid_crop_image, bwperim(worm_mask_disp), [0 1 0]));

title_text = sprintf('Worm size = %s', num2str(worm_size));

text(floor(size(grid_crop_image,1)/2), ...
    5, ...
    title_text, ...
    'BackgroundColor', [.8 .7 .7], ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center' ...
    )

% Print Overlay Summary Figure
screen_DPI = get(0, 'ScreenPixelsPerInch');
set(overlay_figure, 'Units', 'pixels', 'Position', [32, 32, size(grid_crop_image, 2), size(grid_crop_image, 1)]);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gca, 'Units', 'points');
set(overlay_figure, 'Units', 'points', 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
print(overlay_figure, '-dpng', sprintf('-r%d', screen_DPI), image_overlay_filename);

close(overlay_figure);

end
