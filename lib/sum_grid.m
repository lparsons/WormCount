function counts = sum_grid(count_matrix, category_grid)

% category_grid is struct with field for each category
% each field contains nx2 array with x,y grid coordinates for grids in that
% category
% cg = struct('orig', [3,3; 3,4;4,3;4,4], 'left', [1, 1], 'right', [6,6])
%
% count matrix is n x m matrix with counts for each grid location

if nargin < 2
    cg = zeros(size(count_matrix));
    cg(:,1:2) = 1;
    cg(:,3:4) = 2;
    cg(:,5:6) = 3;
    [x, y] = ind2sub(size(cg),find(cg == 1));
    [x2, y2] = ind2sub(size(cg),find(cg == 2));
    [x3, y3] = ind2sub(size(cg),find(cg == 3));
    category_grid = struct('left', [x,y], 'center', [x2,y2], 'right', [x3,y3]);
end
categories = fieldnames(category_grid);
counts = struct();
for c = 1:length(categories)
    category = categories{c};
    ind = sub2ind(size(count_matrix),  category_grid.(category)(:,1),  category_grid.(category)(:,2));
    counts.(category) = sum(count_matrix(ind));
end


% 
% ind = sub2ind(size(worm_counts),  cg.orig(:,1),  cg.orig(:,2))
% sum(worm_counts(ind))

