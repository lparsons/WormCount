% Write mixed data to csv file
% Source: http://francisbarnhart.com/blog/2005/01/19/matlab_pain/
function cellwrite(filename, cellarray, varargin)
switch nargin
    case 2
        delimeter = ',';
        permission = 'at';
    case 3
        delimeter = varargin{1};
        permission = 'at';
    case 4
        delimeter = varargin{1};
        permission = varargin{2};
    otherwise
        error('Unexpected inputs');
end
[rows, cols] = size(cellarray);
fid = fopen(filename, permission);
for i_row = 1:rows
    file_line = '';
    for i_col = 1:cols
        contents = cellarray{i_row, i_col};
        if isnumeric(contents)
            contents = num2str(contents);
        elseif isempty(contents)
            contents = '';
        end
        if i_col < cols
            %file_line = [file_line, contents, ','];
            fprintf(fid, ['%s' delimeter], contents);
        else
            %file_line = [file_line, contents];
            fprintf(fid, '%s\n', contents);
        end
    end
    %count = fprintf(fid, '%s\n', file_line);
end
st = fclose(fid);
