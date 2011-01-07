function [ ci ] = chemotaxis_index( count_data )
%CHEMOTAXIS_INDEX Calculates cheomotaxis index
%   
%   CHEMOTAXIS_INDEX = count_worms_plate(count_data)
%   
%       count_data is a struct with four fields:
%           total, ori, but, eth
%
%       CHEMOTAXIS_INDEX = (but-eth)/(total-ori)

p = inputParser;
p.FunctionName = 'chemotaxis_index';
p.addRequired('count_data', @isstruct);
p.parse(count_data);

% Check that all required fileds are specified in count_data structure
required_fields = {'tot', 'ori', 'but', 'eth'};
for f=1:size(required_fields,2)
    if ~isfield(count_data, required_fields{f})
        error('Missing field: %s', required_fields{f})
    end
end

ci = (count_data.but - count_data.eth) / (count_data.tot - count_data.ori);
end

