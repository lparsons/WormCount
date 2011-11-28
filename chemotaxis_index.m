%   Copyright 2011 Lance R. Parsons
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

function [ ci ] = chemotaxis_index( count_data )
%CHEMOTAXIS_INDEX Calculates cheomotaxis index
%   
%   CHEMOTAXIS_INDEX = chemotaxis_index(count_data)
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

