function mergedData = catData(Data, varargin)
% concatenatedData = catData(Data, varargin)
% CATDATA is a function that concatenates two or more vectors of uneven
% length using vertical (default) or horizontal concatenation. 
%
% REQUIRED INPUT:
%       Data        - A cell array of data vectors in the order that they
%                   should be concatenated.
% OPTIONAL INPUT:
%       Method      - String specifying whether to concatenate data as
%                   column vectors ('col' or 'vert') or row vectors ('horz'
%                   or 'row').
% OUTPUT:
%       catData     - An array containing the concatenated data. If the
%                   arrays are different lengths, nan values fill empty
%                   values.
%       
% Author: Jeffrey J. Nirschl
%         Holzbaur Lab (University of Pennsylvania)
% Date created: 08/17/2015
% Distributable under BSD liscence. 
% Copyright (c) 2015 Jeffrey J. Nirschl
% All rights reserved.
% Last modified: 09/23/2015

% Input checking
low_in = 1; high_in= inf; % Set range of acceptable input arguments
low_out = 0; high_out= 1; % Set range of acceptable output arguments
narginchk(low_in,high_in);
nargoutchk(low_out,high_out);
clear low_in high_in low_out high_out;

% Input validation functions
validateMethod = @(x) any(strcmpi(x,{'col','cols','vert','row','rows','horz'}));

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('Data'); % Data vectors as a cell array of vectors
ip.addParameter('Method','col',validateMethod); % Concatenate data horizontally or vertically
ip.parse(Data, varargin{:});
p = ip.Results;

% Error check
if ~iscell(Data)
    error('Input data must be a cell array of vectors\n');
end

% Validate attributes/ classes
Attrb = {'vector'};
Class = {'numeric','cell'};
for i1 = 1:numel(Data)
    validateattributes(Data{i1},Class,Attrb);
end

% Initialize variables
nData = numel(p.Data); % Find the number of data vectors to concatenate
compiledData = cell(1,nData); % Allocate variable
for j1 = 1:(numel(Data))
    compiledData{j1} = Data{j1}; % Assign Data inputs
end
nSize = cellfun(@(x) numel(x),compiledData,'UniformOutput',true);
nMax  = max(nSize);
        
switch lower(p.Method)
    case {'vert','col'}
        nRow = nMax;
        nCol = nData;
        OUTPUT = nan(nRow,nCol);
        for k1 = 1:numel(compiledData)
            OUTPUT(1:nSize(k1),k1) = compiledData{k1};
        end
    case {'horz','row'}
        nRow = nData;
        nCol =nMax;
        OUTPUT = nan(nRow,nCol);
        for k1 = 1:numel(compiledData)
            OUTPUT(k1,1:nSize(k1)) = compiledData{k1};
        end
    otherwise
        error('Unknown contatenation method.')
end

% Assign output variable
mergedData = OUTPUT;

end