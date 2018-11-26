function [DUP_VAL, DUP_IDX, ORG_IDX, COUNTS] = findDuplicates(X, varargin)
% [DUP_VAL, DUP_IDX, ORG_IDX, COUNTS] = findDuplicates(X, varargin)
% findDuplicates is a function that accepts an input vector, X, and returns
% the values and indices for duplicates and the original (first or last
% element) as well as a count of the number of duplicates for each unique
% value in X. It can handle nan/ inf values.
% 
% Author: Jeffrey J. Nirschl
%         Holzbaur Lab (University of Pennsylvania)
% Date created: 04/12/2016
% Distributable under BSD licence. 
% Copyright (c) 2016 Jeffrey J. Nirschl
% All rights reserved.
% Last modified: 07/23/2016

% Covered under a BSD license
% Copyright (c) 2016 Jeffrey J. Nirschl
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, IM, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Input checking
low_in = 1; high_in= inf; % Set range of acceptable input arguments
low_out = 0; high_out= 4; % Set range of acceptable output arguments
narginchk(low_in,high_in);
nargoutchk(low_out,high_out);
clear low_in high_in low_out high_out;

% Parse input and check for errors
[X, p]      = errorCheck(X, varargin{:});

% Body of function
nanVal      = isnan(X);     % Store  nan indices
infVal      = isinf(X);     % Store  inf indices
X           = X(~isnan(X)); % Remove nan values
X           = X(~isinf(X)); % Remove inf values
[N, EDGES]  = histcounts(X, [(unique(X)); X(end)+1]); % Bin X using edges defined by (unique(X))

DUP         = EDGES(N' > 1);% These are the duplicated values
if ~isempty(DUP)
    IDX     = cell(numel(DUP),1);
    ORG_IDX = cell(numel(DUP),1);
    DUP_IDX = cell(numel(DUP),1);
    for idx = 1:numel(DUP)
        IDX{idx}    = (ismember(X, DUP(idx)));
        ORG_IDX{idx}= find(IDX{idx},1, p.KEEP);
    end
    
    IDX         = any(cell2mat(IDX'),2);
    ORG_IDX     = cell2mat(ORG_IDX);
    IDX(ORG_IDX)= 0;
    DUP_VAL     = unique(X(IDX));
    DUP_IDX     = find(IDX);
else
    IDX         = false(size(X));
    ORG_IDX     = find(ismember(X,unique(X)),1,p.KEEP);
    DUP_VAL     = [];
    DUP_IDX     = [];
end

% Check if final edge exceeds the unique values. Remove final EDGE element 
% if needed.
if and(max(EDGES)>max(unique(X)), numel(EDGES)> numel(unique(X)))
    EDGES(end)= [];         % Remove end index
end

% Assign output

if and(any(nanVal),any(infVal))
    COUNTS      = catData([{[p.NAN_LABEL; EDGES; inf]}, {[sum(nanVal); N'; sum(infVal)]}],  'Method','Col');
elseif any(nanVal)
    COUNTS      = catData([{[p.NAN_LABEL; EDGES]},      {[sum(nanVal); N'             ]}],  'Method','Col');
elseif any(infVal)
    COUNTS      = catData([{[     EDGES; inf]},         {[             N'; sum(infVal)]}],  'Method','Col');
else
    COUNTS      = catData([{[     EDGES     ]},         {[             N'             ]}],  'Method','Col');
end

end


function [X, p] = errorCheck(X, varargin)
% Subfunction to check for input errors

% Input validation functions
validateX           = @(x) all([isvector(x)]);
validateNAN_LABEL   = @(x) all([isscalar(x), isnumeric(x)]);
validateKEEP        = @(x) any(strcmpi(x,{'first','last'}));

% Parse inputs
ip = inputParser;
ip.CaseSensitive= false;
ip.addRequired('X');
ip.addParameter('NAN_LABEL',nan, validateNAN_LABEL);
ip.addParameter('INF_LABEL',inf, validateNAN_LABEL);
ip.addParameter('KEEP','first', validateKEEP);
ip.parse(X, varargin{:});
p = ip.Results;

% Check dependencies
REQ_FUN = {''}; % Check for required M-Files
REQ_FUN_PATH = cell(numel(REQ_FUN),1);
for a1 = 1:numel(REQ_FUN)
    if ~isempty(REQ_FUN{a1})
        REQ_FUN_PATH{a1} = which(REQ_FUN{a1});
        if isempty(REQ_FUN_PATH{a1})
            error([REQ_FUN_PATH{a1} '.m is a required function.']);
        end
    end
end

% Validate attributes
validateattributes(X,{'numeric'},{'vector','nonempty','nonsparse'},...
                   'findDuplicates','X',1)
               
p.X_orig    = X;

if ~iscolumn(X)
    X       = sort(X');
else
    X       = sort(X);
end

end
