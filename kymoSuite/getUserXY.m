function [XY, IM]   = getUserXY(varargin)
% XY = getUserXY(varargin)
% getUserXY accepts an image and allows the user to draw a select points
% for a segmented line or a spline fit to selected points
%
% Author: Jeffrey J. Nirschl
%         Holzbaur Lab (University of Pennsylvania)
% Date created: 11/04/2016
% Distributable under BSD licence. 
% Copyright (c) 2016 Jeffrey J. Nirschl
% All rights reserved.
% Last modified: 11/21/2016

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
low_in = 0; high_in= inf; % Set range of acceptable input arguments
low_out = 0; high_out= 1; % Set range of acceptable output arguments
narginchk(low_in,high_in);
nargoutchk(low_out,high_out);
clear low_in high_in low_out high_out;

% Parse input and check for errors
[p]         = errorCheck(varargin{:});

% Initialize variables
Fig         = p.FIGURE;
if isempty(p.IM)
    IM      = getimage(Fig);
else
    IM      = p.IM;
end

% Get line coordinates from user input
X       = []; Y     = [];
Xtemp   = []; Ytemp = [];

% Set figure and set pointer
hold on;
set(Fig,'Pointer','crosshair');

% Update figure windows and callbacks
figure(Fig);
colormap(p.CMAP)
char = 0;
set(Fig, 'CurrentCharacter','0'); % Set current char

while p.MAX_ITER~=0
    figure(Fig);
    set(gcf, 'CurrentCharacter','0'); % Reset current char
    
    % Key/ button functions
    keydown     = waitforbuttonpress;
    hold on;
    set(Fig,'Pointer','crosshair');

    char        = get(Fig, 'CurrentCharacter'); % Get current char
    
    axes_handle = gca;
    
    pt          = get(axes_handle, 'CurrentPoint');
%     disp(char)
%     fprintf([char '\n'])
    
    switch char
        case {0,'0'}
            % Continue processing on <Left mouse click>
        case {13, '13'} % & how_many ~= 0)
            % Exit the while loop  on <Return>
            break;
       case {'d','','delete'}
            % Delete previous value
            if ~isempty(Xtemp)
                fprintf('Deleting previous entry...\n')
                
                Xtemp   = Xtemp(1:end-1);
                Ytemp   = Ytemp(1:end-1);
                pt      = [];
                IM      = getimage(Fig);
                [xCols, yRows] = size(IM);
                figure(Fig); cla(gca);

                imshow(IM)
                colormap(p.CMAP)
                set(Fig, 'CurrentCharacter','0'); % Set current char
            else
                fprintf('Error: No points detected for deletion!\n')
                continue
            end
        otherwise
            % Continue processing on all other keys
            pt      = [];
    end
    
    % Update X and Y
    if ~isempty(pt)
        Xtemp   = [Xtemp;pt(1,1)];
        Ytemp   = [Ytemp;pt(1,2)];
    end    
    % Check for errors
    if any(diff(Ytemp)<= 0) 
        fprintf('You cannot select Y coordinates from previous time frames!\n')
        
        % Revert XY coords
        Xtemp = Xtemp(1:end-1);
        Ytemp = Ytemp(1:end-1);
    end

    switch lower(p.LINE(1:2))
        case 'se'
            % Plot current point
            plot(Xtemp,Ytemp,'-c.','LineWidth',2); 
            XX_plot = Xtemp;
            YY_plot = Ytemp;

        case 'sp'
            plot(Xtemp,Ytemp,'.s'); % Plot current point
            if numel(Xtemp)> 1
                XX_plot = nan(100*numel(Xtemp)-1,1);
                YY_plot = nan(100*numel(Ytemp)-1,1);
                for i1 = 2:numel(Xtemp)
                    % Spline interpolation
                    XX_spline = linspace(Xtemp(i1-1),Xtemp(i1))';
                    YY_spline = spline(Xtemp,Ytemp,XX_spline);

                    % Values for plotting
                    Idx = (i1-(i1-1))+(100*(i1-2)):(100*(i1-1));
                    XX_plot(Idx,1) = XX_spline;
                    YY_plot(Idx,1) = YY_spline;
                end
                % Plot current point
                plot(XX_plot,YY_plot,'-c','LineWidth',2);
            end
    end            
    
    % Iteration countdown: The iter starts at -1 by default, which means
    % the while loop will continue indefinitely. Set the optional input
    % argument to the desired number of iterations, if desired.
    p.MAX_ITER = p.MAX_ITER-1;

end

if numel(Xtemp)>1
    switch lower(p.LINE(1:2))
        case 'se'
                for i1 = 2:numel(Xtemp)
                    % Linear interpolation
                    XX_linear = linspace(Xtemp(i1-1),Xtemp(i1))';
                    YY_linear = linspace(Ytemp(i1-1),Ytemp(i1))';
                end
                % Plot current point
                plot(XX_plot,YY_plot,'-c','LineWidth',2);

                % Assign XY output
                XY = [XX_plot(:),YY_plot(:)];
        case 'sp'
                XX_plot = nan(100*numel(Xtemp)-1,1);
                YY_plot = nan(100*numel(Ytemp)-1,1);
                for i1 = 2:numel(Xtemp)
                    % Spline interpolation
                    XX_spline = linspace(Xtemp(i1-1),Xtemp(i1))';
                    YY_spline = spline(Xtemp,Ytemp,XX_spline);

                    % Values for plotting
                    Idx = (i1-(i1-1))+(100*(i1-2)):(100*(i1-1));
                    XX_plot(Idx,1) = XX_spline;
                    YY_plot(Idx,1) = YY_spline;
                end
                plot(XX_plot,YY_plot,'-c', 'LineWidth',2); % Plot current point

                % Assign XY output
                XY = [XX_plot(:),YY_plot(:)];
    end
else
    % Assign XY output
    XY = [Xtemp(:),Ytemp(:)];
end

% Remove NaN values
XY(isnan(XY(:,1)),:)    = [];
set(Fig,'Pointer','Arrow');
end


function [p]     = errorCheck(varargin)
% Subfunction to check for input errors

% Input validation functions
validateNUMERIC     = @(x) all([isscalar(x), isnumeric(x), round(x)==x]);
 
ip = inputParser;
ip.CaseSensitive    = false;
validateLINE        = @(x) any(strcmpi(x,{'segmented','seg','segment','se',...
                        'sp','spl','spline','l','li','line'}));
validateCMAP        = @(x) all([isnumeric(x), size(x,2)==3]);
validateINTERP      = @(x) any(strcmpi(x,{'nearest','bicubic','bilinear'}));
validateFIGURE      = @(x) any(strcmpi(class(x),{'matlab.ui.Figure','matlab.graphics.primitive.Image',...
                        'matlab.graphics.axis.Axes'}));
ip.addParameter('IM',[],@isnumeric);
ip.addParameter('XY',[0,0],@ismatrix); % Mx2 array, M rows by 2 col vectors
ip.addParameter('LINE','seg',validateLINE);
ip.addParameter('INTERP','bicubic',validateINTERP);
ip.addParameter('FIGURE',gcf,validateFIGURE);
ip.addParameter('MAX_ITER',-1,@(x) all([isnumeric(x), isscalar(x)])); % Default goes on for infinity
ip.addParameter('CMAP',colormap(gray),validateCMAP);
ip.parse(varargin{:});
p = ip.Results;

% Validate input
validateattributes(p.XY(:,1),{'numeric'},{'numel',numel(p.XY(:,2))},...
    'lineScan','XY'); % Validate XY coordinate attributes

% Check package dependencies
PKG_DEP     = { {'matlab', 'R2015a'}, {'images', '9.3'}, {'stats', '10.1'},...
              };
ERROR_MSG   = { 'MATLAB %s or higher is required.',...
                'Image Processing Toolbox %s or higher is required.',...
                'Statistics and Machine Learning Toolbox %s or higher is required.',...
                'Signal Processing Toolbox %s or higher is required.'};
for idx      = 1:numel(PKG_DEP)
    if verLessThan(PKG_DEP{idx}{:})
        error(ERROR_MSG{idx},PKG_DEP{idx}{2});
    end
end

% Check custom function dependencies
REQ_FUN = {'','',...
            }; % Check for required M-Files
REQ_FUN_PATH = cell(numel(REQ_FUN),1);
for a1 = 1:numel(REQ_FUN)
    if ~isempty(REQ_FUN{a1})
        REQ_FUN_PATH{a1} = which(REQ_FUN{a1});
        if isempty(REQ_FUN_PATH{a1})
            error([REQ_FUN_PATH{a1} '.m is a required function.']);
        end
    end
end

end


% --- Executes on mouse scroll.
function mouse_scroll_Callback(hObject, eventdata, api)
% hObject    handle to mouse_scroll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

scrollInc   = 30;
    
loc = api.getVisibleLocation();
if eventdata.VerticalScrollCount > 0 % scrolled down
    api.setVisibleLocation(loc(1)           , loc(2) +scrollInc);
elseif eventdata.VerticalScrollCount < 0 % scrolled up
    api.setVisibleLocation(loc(1)           , loc(2) -scrollInc);
end

guidata(hObject);
end

