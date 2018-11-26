% Kymograph Analysis function
% A function that returns a structure array containing the parameters
% essential for kymograph analysis. The program then graphs the results by 
% neuron in MATLAB to assess quality of the experiment and inter-neuron
% variability. It then exports data to both a .MAT file and an Excel file.
% See the instruction manual titled Kymo_Suite for more information
%
% 
% All kymographs are made from the retrograde direction to the
% anterograde direction so throughout the analysis all positive
% movements/directions/velocities are anterograde and negative
% movements/directions/velocities are retrograde.
%
% Author: Jeffrey J. Nirschl
%         Holzbaur Lab (University of Pennsylvania)
% Date created: 09/02/2013
% Distributable under BSD licence. 
% Copyright (c) 2017 Jeffrey J. Nirschl
% All rights reserved.
% Last modified: 02/13/2017

% Update list
%   Version 1.5 on 10/30/2013
%       Updated: calculation of direction changes, binning direction changes by
%       net displacement >10um for anterograde etc. Also direction change
%       calculation updated to define a direction change as a change in
%       instantaneous velocity with at least 1um movement in the new direction.
%       Calculations also include run length after direction change.
%   Version 2.0 on 11/15/2013
%       Updated 6/9/2014 to fix GUI bugs
%   Version 2.5 on 12/19/2016
%       Updated to write interpolated pixel values to Excel. Previous
%       versions only output the interpolated instantaneous velocities.
% Distributable under BSD liscence. 
% Copyright (c) 2017 Jeffrey J. Nirschl
% All rights reserved.
%
% Redistribution and use in source and binary forms are permitted provided
% that the above copyright notice and this paragraph are duplicated in all
% such forms and that any documentation, advertising materials, and other
% materials related to such distribution and use acknowledge that the
% software was developed by Jeffrey Nirschl, Holzbaur lab, University of
% Pennsylvania (hereafter referred to as The Author).
% 
% The name of The Author  may not be used to endorse or promote products
% derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED
% WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


function output =kymo_analysis(varargin)
% dbstop if error

% Initialize variables and start timer
tic % Start timer
Dir = ['C:\Users\Holzbaur-Users\Documents\Matlab\',...
    'Select save directory.'];
SaveDir = uigetdir(Dir ,'Select save directory.'); cd(SaveDir);% Creates directory as current working directory
datetime = datestr(now,'dd-mmm-yyyy HH:MM:SS');

% Set options for user input dialog boxes
options.Resize='on'; options.WindowStyle='normal'; options.Interpreter='tex';
inputdlg_text= {'Filename','um/ Pixel','Sec/ Frame','Pause = Inst. Vel < ___ (um/sec)','Direction Change Threshold = __ um',};
try
    inputdlg_user= {varargin{1,1}{1,1}{1,1},varargin{1,1}{1, 1}{1, 2}(12:end),varargin{1,1}{1, 1}{1, 3}(11:end),'0.100', '0.5'};
    start = 1;
catch
    if isempty(varargin{1,1}{2,1})==0;
       inputdlg_user= {varargin{1,1}{2,1}{1,1},varargin{1,1}{2, 1}{1, 2}(12:end),varargin{1,1}{2, 1}{1, 3}(11:end),'0.100', '0.5'};
       start = 2;
    elseif isempty(varargin{1,1}{3,1})==0;
       inputdlg_user= {varargin{1,1}{3,1}{1,1},varargin{1,1}{3, 1}{1, 2}(12:end),varargin{1,1}{3, 1}{1, 3}(11:end),'0.100', '0.5'};
       start = 3;
    else
       inputdlg_user= {varargin{1,1}{4,1}{1,1},varargin{1,1}{4, 1}{1, 2}(12:end),varargin{1,1}{4, 1}{1, 3}(11:end),'0.100', '0.5'};
       start = 4;
    end
end
inputdlg_default= {'Enter title here','0.0715','0.5','0.100', '0.5'};
% ABOVE--This creates a dialog box for userinput parameters with default 
% "0.147um/pixel"AM, "0.0715um/pixel"JN, "0.366s/frame"AM, "0.5s/frame"JN,
% and instantanous velocities < 100nm/sec defined as a pause.
% Update 10/17/2013 --I am changing the threshold for
% defining a pause. Velocities of >84nm/sec travel >10microns in a two
% minute video. Since my threshold for displacement in a 2 minute video is
% 10 microns I feel that pause definition should be consistent with this.
% An instantaneous velocity < 83nm/sec will not travel > 10 microns in a 2
% minute video and I will use this more stringent definition for pauses.

for z=1:size(varargin,2)
    raw = varargin{1,z}; % Import variables into 'raw'
    exp_filename = raw{start,1}{1,1};
    output = cell(0,1);
    try
        param= inputdlg(inputdlg_text,'Input Parameters',1,inputdlg_user,options);
    catch err
        param= inputdlg(inputdlg_text,'Input Parameters',1,inputdlg_default,options);
    end

    fname = [param{1,1} '.mat']; % Filename from user input dialog box
    pixel = str2num(param{2,1}); % Pixel conversion from user input dialog box
    tframe= str2num(param{3,1}); % Sec/ frame from user input dialog box
    pause= str2num(param{4,1}); % Pause length thresh. from user input dialog box
    dir_thresh= str2num(param{5,1}); % To be a direction change the vesicle 
    % must change inst. vel sign AND move more than 1um in the new
    % direction else it is defined as a pause 

    sz= size(raw); % Size of input file
    interp_btn = questdlg('Does your data require interpolation of intermediate X coordinates?');

    % Create Nan matrices
%     output.initial_xy = nan(sz(1),1);
    output.initial_y = nan(sz(1),sz(2)-1);
    output.cumdist= nan(sz(1),sz(2)-1);
    output.netdist= nan(sz(1),sz(2)-1);
    output.cumtime= nan(sz(1),sz(2)-1);
    output.tantero= nan(sz(1),sz(2)-1);
    output.tretro= nan(sz(1),sz(2)-1);
    output.tpause= nan(sz(1),sz(2)-1);
    output.dir_run_num = nan(sz(1),sz(2)-1);
    output.dir_totalswitches= nan(sz(1),sz(2)-1);
    output.dir_antero_switches= nan(sz(1),sz(2)-1);
    output.dir_retro_switches= nan(sz(1),sz(2)-1);
    output.dir_pause_switches= nan(sz(1),sz(2)-1);
    output.dir_less_thresh= nan(sz(1),sz(2)-1);
    output.dir_netmovement= nan(sz(1),sz(2)-1);
    
%     output.instvel{i,idx}= NaN;
%     output.instaccl{i,idx}= NaN;
%     output.dir_indices{i,idx} = NaN;
%     output.dir_sign{i,idx} = NaN;
%     output.dir_run_length{i,idx} = NaN;
%     output.dir_pause{i,idx} = NaN;
%     output.pause_length{i,idx}= NaN;
    
    for i=start:sz(1); % 'For' loop over all rows
        idx= 0;  % Index for saving Output 
        for j=2:sz(2) % 'For' loop over all columns minus the first column
            % Allocate outputs to NaN, they will be overwritten with data
            % below
            idx= idx + 1;
            output.initial_xy{i,idx}= NaN; output.initial_y(i,idx) = NaN;
            output.interpdist{i,idx}= NaN; output.interptime{i,idx}= NaN;
            output.cumdist(i,idx)   = NaN; output.netdist(i,idx)   = NaN;
            output.cumtime(i,idx)   = NaN; output.instvel{i,idx}   = NaN;
            output.instaccl{i,idx}  = NaN; output.tantero(i,idx)   = NaN;
            output.tretro(i,idx)    = NaN; output.tpause(i,idx)    = NaN;
            output.dir_indices{i,idx}= NaN;output.dir_sign{i,idx}  = NaN;
            output.dir_run_num(i,idx)= NaN;output.dir_run_length{i,idx} = NaN;
            output.dir_pause{i,idx} = NaN; output.dir_netmovement(i,idx) = NaN;
            output.dir_totalswitches(i,idx)     = NaN;
            output.dir_antero_switches(i,idx)   = NaN;
            output.dir_retro_switches(i,idx)    = NaN;
            output.dir_pause_switches(i,idx)    = NaN;
            output.dir_less_thresh(i,idx)       = NaN;
            output.pause_length{i,idx}          = NaN;
            
            
            % Process data. Skip empty columns or columns with less than 2 values
            if and(isempty(raw{i,j}) == 0, size(raw{i,j},1) > 2)
                try
                    if strcmp(interp_btn,'Yes') ==1; % interpolation may not be required
                        [~, DUP_IDX]= findDuplicates(raw{i,j}(:,2));
                        if isempty(DUP_IDX)
                            interptime = [raw{i,j}(1,2):1:raw{i,j}(end,2)]'; % Interpolate values of y
                            interpdist = interp1(raw{i,j}(:,2),raw{i,j}(:,1),interptime,'linear','extrap');% Interpolate values of x
                        else
                            fprintf(['Duplicate time values found in Row %d Col %d.', ...
                                ' Averaging values to one position.\n'],i,j)
                            [~,groupCell]   = groupUnique(fliplr(raw{i,j}));
                            groupCell       = cell2mat(cellfun(@(x) [x(1) ,nanmean(x(2:end))], groupCell, 'unif',false));
                            raw{i,j}        = [];
                            raw{i,j}(:,1)   = groupCell(:,2); % Position or X data
                            raw{i,j}(:,2)   = groupCell(:,1); % Time or Y data
                            
                            if and(~isempty(raw{i,j}), size(raw{i,j},1) > 2)
                                interptime  = [raw{i,j}(1,2):1:raw{i,j}(end,2)]'; % Interpolate values of y
                                interpdist  = interp1(raw{i,j}(:,2),raw{i,j}(:,1),interptime,'linear','extrap');% Interpolate values of x
                            else
                                fprintf(['Less than 2 coordinates after averaging duplicate', ...
                                 'time values.\n\tIgnoring vesicle %d in neuron %d. \n'], ...
                                 j-1, i);
                
                                idx = idx +1;
                                continue
                            end
                        end

                    else
                        interptime = raw{i,j}(:,2)'; %interptime = interptime(isnan(interptime)==0);
                        interpdist = raw{i,j}(:,1)'; %interpdist = interpdist(isnan(interpdist)==0);
                    end
                catch err2
                     msgbox(sprintf('Error in Row %s Col %s',num2str(i),num2str(j)));
                     rethrow(err2)
                end

                    % Store Raw Initial values X and Y
                    initial_x = interpdist(1,1);
                    initial_y = interptime(1,1);
                        
                    % Calculate dx and dy parameters from datafile
                    dx= [diff(interpdist);0];
                    dy= [diff(interptime);0];

                    % Convert units
                    dx = dx .* pixel; % Xpixels to micrometers
                    dy = dy .* tframe; % Ypixels to seconds

                    % Cumulative distance
                    cumdist= sum(abs(dx));
                    cumtime= sum(abs(dy));
                    netdist= sum(dx);

                    % Instantaneous velocity and acceleration
                    instvel= dx./dy;
                    instaccl= diff(instvel)./(1*tframe);

                    % Indices of paused/ant/retrograde velocities
                    velpos= find(instvel > pause);%Finds indices of anterograde velocities
                    velneg= find(instvel < pause*-1);%Finds indices of anterograde velocities
                    velzero= find((instvel<=pause) & (instvel>=pause*-1));%Finds indices of pauses

                    % Time of paused/ ant/ retrograde
                    tpause=(length(velzero)*tframe); %total time paused in seconds
                    tantero=(length(velpos)*tframe); %total time moving antero in seconds
                    tretro=(length(velneg)*tframe); %total time moving retro in seconds

                    % Number of direction changes prelim_code: This finds
                    % instances where instantaneous velocity changes,
                    % whether from positive to zero or negative or vice
                    % versa. It then find the indices for later analysis of
                    % directional changes.
                    instvel_sign=zeros(size(instvel));%Paused motiltiy=0
                    instvel_sign(velpos)=1;%Anterograde motility=1
                    instvel_sign(velneg)=-1;%Retrograde motility=-1
                    inst_velA=instvel_sign(1:(end-1));
                    inst_velB=instvel_sign(2:end);
                    dirchanges=(inst_velA~=inst_velB);%directional changes and pauses=1
                    dirchanges_index=find(abs(dirchanges)==1);%Index of directional changes and pauses = 1
                    sz_change = 1; %Use this to change the size dirchanges_index from 1 = row to 2 = column
                    
                    true_dir_change = [];

                    % This creates a matrix 1: size(dirchanges_index) 
                    % called true_dir_change that goes through each dir
                    % change to see if the run after changing direction
                    % passes the direction change threshold
                    for dir_index=1:size(dirchanges_index,sz_change);
                        if dir_index < size(dirchanges_index,sz_change);
                            if sum(abs(dx(dirchanges_index(dir_index)+1:dirchanges_index(dir_index+1)))) >= dir_thresh; % Calculates how far the vesicle moved after changing direction
                                true_dir_change = [true_dir_change, dir_index];
                            end
                        elseif dir_index == size(dirchanges_index,sz_change) && isempty(dirchanges_index)==0;
                            if sum(abs(dx(dirchanges_index(dir_index)+1:end))) >= dir_thresh; % Calculates how far the vesicle moved after changing direction
                                true_dir_change = [true_dir_change, dir_index];
                            end
                        end
                    end
                    
                    % This reassigns the dirchanges_index based on the true
                    % direction changes
                    if isempty(true_dir_change)==0
                        dirchanges_index = dirchanges_index(true_dir_change);
                    else
                        dirchanges_index = [1];
                    end
                    
                    % This uses the indices generated above to calculate
                    % direction/ sign of each run, length of run, total
                    % number of direction changes, and account for pauses.
                    % A vesivle must change direction and move at least X
                    % um (as set by user) in order to be considered a true
                    % direction change.
                    % Initialize variables
                    dir.indices = {};
                    dir.sign= [];
                    dir.run_length= [];
                    dir.switches= [];
                    dir.pause= [];
                    
                    
                    if size(dirchanges_index,sz_change)==1;
                        if isempty(dirchanges_index)==0; % Calculates run length for arrays with only 1 direction change index
                            dir_idx = 1;
                            dir.sign(dir_idx)= 0; % Sign for the current direction change--Antero is positive, Retro is negative
                            dir.run_length(dir_idx) = sum(abs(dx(dirchanges_index(dir_idx):end)));
                        else % Calculates run length for arrays with no direction change indices
                            dir.sign(dir_idx)= sum(dx(1:end)); % Sign for the current direction change--Antero is positive, Retro is negative
                            dir.run_length(dir_idx) = sum(abs(dx(1:end)));
                        end
                    elseif size(dirchanges_index,sz_change)>1;
                        for dir_idx=1:size(dirchanges_index,sz_change)-1;
                            if sum(abs(dx(dirchanges_index(dir_idx)+1:dirchanges_index(dir_idx+1)))) > dir_thresh; % Calculates how far the vesicle moved after changing direction
                                dir.indices{dir_idx} = sprintf('%d:%d',dirchanges_index(dir_idx)+1,dirchanges_index(dir_idx+1)); % Indices for the values in the instvel variable
                                dir.sign(dir_idx) = instvel_sign(dirchanges_index(dir_idx)+1); % Sign for the current direction of an individual run--not for the entire run. Anterograde is positive, Retrograde is negative
                                dir.run_length(dir_idx) = sum(abs(dx(dirchanges_index(dir_idx)+1:dirchanges_index(dir_idx+1)))); % Absolute value of run length
                                dir.switches(dir_idx) = 1;
                                dir.pause(dir_idx) = 0;
                            else
                                dir.indices{dir_idx} = sprintf('%d:%d',dirchanges_index(dir_idx)+1,dirchanges_index(dir_idx+1)); % Indices for the values in the instvel variable
                                dir.sign(dir_idx) = instvel_sign(dirchanges_index(dir_idx)+1); % Sign for the current direction change--Anterograde is positive, Retrograde is negative
                                dir.run_length(dir_idx) = sum(abs(dx(dirchanges_index(dir_idx)+1:dirchanges_index(dir_idx+1)))); % Absolute value of run length
                                dir.switches(dir_idx) = 0;
                                dir.pause(dir_idx) = 1;
                            end
                        end
                    end
                    
                    % Classifies individual vesicles based on Net
                    % displacement >10um as Anterograde (1), < 10um Retrograde (-1),
                    % or -10um< X < 10um Non-motile/ Bidirectional (0)
                    if netdist >=10;
                        dir.netmovement = 1; % anterograde is positive displacement
                    elseif netdist <=-10
                        dir.netmovement = -1;% retrograde is negative displacement
                    else
                        dir.netmovement = 0;% Nonmotile/ bidirectional does not displace >=10um in from start position 
                    end


                    %Determine length of pauses--original-mod            
                    pause_length=[];
                    for k=1:numel(dirchanges_index)
                        if k==1 && (instvel(dirchanges_index(k))>(pause*-1)) && (instvel(dirchanges_index(k))<pause)
                            pause_length=[pause_length,dirchanges_index(k)];
                        elseif k>1 && k<numel(dirchanges_index) && (instvel(dirchanges_index(k))>(pause*-1)) && (instvel(dirchanges_index(k))<pause)
                            pl=(dirchanges_index(k)-dirchanges_index(k-1));
                            pause_length=[pause_length,pl];
                        elseif k==numel(dirchanges_index) && (instvel(end)>(pause*-1)) && (instvel(end)<pause)
                            pl=numel(instvel)-dirchanges_index(k);
                            pause_length=[pause_length,pl];
                        end
                    end
                        
                    % Format data and cell structure for output
                    instvel = instvel';
%                     instvel(end) = 0; 
%                     instaccl(end)= 0;
                    instaccl = instaccl';
                    dir_totalswitches = sum(dir.switches);
                    dir_antero_switches = (dir.netmovement >0).* dir_totalswitches; % Number of switches binned by anterograde displacement >10um
                    dir_retro_switches = (dir.netmovement <0).* dir_totalswitches; % Number of switches binned by retrograde displacement >10um
                    dir_pause_switches = (dir.netmovement==0).* dir_totalswitches; % Number of switches binned by no net movement >10 in either direction
                    dir_less_thresh = sum((dir.pause==1).* dir_totalswitches); % Number of switches that do not meet threshold distance traveled
                    tlength = [tantero, tretro, tpause];

                    % Output
                    output.fname{i,1} = {[fname] [pixel] [tframe] [pause] [dir_thresh] [datetime]};
                    output.initial_xy{i,idx}= [initial_x initial_y];
                    output.initial_y(i,idx) = initial_y;
                    output.interpdist{i,idx}= interpdist';
                    output.interptime{i,idx}= interptime';
                    output.cumdist(i,idx)   = cumdist;
                    output.netdist(i,idx)   = netdist;
                    output.cumtime(i,idx)   = cumtime;
                    output.instvel{i,idx}   = instvel;
                    output.instaccl{i,idx}  = instaccl;
                    output.tantero(i,idx)   = tantero;
                    output.tretro(i,idx)    = tretro;
                    output.tpause(i,idx)    = tpause;
                    output.dir_indices{i,idx} = dir.indices;
                    output.dir_sign{i,idx}  = dir.sign;
                    output.dir_run_num(i,idx) = size(dir.run_length,2);
                    output.dir_run_length{i,idx} = dir.run_length;
                    output.dir_totalswitches(i,idx) = dir_totalswitches;
                    output.dir_antero_switches(i,idx) = dir_antero_switches;
                    % Above calculates the sum of switches binned by net
                    % displacement anterograde >10um. This is repeated
                    % for retrograde/ non-motile below.
                    output.dir_retro_switches(i,idx) = dir_retro_switches;
                    output.dir_pause_switches(i,idx) = dir_pause_switches;
                    output.dir_less_thresh(i,idx) = dir_less_thresh;
                    output.dir_pause{i,idx} = dir.pause;
                    output.dir_netmovement(i,idx) = dir.netmovement;
                    output.pause_length{i,idx}= pause_length;
                    output.compiled{i,idx}= {[fname] [instvel] [cumdist] [cumtime] ...
                        [tlength] [dir_totalswitches] [pause_length]};

            else
                    fprintf(['Less than 2 coordinates found in tracer file ', ...
                             'Row %d, Col %d.\n\tIgnoring vesicle %d in neuron %d. \n'], ...
                             i,j, j-1, i);
                
                    idx = idx+1;
            end % End 'If' statement to only analyze non-empty columns           
        end % End 'For' loop over columns        
    end % End 'For' loop over rows
        
    % Re-organize variables of cell arrays for Instvel and Run Length
    % Instantaneous velocity and Run Length
    [Y, Z] = max(cellfun(@max,(cellfun(@size,output.instvel,'UniformOutput',false))),[],2); % Gets the values (Y) and indices (Z) for the max values of the rows
    output.antero_run_length = {};
    output.retro_run_length = {};
    output.net_velocity = output.netdist./output.cumtime;
    
    for row = 1:size(output.instvel,1); % Calculates the size of the compiled_run_length matrix requried
        size_max(row) = size(cat(2,output.dir_run_length{row,1:end}),2);
    end
    clear row
    size_max = max(size_max);
    output.compiled_run_length = nan(size(output.instvel,1),size_max); % Creates the compiled_run_length matrix filled with Nan values--actual values will replace Nan later
    
    for row = 1:size(output.instvel,1);
        for col = 1:size(output.instvel,2);
            length_instvel = size(output.instvel{row,col},1);
            if length_instvel < Y(row)
                output.instvel{row,col}     = [output.instvel{row,col}(1:end)' ;NaN(Y(row)-length_instvel,1)];
                output.interpdist{row,col}  = [output.interpdist{row,col}(1:end)' ;NaN(Y(row)-length_instvel,1)];
                output.interptime{row,col}  = [output.interptime{row,col}(1:end)' ;NaN(Y(row)-length_instvel,1)];
            end
        end % End 'for loop' over columns
        size1 = size(cat(1,output.instvel{row,1:end}),1);
        
        if isempty(cat(2,output.dir_sign{row,1:end}))==0 && isempty(cat(2,output.dir_run_length{row,1:end}))==0; % Only if variable exists
            try
                size2 = size(cat(2,output.dir_sign{row,1:end}) .* cat(2,output.dir_run_length{row,1:end}),2); 
                output.compiled_run_length(row,1:size2) = cat(2,output.dir_sign{row,1:end}) .* cat(2,output.dir_run_length{row,1:end});
            catch 
                try
                    size2 = size(cat(2,output.dir_run_length{row,1:end}),2);
                    output.compiled_run_length(row,1:size2) = cat(2,output.dir_run_length{row,1:end});
                catch err2
                    rethrow(err2)
                end
            end
            output.antero_run_length(row,1) = {output.compiled_run_length(row,output.compiled_run_length(row,:)>0)};
            output.retro_run_length(row,1) = {output.compiled_run_length(row,output.compiled_run_length(row,:)<0)};
        else
            size2 = size(cat(2,output.dir_run_length{row,1:end}),2);
            output.compiled_run_length(row,1:size2) = cat(2,output.dir_run_length{row,1:end});
            output.antero_run_length(row,1) = {output.compiled_run_length(row,output.compiled_run_length(row,:)>0)};
            output.retro_run_length(row,1) = {output.compiled_run_length(row,output.compiled_run_length(row,:)<0)};
        end
    end % End 'for loop' over rows
    
%     size_compiled_instvel_max = max(sum(cell2mat(cellfun(@(x) size(x,1),output.instvel,'UniformOutput',false)),2));
    output.compiled_instvel     = cell(size(output.instvel,1),1);
    output.compiled_interpdist  = cell(size(output.instvel,1),1);
    output.compiled_interptime  = cell(size(output.instvel,1),1);
    
    for row = 1:size(output.instvel,1);
        % Compiled instvel
        tempInstvelNan  = cellfun(@(x) find(isnan(x),true,'first'),output.instvel(row,1:end),'unif',false);
        tempInstvel     = cellfun(@(x,y) x(1:y+1),output.interptime(row,1:end), tempInstvelNan,'unif',false);
        output.compiled_instvel{row,1}      = cat(1,tempInstvel{1,1:end});
        
        % Compiled interpdist
        tempInterpNan   = cellfun(@(x) find(isnan(x),true,'first'),output.interpdist(row,1:end),'unif',false);
        tempInterpDist  = cellfun(@(x,y) x(1:y+1),output.interpdist(row,1:end), tempInterpNan,'unif',false);
        output.compiled_interpdist{row,1}   = cat(1,tempInterpDist{1,1:end});
        
        % Compiled interptime
        tempInterpNan   = cellfun(@(x) find(isnan(x),true,'first'),output.interptime(row,1:end),'unif',false);
        tempInterpTime  = cellfun(@(x,y) x(1:y+1),output.interptime(row,1:end), tempInterpNan,'unif',false);
        output.compiled_interptime{row,1}   = cat(1,tempInterpTime{1,1:end});
    end
    % Concatenate compiled_interptime and compiled_interpdist
    output.compiled_instvel     = catData(output.compiled_instvel, 'Method', 'Row');
    output.compiled_interptime  = catData(output.compiled_interptime, 'Method', 'Row');
    output.compiled_interpdist  = catData(output.compiled_interpdist, 'Method','Row');
    clear size1 size2
    
    % Compile and save data in output.summary to user input filename
    a =repmat('Neuron',size(output.fname,1),1); output.neuron_label = cellstr(horzcat(a,num2str((1:1:size(a,1))')));
    output.summary{1,1} = output.fname{start,1};
    output.summary{1,2} = [nanmean(output.compiled_instvel,2) nanstd(output.compiled_instvel,[],2)...
        sum(isnan(output.compiled_instvel)==0,2)];
    output.summary{1,3} = [nanmean(output.cumdist,2) nanstd(output.cumdist,[],2)...
        sum(isnan(output.cumdist)==0,2)]; % mean, SD, and N over cumdist
    output.summary{1,4} = [nanmean(output.cumtime,2) nanstd(output.cumtime,[],2)...
        sum(isnan(output.cumtime)==0,2)]; % mean, SD, and N over cumtime
    output.summary{1,5} = [nanmean(output.tantero,2) nanstd(output.tantero,[],2)...
        sum(isnan(output.tantero)==0,2)]; % mean, SD, and N over tantero
    output.summary{1,6} = [nanmean(output.tpause,2) nanstd(output.tpause,[],2)...
        sum(isnan(output.tpause)==0,2)]; % mean, SD, and N over tpause
    output.summary{1,7} = [nanmean(output.tretro,2) nanstd(output.tretro,[],2)...
        sum(isnan(output.tretro)==0,2)]; % mean, SD, and N over tretro

    total_vesicles = sum(isnan(output.netdist)==0,2);
    output.summary{1,8} = [sum(output.netdist>=10,2) total_vesicles ...
        sum(output.netdist>=10,2)./total_vesicles]; % Anterograde vesicles
    output.summary{1,9} = [sum(output.netdist>-10 & output.netdist<10,2) total_vesicles ...
        sum(output.netdist>-10 & output.netdist<10,2)./total_vesicles]; % Paused vesicles
    output.summary{1,10} = [sum(output.netdist<=-10,2) total_vesicles ...
        sum(output.netdist<=-10,2)./total_vesicles]; % Retrograde vesicles
    output.summary{1,11} = total_vesicles;
    output.summary{1,15} = [nanmean(output.dir_totalswitches,2) nanstd(output.dir_totalswitches,[],2)...
        sum(isnan(output.dir_totalswitches)==0,2)]; % mean, SD, and N over total dirchanges
    output.summary{1,16} = [nanmean(output.dir_run_num,2) nanstd(output.dir_run_num,[],2)...
        sum(isnan(output.dir_run_num)==0,2)]; % mean, SD, and N over total run number 
    output.summary{1,17} = [nanmean(abs(output.compiled_run_length),2) nanstd(abs(output.compiled_run_length),[],2)...
        sum(isnan(abs(output.compiled_run_length))==0,2)]; % mean, SD, and N over total run lengths (positive and negative)
    output.summary{1,18} = [nanmean(output.dir_antero_switches,2) nanstd(output.dir_antero_switches,[],2)...
        sum(output.dir_netmovement >0,2)]; % mean, SD, and N over anter dirchanges
    output.summary{1,19} = [nanmean(output.dir_pause_switches,2) nanstd(output.dir_pause_switches,[],2)...
        sum(output.dir_netmovement==0,2)]; % mean, SD, and N over nonmotile/bidirectional dirchanges
    output.summary{1,20} = [nanmean(output.dir_retro_switches,2) nanstd(output.dir_retro_switches,[],2)...
        sum(output.dir_netmovement <0,2)]; % mean, SD, and N over retro dirchanges
    output.summary{1,21} = [nanmean(output.dir_less_thresh,2) nanstd(output.dir_less_thresh,[],2)...
        sum(output.dir_less_thresh>0,2)]; % mean, SD, and N over nonmotile/bidirectional dirchanges
    output.summary{1,22} = [cell2mat(cellfun(@nanmean,output.antero_run_length,'UniformOutput', false)) cell2mat(cellfun(@nanstd,output.antero_run_length,'UniformOutput', false))...
        cell2mat(cellfun(@length,output.antero_run_length,'UniformOutput', false))]; % mean, SD, and N over antero run lengths
    output.summary{1,23} = [cell2mat(cellfun(@nanmean,output.retro_run_length,'UniformOutput', false)) cell2mat(cellfun(@nanstd,output.retro_run_length,'UniformOutput', false))...
        cell2mat(cellfun(@length,output.retro_run_length,'UniformOutput', false))]; % mean, SD, and N over retro run lengths
    output.summary{1,24} = [nanmean(output.net_velocity,2) nanstd(output.net_velocity,[],2)...
        sum(isnan(output.net_velocity)==0,2)]; % mean, SD, and N over anter dirchanges
    
    % Tempsave data
    cd(SaveDir);
    analyzed_data = output;
    fname = [fname(1:end-4)];
    save([fname '_kymo-analysis.mat'],'-mat','analyzed_data'); 
    save([fname '_kymo-analysis.bkp'],'-mat','analyzed_data'); % Saves a backup file .bkp
    
    try
        tempSaveDir = userpath;
        tempSaveDir = tempSaveDir(1:end-1);
        if ~(exist([tempSaveDir filesep 'kymo_suite'],'dir')==7)
            mkdir(tempSaveDir, 'kymo_suite')
        end
        tempSaveDir = [tempSaveDir filesep 'kymo_suite'];
        save([tempSaveDir filesep fname '_kymo-analysis.mat'],'-mat','analyzed_data');
    catch
    end
        

    btn = questdlg('Would you like to graph this data in Matlab?'); % User dialog box to select graph/ not graph
    if strcmp(btn,'Yes') == 1    
        % Initialize variables for plotting graphs
        try
            param= inputdlg({'Graph Title','Color Map','Error sides (one vs two sided)','Width'},...
            'Input Parameters',1,{exp_filename,'jet','1','0.5'},options);
        catch
            param= inputdlg({'Graph Title','Color Map','Error sides (one vs two sided)','Width'},...
            'Input Parameters',1,{'Enter title here','jet','1','0.5'},options);
        end
        sup_title = param{1,1}; % Graph Title from ui dialog box
        bw_color_map = param{2,1}; % Color Map from uui dialog box
        error_sides = str2num(param{3,1}); 
        width= str2num(param{4,1}); 

        bw_title = {'Mean Inst Vel' 'Cum. Dist' 'Cum. Time' 'Time Antero' ...
            'Time Paused' 'Time Retro' 'Antero > 10um' 'Paused' 'Retro > 10um'};
        bw_xlabel = {''};
        bw_ylabel = {'Vel (um/s)' 'Mean Dist (um)' 'Mean Time (s)' 'Mean Time (s)' 'Mean Time (s)' ...
            'Mean Time(s)' 'Frac Events' 'Frac Events' 'Frac Events'};
        gridstatus = 'none';
        legend_type = 'plot';

        figure
        for i=2:10;
            subplot(3,3,i-1); % Create an M by N subplot

            if i<=7 ;
                barvalues = output.summary{1,i}(:,1); % Values to plot
                errors = output.summary{1,i}(:,2);% SD error bars
            else
                barvalues = output.summary{1,i}(:,3); % Values to plot
                errors = zeros(size(output.summary{1,i}(:,1))); %No error bars for Frac Events
            end
            bw_legend = [];
            
            try
                barweb(barvalues, errors, width, [], bw_title{i-1}, bw_xlabel, bw_ylabel{i-1}, bw_color_map, gridstatus, bw_legend, error_sides, legend_type);
            catch
                if all(barvalues==0)==0;
                    errors = zeros(size(output.summary{1,i}(:,1))); % No error bars if there is error in graphing error bars
                    barvalues(isnan(barvalues))=0; % Set barvalues = NaN to zero

%                     barweb(barvalues, errors, width, [], bw_title{i-1}, bw_xlabel, bw_ylabel{i-1}, bw_color_map, gridstatus, bw_legend, error_sides, legend_type);
                end

            end
        end
        [ax,h3]=suplabel(sup_title ,'t');set(h3,'FontSize',18);
        
        % Second Figure graphing direction change data
        bw_title2 = {'Avg Total Switches' 'Avg Runs' 'Avg Run Length' ...
            'Switches Antero > 10um' 'Non-Motile' 'Switches Retro > 10um'...
            'Antero Run Length' 'Retro Run Length' '' };
        bw_xlabel2 = {''};
        bw_ylabel2 = {'Switches/ Track' 'Avg Runs' 'Dist (um)' ' Switches/ Track' 'Switches/ Track' ...
            'Switches/ Track' 'Length (um)' 'Length (um)' '?????'};
        gridstatus = 'none';
        legend_type = 'plot';
        
        figure
        for i=15:20;
            subplot(3,3,i-14); % Create an M by N subplot
            try
                
                if all(output.summary{1,i}(:,1)==0);
                    barvalues = repmat(0.001,size(output.summary{1,i}(:,1)),1);
                    errors = zeros(size(output.summary{1,i}(:,1)));
                else
                    barvalues = output.summary{1,i}(:,1); % Values to plot
                    errors = output.summary{1,i}(:,2)./sqrt(output.summary{1,i}(:,3)); % SD / sqrt(N) = SEM error bars
                    bw_legend = [];
                    barweb(barvalues, errors, width, [], bw_title2{i-14}, bw_xlabel2, bw_ylabel2{i-14}, bw_color_map, gridstatus, bw_legend, error_sides, legend_type);
                end
            catch
                msgbox(sprintf('Kymo_analysis could not plot %s',bw_title2{i-14}));
            end
        end
        [ax,h3]=suplabel(sup_title ,'t');set(h3,'FontSize',18);
        
    end
    
    % Get filename and current directory
    SaveDir = uigetdir(SaveDir,'Select save directory.');
    if ~logical(SaveDir)
        SaveDir = pwd;
    end
    cd(SaveDir);
    analyzed_data = output;
    fname = [fname];
    save([fname '_kymo-analysis.mat'],'-mat','analyzed_data');


    btn2 = questdlg('Output all data to excel?');
    if strcmp(btn2,'Yes') == 1
        warning('off')
        % Initiate Excel ActX server % Format for using---xlswrite2007(file,data,sheet,range)
        Excel = actxserver('Excel.Application');
        File = [SaveDir filesep fname '_kymo-analysis.xlsx'];  %Make sure you put the whole path
        if ~exist(File,'file')
            ExcelWorkbook = Excel.workbooks.Add;
            ExcelWorkbook.SaveAs(File);
            ExcelWorkbook.Close(false);
        else
            msgbox(sprintf('This will overwrite %s to this directory: %s .',fname,[SaveDir filesep]));
        end
        ExcelWorkbook = Excel.workbooks.Open(File);

        % Correct for empty variables
        if isempty(output.antero_run_length)==1;
            output.antero_run_length = 0;
        else
            output.antero_run_length = output.compiled_run_length(output.compiled_run_length>0);
        end
        
        if isempty(output.retro_run_length)==1;
            output.retro_run_length = 0;
        else
            output.retro_run_length = output.compiled_run_length(output.compiled_run_length<0);
        end

        % Write info to Excel
        try    
            File=[SaveDir filesep fname '_kymo-analysis.xlsx'];

            xlswrite2007(File,{'Filename','um/ Pixel','Sec/ Frame','Pause = Inst. Vel < X (um/sec)' 'Direction Change Threshold' 'Time of Analysis'},'Summary','A1');
            xlswrite2007(File,output.summary{1,1},'Summary','A2'); % Filename and info
            xlswrite2007(File,{'Mean','SD','N'},'InstVel','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'InstVel','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','InstVel','A25'); % neuron label
                xlswrite2007(File,output.summary{1,2},'InstVel','B2'); % Instvel mean, SD, and N
                xlswrite2007(File,output.compiled_instvel','InstVel','A26'); % Raw InstVel values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Neuron'},'Interp_X_Pixels','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'Interp_X_Pixels','A2'); % neuron label
                xlswrite2007(File,{'Interp_X_Pixels'},'Interp_X_Pixels','A24'); % neuron label
                xlswrite2007(File,output.neuron_label','Interp_X_Pixels','A25'); % neuron label
                xlswrite2007(File,output.compiled_interpdist','Interp_X_Pixels','A26'); % Raw interpolated X pixel values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Neuron'},'Interp_Y_Pixels','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'Interp_Y_Pixels','A2'); % neuron label
                xlswrite2007(File,{'Interp_Y_Pixels'},'Interp_Y_Pixels','A24'); % neuron label
                xlswrite2007(File,output.neuron_label','Interp_Y_Pixels','A25'); % neuron label
                xlswrite2007(File,output.compiled_interptime','Interp_Y_Pixels','A26'); % Raw interpolated X pixel values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean','SD','N'},'NetVel','B1'); %xlswrite2007(File,interpdist','','A26'); % Raw interpolated pixel values for all vesicles per neuron (neurons separated by column) Titles/ info
                xlswrite2007(File,output.neuron_label,'NetVel','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','NetVel','A25'); % neuron label
                xlswrite2007(File,output.summary{1,24},'NetVel','B2'); % Net vel mean, SD, and N
                xlswrite2007(File,output.net_velocity','NetVel','A26'); % Raw InstVel values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean','SD','N'},'CumDist','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'CumDist','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','CumDist','A25'); % neuron label
                xlswrite2007(File,output.summary{1,3},'CumDist','B2'); % Cumdist mean, SD, and N
                xlswrite2007(File,output.cumdist','CumDist','A26'); % Raw Cumdist values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean','SD','N'},'CumTime','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'CumTime','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','CumTime','A25'); % neuron label
                xlswrite2007(File,output.summary{1,4},'CumTime','B2'); % Cumtime mean, SD, and N
                xlswrite2007(File,output.cumtime','CumTime','A26'); % Raw CumTime values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean Total Switches moving >= threshold','SD','N'},'TotalDirChanges','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'TotalDirChanges','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','TotalDirChanges','A25'); % neuron label
                xlswrite2007(File,output.summary{1,15},'TotalDirChanges','B2'); % Dirchanges mean, SD, and N
                xlswrite2007(File,output.dir_totalswitches','TotalDirChanges','A26'); % Raw Total DirChange values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean Switches binned by net Antero dislacement (>10um)','SD','N'},'AnteroDirChanges','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'AnteroDirChanges','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','AnteroDirChanges','A25'); % neuron label
                xlswrite2007(File,output.summary{1,18},'AnteroDirChanges','B2'); % Dirchanges mean, SD, and N
                xlswrite2007(File,output.dir_antero_switches','AnteroDirChanges','A26'); % Raw Antero DirChange values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean Switches binned by net dislacement <10um','SD','N'},'NonMotileDirChanges','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'NonMotileDirChanges','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','NonMotileDirChanges','A25'); % neuron label
                xlswrite2007(File,output.summary{1,19},'NonMotileDirChanges','B2'); % Dirchanges mean, SD, and N
                xlswrite2007(File,output.dir_pause_switches','NonMotileDirChanges','A26'); % Raw Non-Motile DirChange values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean Switches moving < threshold','SD','N'},'DirChanges<Thresh','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'DirChanges<Thresh','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','DirChanges<Thresh','A25'); % neuron label
                xlswrite2007(File,output.summary{1,21},'DirChanges<Thresh','B2'); % Dirchanges mean, SD, and N
                xlswrite2007(File,output.dir_less_thresh','DirChanges<Thresh','A26'); % Raw DirChange < thresh values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Mean Switches binned by net Retro dislacement (>10um)','SD','N'},'RetroDirChanges','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'RetroDirChanges','A2'); % neuron label
                xlswrite2007(File,output.neuron_label','RetroDirChanges','A25'); % neuron label
                xlswrite2007(File,output.summary{1,20},'RetroDirChanges','B2'); % Dirchanges mean, SD, and N
                xlswrite2007(File,output.dir_retro_switches','RetroDirChanges','A26'); % Raw Retro DirChange values for all vesicles per neuron (neurons separated by column)
            xlswrite2007(File,{'Anterograde Mean (s)','SD','N'},'Time','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'Time','A2'); % neuron label
                xlswrite2007(File,output.summary{1,5},'Time','B2'); % TAntero mean, SD, and N
            xlswrite2007(File,{'Paused Mean (s)','SD','N'},'Time','F1'); % Titles/ info
                xlswrite2007(File,output.summary{1,6},'Time','F2'); % TPause mean, SD, and N
            xlswrite2007(File,{'Retrograde Mean (s)','SD','N'},'Time','J1'); % Titles/ info
                xlswrite2007(File,output.summary{1,7},'Time','J2'); % TRetro mean, SD, and N
            xlswrite2007(File,{'N Antero' 'Total Vesicles' 'Frac-Antero'},'FracEvents','B1'); % Titles/ info
                xlswrite2007(File,output.neuron_label,'FracEvents','A2'); % neuron label
                xlswrite2007(File,output.summary{1,8},'FracEvents','B2'); % Fraction Antero (>10)
            xlswrite2007(File,{'N Non-motile' 'Total Vesicles' 'Frac-Non-motile'},'FracEvents','F1'); % Titles/ info
                xlswrite2007(File,output.summary{1,9},'FracEvents','F2'); % Fraction Non-motile
            xlswrite2007(File,{'N Retro' 'Total Vesicles' 'Frac-Retro'},'FracEvents','J1'); % Titles/ info
                xlswrite2007(File,output.summary{1,10},'FracEvents','J2'); % Fraction Retro (>10)
                xlswrite2007(File,output.neuron_label,'Antero Run Length','A2'); % neuron label
                xlswrite2007(File,output.summary{1,22},'Antero Run Length','B2'); % Avg Run Length all vesicles moving anterograde >= direction change threshold
                xlswrite2007(File,{'All Neurons'},'Antero Run Length','A25'); % neuron label
            if isempty(output.antero_run_length)==0;
                xlswrite2007(File,output.antero_run_length,'Antero Run Length','A26');
            end
                xlswrite2007(File,output.neuron_label,'Retro Run Length','A2'); % neuron label
                xlswrite2007(File,output.summary{1,23},'Retro Run Length','B2'); % Avg Run Length all vesicles moving retrograde >= direction change threshold
                xlswrite2007(File,{'All Neurons'},'Retro Run Length','A25'); % neuron label
            if isempty(output.retro_run_length)==0
                xlswrite2007(File,output.retro_run_length,'Retro Run Length','A26');
            end
                        
        catch err3
            cd(SaveDir);
            analyzed_data = output;
            save([fname '_kymo-analysis.bkp'],'-mat','analyzed_data');% Saves a backup file .bkp only if there is an error. Change extension .bkp to .mat to make a useable MAT file.
            msgbox(['Some variables do not have values and were not written to Excel in the directory' sprintf(' %s %s.',SaveDir,fname)], 'Error');
            fail=1;
%             rethrow(err3);
        end
            msgbox(sprintf('Saved %s to %s',fname,[SaveDir filesep]));
            fprintf(sprintf('Saved %s to \n',fname));
            disp([SaveDir filesep]);

        try
            Excel.ActiveWorkbook.Worksheets.Item('Sheet1').Delete;
            Excel.ActiveWorkbook.Worksheets.Item('Sheet2').Delete;
            Excel.ActiveWorkbook.Worksheets.Item('Sheet3').Delete;
        catch
        end
            
        % Close Excel ActX server
        ExcelWorkbook.Save;
        ExcelWorkbook.Close(false);  % Close Excel workbook.
        Excel.Quit;
        delete(Excel);
        warning('on')
    end
end
toc



