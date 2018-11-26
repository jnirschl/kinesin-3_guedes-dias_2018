function varargout = peakAnalysis(varargin)
% PEAKANALYSIS MATLAB code for peakAnalysis.fig
%      PEAKANALYSIS, by itself, creates a new PEAKANALYSIS or raises the existing
%      singleton*.
%
%      H = PEAKANALYSIS returns the handle to a new PEAKANALYSIS or the handle to
%      the existing singleton*.
%
%      PEAKANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAKANALYSIS.M with the given input arguments.
%
%      PEAKANALYSIS('Property','Value',...) creates a new PEAKANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before peakAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to peakAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help peakAnalysis

% Last Modified by GUIDE v2.5 06-Apr-2017 11:36:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @peakAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @peakAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before peakAnalysis is made visible.
function peakAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to peakAnalysis (see VARARGIN)

% Use OpenGL rendering
 switch lower(computer)
    case {'pcwin','pcwin64'}
        opengl software
     otherwise
 end

% Choose default command line output for peakAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes peakAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Initialize variables
handles.SUCCESS     = true;

% Subfunction to check for input errors and parse input
% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;

% Logical inputs
ip.addParameter('DIR','',@ischar);
ip.addParameter('EMAIL',false,@islogical );
ip.addParameter('EXCEL',false,@islogical ); % Export to Excel
ip.addParameter('PLOT_1',true,@islogical );    % Plot original + overlay
ip.addParameter('PLOT_2',true,@islogical );    % Plot signal density overlap
ip.addParameter('RAND_1',false,@islogical );    % Randomize both SIG_1 (Syn) - default = false
ip.addParameter('RAND_2',true,@islogical );    % Randomize both SIG_2 (PTM) - default = tru
ip.addParameter('NND_BIPLOT',true,@islogical );    % Plot Rand NND above original
ip.addParameter('PARALLEL',false,@islogical );    % Parallel processing
ip.addParameter('CALC_PEAKS',true,@islogical );% Calculate new peaks or use
                                               % existing calculations?

ip.addParameter('EB',false,@islogical ); 
ip.addParameter('EB_dynamics',false,@islogical ); 
ip.addParameter('PAUSE',false,@islogical ); 
ip.addParameter('PLOT_RASTER',false,@islogical );    % Plot original + overlay
ip.addParameter('NORM',false,@islogical ); 
ip.addParameter('PEAK_RAND',true,@islogical ); % Randomize peaks and
                                     % calculate NND?
ip.addParameter('GRAPH_COMPILED',true,@islogical );
ip.addParameter('GRAPH_OUTPUT',true,@islogical );
ip.addParameter('CALC_PTM_SYN_OVERLAP',true,@islogical );
ip.addParameter('CALC_SYN_WIDTH',true,@islogical );
ip.addParameter('CLUSTER_PEAKS',true,@islogical ); % Cluster peaks
ip.addParameter('SET_PEAK_N',false,@islogical );
ip.addParameter('USE_SMOOTH_1',false,@islogical ); % SIG_1 is PTM
ip.addParameter('USE_SMOOTH_2',true,@islogical );  % SIG_2 is SYN
ip.addParameter('CALC_IDX',logical([zeros(1,5),1,zeros(1,6)]),@islogical ); % mean only
% ip.addParameter('BOTH',false,@islogical ); % Randomize both SIG_1 
                                           % and SIG_2 or only SIG_2 %
                                           % delete 20180609
% Numeric inputs
ip.addParameter('px_nm',1,@(x) all([isscalar(x),x>0,isnumeric(x)]));
ip.addParameter('UM_PX',1,@(x) all([isscalar(x),x>0,isnumeric(x)]));
ip.addParameter('S_FRM',1,@(x) all([isscalar(x),x>0,isnumeric(x)]));

ip.addParameter('FILTER',[],@(x) or(isempty(x), all([isscalar(x),x>0,isnumeric(x)])));
ip.addParameter('WIN_SIZE',10, @(x) all([isscalar(x),x>0,isnumeric(x)]));
% ip.addParameter('SIGMA',3, @(x) all([isscalar(x),x>0,isnumeric(x)]));
ip.addParameter('MINPEAKDIST',0.1, @(x) all([isscalar(x),x>0,isnumeric(x)]));
ip.addParameter('KNN',2,@(x) all([isscalar(x),x>0,isnumeric(x)]));
ip.addParameter('NSAMPLE',1e3,@(x) all([isscalar(x),x>100,x<1e5,isnumeric(x)]) );
ip.addParameter('MIN_PEAK_WIDTH',0.05,@(x) all([isscalar(x),x>0,isnumeric(x)]));

% String inputs
ip.addParameter('SAVE_DIR',pwd,@ischar );
ip.addParameter('PEAK_ARGS',{'Annotate','extents', ...
                    'WidthReference','halfheight'},@iscell );
ip.addParameter('SMOOTH','win',@(x) any(strcmpi(x,{'spine','sp', ...
                    'win','window','','na'}) ));
ip.addParameter('EXT','.pdf',@(x) any(strcmpi(x,{'.pdf','.eps', ...
                    '.png','.tiff','.tif'})));
ip.addParameter('INPUT_EXT','.csv',@(x) any(strcmpi(x,{'.csv'})));
ip.addParameter('DISTANCE','sqeuclidean',@(x) any(strcmpi(x,{'sqeuclidean','cityblock','cosine','correlation','hamming'})));
ip.addParameter('EB_ANALYSIS','initiation',@(x) any(strcmpi(x,{'antero','retro','init','initiation','term','termination'})));

ip.addParameter('PEAK_METHOD','MinPeakHeight',@(x) any(strcmpi(x, ...
                    {'MinPeakHeight','Threshold'})));
ip.addParameter('SIG_NAMES',{'PTM','SYN'},@iscell);
ip.addParameter('EMAIL_OPT',{'RECIPIENT','jnirsch@gmail.com'},@iscell);
ip.parse(varargin{:});
p = ip.Results;
% p.MIN_PEAK_WIDTH = p.MIN_PEAK_WIDTH/ power(p.px_nm,-1);

% Update dir if none given
if isempty(p.DIR)
    p.DIR   = uigetdir('Please select the folder containing one subfolder per experimental condition');
end

% Check package dependencies
PKG_DEP     = { {'matlab', 'R2014b'}, {'images', '7.1'}, {'stats', '9.1'},...
              {'signal', '6.1'}};
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
REQ_FUN = {'export_fig','normData','dirList',...
          'slidingWindowFilter_GD','distinguishable_colors',...
          'catData','importCSV','mergeArray','nhist',...
          'sendMATLABMail','firstOrderStats',...
          'xlswrite2007','kCluster','groupUnique',...
          'export2excel','export2csv'}; % Check for required M-Files
REQ_FUN_PATH = cell(numel(REQ_FUN),1);
for a1 = 1:numel(REQ_FUN)
    if ~isempty(REQ_FUN{a1})
        REQ_FUN_PATH{a1} = which(REQ_FUN{a1});
        if isempty(REQ_FUN_PATH{a1})
            error([REQ_FUN_PATH{a1} '.m is a required function.']);
        end
    end
end

ALL_DIR_NAME = {'FILEPATH'};
ALL_DIR = {p.DIR};

for DIR_idx     = 1:numel(ALL_DIR);
    tempPath    = strfind(ALL_DIR{DIR_idx},filesep);
    if ~ischar(ALL_DIR{DIR_idx})
        error('%s must be a character array.',ALL_DIR_NAME{DIR_idx})
    end

    % Ensure that DIR is a directory with subfolders containing .csv files
    if strcmpi(ALL_DIR_NAME{DIR_idx},'FILEPATH')
        % Select new directory for DIR
        if ~(exist(ALL_DIR{DIR_idx},'dir')==7)
            p.DIR   = uigetdir('Please select the folder containing one subfolder per experimental condition');
            p.FILENAME      = ALL_DIR{DIR_idx}((tempPath(end)+1):end);
            p.FILENAME      = strrep(p.FILENAME,p.INPUT_EXT,'');
        else
            p.FILENAME      = ALL_DIR{DIR_idx}((tempPath(end)+1):end);
            p.FILENAME      = strrep(p.FILENAME,p.INPUT_EXT,'');
        end
        
        % Ensure two CSV files are in each subfolder of DIR
        filterSpec          = {'*.csv','CSV files(*.csv)'; '*.*',  'All Files (*.*)'};
        SubFolders          = dirList(ALL_DIR{DIR_idx});
        
        % Go up one directory if there are no subfolders
        if isempty(SubFolders)
            fSepIdx     = strfind(p.DIR,filesep);
            p.DIR     	= p.DIR(1:(fSepIdx(end)));
            ALL_DIR     = {p.DIR};
            SubFolders  = dirList(p.DIR);
            if isempty(SubFolders)
                error(['Could not find any subfolders in the selected directory!',...
                    'Please choose a folder that contains subfolders that you want to analyze.'])
            end
        end
        
        if ~p.EB_dynamics && ~p.PAUSE
            % FOR loop to check for .csv files in each subfolder
            for SubFoldIdx      = 1:numel(SubFolders)
                tempSubFold     = SubFolders(SubFoldIdx).name;

                SubFiles        = dir([ALL_DIR{DIR_idx} filesep tempSubFold filesep filterSpec{1}]);

                if numel(SubFiles)<2
                    errordlg('Please ensure that two ''.csv'' files are in each subfolder of DIR');
                    validateattributes(SubFiles,{'struct'},{'nonempty','numel',2},'peakAnalysis')
                end
            end
        end
        
    end
end     % END FOR loop for DIR_idx
    

% Allocate output variables
handles.p           = p;
handles             = update_params_Callback(hObject, [], handles);
handles.SUCCESS     = true;
handles.ERR.message = 'No error has occured during execution of peakAnalysis.m';
handles.ERR.stack   = [];
handles.EMAIL       = 'peakAnalysis has completed!';

% Get information on EB_ANALYSIS
if logical(get(handles.EB_ANALYSIS_init_term,'Value'));
    handles.p.EB_ANALYSIS   = 'init';%'antero';
else
    handles.p.EB_ANALYSIS   = 'term';%'retro';
end

if handles.p.px_nm==1
    handles.p.XLABEL        = 'Distance (pixels)';
else
    handles.p.XLABEL        = 'Distance (\mum)';
end

% Set CompileData button to off for now
set(handles.graph_compiled,'Enable','off') 

% Set warning to off
warning('off');

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = peakAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Update all parameters in handles.p
function handles = update_params_Callback(hObject, eventdata, handles)
% Update all params

try
    handles.p.EXCEL         = logical(get(handles.export_excel_box,'Value'));

    handles.p.PLOT_1        = logical(get(handles.plot_1,'Value'));
%     handles.p.PLOT_2        = logical(get(handles.plot_2,'Value'));
    handles.p.PEAK_RAND     = logical(get(handles.peak_rand,'Value'));

    handles.p.USE_SMOOTH_1  = logical(get(handles.use_smooth_1,'Value'));
    handles.p.USE_SMOOTH_2  = logical(get(handles.use_smooth_2,'Value'));

    menu_contents           = cellstr(get(handles.smooth_menu,'String'));
    handles.p.SMOOTH        = menu_contents{get(handles.smooth_menu,'Value')};

    handles.p.WIN_SIZE      = str2double(get(handles.win_size_edit,'String'));
    handles.p.SIGMA         = (handles.p.WIN_SIZE/3);

    menu_contents           = cellstr(get(handles.peak_method_menu,'String'));
    handles.p.PEAK_METHOD   = menu_contents{get(handles.peak_method_menu,'Value')};

    SIG_1                   = get(handles.sig_names_1,'String');
    SIG_2                   = get(handles.sig_names_2,'String');
    handles.p.SIG_NAMES     = {SIG_1,SIG_2};

    menu_contents           = cellstr(get(handles.output_format_menu,'String'));
    handles.p.EXT           = menu_contents{get(handles.output_format_menu,'Value')};

    handles.p.NSAMPLES      = str2double(get(handles.nsample,'String'));
    handles.p.KNN           = str2double(get(handles.knn_edit,'String'));
    handles.p.px_nm         = (str2double(get(handles.nm_px,'String')))^-1;
    handles.p.px_um         = (handles.p.px_nm)*1e3;
catch ERR
    disp(ERR.message)
    handles.ERR             = ERR;
    error_Callback(hObject,[],handles);
end

% Update handles structure
guidata(hObject,handles);


% --- Executes on error during processing
function error_Callback(hObject, eventdata, handles)
% Executes on error during processing

if or(handles.SUCCESS,~isempty(handles.ERR.stack))
    if handles.p.EMAIL
        fprintf('Sending email...')
        sendMATLABMail(handles.p.EMAIL_OPT{:});
    end
else
    if handles.p.EMAIL
        sendMATLABMail(handles.p.EMAIL_OPT{:},'SUBJ','ERROR', ...
            'MSG',handles.ERR.message);
    end
end

% Update handles structure
guidata(hObject,handles);


% --- Executes on button press in load_data.
function load_data_Callback(hObject, ~, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update all params
handles = update_params_Callback(hObject, [], handles);

% --- Subfunction to load data from CSV
% Initialize variables
DIR         = handles.p.DIR;
p           = handles.p;
FOLDERS     = dirList(DIR);

% Go up one directory if there are no subfolders
if isempty(FOLDERS)
    fSepIdx = strfind(DIR,filesep);
    DIR     = DIR(1:(fSepIdx(end)));
    handles.p.DIR = DIR;
    FOLDERS     = dirList(DIR);
    if isempty(FOLDERS)
        error(['Could not find any subfolders in the selected directory!',...
            'Please choose a folder that contains subfolders that you want to analyze.'])
    else
        tempPATH                = [DIR filesep FOLDERS(1).name];
        FILES                   = dir([tempPATH filesep '*' p.INPUT_EXT]);
        if isempty(FILES)
            error('No %s files found in:\n%s',p.INPUT_EXT, tempPATH)
        else
            validateattributes(FILES,{'struct'},{'nonempty','numel',2},'Load data')
        end
    end
end

SUCCESS     = true;

% Update signal names
SIG_1       = get(handles.sig_names_1,'String');
SIG_2       = get(handles.sig_names_2,'String');
p.SIG_NAMES = {SIG_1,SIG_2};

% Allocate variables
SAVE_DIR    = cell(numel(FOLDERS),1);

SIGNAL_1_all= cell(numel(FOLDERS),1); % SIGNAL_1 is the PTM
normSIGNAL_1_all= cell(numel(FOLDERS),1); % SIGNAL_1 is the PTM

SIGNAL_2_all= cell(numel(FOLDERS),1); % SIGNAL_2 is synapsin
normSIGNAL_2_all= cell(numel(FOLDERS),1); % SIGNAL_2 is synapsin

% Load EB dynamics data
if p.EB_dynamics
    if ispc
        [~,FILEPATH]= dos('dir *.xls* /b/s');
        
    elseif isunix
        [~,FILEPATH]= unix('find -name ''*.xls*''');
    end
    FILEPATH    = strsplit(FILEPATH,['.' filesep])';
    EMPTY       = cellfun(@isempty,FILEPATH,'unif',true);
    FILEPATH    = FILEPATH(~EMPTY);
    FILEPATH    = cellfun(@(x) fullfile(DIR, x),FILEPATH,'unif',false);
    FILEPATH    = cellfun(@strtrim ,FILEPATH,'unif',false);
    tempSAVE_DIR        = cellfun(@(x) strsplit(x,filesep),FILEPATH,'unif',false);
    
    % Print updates to command window
    fprintf('Loading EB dynamics data ...\n')
    fprintf('\tSignal 1 (stationary) are synapses.\n\tSignal 2 (randomization) are dynamic events.\n\n')

    % Load data
    [LineScan, Init, Term, Info] = load_EB_Dynamics(FILEPATH,'PLOT', p.PLOT_RASTER);

	% Update figure
    set(handles.nm_px,'string',num2str(Info(1).um_px*1e3))
    set(handles.sig_names_1,'string','SYN')
    set(handles.sig_names_2,'string','EB')
    
    % Output
    handles.p.px_um             = power(Info(1).um_px,-1);
    handles.FOLDERS             = FOLDERS;
    handles.SAVE_DIR            = cellfun(@(x) fullfile(x{1:end-1}),tempSAVE_DIR,'unif',false);
    handles.SIGNAL_1_all        = num2cell(catData({LineScan.Int}','method','row'),2)';
        handles.SIGNAL_1_all    = cellfun(@(x) x(:), handles.SIGNAL_1_all,'unif',false);
    handles.normSIGNAL_1_all    = num2cell(catData({LineScan.IntSmooth}','method','row'),2)';
        handles.normSIGNAL_1_all= cellfun(@(x) x(:), handles.normSIGNAL_1_all,'unif',false);
    handles.SIGNAL_2_all        = SIGNAL_2_all;
    handles.normSIGNAL_2_all    = normSIGNAL_2_all;
    handles.EB_Init             = Init;
    handles.EB_Term             = Term;
    handles.EB_Info             = Info;
    handles.EB_line_scan_x      = {LineScan.x}';
elseif p.PAUSE
    if ispc
        [~,FILEPATH]= dos('dir *.xls* /b/s');

    elseif isunix
        [~,FILEPATH]= unix('find -name ''*.xls*''');
    end
    FILEPATH    = strsplit(FILEPATH,['.' filesep])';
    EMPTY       = cellfun(@isempty,FILEPATH,'unif',true);
    FILEPATH    = FILEPATH(~EMPTY);
    FILEPATH    = cellfun(@(x) fullfile(DIR, x),FILEPATH,'unif',false);
    FILEPATH    = cellfun(@strtrim ,FILEPATH,'unif',false);
    tempSAVE_DIR        = cellfun(@(x) strsplit(x,filesep),FILEPATH,'unif',false);

    % Print updates to command window
    fprintf('Loading EB dynamics data ...\n')
    fprintf('\tSignal 1 (stationary) are synapses.\n\tSignal 2 (randomization) are vesicle pauses.\n\n')

    % Load data
    [LineScan, Init, Term, Info] = load_Pause(FILEPATH,'PLOT', p.PLOT_RASTER,...
                                  'UM_PX', p.UM_PX,'S_FRM',p.S_FRM,...
                                  'FILTER',p.FILTER);

    % Update figure
    set(handles.nm_px,'string',num2str(Info(1).um_px*1e3))
    set(handles.sig_names_1,'string','SYN')
    set(handles.sig_names_2,'string','Pause')

    % Output
    handles.p.px_um             = power(Info(1).um_px,-1);
    handles.FOLDERS             = FOLDERS;
    handles.SAVE_DIR            = cellfun(@(x) fullfile(x{1:end-1}),tempSAVE_DIR,'unif',false);
    handles.SIGNAL_1_all        = num2cell(catData({LineScan.Int}','method','row'),2)';
        handles.SIGNAL_1_all    = cellfun(@(x) x(:), handles.SIGNAL_1_all,'unif',false);
    handles.normSIGNAL_1_all    = num2cell(catData({LineScan.IntSmooth}','method','row'),2)';
        handles.normSIGNAL_1_all= cellfun(@(x) x(:), handles.normSIGNAL_1_all,'unif',false);
    handles.SIGNAL_2_all        = SIGNAL_2_all;
    handles.normSIGNAL_2_all    = normSIGNAL_2_all;
    handles.EB_Init             = Init;
    handles.EB_Term             = Term;
    handles.EB_Info             = Info;
    handles.EB_line_scan_x      = {LineScan.x}';
else
    try
        for idxFolder   = 1:numel(FOLDERS)
            tempPATH                = [DIR filesep FOLDERS(idxFolder).name];
            SAVE_DIR{idxFolder}     = [tempPATH filesep 'OUTPUT'];

            % Create SAVE_DIR if it does not exist
            if ~(exist(SAVE_DIR{idxFolder},'dir')==7)
                mkdir(tempPATH,'OUTPUT');
            end

            FILES                   = dir([tempPATH filesep '*' p.INPUT_EXT]);
            if isempty(FILES)
                error('No %s files found in:\n%s',p.INPUT_EXT, tempPATH)
            else
                validateattributes(FILES,{'struct'},{'nonempty','numel',2},'Load data')
            end

            SIG_NAME_COMP           = cellfun(@(x) [x],p.SIG_NAMES,'unif',false);
            SIGNAL_IDX              = ~(cellfun('isempty', ...
                                               regexpi({FILES.name},SIG_NAME_COMP{1})));
            if ~xor(SIGNAL_IDX(1),SIGNAL_IDX(2))
                SIGNAL_IDX          = (cellfun('isempty', ...
                                               regexpi({FILES.name},SIG_NAME_COMP(2))));
            end

            if ~xor(SIGNAL_IDX(1),SIGNAL_IDX(2))
                TITLE               = 'Choose the file to use as Signal 1';
                BTN_CHOICES         = {FILES.name}; BTN_CHOICES = cellfun(@(x) strrep(x,'.csv',''),BTN_CHOICES,'unif',false);
                BTN                 = questdlg('Signal 1',TITLE, BTN_CHOICES{:},BTN_CHOICES{1});
                if ~isempty(BTN)
                    SIGNAL_IDX              = ~(cellfun('isempty', ...
                                               regexpi({FILES.name},BTN)));
                else
                    errordlg(fprintf('Dialog box closed, defaulting to use %s as SIGNAL 1.\n',BTN_CHOICES{1}))
                    fprintf('Dialog box closed, defaulting to use %s as SIGNAL 1.\n',BTN_CHOICES{1})
                    SUCCESS             = false;
                    handles.ERR         = ERR;
                    SIGNAL_IDX(1)       = true;
                    SIGNAL_IDX(2)       = false;
                end
            end

            SIGNAL_1_all{idxFolder} = num2cell(cell2mat(importCSV([tempPATH filesep...
                                FILES(SIGNAL_IDX).name])),1);
            SIGNAL_2_all{idxFolder} = num2cell(cell2mat(importCSV([tempPATH filesep...
                                FILES(~SIGNAL_IDX).name])),1);

            if or(p.EB,numel(SIGNAL_1_all{idxFolder}{1}) ~= numel(SIGNAL_2_all{idxFolder}{1}))
                fprintf(['Assuming EB comet analysis for folder ''%s''. ',...
                        'Generating ''line scan'' for EB based on pixel coordinates...\n'],...
                         FOLDERS(idxFolder).name)
                 fprintf('\tSignal 1 (stationary) is ''%s''.\n\tSignal 2 (randomization) is ''%s''\n',FILES(SIGNAL_IDX).name, FILES(~SIGNAL_IDX).name)
                WITHIN_BOUNDS   = @(x,BOUNDS) (x >= min(BOUNDS) & x <= max(BOUNDS));
                if numel(SIGNAL_1_all{idxFolder}{1}) < numel(SIGNAL_2_all{idxFolder}{1})
                    temp_X_val      = (1:1:numel(SIGNAL_2_all{idxFolder}{1,1}))';
                    X_val_all       = cellfun(@(x) colon(1,1,find(~isnan(x),1,'last'))'./ p.px_um ,SIGNAL_2_all{idxFolder},'unif',false);
                    emptyVec        = cellfun(@(x) zeros(numel(x),1),X_val_all,'unif',false);

                    BOUNDS          = cellfun(@(x) num2cell([x-0.5,x+0.5],2),SIGNAL_1_all{idxFolder},'unif',false); % Find within 1 pixel
                    IDX             = cellfun(@(x) (cellfun(@(y,z) find(WITHIN_BOUNDS(temp_X_val,y),1,'first'),x,'unif',false)),...
                                         BOUNDS,'unif',false);
                    IDX             = cellfun(@(x) cell2mat(x), IDX,'unif',false);

                    SIGNAL_1_all{idxFolder}     = cellfun(@(x,y) assign2one(x,y),emptyVec,IDX,'unif',false);
                    SIGNAL_1_all{idxFolder}     = num2cell(catData(SIGNAL_1_all{idxFolder}),1);
                else
                    temp_X_val      = (1:1:numel(SIGNAL_1_all{idxFolder}{1,1}))';
                    X_val_all       = cellfun(@(x) colon(1,1,find(~isnan(x),1,'last'))'./ p.px_um ,SIGNAL_1_all{idxFolder},'unif',false);
                    emptyVec        = cellfun(@(x) zeros(numel(x),1),X_val_all,'unif',false);

                    BOUNDS          = cellfun(@(x) num2cell([x-0.5,x+0.5],2),SIGNAL_2_all{idxFolder},'unif',false); % Find within 1 pixel
                    IDX             = cellfun(@(x) (cellfun(@(y,z) find(WITHIN_BOUNDS(temp_X_val,y),1,'first'),x,'unif',false)),...
                                         BOUNDS,'unif',false);
                    IDX             = cellfun(@(x) cell2mat(x), IDX,'unif',false);

                    SIGNAL_2_all{idxFolder}     = cellfun(@(x,y) assign2one(x,y),emptyVec,IDX,'unif',false);
                    SIGNAL_2_all{idxFolder}     = num2cell(catData(SIGNAL_2_all{idxFolder}),1);
                end
                fprintf('\tWriting to CSV file in ''OUTPUT'' subfolder...\n');
                csvwrite([SAVE_DIR{idxFolder} filesep strrep(FILES(SIGNAL_IDX).name, '.csv','_SIG_1_pseudoLineScan.csv')], SIGNAL_1_all{idxFolder})
                csvwrite([SAVE_DIR{idxFolder} filesep strrep(FILES(~SIGNAL_IDX).name,'.csv','_SIG_2_pseudoLineScan.csv')], SIGNAL_2_all{idxFolder})
            else
                fprintf('\tSignal 1 (stationary) is ''%s''.\n\tSignal 2 (randomization) is ''%s''\n',FILES(SIGNAL_IDX).name, FILES(~SIGNAL_IDX).name)
            end

            if p.NORM
                normSIGNAL_1_all{idxFolder}= cellfun(@(x) normData(x), ...
                                                     SIGNAL_1_all{idxFolder},'unif',false);
                normSIGNAL_2_all{idxFolder}= cellfun(@(x) normData(x), ...
                                                    SIGNAL_2_all{idxFolder},'unif',false);
            else
                normSIGNAL_1_all{idxFolder}= cellfun(@(x) (x), ...
                                                     SIGNAL_1_all{idxFolder},'unif',false);
                normSIGNAL_2_all{idxFolder}= cellfun(@(x) (x), ...
                                                    SIGNAL_2_all{idxFolder},'unif',false);
            end


            % Error check
            % Check that the line scans per each csv is the same for
            % SIGNAL_1(PTM) and SIGNAL_2(Synapsin)
            validateattributes(normSIGNAL_2_all{idxFolder},{'cell'}, ...
                               {'size',size(normSIGNAL_1_all{idxFolder})});

            % Check that the number of elements for each column is identical
            validateattributes(normSIGNAL_2_all{idxFolder}{1,1}, ...
                               {'double'},{'size',size(normSIGNAL_1_all{idxFolder}{1,1})});
            clear SIGNAL_IDX
        end
    catch ERR
        SUCCESS                 = false;
        handles.ERR             = ERR;
        handles.SUCCESS         = false;
        error_Callback(hObject,[],handles)
    end
    
    % Output
    handles.FOLDERS             = FOLDERS;
    handles.SAVE_DIR            = SAVE_DIR;
    handles.SIGNAL_1_all        = SIGNAL_1_all;
    handles.normSIGNAL_1_all    = normSIGNAL_1_all;
    handles.SIGNAL_2_all        = SIGNAL_2_all;
    handles.normSIGNAL_2_all    = normSIGNAL_2_all;
    handles.p                   = p;
end

% Check success
if SUCCESS
    fprintf('Data loaded successfully.\n');
else
    fprintf('Error loading data!\n');
    disp(ERR.message)
    error_Callback(hObject,[],handles);
end

% Update handles structure
guidata(hObject, handles);


function x = assign2one(x,y)
% Subfunction to assign x to one at the indices of y
x(y)    = 1;


% --- Executes on button press in process_data.
function process_data_Callback(hObject, eventdata, handles)
% hObject    handle to process_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update all params
handles = update_params_Callback(hObject, [], handles);

if ~isfield(handles,'FOLDERS')
    errordlg('Please load data first!')
    return
end

% --- Subfunction to process data
% Initialize variables
DIR             = handles.p.DIR;
FOLDERS         = handles.FOLDERS;
SAVE_DIR        = handles.SAVE_DIR;
SIGNAL_1_all    = handles.SIGNAL_1_all;
normSIGNAL_1_all= handles.normSIGNAL_1_all;
SIGNAL_2_all    = handles.SIGNAL_2_all;
normSIGNAL_2_all= handles.normSIGNAL_2_all;
p               = handles.p;
SUCCESS         = handles.SUCCESS;
ERR             = [];
p.EXCEL         = false;
p.EB_ANALYSIS   = handles.p.EB_ANALYSIS;
warning off

try    
    for idxFolder   = 1:numel(FOLDERS) % FOR loop over all folders
        
        % Search for compiled data
        if exist([SAVE_DIR{idxFolder} filesep 'compiled_Data.mat'], ...
                 'file')==2
            fprintf('Loading pre-calculated compiled_Data file...\n')
            load([SAVE_DIR{idxFolder} filesep 'compiled_Data.mat']);
        else
            fprintf('Processing folder %s...\n',FOLDERS(idxFolder).name)
        end

        % Check whether peaks have already been detected
        if p.EB_dynamics || p.PAUSE
            p.CALC_PEAKS    = false;
            p.PLOT_1        = false;
            p.CALC_SYN_WIDTH= false;
            p.CALC_PTM_SYN_OVERLAP = false;
            normSIGNAL_1    = normSIGNAL_1_all;
            normSIGNAL_2    = normSIGNAL_1_all;
            X_val_all       = cellfun(@(x) colon(1,1,find(~isnan(x),1,'last'))'./ p.px_um ,normSIGNAL_1,'unif',false);
            X_val_all       = catData(X_val_all); X_val_all = num2cell(X_val_all,1);
            tempSAVE_DIR    = [filesep SAVE_DIR{idxFolder}];
        else
            % For each subfolder there should be a cell array of column vectors
            normSIGNAL_1    = normSIGNAL_1_all{idxFolder};
            normSIGNAL_2    = normSIGNAL_2_all{idxFolder};
            X_val_all       = cellfun(@(x) colon(1,1,find(~isnan(x),1,'last'))'./ p.px_um ,normSIGNAL_1,'unif',false);
            X_val_all       = catData(X_val_all); X_val_all = num2cell(X_val_all,1);

            tempSAVE_DIR    = SAVE_DIR{idxFolder};

            [orig_normSIGNAL_1, orig_normSIGNAL_2, normSIGNAL_1, normSIGNAL_2, ...
              sm_normSIGNAL_1,   sm_normSIGNAL_2] = smooothPeaks(normSIGNAL_1, normSIGNAL_2, p);
        end
        

        %% Find and plot peaks in signals
        % This code keeps one set of peaks static and randomizes position of the 
        % second set of peaks. It addresses the question: would the distances
        % between the peaks of Signal_1 and Signal_2 be expected if Signal_2 were
        % randomly distributed.

        % Plot original + overlay of the two signals for all line scans
        %   20160612 Updated this section to use the appropriate length X_val
        %   vector
        if p.PLOT_1
            fprintf('Plotting and exporting original signal plus overlay...\n')
            for idx     = 1:size(normSIGNAL_1,2);
                [X_val , temp_orig_normSIGNAL_1, nonNAN_IDX]  = removeNan(X_val_all{1,idx},orig_normSIGNAL_1{1,idx});
                
                validateattributes(temp_orig_normSIGNAL_1,{'numeric'},{'nonempty','nonnan','numel',numel(X_val)});
                
                temp_orig_normSIGNAL_2  = orig_normSIGNAL_2{1,idx}(nonNAN_IDX);
                temp_sm_normSIGNAL_1    = sm_normSIGNAL_1{1,idx}(nonNAN_IDX);
                temp_sm_normSIGNAL_2    = sm_normSIGNAL_2{1,idx}(nonNAN_IDX);
                temp_normSIGNAL_1       = normSIGNAL_2{1,idx}(nonNAN_IDX);
                temp_normSIGNAL_2       = normSIGNAL_2{1,idx}(nonNAN_IDX);
                
                try
                    % Plot the overlay of the original signals
                    figure(1),  plot(X_val,  temp_orig_normSIGNAL_1); hold on;
                                plot(X_val,  temp_orig_normSIGNAL_2);
                    xlabel(p.XLABEL,'interpreter','tex');
                    ylabel('Normalized Signal Amplitude (AU)','interpreter','tex');
                    legend(p.SIG_NAMES);
                    ylim([-0.2,1.2]);
                    xlim('auto')
                    try
                        error('!') %JN 2017
                        export_fig([tempSAVE_DIR filesep ...
                                    sprintf(['Original_Line_Scans_overlay_%d' p.EXT],idx)]); % Save overlay of original signals
                    catch
                        p.EXT = strrep(p.EXT,'.','');
                        handles.p.EXT = p.EXT;
                        savefig([tempSAVE_DIR filesep ...
                                    sprintf(['Original_Line_Scans_overlay_%d' '.fig'],idx)]); % Save overlay of original signals
%                         saveas(gcf, [tempSAVE_DIR filesep ...
%                                     sprintf(['Original_Line_Scans_overlay_%d'],idx)], p.EXT); % Save overlay of original signals
                    end


                    figure(2),  plot(X_val, temp_sm_normSIGNAL_1); hold on;
                                plot(X_val,  temp_sm_normSIGNAL_2);
                    xlabel(p.XLABEL,'interpreter','tex');
                    ylabel('Smoothed Normalized Signal Amplitude (AU)','interpreter','tex');
                    legend(p.SIG_NAMES);
                    ylim([-0.2,1.2]);
                    xlim('auto')
                    try
                        error('!') %JN 2017
                        export_fig([tempSAVE_DIR filesep ...
                                    sprintf(['Smoothed_SIG-2_Line_Scans_overlay_%d' p.EXT],idx)]); % Save overlay of original signals
                    catch
                        savefig([tempSAVE_DIR filesep ...
                                    sprintf(['Smoothed_SIG-2_Line_Scans_overlay_%d.fig'],idx)]); % Save overlay of original signals
%                         saveas(gcf, [tempSAVE_DIR filesep ...
%                                     sprintf(['Smoothed_SIG-2_Line_Scans_overlay_%d'],idx)],p.EXT); % Save overlay of original signals
                    end

                    % Plot the identified peaks
                    figure(3);
                    if p.USE_SMOOTH_1
                        findpeaks(temp_normSIGNAL_1     ,X_val, p.PEAK_METHOD,0.15,p.PEAK_ARGS{:},...
                                    'MinPeakDistance',p.MINPEAKDIST); % Peaks 1 (static)
                    else
                        findpeaks(temp_orig_normSIGNAL_1,X_val, p.PEAK_METHOD,0.15,p.PEAK_ARGS{:},...
                                    'MinPeakDistance',p.MINPEAKDIST); % Peaks 1 (static)
                    end

                    title(['SIG 1 Peaks ' p.SIG_NAMES{1}])
                    xlabel(p.XLABEL,'interpreter','tex');
                    ylabel('Normalized Signal Amplitude (AU)','interpreter','tex');
                    ylim([-0.2,1.2]);
                    xlim('auto')
                    try
                        error('!') %JN 2017
                        export_fig([tempSAVE_DIR filesep sprintf([p.SIG_NAMES{1} '_Peaks_%d' ...
                                            p.EXT],idx)]); % Save overlay of original signals
                    catch
                        savefig([tempSAVE_DIR filesep sprintf([p.SIG_NAMES{1} '_Peaks_%d' ...
                                            '.fig'],idx)]); % Save overlay of original signals
%                         saveas(gcf, [tempSAVE_DIR filesep sprintf([p.SIG_NAMES{1} '_Peaks_%d' ...
%                                             ],idx)], p.EXT); % Save overlay of original signals
                    end

                    figure(4); 
                    if p.USE_SMOOTH_2
                        findpeaks(temp_normSIGNAL_2       ,X_val ,p.PEAK_METHOD,0.05,p.PEAK_ARGS{:},...
                                    'MinPeakDistance',p.MINPEAKDIST); % Peaks 2 (to-randomize)
                    else
                        findpeaks(temp_orig_normSIGNAL_2  ,X_val,p.PEAK_METHOD,0.05,p.PEAK_ARGS{:},...
                                    'MinPeakDistance',p.MINPEAKDIST); % Peaks 2 (to-randomize)
                    end

                    title(['SIG 2 Peaks ' p.SIG_NAMES{2}])
                    xlabel(p.XLABEL,'interpreter','tex');
                    ylabel('Normalized Signal Amplitude (AU)','interpreter','tex');
                    ylim([-0.2,1.2]);
                    xlim('auto')
                    try
                        error('!') %JN 2017
                        export_fig([tempSAVE_DIR filesep sprintf([p.SIG_NAMES{2} '_Peaks_%d' ...
                                            p.EXT],idx)]); % Save overlay of original signals
                    catch
                        savefig([tempSAVE_DIR filesep sprintf([p.SIG_NAMES{2} '_Peaks_%d' ...
                                            '.fig'],idx)]); % Save overlay of original signals
%                         saveas(gcf, [tempSAVE_DIR filesep sprintf([p.SIG_NAMES{2} '_Peaks_%d' ...
%                                             ],idx)], p.EXT); % Save overlay of original signals
                    end
                    close 1 2 3 4
                catch
                    close 1 2 3 4
                end
            end
            
            clear temp_orig_normSIGNAL_1    temp_orig_normSIGNAL_2 ...
                      temp_sm_normSIGNAL_1      temp_sm_normSIGNAL_2 ...
                      temp_normSIGNAL_1         temp_normSIGNAL_2 X_val
        end



        %% Save the locations of the peaks
        if p.CALC_PEAKS
            [X_val_1 , temp_normSIGNAL_1, ~]  = ...
                        cellfun(@(x,y) removeNan(x,y), X_val_all, normSIGNAL_1,'unif',false);
            [X_val_2 , temp_normSIGNAL_2, ~]  = ...
                        cellfun(@(x,y) removeNan(x,y), X_val_all, normSIGNAL_2,'unif',false);
            
            fprintf('Calculating peaks...\n')
            [PKS_1, LOC_1, W_1]= cellfun(@(x,y) findpeaks(x,y,p.PEAK_METHOD,0.15,p.PEAK_ARGS{:},...
                                                    'MinPeakDistance',p.MINPEAKDIST),...
                                                temp_normSIGNAL_1, X_val_1,'unif',false);
            [PKS_2, LOC_2, W_2]= cellfun(@(x,y) findpeaks(x,y, p.PEAK_METHOD,0.05,p.PEAK_ARGS{:},...
                                                    'MinPeakDistance',p.MINPEAKDIST),...
                                                temp_normSIGNAL_2, X_val_2,'unif',false);

            numPK_1     = cellfun(@(x) numel(x), LOC_1,'unif',false);
            numPK_2     = cellfun(@(x) numel(x), LOC_2,'unif',false);
            
            clear temp_normSIGNAL_1 temp_normSIGNAL_2 X_val_1 X_val_2
        else
            PKS_1       = cellfun(@(x,y) x(y), normSIGNAL_1, {handles.EB_Info.syn_pk},'unif',false);
            LOC_1       = {handles.EB_Info.syn_pk_um};
            W_1         = {handles.EB_Info.syn_width};
            
            fprintf('Processing %s data\n',lower(p.EB_ANALYSIS))
            p.EB_ANALYSIS = 'init'; % Edit locally JN 20180626
            switch lower(p.EB_ANALYSIS)
                case {'initiation', 'init','antero','anterograde'}
                    PKS_2       = cellfun(@(x) ones(numel(x),1),{handles.EB_Init.x},'unif',false);
                    LOC_2       = {handles.EB_Init.x_um};
                    LOC_2       = cellfun(@(x) x',LOC_2,'unif',false);
                    W_2         = PKS_2;
                case {'termination','term','retro','retrograde'}
                    PKS_2       = cellfun(@(x) ones(numel(x),1),{handles.EB_Term.x},'unif',false);
                    LOC_2       = {handles.EB_Term.x_um};
                    LOC_2       = cellfun(@(x) x',LOC_2,'unif',false);
                    W_2         = PKS_2;
                otherwise
                    error('Unknown analysis method for EB/ Pause analysis.')
            end
            
            numPK_1     = cellfun(@(x) numel(x), LOC_1,'unif',false);
            numPK_2     = cellfun(@(x) numel(x), LOC_2,'unif',false);
        end

        if p.CALC_SYN_WIDTH
            fprintf('Calculating synapsin peak width...\n')
            % Calculate the Synapsin peak width
            LOC_2_geq_2     = cellfun(@(x) numel(x)> 1,LOC_2,'unif',true);
            idxbest         = cell(size(LOC_2)); idxbest(~LOC_2_geq_2)  = {1};
            cBest           = cell(size(LOC_2)); cBest(~LOC_2_geq_2)    = LOC_2(~LOC_2_geq_2);
            OptK            = cell(size(LOC_2)); OptK(~LOC_2_geq_2)     = {1};

            if p.CLUSTER_PEAKS
                fprintf('K-means clustering of peaks temporarily disabled...\n')
%                 fprintf('K-means clustering of synapsin peaks...\n\tThis may take some time...\n')
% 
%                     [idxbest(LOC_2_geq_2), cBest(LOC_2_geq_2), OptK(LOC_2_geq_2), ~, ~]= ...
%                         cellfun(@(x) kCluster(x'.*(p.px_um),'K',1:1:max([1,numel(x)]),...
%                                             'PARALLEL',p.PARALLEL,...
%                                             'DISTANCE',p.DISTANCE), LOC_2(LOC_2_geq_2),'unif',false);
%                 
%                 cBest(LOC_2_geq_2)= cellfun(@(x) x./(p.px_um),cBest(LOC_2_geq_2),'unif',false);
%                 groupedData_2   = cell(size(LOC_2));
%                     groupedData_2(~LOC_2_geq_2)     = W_2(~LOC_2_geq_2);
%                     groupedData_2_sum(~LOC_2_geq_2) = W_2(~LOC_2_geq_2);
%                 groupedData_2(LOC_2_geq_2)   = cellfun(@(x,y) groupUnique([x(:),y(:)]),idxbest(LOC_2_geq_2),W_2(LOC_2_geq_2),'unif',false);
%                 groupedData_2_sum(LOC_2_geq_2)= cellfun(@(x) nansum(x(:,2:end),2),groupedData_2(LOC_2_geq_2),'unif',false);
%                 numPK_2         = cellfun(@(x) numel(unique(x)),idxbest,'unif',false);
%                 
%                 fprintf('Re-calculating synapsin peaks given peak clustering...\n')
%                 PKS_2           = numPK_2;
%                 W_2             = cellfun(@(x) x',groupedData_2_sum,'unif',false);
%                 LOC_2           = cellfun(@(x) x',cBest,'unif',false);
            end

            % Re-calculate synapsin peaks based on new peak number
            if false
                if and(p.SET_PEAK_N, p.CLUSTER_PEAKS)
                    fprintf('Re-calculating synapsin peaks given peak clustering...\n\tThis may take some time...\n')

                    [X_val , temp_sm_normSIGNAL_2, ~]  = ...
                            cellfun(@(x,y) removeNan(x,y), X_val_all, sm_normSIGNAL_2,'unif',false);

                    [PKS_2, LOC_2, W_2]= cellfun(@(x,y,z) findpeaks(x,z,p.PEAK_METHOD,0.05,p.PEAK_ARGS{:},'NPeaks',y+1,...
                                                        'MinPeakDistance',p.MINPEAKDIST),...
                                                        temp_sm_normSIGNAL_2, numPK_2, X_val,'unif',false);

                    clear X_val temp_sm_normSIGNAL_2

                    if p.PLOT_1
                        for cellIdx = 1:numel(normSIGNAL_2)
                            [X_val , temp_sm_normSIGNAL_2]  = ...
                                removeNan(X_val_all{1,idx},sm_normSIGNAL_2{1,idx});

                            figure(1)
                            findpeaks(temp_sm_normSIGNAL_2,X_val,p.PEAK_METHOD,0.05,p.PEAK_ARGS{:},...
                                'NPeaks',numPK_2{1,cellIdx}+1,'MinPeakDistance',p.MINPEAKDIST);
                            title('Smoothed Synapsin Peaks')
                            xlabel(p.XLABEL,'interpreter','tex');
                            ylabel('Normalized Signal Amplitude (AU)','interpreter','tex');
                            ylim([-0.2,1.2]);
                            xlim('auto')
                            try
                                error('!') %JN 2017
                                export_fig([tempSAVE_DIR filesep sprintf(['smooth_Synapsin_Peaks_%d' p.EXT],cellIdx)]); % Save overlay of original signals
                            catch
                                savefig([tempSAVE_DIR filesep sprintf(['smooth_Synapsin_Peaks_%d.fig'],cellIdx)]); % Save overlay of original signals
%                                 saveas(gcf, [tempSAVE_DIR filesep sprintf(['smooth_Synapsin_Peaks_%d'],cellIdx)], p.EXT); % Save overlay of original signals
                            end
                            close(1)
                        end
                    end
                elseif and(~p.SET_PEAK_N, p.CLUSTER_PEAKS)
                    fprintf('Re-calculating synapsin peaks given peak clustering...\n')
                    [X_val , temp_sm_normSIGNAL_2]  = ...
                            cellfun(@(x,y) removeNan(x,y), X_val_all, sm_normSIGNAL_2,'unif',false);

                    [PKS_2, LOC_2, W_2]= cellfun(@(x,z) findpeaks(x,z, p.PEAK_METHOD,0.05,p.PEAK_ARGS{:}),...
                                                        temp_sm_normSIGNAL_2, X_val, 'unif',false);

                    clear X_val temp_sm_normSIGNAL_2                                
                    if p.PLOT_1
                        for cellIdx = 1:numel(normSIGNAL_2)
                            [X_val , temp_sm_normSIGNAL_2]  = ...
                                removeNan(X_val_all{1,cellIdx},sm_normSIGNAL_2{1,cellIdx});

                            figure(1)
                            findpeaks(temp_sm_normSIGNAL_2, X_val, p.PEAK_METHOD,0.05,p.PEAK_ARGS{:},...
                                'MinPeakDistance',p.MINPEAKDIST)
                            title('Smoothed Synapsin Peaks')
                            xlabel(p.XLABEL,'interpreter','tex');
                            ylabel('Normalized Signal Amplitude (AU)','interpreter','tex');
                            ylim([-0.2,1.2]);
                            xlim('auto')
                            try
                                error('!') %JN 2017
                                export_fig([tempSAVE_DIR filesep sprintf(['smooth_Synapsin_Peaks_%d' p.EXT],cellIdx)]); % Save overlay of original signals
                            catch
                                savefig([tempSAVE_DIR filesep sprintf(['smooth_Synapsin_Peaks_%d.fig'],cellIdx)]); % Save overlay of original signals
%                                 saveas(gcf,[tempSAVE_DIR filesep sprintf(['smooth_Synapsin_Peaks_%d'],cellIdx)], p.EXT); % Save overlay of original signals
                            end
                            close(1)
                        end
                    end
                end
            end
           
            numPK_1     = cellfun(@(x) numel(x), LOC_1,'unif',false);
            numPK_2     = cellfun(@(x) numel(x), LOC_2,'unif',false);
        end


        if p.PEAK_RAND
            fprintf('Randomizing peaks and calculating nearest neighbor distances...\n')
            [X_val_1 , temp_normSIGNAL_1, ~]  = ...
                        cellfun(@(x,y) removeNan(x,y), X_val_all, normSIGNAL_1,'unif',false);
            [X_val_2 , temp_normSIGNAL_2, ~]  = ...
                        cellfun(@(x,y) removeNan(x,y), X_val_all, normSIGNAL_2,'unif',false);
            
            
            % Allocate output variables
            randPK_1    = cell(p.NSAMPLE, size(normSIGNAL_1,2));
            randPK_2    = cell(p.NSAMPLE, size(normSIGNAL_2,2));

            % Print to command window
            % JN 20180609
            printTextRand   = {'Signal 1','Signal 2'};
            fprintf('Randomizing %s\n', printTextRand{[p.RAND_1, p.RAND_1]})
            
            for repIdx = 1:p.NSAMPLE
                if and(p.RAND_1, p.RAND_2)
                    % Randomize both
                    error('peakAnalysis has temporarily disabled randomizing both signals!')
                    randPK_1(repIdx,:)  = cellfun(@(x,y,z) datasample(z,y,'replace',true)',...
                                                temp_normSIGNAL_1, numPK_1, X_val_1,'unif',false);
                    randPK_2(repIdx,:)  = cellfun(@(x,y,z) datasample(z,y,'replace',true)',...
                                                temp_normSIGNAL_2, numPK_2, X_val_2,'unif',false);
                
                elseif p.RAND_1
                    % Randomize SIG1
                    randPK_1(repIdx,:)  = cellfun(@(x,y,z) datasample(z,y,'replace',true)',...
                                                temp_normSIGNAL_1, numPK_1, X_val_1,'unif',false);
                elseif p.RAND_2
                    % Randomize SIG_2
                    randPK_2(repIdx,:)  = cellfun(@(x,y,z) datasample(z,y,'replace',true)',...
                                                temp_normSIGNAL_2, numPK_2, X_val_2,'unif',false);
                end
            end

            % Calculate distance between both sets of randomized peaks
%             Idx_Dist_2  = cell(1, size(normSIGNAL_1,2)); % Why do I not need this JN 20180609?
            Dist_2      = cell(1, size(normSIGNAL_1,2));

            %% knnsearch finds the nearest K neighbors in X, for each value in Y
            for colIdx  = 1:size(normSIGNAL_1,2)
                % JN 20180609
                if and(p.RAND_1, p.RAND_2)
                    error('peakAnalysis has temporarily disabled randomizing both signals!')
                    [~, Dist_2{colIdx}] = cellfun(@(x,y) knnsearch(x(:),y(:),'k',p.KNN), randPK_1(:,colIdx),...
                                                         randPK_2(:,colIdx),'unif',false);
                elseif p.RAND_1
                    % Random Sig_1 with Original Sig_2
                    [Idx_Dist_2, Dist_2{colIdx}] = cellfun(@(x,y) knnsearch(x(:),y(:),'k',p.KNN), randPK_1(:,colIdx),...
                                     repmat(LOC_2(1,colIdx),p.NSAMPLE,1),'unif',false);                
                elseif p.RAND_2
                    % Original Sig_1 with Random Sig_2 (Default)
                    [Idx_Dist_2, Dist_2{colIdx}] = cellfun(@(x,y) knnsearch(x(:),y(:),'k',p.KNN),repmat(LOC_1(1,colIdx),p.NSAMPLE,1),...
                                                         randPK_2(:,colIdx),'unif',false);
                end
                
                % JN 20180609 - need to update randpk2 in signum
                Idx_Dist_2      = cellfun(@(x) num2cell(x,1),Idx_Dist_2,'unif',false);
                Idx_Dist_2      = vertcat(Idx_Dist_2{:});
                
                % Added check for numel==1 JN 20170529
                if numel(LOC_1{1,colIdx})==1 % include numel(LOC_2{1,colIdx})==1
                    % Z is not right here
                    if p.RAND_1
                        signDist        = cell2mat(cellfun(@(Syn,y,z) sign(z(:)-Syn(y)),...
                                            repmat(randPK_1(:,colIdx),1,p.KNN), ...
                                            Idx_Dist_2, ...
                                            repmat(LOC_2(1,colIdx),p.NSAMPLE,p.KNN), 'unif',false));
                    else
                        signDist        = cell2mat(cellfun(@(Syn,y,z) sign(z(:)-Syn(y)),...
                                            repmat(LOC_1(1,colIdx),p.NSAMPLE,p.KNN), ...
                                            Idx_Dist_2, ...
                                            repmat(randPK_2(:,colIdx),1,p.KNN), 'unif',false));
                    end
                else
                    if p.RAND_1
                        %% Update
                        signDist        = cell2mat(cellfun(@(Syn,y,z) sign(z-Syn(y)),...
                                            repmat(randPK_1(:,colIdx),1,p.KNN), ...
                                            Idx_Dist_2, ...
                                            repmat(LOC_2(1,colIdx),p.NSAMPLE,p.KNN), 'unif',false));
                    else
                        signDist        = cell2mat(cellfun(@(Syn,y,z) sign(z-Syn(y)),... % JN  20170529
                                        repmat(LOC_1(1,colIdx),p.NSAMPLE,p.KNN), ...
                                        Idx_Dist_2, ...
                                        repmat(randPK_2(:,colIdx),1,p.KNN), 'unif',false));
                    end
                    
                end
                Dist_2{colIdx}  = cell2mat(Dist_2{colIdx});
                Dist_2{colIdx}(:,1:p.KNN) = Dist_2{colIdx}(:,1:p.KNN).* (signDist(:));
            end
            nColSmall           = cellfun(@(x) size(x,2),Dist_2,'unif',true)<2;
            Dist_2(nColSmall)   = cellfun(@(x) [x nan(numel(x),1)], Dist_2(nColSmall),'unif',false);

%             Idx_orig   = cell(1, size(normSIGNAL_1,2));
            Dist_orig   = cell(1, size(normSIGNAL_1,2));
            for colIdx  = 1:size(normSIGNAL_1,2)
                % Find the distances in the original data
                [Idx_orig, Dist_orig{colIdx}] = cellfun(@(x,y) knnsearch(x(:),y(:),'k',p.KNN),LOC_1(:,colIdx),...
                                                         LOC_2(:,colIdx),'unif',false);
                    
                Idx_orig                = cellfun(@(x) num2cell(x,1),Idx_orig,'unif',false);
                Idx_orig                = vertcat(Idx_orig{:});
                if numel(LOC_1{1,colIdx})==1
                    signDist                = cell2mat(cellfun(@(Syn,y,z) sign(z(:)-Syn(y(:,1)))',... % Added JN  20170529
                                                repmat(LOC_1(:,colIdx),1,p.KNN), ...
                                                Idx_orig, ...
                                                repmat(LOC_2(:,colIdx),1,p.KNN), 'unif',false));
                else
                    signDist                = cell2mat(cellfun(@(Syn,y,z) sign(z-Syn(y(:,1)))',... % Added JN  20170529
                                                repmat(LOC_1(:,colIdx),1,p.KNN), ...
                                                Idx_orig, ...
                                                repmat(LOC_2(:,colIdx),1,p.KNN), 'unif',false));
                end
                Dist_orig{colIdx}       = cell2mat(Dist_orig{colIdx});
                Dist_orig{colIdx}(:,1:p.KNN)  = Dist_orig{colIdx}(:,1:p.KNN).* (signDist(:));
            end
            Dist_orig(nColSmall)   = cellfun(@(x) [x nan(numel(x),1)], Dist_orig(nColSmall),'unif',false);
            clear temp_normSIGNAL_1 temp_normSIGNAL_2 X_val_1 X_val_2
        else
            Dist_2          = cellfun(@(x) repmat(nan,size(x,1),size(x,2)),...
                                PKS_1,'unif',false);
            Dist_orig       = cellfun(@(x) repmat(nan,size(x,1),size(x,2)),...
                                PKS_2,'unif',false);
        end



        if p.CALC_PTM_SYN_OVERLAP 
            fprintf('Calculating PTM signal density within synapsin FWHM\n')
            temp_X_val          = (1:1:numel(normSIGNAL_1{1}))./p.px_um;
            [X_val_1 , temp_orig_normSIGNAL_1, ~]  = ...
                        cellfun(@(x,y) removeNan(x,y), X_val_all, orig_normSIGNAL_1,'unif',false);
            [X_val_2 , temp_orig_normSIGNAL_2, ~]  = ...
                        cellfun(@(x,y) removeNan(x,y), X_val_all, orig_normSIGNAL_2,'unif',false);

            % Check vector length
            errIdx              = cellfun(@(x,y) length(x)~=length(y), temp_orig_normSIGNAL_1, temp_orig_normSIGNAL_2,'unif',true); 
            if any(errIdx)
                fprintf('Signal 1 and Signal 2 have different lengths: \t%d\n',find(errIdx))
            end
            
            LOC_2_OVERLAP       = cellfun(@(x) x(:)',LOC_2,'unif',false);
            W_2_OVERLAP         = cellfun(@(x) x(:)',W_2,'unif',false);
            
            % Compute the PTM signal density within the Synapsin FWHM
            W_2_OVERLAP         = cellfun(@(x,y,z) [x-(y./2); min([x+(y./2); repmat(max(z),1,numel(x))],[],1)],...
                              LOC_2_OVERLAP,W_2_OVERLAP,X_val_2,'unif',false); % Added 20160603
            W_2_OVERLAP         = cellfun(@(x) num2cell(x,1),W_2_OVERLAP,'unif',false);
            W_2_OVERLAP         = cellfun(@(x) cellfun(@(y) y',x,'unif',0),W_2_OVERLAP,'unif',false);
                                   
            WITHIN_BOUNDS       = @(x,BOUNDS) (x >= min(BOUNDS) & x <= max(BOUNDS));
            [PK_W_LOC_2]= cellfun(@(x) cellfun(@(y) find(WITHIN_BOUNDS(temp_X_val,y)),x,'unif',false),...
                            W_2_OVERLAP,'unif',false); % Added 20160603
            
            PK_W_LOC_2  = cellfun(@(x) unique(horzcat(x{:})), PK_W_LOC_2 ,'unif',false);
            
            PTM_SYN_OVERLAP     = cellfun(@(x,y) x(y(:)), temp_orig_normSIGNAL_1, PK_W_LOC_2, 'unif', false);
            emptyCell           = cellfun(@isempty, PK_W_LOC_2,'unif',true);
            
            warning('off')
            stats_PTM_SYN_OVERLAP= cell(1,numel(PTM_SYN_OVERLAP));
            stats_PTM_SYN_OVERLAP(~emptyCell)= cellfun(@(x) firstOrderStats(x,'CALC_IDX',p.CALC_IDX),...
                PTM_SYN_OVERLAP(~emptyCell),'unif',false);
            [~,STATS]            = firstOrderStats(rand(50),'CALC_IDX',p.CALC_IDX);
       
            ColorVec = distinguishable_colors(25);
            ColorVec = repmat({num2cell(ColorVec,2)},1,numel(PKS_2));
            
            clear temp_normSIGNAL_1 temp_normSIGNAL_2 X_val_1 X_val_2
            if p.PLOT_2
                
                fprintf('Plotting synapspin FWHM on PTM signal...\n')
                for idx = 1:size(normSIGNAL_1,2);
                    % Updated 20160612 to use X_val of length equal to the
                    % line scan length
                    try
                        [X_val_1 , temp_normSIGNAL_1, nonNAN_IDX]  = ...
                            removeNan(X_val_all{1,idx},normSIGNAL_1{1,idx});
                        [X_val_2 , temp_normSIGNAL_2, nonNAN_IDX]  = ...
                            removeNan(X_val_all{1,idx},normSIGNAL_2{1,idx});

                        tempColorVec = cellfun(@(x) x,ColorVec,'unif',false); % (1:numel(W_2_OVERLAP{idx}))
                        figure(1),  plot(X_val_1,  temp_normSIGNAL_1);
                        hold on,    plot(X_val_2,  temp_normSIGNAL_2);

                        cellfun(@(x,y) line(x,[1.1,1.1],'LineWidth',5,'Color',y),...
                            W_2_OVERLAP{idx}',tempColorVec{idx}(1:numel(W_2_OVERLAP{idx})),'unif',false);

                        xlabel(p.XLABEL,'interpreter','tex');
                        ylabel('Normalized Signal Amplitude (AU)','interpreter','tex');
                        legend({p.SIG_NAMES{1},p.SIG_NAMES{2},'SIG2 FWHM'});
                        ylim([-0.2,1.2]);
                        try
                            error('!') %JN 2017
                            export_fig([tempSAVE_DIR filesep ...
                                    sprintf(['Original_Line_Scans_overlay_SIG_DENSITY_%d' p.EXT],idx)]); % Save overlay of original signals
                        catch
                            p.EXT = strrep(p.EXT,'.','');
                            handles.p.EXT = p.EXT;
                            savefig([tempSAVE_DIR filesep ...
                                    sprintf(['Original_Line_Scans_overlay_SIG_DENSITY_%d.fig'],idx)]); % Save overlay of original signals
%                             saveas(gcf, [tempSAVE_DIR filesep ...
%                                     sprintf(['Original_Line_Scans_overlay_SIG_DENSITY_%d'],idx)],p.EXT); % Save overlay of original signals
                        end
                        close(1)
                    catch
                        close(1)
                    end
                end    
            end
        else
            PTM_SYN_OVERLAP = [];
            stats_PTM_SYN_OVERLAP = [];
        end

        
        if and(p.GRAPH_OUTPUT,p.PEAK_RAND)
            fprintf('Graphing nearest neighbor distances for randomized data...\n')
            % Graph the output
            color_vec   = distinguishable_colors(2);
            for colIdx  = 1:size(normSIGNAL_1,2)
                try
                    figure(1);
                    title(sprintf('Randomized distribution (%s)',p.EB_ANALYSIS))
                    printText   = mergeArray(1:p.KNN,nanmean(Dist_2{colIdx}(:,1:p.KNN),1));
                    legendText  = sprintf('knn-%d mu=%0.2f,',printText);
                    legendText(end) = []; % Remove last comma
                    nhist(num2cell(Dist_2{colIdx},1)','legend',strsplit(legendText,','),...
                        'box','xlabel',['KNN' p.XLABEL],'color',color_vec);
                    Xplot   = nanmean(Dist_orig{colIdx}(:,1:p.KNN),1);
                    if numel(Xplot)< 2
                        tempX= nan(1,2); tempX(1:numel(Xplot)) = Xplot;
                        Xplot= tempX;
                    end
                    Xplot   = repmat(Xplot,2,1);
                    Yplot   = repmat([0;max(ylim)],1,numel(Xplot)/2);
                    cellfun(@(x,y,c) line('XData',x,'YData',y,'Color',c), ...
                            num2cell(Xplot,1), num2cell(Yplot,1), num2cell(color_vec,2)','unif',false);
                    try
                        error('!') %JN 2017
                        export_fig([tempSAVE_DIR  filesep ...
                                sprintf(['NN_Dist_Line_Scans_%d_%s' p.EXT],colIdx,p.EB_ANALYSIS)]); % Save figure
                    catch
                        savefig([tempSAVE_DIR  filesep ...
                                sprintf(['NN_Dist_Line_Scans_%d_%s.fig'],colIdx,p.EB_ANALYSIS)]); % Save figure
%                         saveas(gcf,[tempSAVE_DIR  filesep ...
%                                 sprintf(['NN_Dist_Line_Scans_%d_%s' p.EXT],colIdx,p.EB_ANALYSIS)]); % Save figure
                    end
                catch
                    fprintf('Error idx %d\n',colIdx)
                end
            end
        end

        if and(p.GRAPH_COMPILED,p.PEAK_RAND)
            % Graph the compiled output
            figure(1)
            clf
            color_vec       = distinguishable_colors(p.KNN);
            compiled        = num2cell(cell2mat(vertcat(Dist_2')),1);
            compiled_D_orig = num2cell(cell2mat(vertcat(Dist_orig')),1);
            
            printText   = mergeArray(1:p.KNN,nanmean(cell2mat(compiled(:,1:p.KNN)),1));
            legendText  = sprintf('knn-%d mu=%0.2f,',printText);
            legendText(end) = []; % Remove last comma
            if p.NND_BIPLOT
                AX(1)= subplot(2,1,1);
            end
            title(sprintf('Randomized distribution (%s)',p.EB_ANALYSIS))
            nhist(compiled(:,1:p.KNN),'legend',strsplit(legendText,','),...
                    'box','color',color_vec);
            
            % Plot line for mean of original data    
            Xplot           = nanmean(cell2mat(compiled_D_orig(1:p.KNN)),1);
            Xplot           = repmat(Xplot,p.KNN,1);
            Yplot           = repmat([0;max(ylim)],1,numel(Xplot)/p.KNN);
            cellfun(@(x,y,c) line('XData',x,'YData',y,'Color',c), num2cell(Xplot,1), num2cell(Yplot,1), num2cell(color_vec,2)','unif',false);
            X_LIM_1         = xlim(gca);
            Y_LIM_1         = ylim(gca);
            
            %
            if p.NND_BIPLOT
                AX(2)= subplot(2,1,2);
                title(sprintf('Original distribution (%s)',p.EB_ANALYSIS))
                printText   = mergeArray(1:p.KNN,nanmean(cell2mat(compiled_D_orig(:,1:p.KNN)),1));
                legendText  = sprintf('knn-%d mu=%0.2f,',printText);
                legendText(end) = []; % Remove last comma
                nhist(compiled_D_orig(:,1:p.KNN),'legend',strsplit(legendText,','),...
                        'box','xlabel',['KNN ' p.XLABEL],'color',color_vec);
                X_LIM_2         = xlim(gca);
                Y_LIM_2         = ylim(gca);

                XLIM_final      =   [min([X_LIM_1(:);X_LIM_2(:)]), max([X_LIM_1(:);X_LIM_2(:)])];
                YLIM_final      =   [min([Y_LIM_1(:);Y_LIM_2(:)]), max([Y_LIM_1(:);Y_LIM_2(:)])];

                xlim(XLIM_final); ylim(YLIM_final);
                axes(AX(1)); xlim(XLIM_final); ylim(YLIM_final);
                linkaxes(AX,'xy');
            end


            try
                error('!') %JN 2017
                export_fig([filesep tempSAVE_DIR  filesep [sprintf('NN_Dist_Line_Scans_compiled_%s',p.EB_ANALYSIS) p.EXT]]); % Save figure
            catch
                savefig([tempSAVE_DIR  filesep ['NN_Dist_Line_Scans_compiled.fig']]); % Save figure
%                 saveas(gcf, [tempSAVE_DIR  filesep ['NN_Dist_Line_Scans_compiled' p.EXT]]); % Save figure
            end

            close(1);
            clear temp_normSIGNAL_1     temp_normSIGNAL_2   X_val
        end

        %% Save data
        SaveFilename    = [tempSAVE_DIR  filesep sprintf('compiled_Data_%s.mat',p.EB_ANALYSIS)];

        if exist('compiled_data','var')
            fprintf('Over-writing previous compiled_Data...')
%             compiled_Data   = struct('Dist_2',{compiled_Data.Dist_2},'Dist_orig',{compiled_Data.Dist_orig},...
%                                      'RandPK_2', num2cell(randPK_2,1), ...
%                                      'PKS_1',PKS_1,'PKS_2',PKS_2,...
%                                      'LOC_1',LOC_1,'LOC_2',LOC_2',...
%                                      'W_1',W_1,'W_2',W_2, 'numPK_1',numPK_1,'numPK_2',numPK_2,...
%                                      'PTM_SYN_OVERLAP', {compiled_Data.PTM_SYN_OVERLAP},...
%                                      'stats_PTM_SYN_OVERLAP',{compiled_Data.stats_PTM_SYN_OVERLAP});
            compiled_Data   = struct('Dist_2',Dist_2,'Dist_orig',Dist_orig,...
                                     'RandPK_1', num2cell(randPK_1,1), ...
                                     'RandPK_2', num2cell(randPK_2,1), ...
                                     'PKS_1',PKS_1,'PKS_2',PKS_2,...
                                     'LOC_1',LOC_1,'LOC_2',LOC_2,...
                                     'W_1',W_1,'W_2',W_2, 'numPK_1',numPK_1,'numPK_2',numPK_2,...
                                     'PTM_SYN_OVERLAP', PTM_SYN_OVERLAP,...
                                     'stats_PTM_SYN_OVERLAP',stats_PTM_SYN_OVERLAP,...
                                     'Analysis',p.EB_ANALYSIS);
            save(SaveFilename,'compiled_Data');
        else
            compiled_Data   = struct('Dist_2',Dist_2,'Dist_orig',Dist_orig,...
                                     'RandPK_1', num2cell(randPK_1,1), ...
                                     'RandPK_2', num2cell(randPK_2,1), ...
                                     'PKS_1',PKS_1,'PKS_2',PKS_2,...
                                     'LOC_1',LOC_1,'LOC_2',LOC_2,...
                                     'W_1',W_1,'W_2',W_2, 'numPK_1',numPK_1,'numPK_2',numPK_2,...
                                     'PTM_SYN_OVERLAP', PTM_SYN_OVERLAP,...
                                     'stats_PTM_SYN_OVERLAP',stats_PTM_SYN_OVERLAP,...
                                     'Analysis',p.EB_ANALYSIS);
            save(SaveFilename,'compiled_Data');
        end

        % Export to Excel
        if p.EXCEL
            error('RandPk_1 does not export 20180609')
            SHEETS = {'OriginalDistance','RandomizationDistance','Rand_PK_2',...
                      ['PK_1_AMP_' p.SIG_NAMES{1}], ['PK_1_Loc_' p.SIG_NAMES{1}], ...
                      ['PK_1_FWHM_' p.SIG_NAMES{1}], ['N_PK_1_' p.SIG_NAMES{1}], ...
                      ['PK_2_AMP_' p.SIG_NAMES{2}], ['PK_2_Loc_' p.SIG_NAMES{2}], ...
                      ['PK_2_FWHM_' p.SIG_NAMES{2}], ['N_PK_2_' p.SIG_NAMES{2}], ...
                      [p.SIG_NAMES{1},'_',p.SIG_NAMES{2},'_','_overlap'],...
                      [p.SIG_NAMES{1},'_',p.SIG_NAMES{2},'_','_overlap_stats']};

            DATA   = { {compiled_Data.Dist_orig}' , {compiled_Data.Dist_2}', randPK_2,...
                      PKS_1(:), LOC_1(:),...
                        W_1(:), numPK_1(:),...
                      PKS_2(:), LOC_2(:),...
                        W_2(:), numPK_2(:),...
                      {compiled_Data.PTM_SYN_OVERLAP}', {compiled_Data.stats_PTM_SYN_OVERLAP}' };
            DATA_INFO = {};%{{'NND_1','NND_2'},{'NND_1_Rand','NND_2_Rand'},...
                        %  strcat(p.SIG_NAMES{1},{'_Peak_AMP'}),...
                        % strcat(p.SIG_NAMES{1},{'_Peak_N'}),...
                        %  strcat(p.SIG_NAMES{1},{'_Peak_N'}),...
                        % };};
            try

                switch lower(computer)
                    case {'pcwin','pcwin64'}
                        [EXPORT_SUCCESS] = export2excel(SHEETS,DATA,...
                            'DATA_INFO',DATA_INFO,'FNAME',FOLDERS(idxFolder).name,...
                            'SAVEDIR',tempSAVE_DIR);
                    otherwise
                        fprintf('\nThe function export2Excel is not supported on your operating sustem.\n');
                        fprintf('Using export2csv instead\n')
                        [EXPORT_SUCCESS] = export2csv(SHEETS,DATA,...
                            'DATA_INFO',DATA_INFO,'FNAME',FOLDERS(idxFolder).name,...
                            'SAVEDIR',tempSAVE_DIR);
                end


            catch ERR
                disp(ERR.message)
                fprintf('Error writing to Excel!\nUsing export2csv instead\n')
                [EXPORT_SUCCESS] = export2csv(SHEETS,DATA,...
                    'DATA_INFO',DATA_INFO,'FNAME',FOLDERS(idxFolder).name,...
                    'SAVEDIR',tempSAVE_DIR);
            end
        else
            EXPORT_SUCCESS      = true;
        end

    end % END of FOR loop over all folders

catch ERR
    fprintf(ERR.message)
    rethrow(ERR)
    SUCCESS         = false;
    EXPORT_SUCCESS  = nan;
    error_Callback(hObject,[],handles)
end

% Allocate output variables
handles.SUCCESS         = SUCCESS;
handles.EXPORT_SUCCESS  = EXPORT_SUCCESS;
handles.ERR             = ERR;
handles.p               = p;

warning('on')
fprintf('Finished processing data!\n')
if p.EMAIL
    fprintf('Sending email...')
    sendMATLABMail(handles.p.EMAIL_OPT{:});
end

% Update handles structure
guidata(hObject, handles);


function Y = gaussFilter(X,WIN_SIZE, SIGMA)
% Subfunction for Gaussian filtering

nPixels     = length(X);

xArray      = linspace(-WIN_SIZE/ 2, WIN_SIZE/ 2, WIN_SIZE);
gFilt       = exp(-xArray .^ 2 / (2 * SIGMA ^ 2));
gFilt       = gFilt/ sum (gFilt); % normalize

Y           = conv (X, gFilt, 'same');


function [orig_normSIGNAL_1, orig_normSIGNAL_2, normSIGNAL_1, normSIGNAL_2, ...
          sm_normSIGNAL_1,   sm_normSIGNAL_2] = smooothPeaks(normSIGNAL_1, normSIGNAL_2, p)
% SubFunction to smooth peaks using specified method

switch lower(p.SMOOTH)
    case {'gauss','g','gaussian'}
        fprintf('Gaussian smoothing\n');
        % SIG 1
        if p.USE_SMOOTH_1
            sm_normSIGNAL_1     = cellfun(@(x) gaussFilter(x,p.WIN_SIZE, p.SIGMA), ...
                                          normSIGNAL_1,'unif',false);
        else
            sm_normSIGNAL_1     = normSIGNAL_1;
        end

        % SIG 2
        if p.USE_SMOOTH_2
            sm_normSIGNAL_2     = cellfun(@(x) gaussFilter(x,p.WIN_SIZE, p.SIGMA), ...
                                          normSIGNAL_2,'unif',false);
        else
            sm_normSIGNAL_2     = normSIGNAL_2;
        end
    case {'spline','sp'}
        fprintf('Spline smoothing\n');
        % SIG 1
        if p.USE_SMOOTH_1
            sm_normSIGNAL_1     = cellfun(@(x) smooth1q(x,[]), ...
                                          normSIGNAL_1,'unif',false);
        else
            sm_normSIGNAL_1     = normSIGNAL_1;
        end

        % SIG 2
        if p.USE_SMOOTH_2
            sm_normSIGNAL_2     = smooth1q(normSIGNAL_2,3);
        else
            sm_normSIGNAL_2     = normSIGNAL_2;
        end
    case {'win','window','sliding win mean'}
        fprintf('Sliding window mean smoothing\n');
        if p.USE_SMOOTH_1
            sm_normSIGNAL_1     = ...
                cellfun(@(x) slidingWindowFilter_GD(x,p.WIN_SIZE,@(x) nanmean(x,2)),normSIGNAL_1,'unif',false);
        else
            sm_normSIGNAL_1     = normSIGNAL_1;
        end

        if p.USE_SMOOTH_2
            sm_normSIGNAL_2     = ...
                cellfun(@(x) slidingWindowFilter_GD(x,p.WIN_SIZE,@(x) nanmean(x,2)),normSIGNAL_2,'unif',false);
        else
            sm_normSIGNAL_2     = normSIGNAL_2;
        end
    case {'','na'}
        fprintf('No smoothing\n')
        sm_normSIGNAL_1     = normSIGNAL_1;
        sm_normSIGNAL_2     = normSIGNAL_2;
end

% Update normSIGNAL
if and(p.USE_SMOOTH_1, exist('sm_normSIGNAL_1','var'))
    fprintf('Using smoothed signal 1\n')
    orig_normSIGNAL_1= normSIGNAL_1;
    normSIGNAL_1     = sm_normSIGNAL_1;
%         p.SIG_NAMES{1} = ['sm_', p.SIG_NAMES{1}];
else
    orig_normSIGNAL_1= normSIGNAL_1;
end

if and(p.USE_SMOOTH_2, exist('sm_normSIGNAL_2','var'))
    fprintf('Using smoothed signal 2\n')
    orig_normSIGNAL_2= normSIGNAL_2;
    normSIGNAL_2     = sm_normSIGNAL_2;
%         p.SIG_NAMES{2} = ['sm_', p.SIG_NAMES{2}];
else
    orig_normSIGNAL_2= normSIGNAL_2;
end



function [X,Y,nonNAN_IDX]  = removeNan(X,Y)
nonNAN_IDX  = (1:find(~isnan(X),1,'last'));
X           = X(nonNAN_IDX);
Y           = Y(nonNAN_IDX);


% --- Executes on button press in graph_compiled.
function graph_compiled_Callback(hObject, eventdata, handles)
% hObject    handle to graph_compiled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update all params
handles = update_params_Callback(hObject, eventData, handles);

% --- Subfunction to graph compiled data
errordlg('Not yet completed')
error('Not yet completed')

% Initialize variables
FOLDERS         = handles.FOLDERS;
SAVE_DIR        = handles.SAVE_DIR;
SIGNAL_1_all    = handles.SIGNAL_1_all;
normSIGNAL_1_all= handles.normSIGNAL_1_all;
SIGNAL_2_all    = handles.SIGNAL_2_all;
normSIGNAL_2_all= handles.normSIGNAL_2_all;
p               = handles.p;

try
    % Compile results
    um_px = (p.px_nm^-1);
    FNAME_1 = inputdlg('Name:','Name for data 1')
    [FILENAME_1, PATH_1] = getfile(p.EXT,['Select the compiled_Data.mat file ' ...
                        'to load'],'MultiSelect','off');
    load([PATH_1 FILE_1])
    FILE_1 = compiled_data;
    clear compiled_data

    FNAME_2 = inputdlg('Name:','Name for data 2')
    [FILENAME_2, PATH_2] = getfile(p.EXT,['Select the compiled_Data.mat file ' ...
                        'to load'],'MultiSelect','off');
    load([PATH_2 FILE_2])
    FILE_2 = compiled_data;
    clear compiled_data

    DIST_1          = num2cell(vertcat(FILE_1.Dist_orig),1);
    DIST_2          = num2cell(vertcat(FILE_2.Dist_orig),1);

    GLUT_DIST_2     = catData([DIST_1, DIST_2]);
    GLUT_DIST_2     = mergeArray(GLUT_DIST_2(:,3:4),GLUT_DIST_2(:, ...
                                                      1:2));

    boxplot(GLUT_DIST_2.*um_px,'symbol','r+')
    set(gca,'XTickLabel',{[FNAME_1 'nn1'],[FNAME_2 'nn1'],[FNAME_1 ...
                        'nn2'],[FNAME_2 '' nn2']});
    ylabel(p.XLABEL)
    xlabel('Nearest K PTM peaks to Synapsin')

    %[P,ABOVATAB,STATS]= kruskalwallis(GLUT_DIST_2);
    %[C,MU,H,~] = multcompare(STATS);

    DIST_2_ALL       = {vertcat(DIST_2{:}) ./ (p.px_um)};
    DIST_1_ALL      = {vertcat(DIST_1{:})./ (p.px_um)};

    boxplot(catData([DIST_2_ALL,DIST_1_ALL]))
    set(gca,'XTickLabel',{'Tyr','Glut'})
    ylabel('Distance (\mum)')
    xlabel(sprintf('Nearest 1-%d PTM peaks to Synapsin',p.KNN))

    % Plot signal overlap (raw signal overlap)
    GLUT_OVERLAP    = vertcat(GLUT.PTM_SYN_OVERLAP);
    TYR_OVERLAP     = vertcat(TYR.PTM_SYN_OVERLAP);

    boxplot(catData({TYR_OVERLAP,GLUT_OVERLAP}))
    set(gca,'XTickLabel',{'Tyr','Glut'})
    ylabel('Normalized intensity (AU)'); ylim([0,1]);
    xlabel('PTM Signal Density within FWHM of Synapsin Peak')
%     export_fig('synapsin_PTM_overlap_all_intensity.png')
    savefig('synapsin_PTM_overlap_all_intensity.fig')
    saveas(gcf, 'synapsin_PTM_overlap_all_intensity.png')
    close all

    % Plot integrated signal intensity
    stats_GLUT_OVERLAP  = horzcat(GLUT.stats_PTM_SYN_OVERLAP);
    stats_TYR_OVERLAP   = horzcat(TYR.stats_PTM_SYN_OVERLAP);

    GLUT_IntegratedInt  = stats_GLUT_OVERLAP(1,:); % Summed intensity
    TYR_IntegratedInt   = stats_TYR_OVERLAP(1,:); % Summed intensity

    boxplot(catData({TYR_IntegratedInt,GLUT_IntegratedInt}))
    set(gca,'XTickLabel',{'Tyr','Glut'})
    ylabel('Integrated signal intensity (AU)');
    xlabel('PTM Signal Density within FWHM of Synapsin Peak')
%     export_fig('synapsin_PTM_overlap_integrated_intensity.png')
    savefig('synapsin_PTM_overlap_integrated_intensity.fig')
    saveas(gcf, 'synapsin_PTM_overlap_integrated_intensity.png')
    close all

    % Plot mean signal intensity
    GLUT_MeanInt        = stats_GLUT_OVERLAP(6,:); % Mean intensity
    TYR_MeanInt         = stats_TYR_OVERLAP(6,:); % Mean intensity

    boxplot(catData({TYR_MeanInt,GLUT_MeanInt}))
    set(gca,'XTickLabel',{'Tyr','Glut'})
    ylabel('Mean signal intensity (AU)'); ylim([0,1]);
    xlabel('PTM Signal Density within FWHM of Synapsin Peak')
%     export_fig('synapsin_PTM_overlap_mean_intensity.png')
    savefig('synapsin_PTM_overlap_mean_intensity.fig')
    saveas(gcf, 'synapsin_PTM_overlap_mean_intensity.png')
    close all

    % Plot median signal intensity
    GLUT_MedianInt      = stats_GLUT_OVERLAP(7,:); % Median intensity
    TYR_MedianInt       = stats_TYR_OVERLAP(7,:); % Median intensity

    boxplot(catData({TYR_MedianInt,GLUT_MedianInt}))
    set(gca,'XTickLabel',{'Tyr','Glut'})
    ylabel('Median signal intensity (AU)'); ylim([0,1]);
    xlabel('PTM Signal Density within FWHM of Synapsin Peak')
%     export_fig('synapsin_PTM_overlap_median_intensity.png')
    savefig('synapsin_PTM_overlap_median_intensity.fig')
    saveas(gcf, 'synapsin_PTM_overlap_median_intensity.png')
    close all

    % Graph the randomization data
    color_vec           = distinguishable_colors(p.KNN);

    % GLUT
    GLUT_compiled        = num2cell((vertcat(GLUT.Dist_2)).*p.px_nm,1);
    GLUT_compiled_D_orig = num2cell((vertcat(GLUT.Dist_orig)).*p.px_nm,1);

    printText   = mergeArray(1:p.KNN,nanmean(cell2mat(GLUT_compiled(:,1:p.KNN)),1));
    legendText  = sprintf('knn-%d mu=%0.2f,',printText);
    legendText(end) = []; % Remove last comma
    nhist(GLUT_compiled(:,1:p.KNN),'legend',strsplit(legendText,','),...
            'box','xlabel','KNN Distance (microns)','color',color_vec);

    Xplot           = nanmean(cell2mat(GLUT_compiled_D_orig(1:p.KNN)),1);
    Xplot           = repmat(Xplot,p.KNN,1);
    Yplot           = repmat([0;max(ylim)],1,numel(Xplot)/2);
    cellfun(@(x,y,c) line('XData',x,'YData',y,'Color',c), num2cell(Xplot,1), num2cell(Yplot,1), num2cell(color_vec,p.KNN)','unif',false);
    %     export_fig([tempSAVE_DIR  filesep ['NN_Dist_Line_Scans_compiled_all' EXT]]); % Save figure


    % TYR
    TYR_compiled        = num2cell((vertcat(TYR.Dist_2)).*p.px_nm,1);
    TYR_compiled_D_orig = num2cell((vertcat(TYR.Dist_orig)).*p.px_nm,1);

    printText   = mergeArray(1:p.KNN,nanmean(cell2mat(TYR_compiled(:,1:p.KNN)),1));
    legendText  = sprintf('knn-%d mu=%0.2f,',printText);
    legendText(end) = []; % Remove last comma
    nhist(TYR_compiled(:,1:p.KNN),'legend',strsplit(legendText,','),...
            'box','xlabel','KNN Distance (microns)','color',color_vec);

    Xplot           = nanmean(cell2mat(TYR_compiled_D_orig(1:2)),1);
    Xplot           = repmat(Xplot,2,1);
    Yplot           = repmat([0;max(ylim)],1,numel(Xplot)/2);
    cellfun(@(x,y,c) line('XData',x,'YData',y,'Color',c), num2cell(Xplot,1), num2cell(Yplot,1), num2cell(color_vec,2)','unif',false);
    %     export_fig([tempSAVE_DIR  filesep ['NN_Dist_Line_Scans_compiled_all' EXT]]); % Save figure


    % TYR vs GLUT
    TYR_NN2_orig    = vertcat(TYR_compiled_D_orig{:});
    GLUT_NN2_orig   = vertcat(GLUT_compiled_D_orig{:});
    compiled_TYR_GLUT_NN2_orig = catData({TYR_NN2_orig,GLUT_NN2_orig});
    
    printText   = mergeArray(1:p.KNN,nanmean(compiled_TYR_GLUT_NN2_orig,1));
    legendText  = sprintf('knn1-2-%d mu=%0.2f,',printText);
    legendText(end) = []; % Remove last comma
    nhist({TYR_NN2_orig,GLUT_NN2_orig},'legend',strsplit(legendText,','),...
            'box','xlabel','KNN Distance (microns)','color',color_vec);

    Xplot           = nanmean(compiled_TYR_GLUT_NN2_orig,1);
    Xplot           = repmat(Xplot,2,1);
    Yplot           = repmat([0;max(ylim)],1,numel(Xplot)/2);
    cellfun(@(x,y,c) line('XData',x,'YData',y,'Color',c), num2cell(Xplot,1), num2cell(Yplot,1), num2cell(color_vec,2)','unif',false);
    %     export_fig([tempSAVE_DIR  filesep ['NN_Dist_Line_Scans_compiled_TYR_v_GLUT_orig_NN2' EXT]]); % Save figure

    %     figure, boxplot(compiled_TYR_GLUT_NN2_orig,Opt.Box{:},'symbol','r+')
    [P,ABOVATAB,STATS]= ranksum(compiled_TYR_GLUT_NN2_orig);
    [C,MU,H,~] = multcompare(STATS);
    
    close(1);
catch ERR
    SUCCESS         = false;
    handles.SUCCESS = SUCCESS;
    handles.ERR     = ERR;
    error_Callback(hObject,[],handles)
end

% Allocate output variables
handles.SUCCESS     = SUCCESS;
handles.ERR         = ERR;
handles.p           = p;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in export_excel_box.
function export_excel_box_Callback(hObject, eventdata, handles)
% hObject    handle to export_excel_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of export_excel_box
handles.p.EXCEL         = logical(get(handles.export_excel_box,'Value'));

% Update all params
update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in plot_1.
function plot_1_Callback(hObject, eventdata, handles)
% hObject    handle to plot_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_1
handles.p.PLOT_1        = logical(get(handles.plot_1,'Value'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in peak_rand.
function peak_rand_Callback(hObject, eventdata, handles)
% hObject    handle to peak_rand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of peak_rand
handles.p.PEAK_RAND     = logical(get(handles.peak_rand,'Value'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in use_smooth_1.
function use_smooth_1_Callback(hObject, eventdata, handles)
% hObject    handle to use_smooth_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_smooth_1
handles.p.SMOOTH_1      = logical(get(handles.use_smooth_1,'Value'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in use_smooth_2.
function use_smooth_2_Callback(hObject, eventdata, handles)
% hObject    handle to use_smooth_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_smooth_2
handles.p.SMOOTH_2      = logical(get(handles.use_smooth_2,'Value'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rand_both.
function rand_both_Callback(hObject, eventdata, handles)
% hObject    handle to rand_both (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rand_both
handles.p.PEAK_RAND      = logical(get(handles.rand_both,'Value'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in smooth_menu.
function smooth_menu_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns smooth_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from smooth_menu
menu_contents           = cellstr(get(handles.smooth_menu,'String'));
handles.p.PEAK_RAND     = menu_contents{get(handles.smooth_menu,'Value')};

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function smooth_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);


function win_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to win_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of win_size_edit as text
%        str2double(get(hObject,'String')) returns contents of win_size_edit as a double
handles.p.WIN_SIZE      = str2double(get(handles.win_size_edit,'String'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function win_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to win_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in peak_method_menu.
function peak_method_menu_Callback(hObject, eventdata, handles)
% hObject    handle to peak_method_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns peak_method_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from peak_method_menu
menu_contents           = cellstr(get(handles.peak_method_menu,'String'));
handles.p.PEAK_METHOD   = menu_contents{get(handles.peak_method_menu,'Value')};

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function peak_method_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peak_method_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);


function sig_names_1_Callback(hObject, eventdata, handles)
% hObject    handle to sig_names_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sig_names_1 as text
%        str2double(get(hObject,'String')) returns contents of sig_names_1 as a double
SIG_1               = get(handles.sig_names_1,'String');
SIG_2               = get(handles.sig_names_2,'String');

handles.p.SIG_NAMES = {SIG_1,SIG_2};

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function sig_names_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sig_names_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);



function sig_names_2_Callback(hObject, eventdata, handles)
% hObject    handle to sig_names_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sig_names_2 as text
%        str2double(get(hObject,'String')) returns contents of sig_names_2 as a double
SIG_1                   = get(handles.sig_names_1,'String');
SIG_2                   = get(handles.sig_names_2,'String');

handles.p.SIG_NAMES     = {SIG_1,SIG_2};

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sig_names_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sig_names_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in output_format_menu.
function output_format_menu_Callback(hObject, eventdata, handles)
% hObject    handle to output_format_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns output_format_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from output_format_menu
menu_contents           = cellstr(get(handles.output_format_menu,'String'));
handles.p.EXT           = menu_contents{get(handles.output_format_menu,'Value')};

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function output_format_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_format_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);



function nsample_Callback(hObject, eventdata, handles)
% hObject    handle to nsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nsample as text
%        str2double(get(hObject,'String')) returns contents of nsample as a double
handles.p.NSAMPLES      = str2double(get(handles.nsample,'String'));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nsample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);



function nm_px_Callback(hObject, eventdata, handles)
% hObject    handle to nm_px (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nm_px as text
%        str2double(get(hObject,'String')) returns contents of nm_px as a double

handles.p.px_nm      = (str2double(get(handles.nm_px,'String')))^-1;

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nm_px_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nm_px (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Update handles structure
guidata(hObject, handles);



function knn_edit_Callback(hObject, eventdata, handles)
% hObject    handle to knn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of knn_edit as text
%        str2double(get(hObject,'String')) returns contents of knn_edit as a double

handles.p.KNN      = (str2double(get(handles.knn_edit,'String')));

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function knn_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to knn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EB_ANALYSIS_init_term.
function EB_ANALYSIS_init_term_Callback(hObject, eventdata, handles)
% hObject    handle to EB_ANALYSIS_init_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EB_ANALYSIS_init_term

if logical(get(handles.EB_ANALYSIS_init_term,'Value'));
    handles.p.EB_ANALYSIS   = 'antero';
else
    handles.p.EB_ANALYSIS   = 'retro';
end

% Update all params
handles = update_params_Callback(hObject, [], handles);

% Update handles structure
guidata(hObject, handles);