% Kymograph Manual Vesicle Identification function_version 2
%  
% def kymo_vesicle:
% ''' (Tiff image) -> <Nrow x 5 cell>
% A function that returns a Nrow by X Ncol cell containing the file information
%  and the parameters xy coordinates of Ncol kymographs of Nrow number
% of neurons  for analysis in  kymo_analysis.m.
%
% The structure of the output file is as follows:
% Each row is 1 neuron. Columns define information associated with that
% neuron. Date/time, researcher, experiment name/number, and other comments 
% are stored in the first cell. Each remaining cell for Nkymo number of
% columns loads 1 kymograph trace associated with that neuron and stores all of 
% the xy coordinates for all vesicles drawn on each kymograph image.
% 
% Author: Jeffrey J. Nirschl, Holzbaur Lab (University of Pennsylvania)
% Date created: 09/03/2013
% Updated 09/16/2013 to include Excel export and to save reference image of
% kymographs with traces
% Updated 09/25/2013 to add 10um reference bar "view kymograph"
% Updated 06/09/2014 to fix GUI bugs
% Updated 11/21/2016 to add new features--iterative plotting with getUserXY
% Updated 05/28/2018 
% Distributable under BSD liscence. 
% Copyright (c) 2016 Jeffrey J. Nirschl
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


function varargout = kymo_tracer(varargin)
% KYMO_TRACER MATLAB code for kymo_tracer.fig
%      KYMO_TRACER, by itself, creates a new KYMO_TRACER or raises the existing
%      singleton*.
%
%      H = KYMO_TRACER returns the handle to a new KYMO_TRACER or the handle to
%      the existing singleton*.
%
%      KYMO_TRACER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KYMO_TRACER.M with the given input arguments.
%
%      KYMO_TRACER('Property','Value',...) creates a new KYMO_TRACER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kymo_tracer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kymo_tracer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kymo_tracer

% Last Modified by GUIDE v2.5 19-Jan-2014 17:48:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kymo_tracer_OpeningFcn, ...
                   'gui_OutputFcn',  @kymo_tracer_OutputFcn, ...
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


% --- Executes just before kymo_tracer is made visible.
function kymo_tracer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kymo_tracer (see VARARGIN)

% Choose default command line output for kymo_tracer
handles.output = hObject;

validateCMAP    = @(x) any(strcmpi(x,{'gray','red hot', 'red-hot', 'red',...
                    'gfb','green-fire-blue','green fire blue',...
                    'fire'}));
                
% set home dir
if ispc
    userDir = winqueryreg('HKEY_CURRENT_USER',...
        ['Software\Microsoft\Windows\CurrentVersion\' ...
         'Explorer\Shell Folders'],'Personal');
else
    userDir = char(java.lang.System.getProperty('user.home'));
end

% Load data from varargin, if available(REVIEW)
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('Data',cell(0,1),@iscell); 
ip.addParameter('Date',datestr(now, 'mm-DD-YYYY'),@ischar);
ip.addParameter('MainDir',userDir,@isdir);
ip.addParameter('SaveDir',userDir,@isdir);
ip.addParameter('SaveDir_Default',pwd,@isdir);
ip.addParameter('Colormap','red hot',@ischar);
ip.parse(varargin{:});
p   = ip.Results;

if ~isempty(p.Data)
    MainData        = p;
    handles.Date    =MainData.Date;
    handles.MainDir = MainData.MainDir;
    handles.SaveDir = MainData.SaveDir;
    handles.SaveDir_Default = MainData.SaveDir_Default;
    handles.output  = p.Data;
    handles.neuroncounter= size(p.Data,1);
    handles.kymocounter = max([find(cellfun(@isempty, p.Data(end,:), 'unif',true),...
                                1, 'first')-1,1]);
else
    % Set date/dir
    handles.Date = p.Date;
    handles.MainDir = p.MainDir;
    handles.SaveDir = p.SaveDir;
    handles.SaveDir_Default = p.SaveDir_Default ;
    handles.output      = cell(0,1);
    handles.kymocounter = 1;
    handles.neuroncounter= 1;
end

% Initialize variables
set(handles.pixel_edit,'String','0.0715'); %Input the um/pixel
set(handles.frame_edit,'String','0.5'); %Input the time (secs) per frame
set(handles.experiment_details,'String','GD001_XR001- Temp');
% set(handles.photobleach_checkbox,'Value',0);
handles.files       = cell(0,1);
handles.filespath   = cell(0,1);
handles.excel       = cell(0,1);
handles.pixel_conv  = str2num(get(handles.pixel_edit,'String'));
handles.SaveDir     = uigetdir(userDir, 'Select save directory.');
handles.Colormap    = p.Colormap;
handles.cmap        = getLUT(handles.Colormap);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kymo_tracer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kymo_tracer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in filelist.
function filelist_Callback(hObject, eventdata, handles)
% hObject    handle to filelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filelist
% contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filelist

handles.files = get(hObject,'String');
handles.Selectedfiles=get(hObject,'Value');
cd(handles.MainDir);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function filelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in photobleach_checkbox.
function photobleach_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to photobleach_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of photobleach_checkbox

handles.checkbox_value = get(handles.photobleach_checkbox,'Value');
guidata(hObject,handles)

function pixel_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_edit as text
%        str2double(get(hObject,'String')) returns contents of pixel_edit as a double

% handles.pixel_edit = get(hObject,'String');
handles.pixel_conv = str2num(get(hObject,'String'));
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function pixel_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_edit as text
%        str2double(get(hObject,'String')) returns contents of frame_edit as a double

handles.frame_conv = str2num(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function frame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_checkbox.
function plot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to plot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_checkbox


% --- Executes on button press in loadtiff.
function loadtiff_Callback(hObject, eventdata, handles)
% hObject    handle to loadtiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TIFF load
cd(handles.MainDir);
% handles.files = get(handles.filelist,'String');
% handles.Selectedfiles=get(handles.filelist,'Value');
[FILENAME, PATHNAME]=uigetfile({'*.tif';'*'},'MultiSelect','on');
handles.filedir = PATHNAME;

if(PATHNAME~=0)
    
    if(ischar(FILENAME))
        FNAME=cell(1);
        FNAME{1}=FILENAME;
        clear(FILENAME);
        FILENAME=FNAME;
    end
    PATHCELL=cell(size(FILENAME));
    for i=1:size(FILENAME,2)
        PATHCELL{i}=PATHNAME;
    end
    handles.files(size(handles.files,1)+1:size(handles.files,1)+size(FILENAME,2),1)=FILENAME';
    handles.filespath(size(handles.filespath,1)+1:size(handles.filespath,1)+size(FILENAME,2),1)=PATHCELL';
    set(handles.filelist,'String',handles.files);
    cd(handles.MainDir);
end

handles.files = get(handles.filelist,'String');
handles.Selectedfiles=get(handles.filelist,'Value');

cd(handles.MainDir);

guidata(hObject,handles)


% --- Executes on button press in new_tracing.
function new_tracing_Callback(hObject, eventdata, handles)
% hObject    handle to new_tracing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.checkbox_value = get(handles.photobleach_checkbox,'Value');
handles.flist = get(handles.filelist,'String');
handles.Selectedfiles=get(handles.filelist,'Value');
filename = handles.flist(handles.Selectedfiles);
width = handles.info.Width;
height = handles.info.Height;
fail = 0;

if handles.kymocounter == 1
    handles.output{handles.neuroncounter,1} = {[get(handles.experiment_details,'String')]...
    ['um/ pixel= ' get(handles.pixel_edit,'String')] ['s/ frame= ' get(handles.frame_edit,'String')];... 
    [filename] [get(handles.neuron_text,'String')] [get(handles.genotype_text,'String')...
    '_' get(handles.condition_text,'String')]; ['x pixels = ' num2str(width)] ['y pixels = ' num2str(height)]...
    []};
end

% Set new figure to 99
tempIM  = getimage(handles.kymograph);
[xCols, yRows]      = size(tempIM);
H_NEWFIG            = figure(99);
set(H_NEWFIG, 'Units','Pixels','Position',[10, 10, yRows+10, xCols+10]);
AxesH               = axes('Units','Pixels',...
                        'Position',[10, 10, yRows+10, xCols+10]);

hIM                 = imagesc(tempIM); % Create new figure w/ image
H_SCROLL_PANEL      = imscrollpanel(H_NEWFIG, hIM); % Create new figure w/ image
handles.temp_api    = iptgetapi(H_SCROLL_PANEL);
% set(H_NEWFIG,'WindowScrollWheelFcn',{@mouse_scroll_Callback handles.temp_api})
axis image, axis off
grid off

colormap(handles.cmap)
ctrace  = getUserXY('FIGURE',H_NEWFIG,'CMAP',handles.cmap);
close(H_NEWFIG)
linewid = 2;

if fail==0;
    handles.kymocounter = handles.kymocounter + 1;
    if handles.checkbox_value ==1; % Photobleach kymographs are drawn differently--starting
        if sum([0; diff(ctrace(:,1).*handles.pixel_conv)]) >= 10; % Retrograde is positive for FRAP assay
            color =  '-rx'; % Red color for retrograde
        elseif sum([0; diff(ctrace(:,1).*handles.pixel_conv)]) <= -10 % Anterograde is negative for FRAP assay
            color = '-bx'; % Blue/ Azul for anterograde
        else
            color = '-gx';
        end
    else
        if sum([0; diff(ctrace(:,1).*handles.pixel_conv)]) >= 10; % Anterograde is positive
            color =  '-bx'; % Blue/ Azul for Anterograde
        elseif sum([0; diff(ctrace(:,1).*handles.pixel_conv)]) <= -10 % Retrograde is negative
            color = '-rx'; % Red color for retrograde
        else
            color = '-gx';
        end
    end

    figure(handles.kymograph);
    axis off;
    hold on;
    if ~isempty(ctrace)
        plot(ctrace(:,1),ctrace(:,2),...
            color,'LineWidth',linewid,'MarkerSize',8);
        text(ctrace(1,1),ctrace(1,2),['\color{red}' num2str(handles.kymocounter-1)],...
            'VerticalAlignment','top',...
            'HorizontalAlignment','left',...
            'FontSize',14);
        handles.output{handles.neuroncounter,handles.kymocounter}   = [ctrace];
        handles.excel{handles.neuroncounter,handles.kymocounter}    = {ctrace};
    else
        errordlg('Track is empty!')
        handles.kymocounter = handles.kymocounter - 1;
    end
end

% Save .Mat and .bkp files
data = handles.output;
fname = [get(handles.experiment_details,'String') '_TEMPSAVE_kymo-tracer.mat'];
save(fullfile(handles.SaveDir, fname),'-mat','data');
save(fullfile(handles.SaveDir, strrep(fname, '.mat', '.bkp')),'-mat','data');


guidata(hObject,handles)


% --- Executes on button press in markneuron.
function markneuron_Callback(hObject, eventdata, handles)
% hObject    handle to markneuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


number = get(handles.neuron_text,'String');
neuroncount = str2num(number(strfind(number,' ')+1:end));

if (neuroncount == handles.neuroncounter)==1;
    handles.neuroncounter = handles.neuroncounter + 1;
else
    handles.neuroncounter = neuroncount;
end
set(handles.neuron_text,'String',['Neuron ' num2str(handles.neuroncounter)]);
handles.kymocounter = 1;

% Save .Mat and .bkp files
data = handles.output;
fname = [get(handles.experiment_details,'String') '_TEMPSAVE_kymo-tracer.mat'];
save(fullfile(handles.SaveDir, fname),'-mat','data');
save(fullfile(handles.SaveDir, strrep(fname, '.mat', '.bkp')),'-mat','data');

guidata(hObject,handles)



% --- Executes on button press in viewkymo.
function viewkymo_Callback(hObject, eventdata, handles)
% hObject    handle to viewkymo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get filename and path
handles.checkbox_value = get(handles.photobleach_checkbox,'Value');
handles.files = get(handles.filelist,'String');
handles.Selectedfiles=get(handles.filelist,'Value');

a= handles.Selectedfiles;
file= handles.files{a};
fname= file;
dirnew= handles.filespath{a};
cd(dirnew);

% Load image, convert to grayscale, and show image in current axes
try
    if handles.checkbox_value ==1;
        im = imread([dirnew file]); 
        im = imrotate(im,-90);
        handles.info = imfinfo([dirnew file]);

        if strcmpi(handles.info(1,1).ColorType, 'grayscale');
            im = mat2gray(im);
        elseif strcmpi(handles.info(1,1).ColorType, 'truecolor');
            im = rgb2gray(im); 
        end
        
        [xCols, yRows]      = size(im);
        handles.kymograph   = figure('Units','Pixels',...
                                'Position',[10, 10, yRows+10, xCols+10]);

        AxesH               = axes('Units','Pixels',...
                                'Position',[10, 10, yRows+10, xCols+10]);
        
        hIM                 = imagesc(im, 'Parent',AxesH); % Create new figure w/ image
        H_SCROLL_PANEL      = imscrollpanel(handles.kymograph, hIM); % Create new figure w/ image
        handles.main_api    = iptgetapi(H_SCROLL_PANEL);
        axis image, axis off
        grid off

        colormap(handles.cmap)
        set(handles.kymograph,'Name',fname,'NumberTitle','off')
        title('\leftarrow Anterograde                  Retrograde  \rightarrow','Color','k','FontSize',18)
        line([50 50],[0 size(im,1)],'Color',[1 1 0]);
        
        pixel_edit = (str2num(get(handles.pixel_edit,'String')));
        handles.imdist = imdistline(gca,[1,round(10/pixel_edit)+1],[50 50]);
        
        api = handles.imdist;
        api.setLabelTextFormatter(sprintf('%02.0f um',10));
        
        tempFIG     = imcontrast(handles.kymograph); % Adjust contrast of image
        
    else
        im = imread([dirnew file]); 
        handles.info = imfinfo([dirnew file]);

        if strcmpi(handles.info(1,1).ColorType, 'grayscale');
            im = mat2gray(im);
        elseif strcmpi(handles.info(1,1).ColorType, 'truecolor');
            im = rgb2gray(im); 
        end
        
        [xCols, yRows]      = size(im);
        handles.kymograph   = figure('Units','Pixels',...
                                'Position',[10, 10, yRows+10, xCols+10]);

        AxesH               = axes('Units','Pixels',...
                                'Position',[10, 10, yRows+10, xCols+10]);
        
        hIM                 = imagesc(im, 'Parent',AxesH); % Create new figure w/ image
        H_SCROLL_PANEL      = imscrollpanel(handles.kymograph, hIM); % Create new figure w/ image
        handles.main_api    = iptgetapi(H_SCROLL_PANEL);
%         set(handles.kymograph,'WindowScrollWheelFcn',{@mouse_scroll_Callback handles.main_api})
        axis image, axis off
        grid off
        
        
        colormap(handles.cmap)
        set(handles.kymograph,'Name',fname,'NumberTitle','off')
        title('\leftarrow Retrograde                  Anterograde  \rightarrow','Color','k','FontSize',18)
        
        pixel_edit = (str2num(get(handles.pixel_edit,'String')));
        handles.imdist = imdistline(gca,[1,round(10/pixel_edit)+1],[50 50]);
        api = handles.imdist;
        api.setLabelTextFormatter(sprintf('%02.0f um',10));
        
        tempFIG     = imcontrast(handles.kymograph); % Adjust contrast of image 
    end
catch err
    msgbox('Could not load image. Re-select an image and try again.')
    rethrow(err)
end

% Update text in figure
set(handles.filename,'String',[dirnew fname]); % Update filename displayed
guidata(hObject,handles)


% --- Executes on button press in save_refim.
function save_refim_Callback(hObject, eventdata, handles)
% hObject    handle to save_refim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fname = get(handles.kymograph,'Name');
saveas(handles.kymograph, fullfile(handles.SaveDir, [fname(1:end-4) '_kymo-tracer.tif']),'tif'); % Saves Matlabfigure with vesicle tracks as tif
imwrite(getimage(handles.kymograph), fullfile(handles.SaveDir, [fname(1:end-4) '_contrast.tif']),'tif'); % Saves contrast-adjusted reference image

% Save .Mat and .bkp files
data = handles.output;
fname = [get(handles.experiment_details,'String') '_TEMPSAVE_kymo-tracer.mat'];
save(fullfile(handles.SaveDir, fname),'-mat','data');
save(fullfile(handles.SaveDir, strrep(fname, '.mat', '.bkp')),'-mat','data');

guidata(hObject,handles)

% --- Executes on button press in save_close.
function save_close_Callback(hObject, eventdata, handles)
% hObject    handle to save_refim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save .Mat and .bkp files
data = handles.output;
fname = [get(handles.experiment_details,'String') '_kymo-tracer.mat'];
save(fullfile(handles.SaveDir, fname),'-mat','data');
save(fullfile(handles.SaveDir, strrep(fname, '.mat', '.bkp')),'-mat','data');


clear all; close all

% --- Executes on button press in excel_export.
% Pedro export
function excel_export_Callback(hObject, eventdata, handles)
% hObject    handle to excel_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.files = get(handles.filelist,'String');
handles.Selectedfiles=get(handles.filelist,'Value');
% handles.experiment_details = get(handles.experiment_details,'String');
try
    excelname = handles.output{1,1}{1,1};
catch
    if isempty(handles.output{2,1})==0;
        excelname = handles.output{2,1}{1,1};
    elseif isempty(handles.output{3,1})==0;
        excelname = handles.output{3,1}{1,1};
    elseif isempty(handles.output{4,1})==0;
        excelname = handles.output{4,1}{1,1};
    elseif isempty(handles.output{5,1})==0;
        excelname = handles.output{5,1}{1,1};
    elseif isempty(handles.output{6,1})==0;
        excelname = handles.output{6,1}{1,1};
    elseif isempty(handles.output{7,1})==0;
        excelname = handles.output{7,1}{1,1};
    elseif isempty(handles.output{8,1})==0;
        excelname = handles.output{8,1}{1,1};
    elseif isempty(handles.output{9,1})==0;
        excelname = handles.output{9,1}{1,1};
    else
        excelname = handles.output{10,1}{1,1};
    end
end

handles.SaveDir = uigetdir(handles.SaveDir, 'Select a folder to save Excel file.');


% Identify currently selected file
cd(handles.MainDir);
try
    file=handles.files{handles.Selectedfiles};
    fname=file;
    dirnew=handles.filespath{handles.Selectedfiles};
catch
    file=handles.output{1,1}{1,1};
    dirnew= handles.SaveDir;
end


% Save .Mat and .bkp files
data = handles.output;
fname = [get(handles.experiment_details,'String') '_kymo-tracer.mat'];
save(fullfile(handles.SaveDir, fname),'-mat','data');
save(fullfile(handles.SaveDir, strrep(fname, '.mat', '.bkp')),'-mat','data');

try
    tempSaveDir = userpath;
    tempSaveDir = tempSaveDir(1:end-1);
    if ~(exist([tempSaveDir filesep 'kymo_suite'],'dir')==7)
        mkdir(tempSaveDir, 'kymo_suite')
    end
    tempSaveDir = [tempSaveDir filesep 'kymo_suite'];
    save([tempSaveDir filesep fname '_kymo-analysis.mat'],'data');
catch
end


% Initiate Excel ActX server % Format for using---xlswrite2007(file,data,sheet,range)
Excel = actxserver('Excel.Application');
File = [handles.SaveDir '\' excelname '_kymo-tracer.xlsx'];  %Make sure you put the whole path
if ~exist(File,'file')
    ExcelWorkbook = Excel.workbooks.Add;
    ExcelWorkbook.SaveAs(File);
    ExcelWorkbook.Close(false);
end
ExcelWorkbook = Excel.workbooks.Open(File);
range=Excel.ActiveSheet.UsedRange.Address; % Finds used range of active excel sheet

% Write info to Excel
warning off
try    
    File=[handles.SaveDir '\' excelname '_kymo-tracer.xlsx'];
    sz =  size(handles.output);
     for i=1:sz(1);
         for j=1:sz(2);
              if isempty(handles.output{i,j}) == 0;
                  if j ==1;
                      idx = 4;
                      xlswrite2007(File,handles.output{i,j},handles.output{i,1}{2,2},'A1');
                      xlswrite2007(File,['x' 'y'],handles.output{i,1}{2,2},['B' num2str(idx)]);
                      idx = idx+1;
                  else 
                      xlswrite2007(File,handles.output{i,j},handles.output{i,1}{2,2},['B' num2str(idx)]); % Write XY coordinates
                      sz2 = size(handles.output{i,j});
                      idx = idx + sz2(1) +1;
                  end
              end
         end
     end
    
    msgbox(sprintf('Saved %s to %s',fname,[handles.SaveDir '\']));
    fprintf(sprintf('Saved %s to \n',fname));
    disp([handles.SaveDir '\']);
    
catch err2
    msgbox(sprintf('Failed to write to excel to %s\%s.',handles.SaveDir,fname), 'Error')
    fail=1;
    rethrow(err2)
end
warning on

try % Delete the empty sheets 1, 2, 3 in the newly created Excel file
    Excel.ActiveWorkbook.Worksheets.Item('Sheet1').Delete;
    Excel.ActiveWorkbook.Worksheets.Item('Sheet2').Delete;
    Excel.ActiveWorkbook.Worksheets.Item('Sheet3').Delete;
catch
    fprintf('There are no excel sheets to delete.\n');
end

    % Close Excel ActX server
ExcelWorkbook.Save
ExcelWorkbook.Close(false)  % Close Excel workbook.
Excel.Quit;
delete(Excel);
    

guidata(hObject,handles);


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


function LUT    = getLUT(STR)
% Subfunction to get LUT

switch lower(STR)
    case {'gfb','green-fire-blue','green fire blue'}
        LUT     = [0,0,0;0,1,2;0,2,5;0,3,7;0,5,10;0,6,13;0,7,15;0,9,18;0,10,21;0,11,23;0,13,26;0,14,29;0,15,31;0,17,34;0,18,37;0,19,39;0,21,42;0,22,45;0,23,47;0,25,50;0,26,53;0,27,55;0,29,58;0,30,61;0,31,63;0,33,66;0,34,69;0,35,71;0,37,74;0,38,77;0,39,79;0,41,82;0,42,85;0,43,87;0,45,90;0,46,92;0,47,95;0,49,98;0,50,100;0,51,103;0,53,106;0,54,108;0,55,111;0,57,114;0,58,116;0,59,119;0,61,122;0,62,124;0,63,127;0,65,130;0,66,132;0,67,135;0,69,138;0,70,140;0,71,143;0,73,146;0,74,148;0,75,151;0,77,154;0,78,156;0,79,159;0,81,162;0,82,164;0,85,170;0,85,170;0,87,167;0,90,164;0,92,162;0,95,159;0,98,156;0,100,154;0,103,151;0,106,148;0,108,146;0,111,143;0,114,140;0,116,138;0,119,135;0,122,132;0,124,130;0,127,127;0,130,124;0,132,122;0,135,119;0,138,116;0,140,114;0,143,111;0,146,108;0,148,106;0,151,103;0,154,100;0,156,98;0,159,95;0,162,92;0,164,90;0,167,87;0,170,85;0,172,82;0,175,79;0,177,77;0,180,74;0,183,71;0,185,69;0,188,66;0,191,63;0,193,61;0,196,58;0,199,55;0,201,53;0,204,50;0,207,47;0,209,45;0,212,42;0,215,39;0,217,37;0,220,34;0,223,31;0,225,29;0,228,26;0,231,23;0,233,21;0,236,18;0,239,15;0,241,13;0,244,10;0,247,7;0,249,5;0,255,0;0,255,0;3,255,0;7,255,0;11,255,0;15,255,0;19,255,0;23,255,0;27,255,0;31,255,0;35,255,0;39,255,0;43,255,0;47,255,0;51,255,0;55,255,0;59,255,0;63,255,0;67,255,0;71,255,0;75,255,0;79,255,0;83,255,0;87,255,0;91,255,0;95,255,0;99,255,0;103,255,0;107,255,0;111,255,0;115,255,0;119,255,0;123,255,0;127,255,0;131,255,0;135,255,0;139,255,0;143,255,0;147,255,0;151,255,0;155,255,0;159,255,0;163,255,0;167,255,0;171,255,0;175,255,0;179,255,0;183,255,0;187,255,0;191,255,0;195,255,0;199,255,0;203,255,0;207,255,0;211,255,0;215,255,0;219,255,0;223,255,0;227,255,0;231,255,0;235,255,0;239,255,0;243,255,0;247,255,0;255,255,0;255,255,0;255,255,3;255,255,7;255,255,11;255,255,15;255,255,19;255,255,23;255,255,27;255,255,31;255,255,35;255,255,39;255,255,43;255,255,47;255,255,51;255,255,55;255,255,59;255,255,63;255,255,67;255,255,71;255,255,75;255,255,79;255,255,83;255,255,87;255,255,91;255,255,95;255,255,99;255,255,103;255,255,107;255,255,111;255,255,115;255,255,119;255,255,123;255,255,127;255,255,131;255,255,135;255,255,139;255,255,143;255,255,147;255,255,151;255,255,155;255,255,159;255,255,163;255,255,167;255,255,171;255,255,175;255,255,179;255,255,183;255,255,187;255,255,191;255,255,195;255,255,199;255,255,203;255,255,207;255,255,211;255,255,215;255,255,219;255,255,223;255,255,227;255,255,231;255,255,235;255,255,239;255,255,243;255,255,247;255,255,255]./255;
    case {'red hot', 'red-hot', 'red'}
        LUT     = [0,0,0;0,0,0;0,0,0;0,0,0;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;4,0,0;7,0,0;10,0,0;13,0,0;16,0,0;20,0,0;23,0,0;27,0,0;30,0,0;33,0,0;36,0,0;39,0,0;42,0,0;45,0,0;48,0,0;52,0,0;55,0,0;58,0,0;61,0,0;65,0,0;68,0,0;71,0,0;74,0,0;78,0,0;81,0,0;84,0,0;87,0,0;90,0,0;93,0,0;96,0,0;99,0,0;103,0,0;106,0,0;110,0,0;113,0,0;117,0,0;120,0,0;123,0,0;126,0,0;130,0,0;133,0,0;136,0,0;139,0,0;142,0,0;145,0,0;148,0,0;151,0,0;155,0,0;158,0,0;161,0,0;164,0,0;168,0,0;171,0,0;174,0,0;177,0,0;181,0,0;184,0,0;187,0,0;190,0,0;193,0,0;196,0,0;200,0,0;203,0,0;207,0,0;210,0,0;213,0,0;216,0,0;220,0,0;223,0,0;226,0,0;229,0,0;233,0,0;236,0,0;239,0,0;242,0,0;245,0,0;247,0,0;250,1,0;252,2,0;255,3,0;255,6,0;255,9,0;255,12,0;255,16,0;255,19,0;255,22,0;255,25,0;255,29,0;255,32,0;255,35,0;255,38,0;255,42,0;255,45,0;255,48,0;255,51,0;255,55,0;255,58,0;255,61,0;255,64,0;255,68,0;255,71,0;255,74,0;255,77,0;255,81,0;255,84,0;255,87,0;255,90,0;255,94,0;255,97,0;255,100,0;255,103,0;255,106,0;255,109,0;255,112,0;255,115,0;255,119,0;255,122,0;255,126,0;255,129,0;255,133,0;255,136,0;255,139,0;255,142,0;255,145,0;255,148,0;255,151,0;255,154,0;255,158,0;255,161,0;255,164,0;255,167,0;255,171,0;255,174,0;255,177,0;255,180,0;255,184,0;255,187,0;255,190,0;255,193,0;255,197,0;255,200,0;255,203,0;255,206,0;255,209,0;255,212,0;255,216,0;255,219,0;255,223,0;255,226,0;255,229,0;255,232,0;255,236,0;255,239,0;255,242,0;255,245,0;255,249,0;255,250,1;255,252,3;255,253,4;255,255,6;255,255,9;255,255,12;255,255,15;255,255,19;255,255,22;255,255,25;255,255,28;255,255,32;255,255,35;255,255,38;255,255,41;255,255,45;255,255,48;255,255,51;255,255,54;255,255,58;255,255,61;255,255,64;255,255,67;255,255,71;255,255,74;255,255,77;255,255,80;255,255,84;255,255,87;255,255,90;255,255,93;255,255,97;255,255,100;255,255,103;255,255,106;255,255,110;255,255,113;255,255,116;255,255,119;255,255,122;255,255,125;255,255,128;255,255,131;255,255,135;255,255,138;255,255,142;255,255,145;255,255,149;255,255,152;255,255,155;255,255,158;255,255,161;255,255,164;255,255,167;255,255,170;255,255,174;255,255,177;255,255,180;255,255,183;255,255,187;255,255,190;255,255,193;255,255,196;255,255,200;255,255,203;255,255,206;255,255,209;255,255,213;255,255,216;255,255,219;255,255,222;255,255,225;255,255,228;255,255,232;255,255,235;255,255,239;255,255,242;255,255,245;255,255,248;255,255,252;255,255,252;255,255,253;255,255,254;255,255,255;255,255,255;255,255,255;255,255,255;255,255,255;255,255,255;255,255,255;255,255,255]./255;
    case {'fire'}
        LUT     = [0,0,0;0,0,7;0,0,15;0,0,22;0,0,30;0,0,38;0,0,45;0,0,53;0,0,61;0,0,65;0,0,69;0,0,74;0,0,78;0,0,82;0,0,87;0,0,91;1,0,96;4,0,100;7,0,104;10,0,108;13,0,113;16,0,117;19,0,121;22,0,125;25,0,130;28,0,134;31,0,138;34,0,143;37,0,147;40,0,151;43,0,156;46,0,160;49,0,165;52,0,168;55,0,171;58,0,175;61,0,178;64,0,181;67,0,185;70,0,188;73,0,192;76,0,195;79,0,199;82,0,202;85,0,206;88,0,209;91,0,213;94,0,216;98,0,220;101,0,220;104,0,221;107,0,222;110,0,223;113,0,224;116,0,225;119,0,226;122,0,227;125,0,224;128,0,222;131,0,220;134,0,218;137,0,216;140,0,214;143,0,212;146,0,210;148,0,206;150,0,202;152,0,199;154,0,195;156,0,191;158,0,188;160,0,184;162,0,181;163,0,177;164,0,173;166,0,169;167,0,166;168,0,162;170,0,158;171,0,154;173,0,151;174,0,147;175,0,143;177,0,140;178,0,136;179,0,132;181,0,129;182,0,125;184,0,122;185,0,118;186,0,114;188,0,111;189,0,107;190,0,103;192,0,100;193,0,96;195,0,93;196,1,89;198,3,85;199,5,82;201,7,78;202,8,74;204,10,71;205,12,67;207,14,64;208,16,60;209,19,56;210,21,53;212,24,49;213,27,45;214,29,42;215,32,38;217,35,35;218,37,31;220,40,27;221,43,23;223,46,20;224,48,16;226,51,12;227,54,8;229,57,5;230,59,4;231,62,3;233,65,3;234,68,2;235,70,1;237,73,1;238,76,0;240,79,0;241,81,0;243,84,0;244,87,0;246,90,0;247,92,0;249,95,0;250,98,0;252,101,0;252,103,0;252,105,0;253,107,0;253,109,0;253,111,0;254,113,0;254,115,0;255,117,0;255,119,0;255,121,0;255,123,0;255,125,0;255,127,0;255,129,0;255,131,0;255,133,0;255,134,0;255,136,0;255,138,0;255,140,0;255,141,0;255,143,0;255,145,0;255,147,0;255,148,0;255,150,0;255,152,0;255,154,0;255,155,0;255,157,0;255,159,0;255,161,0;255,162,0;255,164,0;255,166,0;255,168,0;255,169,0;255,171,0;255,173,0;255,175,0;255,176,0;255,178,0;255,180,0;255,182,0;255,184,0;255,186,0;255,188,0;255,190,0;255,191,0;255,193,0;255,195,0;255,197,0;255,199,0;255,201,0;255,203,0;255,205,0;255,206,0;255,208,0;255,210,0;255,212,0;255,213,0;255,215,0;255,217,0;255,219,0;255,220,0;255,222,0;255,224,0;255,226,0;255,228,0;255,230,0;255,232,0;255,234,0;255,235,4;255,237,8;255,239,13;255,241,17;255,242,21;255,244,26;255,246,30;255,248,35;255,248,42;255,249,50;255,250,58;255,251,66;255,252,74;255,253,82;255,254,90;255,255,98;255,255,105;255,255,113;255,255,121;255,255,129;255,255,136;255,255,144;255,255,152;255,255,160;255,255,167;255,255,175;255,255,183;255,255,191;255,255,199;255,255,207;255,255,215;255,255,223;255,255,227;255,255,231;255,255,235;255,255,239;255,255,243;255,255,247;255,255,251;255,255,255;255,255,255;255,255,255;255,255,255]./255;
    otherwise
        LUT     = colormap('gray');
end


