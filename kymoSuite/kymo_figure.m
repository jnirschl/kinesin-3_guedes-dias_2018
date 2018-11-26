function varargout = kymo_figure(varargin)
% KYMO_FIGURE MATLAB code for kymo_figure.fig
%      KYMO_FIGURE, by itself, creates a new KYMO_FIGURE or raises the existing
%      singleton*.
%
%      H = KYMO_FIGURE returns the handle to a new KYMO_FIGURE or the handle to
%      the existing singleton*.
%
%      KYMO_FIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KYMO_FIGURE.M with the given input arguments.
%
%      KYMO_FIGURE('Property','Value',...) creates a new KYMO_FIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kymo_figure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kymo_figure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Author: Jeffrey J. Nirschl, Holzbaur Lab (University of Pennsylvania)
% Last Updated: 08/30/2017
% Distributable under BSD liscence. 
% Copyright (c) 2014-2017 Jeffrey J. Nirschl
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



% Edit the above text to modify the response to help kymo_figure

% Last Modified by GUIDE v2.5 11-Aug-2014 17:36:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kymo_figure_OpeningFcn, ...
                   'gui_OutputFcn',  @kymo_figure_OutputFcn, ...
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


% --- Executes just before kymo_figure is made visible.
function kymo_figure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kymo_figure (see VARARGIN)

% Choose default command line output for kymo_figure
handles.output = hObject;

% Determine where to start reading datafile
try
    inputdlg_user= {varargin{1,1}{1, 1}{1, 2}(12:end),varargin{1,1}{1, 1}{1, 3}(11:end)};
    handles.start = 1;
catch
    if isempty(varargin{1,1}{2,1})==0;
       inputdlg_user= {varargin{1,1}{2, 1}{1, 2}(12:end),varargin{1,1}{2, 1}{1, 3}(11:end)};
       handles.start = 2;
    elseif isempty(varargin{1,1}{3,1})==0;
       inputdlg_user= {varargin{1,1}{3, 1}{1, 2}(12:end),varargin{1,1}{3, 1}{1, 3}(11:end)};
       handles.start = 3;
    else
       inputdlg_user= {varargin{1,1}{4, 1}{1, 2}(12:end),varargin{1,1}{4, 1}{1, 3}(11:end)};
       handles.start = 4;
    end
end


% Initialize variables
if isempty(varargin)==1
    msgbox('Please restart and enter a kymo_tracer file!','Error');
    error('Please restart and enter a kymo_tracer file!','Error');
%    datafile = uigetfile({'*.mat';'*'},'Multiselect','off','Please select kymo_tracer file to load.');
%    datafile= load(datafile);
%    handles.data = datafile(:,:);
else
   handles.data = varargin{1,1};
end
handles.LineThick = 1;
handles.Density_height = 25;
handles.SchematicBox = 1;
handles.DistThreshBox = 1;
handles.DistThresh = 10;
handles.NeuronNum = 1;
handles.SpecifyNeuron = false;
handles.TraceStart = 1;
handles.TraceEnd = size(handles.data,2);
set(handles.ColorCodeCheckbox,'Value',0);
% set(handles.KymoSchematicCheck,'Value',1);

handles.FileDir = uigetdir('','Select the directory where image files are located.'); %Create working dir
handles.SaveDir = handles.FileDir;

handles.Files=dir([handles.FileDir filesep '*_contrast.tif']);
x = struct2cell(handles.Files); handles.Files=x(1,:)';
handles.filespath=handles.FileDir;

set(handles.listbox1,'String',handles.Files);

handles.Selectedfiles=get(handles.listbox1,'Value');

% Set options for user input dialog boxes
options.Resize='on'; options.WindowStyle='normal'; options.Interpreter='tex';
inputdlg_text= {'um/ Pixel','Sec/ Frame'};
% inputdlg_user= {varargin{1,1}{1, 1}{1, 2}(12:end),varargin{1,1}{1, 1}{1, 3}(11:end)};
% inputdlg_default= {varargin{1,1}{1, 1}{1, 2}(12:end),varargin{1,1}{1, 1}{1, 3}(11:end)};


    try
        param= inputdlg(inputdlg_text,'Input Parameters',1,inputdlg_user,options);
    catch err
        param= inputdlg(inputdlg_text,'Input Parameters',1,inputdlg_default,options);
    end

handles.pixel = str2num(param{1,1}); % Pixel conversion from user input dialog box
handles.tframe= str2num(param{2,1}); % Sec/ frame from user input dialog box


% Update Neuron Number
file= handles.Files{1};
fname= file;
dirnew= handles.filespath;

Nstop= strfind(handles.Files(handles.Selectedfiles),'_contrast'); Nstop =cell2mat(Nstop);
fname = handles.Files(handles.Selectedfiles); fname = cell2mat(fname);
fname = [fname(1:Nstop-1) '.tif'];

try
    for ii = handles.start:size(handles.data,1)
        if strcmp(handles.data{ii, 1}{2, 1}{1,1},fname) ==1;
            row = ii;
            handles.NeuronNum = row;
            set(handles.neuron_num,'String',num2str(row));
        end
    end

catch
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kymo_figure wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kymo_figure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

handles.Selectedfiles=get(handles.listbox1,'Value');
a= handles.Selectedfiles;
file= handles.Files{a};
fname= file;
dirnew= handles.filespath;

Nstop= strfind(handles.Files(handles.Selectedfiles),'_contrast'); Nstop =cell2mat(Nstop);
fname = handles.Files(handles.Selectedfiles); fname = cell2mat(fname);
fname = [fname(1:Nstop-1) '.tif'];

try
    for ii = handles.start:size(handles.data,1)
        if strcmp(handles.data{ii, 1}{2, 1}{1,1},fname) ==1;
            row = ii;
            handles.NeuronNum = row;
            set(handles.neuron_num,'String',num2str(row));
        end
    end

catch
end


% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in LoadImageBtn.
function LoadImageBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImageBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load tiffs
cd(handles.FileDir);
% handles.files = get(handles.filelist,'String');
handles.Selectedfiles=get(handles.listbox1,'Value');

a= handles.Selectedfiles;
file= handles.Files{a};
fname= file;
dirnew= handles.filespath;
cd(dirnew);

iptsetpref('ImshowBorder','tight');
im = imread([dirnew filesep file]); 

handles.kymograph = figure; imshow(im); % Create new figure w/ image
set(handles.kymograph,'Name',fname,'NumberTitle','off')

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in RefreshBtn.
function RefreshBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RefreshBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Selectedfiles=get(handles.listbox1,'Value');

a= handles.Selectedfiles;
file= handles.Files{a};
fname= file;
dirnew= handles.filespath;
cd(dirnew);

iptsetpref('ImshowBorder','tight');
im = imread([dirnew filesep file]);

if handles.SchematicBox ==1
    if isfield(handles, 'Schematic')
        if isvalid(handles.Schematic)
            figure(handles.Schematic)
        else
            handles.Schematic = figure;
        end
    else
        handles.Schematic = figure;
    end
    imshow(ones(size(im))); % Create new figure w/ image
    set(handles.Schematic,'Name',[fname(1:end-4) '_schematic.tif'],'NumberTitle','off');
    
    if isfield(handles, 'Density')
        if isvalid(handles.Density)
            figure(handles.Density)
        else
            handles.Density= figure;
        end
    else
        handles.Density= figure;
    end
    imshow(ones(handles.Density_height, size(im,2))); % Create new figure w/ image
    set(handles.Density,'Name',[fname(1:end-4) '_density.tif'],'NumberTitle','off');
end

Nstop= strfind(handles.Files(handles.Selectedfiles),'_contrast'); Nstop =cell2mat(Nstop);
fname = handles.Files(handles.Selectedfiles); fname = cell2mat(fname);
fname = [fname(1:Nstop-1) '.tif'];

row = str2double(get(handles.neuron_num,'String'));
if handles.SpecifyNeuron == 1
    row = str2double(get(handles.neuron_num,'String'));
else
    try
        for ii = handles.start:size(handles.data,1)
            if isempty(handles.data{ii, 1})==0
                if strcmp(handles.data{ii, 1}{2, 1}{1,1},fname) ==1;
                    row = ii;
                end
            end
        end
    catch
        msgbox('Unable to match kymograph image and neuron number. Please try specifying the neuron number using the box provided. See the programmer if problem persists.','Error');
        error('Unable to match kymograph image and neuron number. Please try specifying the neuron number using the box provided.','Error');
    end
end
    
handles.DistThreshBox = get(handles.ColorCodeCheckbox, 'value');
for jj = handles.TraceStart+1:handles.TraceEnd;
    if isempty(handles.data{row,jj})==0
        if handles.DistThreshBox ==1      
            if sum(diff(handles.data{row,jj}(:,1)).*handles.pixel) >handles.DistThresh; % Anterograde is positive
                color =  '-b'; % Blue/ Azul for Anterograde
            elseif sum(diff(handles.data{row,jj}(:,1)).*handles.pixel) < -handles.DistThresh % Retrograde is negative
                color = '-r'; % Red color for retrograde
            else 
                color = '-g';
            end
        else
            color = '-r';
        end

        x = handles.data{row,jj}(:,1);
        y = handles.data{row,jj}(:,2);

        try
            figure(handles.kymograph);
        catch
            handles.kymograph = figure; imshow(im); % Create new figure w/ image
            set(handles.kymograph,'Name',fname,'NumberTitle','off')
        end
        
        hold on; axis off;
        plot(x,y,color,'LineWidth',handles.LineThick);

        if handles.SchematicBox ==1;     
            figure(handles.Schematic); hold on; axis off;
            plot(x,y,'k-','LineWidth',handles.LineThick);
        end
    end
end

if handles.SchematicBox ==1
    if handles.LineThick > 1
        errordlg('Line thickness must be set to 1 for accurate density estimation!')
        warning('Line thickness must be set to 1 for accurate density estimation!')
        disp('Density estimation will not be saved. Please use a line thickness of 1 to save a density plot.')
    else
        % Plot density images
        ax              = handles.Schematic;
        schematic_im    = getframe();
        schematic_im    = logical(imcomplement(schematic_im.cdata));
        schematic_im    = bwmorph(schematic_im(:,:,1),'thin',inf);
        density_line    = nansum(schematic_im,1);
        csv_save_name   = strrep(get(handles.Density,'Name'),'contrast_density.tif',...
                        'line_profile_density.csv');
                    
        if exist(csv_save_name,'file')
            FID         = fopen(csv_save_name,'r');
            dataArray   = textscan(FID, '%f%f%f%[^\n\r]', 'Delimiter', ',',...
                                'EmptyValue', NaN, ...
                                'ReturnOnError', false);
            fclose(FID);
            csvwrite(csv_save_name, [1:size(im,2); density_line; dataArray{1,3}']');
        else
            csvwrite(csv_save_name, [1:size(im,2); density_line]',0,0);
        end
        density_im      = repmat(density_line,handles.Density_height,1);

        figure(handles.Density); hold on; axis off;
        iptsetpref('ImshowBorder','tight');
        imshow(normData(density_im))
        colormap(gca, pk_cmap('heat'))
        % colormap(gray)
        if max(density_line) <= 32
            caxis([0,32])
        elseif max(density_line) <= 64
            caxis([0,64])
        elseif max(density_line) <= 128
            caxis([0,128])
        elseif max(density_line) <= 255
            caxis([0,255])
        else
            caxis([0, 2^16])
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function LineThickSlider_Callback(hObject, eventdata, handles)
% hObject    handle to LineThickSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.LineThick = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LineThickSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineThickSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in ColorCodeCheckbox.
function ColorCodeCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to ColorCodeCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ColorCodeCheckbox
handles.DistThreshBox = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);


function DistThreshEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DistThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DistThreshEdit as text
%        str2double(get(hObject,'String')) returns contents of DistThreshEdit as a double
handles.DistThresh = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DistThreshEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DistThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in KymoSchematicCheck.
function KymoSchematicCheck_Callback(hObject, eventdata, handles)
% hObject    handle to KymoSchematicCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of KymoSchematicCheck
handles.SchematicBox = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in SaveBtn.
function SaveBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd(handles.SaveDir);
fname = get(handles.kymograph,'Name');

ax1         = gca(handles.kymograph);
kymo_save   = getframe(ax1);
fname1      = [fname(1:end-4) '_traced.tif' ];
imwrite(kymo_save.cdata,fname1,'tif');

if handles.SchematicBox ==1
    fname2          = get(handles.Schematic,'Name');
    ax2             = gca(handles.Schematic);
    schem_save      = getframe(ax2);
    imwrite(schem_save.cdata,fname2,'tif')
    
    fname3          = get(handles.Density,'Name');
    ax3             = gca(handles.Density);
    density_rgb     = getframe(ax3);
    
    density_gray    = getimage(handles.Density);
    
    figure(handles.Density);
    cbh = colorbar;
    
    set(cbh,'YTick', ...
            colon(min(caxis), ...
                  ceil(range(caxis))/2, ...
                  max(caxis)) )
    density_rgb_bar = getframe(handles.Density);
    savefig(handles.Density, strrep(fname3,'.tif','_colorbar.fig'))
    
    imwrite(density_rgb.cdata,fname3,'tif')
    imwrite(density_rgb_bar.cdata,strrep(fname3,'.tif','_colorbar.tif'),'tif')
    imwrite(uint8(density_gray),strrep(fname3,'.tif','_gray.tif'),'tif')
    imwrite(uint8(density_gray),strrep(fname3,'.tif','_gray.tif'),'tif')
    try
        figure(handles.kymograph)
        export_fig(strrep(fname1,'.tif','.eps'));
        
        figure(handles.Schematic)
        export_fig(strrep(fname1,'.tif','.eps'));
        
        figure(handles.Density)
        export_fig(strrep(fname1,'.tif','.eps'));
    catch
        saveas(ax1, strrep(fname1,'.tif','.svg'));
        saveas(ax2, strrep(fname2,'.tif','.svg'));
        saveas(ax3, strrep(fname3,'.tif','.svg'));
    end
end

% Update handles structure
guidata(hObject, handles);



function neuron_num_Callback(hObject, eventdata, handles)
% hObject    handle to neuron_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron_num as text
%        str2double(get(hObject,'String')) returns contents of neuron_num as a double
handles.NeuronNum       = str2double(get(hObject,'String'));
handles.SpecifyNeuron   = true

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function neuron_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function trace_start_Callback(hObject, eventdata, handles)
% hObject    handle to trace_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trace_start as text
%        str2double(get(hObject,'String')) returns contents of trace_start as a double
handles.TraceStart = floor(str2double(get(hObject,'String')));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function trace_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trace_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trace_end_Callback(hObject, eventdata, handles)
% hObject    handle to trace_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trace_end as text
%        str2double(get(hObject,'String')) returns contents of trace_end as a double
temp = (get(hObject,'String'));
if strcmpi(temp,'end')==1
    handles.TraceEnd = size(handles.data,2);
else
    handles.TraceEnd = str2double(get(hObject,'String'))+1;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function trace_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trace_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
