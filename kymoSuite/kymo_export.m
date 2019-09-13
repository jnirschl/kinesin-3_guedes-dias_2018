function varargout = kymo_export(varargin)
% KYMO_EXPORT MATLAB code for kymo_export.fig
%      KYMO_EXPORT, by itself, creates a new KYMO_EXPORT or raises the existing
%      singleton*.
%
%      H = KYMO_EXPORT returns the handle to a new KYMO_EXPORT or the handle to
%      the existing singleton*.
%
%      KYMO_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KYMO_EXPORT.M with the given input arguments.
%
%      KYMO_EXPORT('Property','Value',...) creates a new KYMO_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kymo_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kymo_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kymo_export

% Last Modified by GUIDE v2.5 03-Nov-2014 18:33:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kymo_export_OpeningFcn, ...
                   'gui_OutputFcn',  @kymo_export_OutputFcn, ...
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


% --- Executes just before kymo_export is made visible.
function kymo_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kymo_export (see VARARGIN)

% Choose default command line output for kymo_export

handles.output = hObject;
try
   handles.analyzed_data = varargin{1,1};
catch err
   msgbox('Please load your analyzed datafile and enter kymo_expot(analyzed_data).','Error'); 
end
handles.fileID = handles.analyzed_data.fname{1,1}{1,1}(1:end-4);
handles.saveDir = uigetdir('','Select save directory.');
handles.date = datestr(now(),'yyyy-mm-dd');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kymo_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kymo_export_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in non_pause_instvel.
function non_pause_instvel_Callback(hObject, eventdata, handles)
% hObject    handle to non_pause_instvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     non_pause_instvel = handles.analyzed_data.compiled_instvel(abs(handles.analyzed_data.compiled_instvel)>0.1);
     x = nanmean(non_pause_instvel); sd = nanstd(non_pause_instvel,0,1); n = sum(isnan(non_pause_instvel)==0);
     fprintf('Non-Paused Instantaneous velocity for all vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
     dlmwrite([handles.fileID '_non_pause_instvel_' handles.date '.txt'],non_pause_instvel,'\n');
     mat2clip(non_pause_instvel);
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end


% --- Executes on button press in retro_instvel.
function retro_instvel_Callback(hObject, eventdata, handles)
% hObject    handle to retro_instvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     retro_instvel = handles.analyzed_data.compiled_instvel(handles.analyzed_data.compiled_instvel<-0.1);
     if ~isempty(retro_instvel)
        x = nanmean(retro_instvel); sd = nanstd(retro_instvel,0,1); n = sum(isnan(retro_instvel)==0);
        fprintf('Non-Paused Instantaneous velocity in Retrograde Direction for all vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_retro_instvel_' handles.date '.txt'],retro_instvel,'\n');
        mat2clip(retro_instvel);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow(err)
end

% --- Executes on button press in antero_instvel.
function antero_instvel_Callback(hObject, eventdata, handles)
% hObject    handle to antero_instvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
    antero_instvel = handles.analyzed_data.compiled_instvel(handles.analyzed_data.compiled_instvel>0.1);
    if ~isempty(antero_instvel)
        x = nanmean(antero_instvel); sd = nanstd(antero_instvel,0,1); n = sum(isnan(antero_instvel)==0);
        fprintf('Non-Paused Instantaneous velocity in Anterograde Direction for all vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_antero_instvel_' handles.date '.txt'],antero_instvel,'\n');
        mat2clip(antero_instvel);
    end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow(err)
end

% --- Executes on button press in cum_dist.
function cum_dist_Callback(hObject, eventdata, handles)
% hObject    handle to cum_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     cum_dist = handles.analyzed_data.cumdist(isnan(handles.analyzed_data.cumdist)==0);
     if ~isempty(cum_dist)
        x = nanmean(cum_dist); sd = nanstd(cum_dist,0,1); n = sum(isnan(cum_dist)==0);
        fprintf('CumDist All Vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_cum_dist_' handles.date '.txt'],cum_dist,'\n');
        mat2clip(cum_dist);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end


% --- Executes on button press in net_dist.
function net_dist_Callback(hObject, eventdata, handles)
% hObject    handle to net_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     net_dist = handles.analyzed_data.netdist(isnan(handles.analyzed_data.netdist)==0);
     if ~isempty(net_dist)
        x = nanmean(net_dist); sd = nanstd(net_dist,0,1); n = sum(isnan(net_dist)==0);
        fprintf('Net Displacement All Vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_net_dist_' handles.date '.txt'],net_dist,'\n');
        mat2clip(net_dist);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end

% --- Executes on button press in frac_antero.
function frac_antero_Callback(hObject, eventdata, handles)
% hObject    handle to frac_antero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     frac_antero = sum(handles.analyzed_data.netdist>=10,2)./sum(isnan(handles.analyzed_data.netdist)==0,2);
     if ~isempty(frac_antero)
         x = nanmean(frac_antero); sd = nanstd(frac_antero,0,1); n = sum(isnan(frac_antero)==0);
         fprintf('Fraction of Events with Net Retrograde Displacement >10um) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
         dlmwrite([handles.fileID '_frac_antero_' handles.date '.txt'],frac_antero,'\n');
         mat2clip(frac_antero);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end

% --- Executes on button press in frac_pause.
function frac_pause_Callback(hObject, eventdata, handles)
% hObject    handle to frac_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     frac_pause = sum(handles.analyzed_data.netdist>-10 & handles.analyzed_data.netdist<10,2)./sum(isnan(handles.analyzed_data.netdist)==0,2);
     if ~isempty(frac_pause)
         x = nanmean(frac_pause); sd = nanstd(frac_pause,0,1); n = sum(isnan(frac_pause)==0);
         fprintf('Fraction of Events -10um < Net Displacement <10um) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
         dlmwrite([handles.fileID '_frac_pause_' handles.date '.txt'],frac_pause,'\n');
         mat2clip(frac_pause);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end

% --- Executes on button press in frac_retro.
function frac_retro_Callback(hObject, eventdata, handles)
% hObject    handle to frac_retro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
    frac_retro = sum(handles.analyzed_data.netdist<=-10,2)./sum(isnan(handles.analyzed_data.netdist)==0,2);
    if isempty(frac_retro)
        x = nanmean(frac_retro); sd = nanstd(frac_retro,0,1); n = sum(isnan(frac_retro)==0);
        fprintf('Fraction of Events with Net Retrograde Displacement >10um) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_frac_retro_' handles.date '.txt'],frac_retro,'\n');
        mat2clip(frac_retro);
    end 
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end

% --- Executes on button press in antero_switches.
function antero_switches_Callback(hObject, eventdata, handles)
% hObject    handle to antero_switches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     antero_switches = handles.analyzed_data.dir_totalswitches(handles.analyzed_data.netdist>=10);
     if ~isempty(antero_switches)
        x = nanmean(antero_switches); sd = nanstd(antero_switches,0,1); n = sum(isnan(antero_switches)==0);
        fprintf('Switches > Threshold in Vesicles with Net Anterograde Displacement >10um) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_antero_switches_' handles.date '.txt'],antero_switches,'\n');
        mat2clip(antero_switches);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end

% --- Executes on button press in retro_switches.
function retro_switches_Callback(hObject, eventdata, handles)
% hObject    handle to retro_switches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     retro_switches = handles.analyzed_data.dir_totalswitches(handles.analyzed_data.netdist<=-10);
     if ~isempty(retro_switches)
         x = nanmean(retro_switches); sd = nanstd(retro_switches,0,1); n = sum(isnan(retro_switches)==0);
        fprintf('Switches > Threshold in Vesicles with Net Retrograde Displacement >10um) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_retro_switches_' handles.date '.txt'],retro_switches,'\n');
        mat2clip(retro_switches);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end

% --- Executes on button press in pause_switches.
function pause_switches_Callback(hObject, eventdata, handles)
% hObject    handle to pause_switches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     pause_switches = handles.analyzed_data.dir_totalswitches(handles.analyzed_data.netdist>-10 & handles.analyzed_data.netdist<10);
     if ~isempty(pause_switches)
        x = nanmean(pause_switches); sd = nanstd(pause_switches,0,1); n = sum(isnan(pause_switches)==0);
        fprintf('Switches > Threshold in Vesicles with Net Displacement in Either Direction <10um) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_pause_switches_' handles.date '.txt'],pause_switches,'\n');
        mat2clip(pause_switches);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end


% --- Executes on button press in total_switches.
function total_switches_Callback(hObject, eventdata, handles)
% hObject    handle to total_switches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     total_switches = handles.analyzed_data.dir_totalswitches(isnan(handles.analyzed_data.netdist)==0);
     if ~isempty(total_switches)
        x = nanmean(total_switches); sd = nanstd(total_switches,0,1); n = sum(isnan(total_switches)==0);
        fprintf('Switches > Threshold in All Vesicles) \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_total_switches_' handles.date '.txt'],total_switches,'\n');
        mat2clip(total_switches);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end


% --- Executes on button press in all_instvel.
function all_instvel_Callback(hObject, eventdata, handles)
% hObject    handle to all_instvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     all_instvel = handles.analyzed_data.compiled_instvel(isnan(handles.analyzed_data.compiled_instvel)==0);
     if ~isempty(all_instvel)
         x = nanmean(all_instvel); sd = nanstd(all_instvel,0,1); n = sum(isnan(all_instvel)==0);
        fprintf('All Instantaneous velocity for all vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_all_instvel_' handles.date '.txt'],all_instvel,'\n');
        mat2clip(all_instvel);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end


% --- Executes on button press in net_vel.
function net_vel_Callback(hObject, eventdata, handles)
% hObject    handle to net_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.saveDir);
try
     netvel = handles.analyzed_data.net_velocity(isnan(handles.analyzed_data.net_velocity)==0);
     if ~isempty(netvel)
        x = nanmean(netvel); sd = nanstd(netvel,0,1); n = sum(isnan(netvel)==0);
        fprintf('Net Velocity All Vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_net_vel_' handles.date '.txt'],netvel,'\n');
        mat2clip(netvel);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end


% --- Executes on button press in CumTime.
function CumTime_Callback(hObject, eventdata, handles)
% hObject    handle to CumTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd(handles.saveDir);
try
     cum_time = handles.analyzed_data.cumtime(isnan(handles.analyzed_data.cumtime)==0);
     if ~isempty(cum_time)
         x = nanmean(cum_time); sd = nanstd(cum_time,0,1); n = sum(isnan(cum_time)==0);
        fprintf('CumTime All Vesicles \n mean = %0-6.4f \n std = %0-6.4f \n n = %0-6.0f \n',x,sd,n)
        dlmwrite([handles.fileID '_cum_time_' handles.date '.txt'],cum_time,'\n');
        mat2clip(cum_time);
     end
catch err
    msgbox('Unable to  export to .txt or copy to clipboard. Make sure you have loaded a file and that Mat2Clip function is added to the Matlab Path.','Error');
    rethrow err
end
