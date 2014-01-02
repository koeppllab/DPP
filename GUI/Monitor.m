function varargout = Monitor(varargin)
% MONITOR MATLAB code for Monitor.fig
%      MONITOR, by itself, creates a new MONITOR or raises the existing
%      singleton*.
%
%      H = MONITOR returns the handle to a new MONITOR or the handle to
%      the existing singleton*.
%
%      MONITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MONITOR.M with the given input arguments.
%
%      MONITOR('Property','Value',...) creates a new MONITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Monitor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Monitor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Monitor

% Last Modified by GUIDE v2.5 07-Apr-2013 16:04:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Monitor_OpeningFcn, ...
                   'gui_OutputFcn',  @Monitor_OutputFcn, ...
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


% --- Executes just before Monitor is made visible.
function Monitor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Monitor (see VARARGIN)

% Choose default command line output for Monitor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Monitor wait for user response (see UIRESUME)
% uiwait(handles.Monitor);


Data = get(hObject, 'UserData');


% Initialize GUI
cellsIdx = Data.CellsIdx;
set(handles.CellMenu, 'String', cellsIdx);

kineticModel = Data.Model.KineticModel;
NumReactions = kineticModel.NumReactions;

parameterNames = kineticModel.ReactionNames;
parameterNames{end+1} = 'Measurement Error';

% if heterogeneity present
if (~isempty(Data.Model.HyperPriorMu))
    parameterNames{end+1} = 'a vs. b (extr)';
end

set(handles.ParameterMenu, 'String', parameterNames);

if (Data.UpdateStateDistribution == 1)
    DrawStateStatistics(handles);
end

if (Data.UpdateParameterDistribution == 1)
    DrawParameterStatistics(handles);
end

Data.UpdateStateDistribution = 0;
Data.UpdateParameterDistribution = 0;

set(handles.Monitor, 'UserData', Data);


function DrawStateStatistics(handles)

Data = get(handles.Monitor, 'UserData');

timeIdx = Data.TimeIndex;

if (timeIdx > 0)

    
    % set parameters from UserData
    Posterior = Data.Posterior;

    cells = Data.Cells;
    cellsIdx = Data.CellsIdx;
    options = Data.qPlotOptions;

    % read cell-to-plot from GUI
    cellToPlot = get(handles.CellMenu, 'Value');
    
    options.cellIdx = cellToPlot;
    options.numPoints = timeIdx * 10;
    options.FigureHandlesCells = handles.StateAxes;
    
    %Call plotting function
    PlotStateStatistics(Posterior, options, cells);
    
end

function DrawParameterStatistics(handles)

Data = get(handles.Monitor, 'UserData');     

Posterior = Data.Posterior;
model = Data.Model;

NumReactions = model.NumReactions;

selParameterIdx = get(handles.ParameterMenu, 'Value');

%% Hardcoded up to now: CLEAN UP.

numBins = 20;

if (selParameterIdx <= NumReactions) % kinetic parameter
    aPrior = model.RandomEffectA(selParameterIdx);
    bPrior = model.RandomEffectB(selParameterIdx);
    Samples = ...
        DrawParameterPosteriorSamples(Posterior, aPrior, bPrior, ...
        selParameterIdx, 3);
    
    
    [h, b] = EvaluateMarkovChain(Samples, 1, numBins);

    p = area(handles.ParameterAxes, b, h);
    set(p, 'FaceColor', [0.8, 0.8, 0.8]);
    
    xlabel(handles.ParameterAxes, ...
        [model.KineticModel.ReactionNames{selParameterIdx} 'in 1/s']);
    ylabel('Posterior');
elseif (selParameterIdx == NumReactions + 1) % measurement error
    a = model.NoiseA;
    b = model.NoiseB;

    Samples = 1 ./ sqrt(gamrnd(a, 1/b, 10000, 1));
    
    [h, b] = EvaluateMarkovChain(Samples, 1, numBins);

    p = area(handles.ParameterAxes, b, h);
    set(p, 'FaceColor', [0.8, 0.8, 0.8]);
    
    xlabel(handles.ParameterAxes, 'Measurement Error (\sigma)');
    ylabel('Posterior');
elseif (selParameterIdx == NumReactions + 2) % a (extrinsic)
    aPrior = model.RandomEffectA(1); %doesn't matter!
    bPrior = model.RandomEffectB(1);
    
    [~, data] = DrawParameterPosteriorSamples(Posterior, aPrior, bPrior, ...
        1, 3);
    
    data = log(data([2, 1], :)');
    
    num2DBins = numBins;
    [counts, centers] = hist3(data, [num2DBins, num2DBins]);
    contourf(handles.ParameterAxes, ...
        centers{2}, centers{1}, counts, 200, 'LineColor', 'none'); 
    CM = colormap('gray');
    CM = 1-CM;
    colormap(handles.ParameterAxes, CM);
    xlabel(handles.ParameterAxes, 'ln \alpha');
    ylabel(handles.ParameterAxes, 'ln \beta');
end

box(handles.ParameterAxes, 'off');


% --- Outputs from this function are returned to the command line.
function varargout = Monitor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ParameterMenu.
function ParameterMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ParameterMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ParameterMenu

Data = get(handles.Monitor, 'UserData');
Data.UpdateParameterDistribution = 1;
set(handles.Monitor, 'UserData', Data);

if (Data.Running == 0)
   DrawParameterStatistics(handles); 
end

% --- Executes during object creation, after setting all properties.
function ParameterMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CellMenu.
function CellMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CellMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CellMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CellMenu
Data = get(handles.Monitor, 'UserData');
Data.UpdateStateDistribution = 1;
set(handles.Monitor, 'UserData', Data);

% If not running, update plot
if (Data.Running == 0)
   DrawStateStatistics(handles); 
end

% --- Executes during object creation, after setting all properties.
function CellMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close Monitor.
function Monitor_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Monitor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%set(handles.Monitor, 'visible', 'off');
%drawnow;
%delete(hObject);


% --- Executes on selection change in ParameterMenu.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to ParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ParameterMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ParameterMenu


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CellMenu.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to CellMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CellMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CellMenu


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
