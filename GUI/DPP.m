function varargout = DPP(varargin)
% DPP MATLAB code for DPP.fig
%      DPP, by itself, creates a new DPP or raises the existing
%      singleton*.
%
%      H = DPP returns the handle to a new DPP or the handle to
%      the existing singleton*.
%
%      DPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DPP.M with the given input arguments.
%
%      DPP('Property','Value',...) creates a new DPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DPP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DPP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DPP

% Last Modified by GUIDE v2.5 09-Apr-2013 10:22:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DPP_OpeningFcn, ...
    'gui_OutputFcn',  @DPP_OutputFcn, ...
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


% --- Executes just before DPP is made visible.
function DPP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DPP (see VARARGIN)




% Choose default command line output for DPP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DPP wait for user response (see UIRESUME)
% uiwait(handles.DPP);


%% Initialize GUI
Data = get(handles.DPP, 'UserData');

% Load parameter table first!
parameterTable = handles.ParameterTable;
set(parameterTable, 'Data', {}, 'ColumnName', {'Name', 'Type', 'a', 'b'});


modelList = dir('../Models/*.mat');

if (isempty(modelList))
    fprintf('No kinetic models found!\n');
    set(handles.KineticModelMenu, 'String', 'No Model Found!');
    return;
end

modelNames = cell(length(modelList), 1);

for i=1:length(modelList)
    modelData = load(['../Models/' modelList(i).name]);
    Data.Models{i} = modelData.kineticModel;
    modelNames{i} = Data.Models{i}.Name;
end

set(handles.KineticModelMenu, 'String', modelNames);
LoadModelToGUI(Data.Models{1}, handles);

% Store kineticModel to UserData
Data.kineticModel = Data.Models{1};

muHyperTable = handles.MuHyperTable;
set(muHyperTable, 'Data', {1; 1000},...
    'RowName', {'a', 'b'});

sigmaHyperTable = handles.SigmaHyperTable;
set(sigmaHyperTable, 'Data', {2, 0; 0, 2}, 'ColumnName', {'a', 'b'},...
    'RowName', {'a', 'b'}, 'ColumnEditable', [true, true]);


%% Initialze Input Data

inputParams.inputTimes = 100*60;
inputParams.inputLevels = 1;

Data.InputParams = inputParams;

Data.Running = 0; % algorithm not running at the beginning
set(handles.DPP, 'UserData', Data);
DrawInputFunction({}, handles);

if (matlabpool('size')>0)
    set(handles.ParallelButton, 'Value', 1);
end





% --- Outputs from this function are returned to the command line.
function varargout = DPP_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ModelName_Callback(hObject, eventdata, handles)
% hObject    handle to ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ModelName as text
%        str2double(get(hObject,'String')) returns contents of ModelName as a double


% --- Executes during object creation, after setting all properties.
function ModelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LoadModelToGUI(kineticModel, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('Loading kinetic model...\n');

NumReactions = kineticModel.NumReactions;
table = handles.ParameterTable;

for i=1:NumReactions
    Data{i, 1} = kineticModel.ReactionNames{i};
    Data{i, 2} = 'Homogeneous';
    Data{i, 3} = 3;
    Data{i, 4} = 3/kineticModel.c(i);
end


set(table, 'Data', Data);

speciesMenu = handles.SpeciesMenu;
set(speciesMenu, 'String', kineticModel.SpeciesNames);


proteinIdx =strcmp(kineticModel.SpeciesNames, 'Protein');
if (sum(proteinIdx)==1)
    set(speciesMenu, 'Value', find(proteinIdx));
end

set(handles.MuHyperText, 'visible', 'off');
set(handles.SigmaHyperText, 'visible', 'off');
set(handles.MuHyperTable, 'visible', 'off');
set(handles.SigmaHyperTable, 'visible', 'off');

set(handles.NoHeterogeneousReactionText, 'visible', 'on');

drawnow;


% --- Executes when selected cell(s) is changed in ParameterTable.
function ParameterTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to ParameterTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in ParameterTable.
function ParameterTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ParameterTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


row = eventdata.Indices(1);
col = eventdata.Indices(2);
data = eventdata.EditData;

AllData = get(hObject, 'Data');

if (col==2)
    
    hetData = AllData(:, 2);
    hetIdx = find(strcmp('Heterogeneous', hetData));
    prevHetIdx = setdiff(hetIdx, row);
    
    if (strcmp(data, 'Heterogeneous'))
        if (~isempty(prevHetIdx))
            AllData{prevHetIdx, 2} = 'Homogeneous';
            set(hObject, 'Data', AllData);
        else
            set(handles.MuHyperText, 'visible', 'on');
            set(handles.SigmaHyperText, 'visible', 'on');
            set(handles.MuHyperTable, 'visible', 'on');
            set(handles.SigmaHyperTable, 'visible', 'on');
            
            set(handles.NoHeterogeneousReactionText, 'visible', 'off');
        end
        
    else
        if (isempty(prevHetIdx))
            set(handles.MuHyperText, 'visible', 'off');
            set(handles.SigmaHyperText, 'visible', 'off');
            set(handles.MuHyperTable, 'visible', 'off');
            set(handles.SigmaHyperTable, 'visible', 'off');
            
            set(handles.NoHeterogeneousReactionText, 'visible', 'on');
        end
    end
end


% --- Executes on selection change in SpeciesMenu.
function SpeciesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SpeciesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SpeciesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SpeciesMenu


% --- Executes during object creation, after setting all properties.
function SpeciesMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpeciesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NoiseModelMenu.
function NoiseModelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NoiseModelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NoiseModelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NoiseModelMenu


% --- Executes during object creation, after setting all properties.
function NoiseModelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NoiseModelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on selection change in ActivationTypeMenu.
function ActivationTypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ActivationTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ActivationTypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ActivationTypeMenu

contents = cellstr(get(hObject,'String'));
selectedItem = contents{get(hObject,'Value')};

if (strcmp(selectedItem, 'Temporal'))
    set(handles.InputLoadButton, 'visible', 'on');
    set(handles.InputPathText, 'visible', 'on');
    set(handles.InputPathEdit, 'visible', 'on');
else
    set(handles.InputLoadButton, 'visible', 'off');
    set(handles.InputPathText, 'visible', 'off');
    set(handles.InputPathEdit, 'visible', 'off');
    
    inputParams.inputTimes = 100*60;
    inputParams.inputLevels = 1;
    
    Data = get(handles.DPP, 'UserData');
    Data.InputParams = inputParams;
    set(handles.DPP, 'UserData', Data);
    
    DrawInputFunction(inputParams, handles);
end



% --- Executes during object creation, after setting all properties.
function ActivationTypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ActivationTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InputPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InputPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputPathEdit as text
%        str2double(get(hObject,'String')) returns contents of InputPathEdit as a double


% --- Executes during object creation, after setting all properties.
function InputPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in InputLoadButton.
function InputLoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to InputLoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [FileName,PathName,FilterIndex] = uigetfile('*.mat', 'Load input file');
    
    if (FilterIndex == 0)
        return;
    end
    inputData = load([PathName FileName]);
    inputParams = inputData.inputParams;
    
    Data = get(handles.DPP, 'UserData');
    Data.InputParams = inputParams;
    set(handles.DPP, 'UserData', Data);
    set(handles.InputPathEdit, 'String', FileName)
    
    DrawInputFunction(inputParams, handles);
catch err
    fprintf('Invalid input file!\n');
end




function DrawInputFunction(inputParams, handles)

if (~isempty(inputParams))
    times = [0, inputParams.inputTimes];
    levels = [inputParams.inputLevels, inputParams.inputLevels(end)];
    
    stairs(handles.InputAxes, times/60, levels,...
        'g', 'LineWidth', 2);
else
    inputAxes = handles.InputAxes;
    plot(handles.InputAxes, [0, 100], [1, 1], 'g', 'LineWidth', 2);
    
end

box(handles.InputAxes, 'off');
xlabel(handles.InputAxes, 'Time in min');
ylabel(handles.InputAxes, 'Input');



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ActivationTypeMenu.
function ActivationTypeMenu_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ActivationTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function MeasurementPriorShapeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MeasurementPriorShapeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeasurementPriorShapeEdit as text
%        str2double(get(hObject,'String')) returns contents of MeasurementPriorShapeEdit as a double


% --- Executes during object creation, after setting all properties.
function MeasurementPriorShapeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasurementPriorShapeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MeasurementPriorScalingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MeasurementPriorScalingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeasurementPriorScalingEdit as text
%        str2double(get(hObject,'String')) returns contents of MeasurementPriorScalingEdit as a double


% --- Executes during object creation, after setting all properties.
function MeasurementPriorScalingEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasurementPriorScalingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MeasurementPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MeasurementPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeasurementPathEdit as text
%        str2double(get(hObject,'String')) returns contents of MeasurementPathEdit as a double


% --- Executes during object creation, after setting all properties.
function MeasurementPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasurementPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadMeasurementButton.
function LoadMeasurementButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMeasurementButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [FileName,PathName,FilterIndex] = uigetfile('*.mat', 'Load input file');
    
    if (FilterIndex == 0)
        return;
    end
    
    data = load([PathName FileName]);
    cells = data.cells;
    
    % Store measurements in UserData structure
    Data = get(handles.DPP, 'UserData');
    Data.cells = cells;
    set(handles.DPP, 'UserData', Data);
    
    NumCells = length(cells);
    NumCellsStr = ['Number of Cells: ' num2str(NumCells)];
    set(handles.NumberOfCellsText, 'String', NumCellsStr);
    set(handles.MeasurementPathEdit, 'String', FileName);
    
    cellIdx = 1:NumCells;
    str = num2str(cellIdx(1));
    for k=2:NumCells
        str = [str ',' num2str(cellIdx(k))];
    end
    set(handles.CellsToMonitorEdit, 'String', str);
    
    DrawMeasurements(cells, handles);
catch
    fprintf('Invalid Measurement File!\n');
end

function DrawMeasurements(cells, handles)

proteinAxes = handles.ProteinAxes;
numCells = length(cells);

for k=1:numCells
    
    time = [0, cells{k}.MeasurementTime];
    protein = [0, cells{k}.Measurement];
    
    p = plot(proteinAxes, time/60, protein, 'r-o'); hold on;
    set(p, 'LineWidth', 1);
    set(p, 'MarkerFaceColor', 'b');
    set(p, 'MarkerEdgeColor', 'b');
    set(p, 'MarkerSize', 3);
    
    xlabel(proteinAxes, 'Time in min');
    ylabel(proteinAxes, 'Abundance');
    box(proteinAxes, 'off');
    
    drawnow;
end

hold off;



function NumberOfParticlesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to NumberOfParticlesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumberOfParticlesEdit as text
%        str2double(get(hObject,'String')) returns contents of NumberOfParticlesEdit as a double


% --- Executes during object creation, after setting all properties.
function NumberOfParticlesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberOfParticlesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BurnInEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BurnInEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BurnInEdit as text
%        str2double(get(hObject,'String')) returns contents of BurnInEdit as a double


% --- Executes during object creation, after setting all properties.
function BurnInEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BurnInEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function callbackIteration(handles)
%drawnow;

Data = get(handles.DPP, 'UserData');

if (Data.Running == 0)
    error('DPP:STOP', '');
end

%% update Monitor GUI
Monitor;

function callbackTime(handles, Posterior, timeIndex)

Data = get(handles.DPP, 'UserData');
MonitorData = get(Data.MonitorHandle, 'UserData');
MonitorData.Posterior = Posterior;
MonitorData.TimeIndex = timeIndex;
MonitorData.UpdateStateDistribution = 1;
MonitorData.UpdateParameterDistribution = 1;
Monitor('UserData', MonitorData);


% --- Executes on button press in StopButton.
function StopButton_Callback(hObject, eventdata, handles)
% hObject    handle to StopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = get(handles.DPP, 'UserData');
Data.Running = 0;
set(handles.DPP, 'UserData', Data);


% --- Executes when user attempts to close DPP.
function DPP_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to DPP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

Data = get(hObject, 'UserData');
if (isfield(Data, 'MonitorHandle'));
    delete(findobj(Data.MonitorHandle));
end
delete(hObject);



function CellsToMonitorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CellsToMonitorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellsToMonitorEdit as text
%        str2double(get(hObject,'String')) returns contents of CellsToMonitorEdit as a double


% --- Executes during object creation, after setting all properties.
function CellsToMonitorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellsToMonitorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SwapCheckbox.
function SwapCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to SwapCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SwapCheckbox


% --- Executes on selection change in KineticModelMenu.
function KineticModelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to KineticModelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns KineticModelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from KineticModelMenu

Data = get(handles.DPP, 'UserData');

modelIdx = get(hObject, 'Value');
kineticModel = Data.Models{modelIdx};
LoadModelToGUI(kineticModel, handles);

% Store kinetic model to UserData
Data.kineticModel = kineticModel;
set(handles.DPP, 'UserData', Data);


% --- Executes during object creation, after setting all properties.
function KineticModelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KineticModelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ParallelButton.
function ParallelButton_Callback(hObject, eventdata, handles)
% hObject    handle to ParallelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ParallelButton

StartPool = get(hObject, 'Value');

if (StartPool == 1)
    numCurrWorkers = matlabpool('size');
    
    if (numCurrWorkers>0)
        fprintf('Matlabpool already running (%d workers)\n', numCurrWorkers);
    else
        matlabpool;
    end
else
    matlabpool close;
end


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RunButton

Data = get(handles.DPP, 'UserData');


if (Data.Running==0) % "On-Button"
    %% Get all relevant data from the UserData structure
    try
        fprintf('Starting DPP...\n');
        
        H = 0; %initialize Monitor handle to zero
        
        Data.Running = 1;
        set(handles.DPP, 'UserData', Data);
        set(hObject, 'String', 'Stop!');
        
        try
            cells = Data.cells;
        catch err
            error('DPP:NoDataLoaded', 'No protein measurements loaded!\n');
        end
        
        NumCellsAll = length(cells);
        
        % Read monitored/stored cells from GUI
        try
            cellIdx = str2num(get(handles.CellsToMonitorEdit, 'String'));
        catch err
            error('DPP:InvalidCellEnumeration', 'Invalid enumeration of cells to monitor!');
        end
        
        if (max(cellIdx)>NumCellsAll)
            error('DPP:InvalidCellEnumeration',...
                'Indices of cells to monitor are limited to the number of cells available (%d)\n', ...
                NumCells);
        end

        cells = cells(cellIdx);
        NumCells = length(cells);
        
        % Kinetic parameters and input
        kineticModel = Data.kineticModel;
        inputParams = Data.InputParams;
        inputParams.inputLevels(end+1) = 0;
        inputParams.inputRateIndex = int32(kineticModel.InputRateIdx);
        NumReactions = kineticModel.NumReactions;
        NumSpecies = kineticModel.NumSpecies;
        
        parData = get(handles.ParameterTable, 'Data');
        
        aPrior = cell2mat(parData(:, 3));
        bPrior = cell2mat(parData(:, 4));
        
        hetReactionIdx = strcmp(parData(:, 2), 'Heterogeneous');
        targetReactionIdx = (1:NumReactions)';
        
        % Hyperparameters
        if (sum(hetReactionIdx)==0) %no heterogeneous reaction
            MuHyper = [];
            SigmaHyper = [];
        else
            data = cell2mat(get(handles.MuHyperTable, 'Data'));
            MuHyper = log(data);
            
            data = cell2mat(get(handles.SigmaHyperTable, 'Data'));
            SigmaHyper = data;
        end
        
        % Morphological features not supported yet in GUI
        MuMorph = [];
        SigmaMorph = [];
        
        
        % Measurement noise parameters
        measurementDensityIdx = get(handles.NoiseModelMenu, 'Value');
        if (measurementDensityIdx == 1)
            mDens = 'logn';
        else
            mDens = 'normal';
        end
        
        measurementA = str2double(get(handles.MeasurementPriorShapeEdit, 'String'));
        measurementB = str2double(get(handles.MeasurementPriorScalingEdit, 'String'));
        
        % first parameter (0.15) is ignored if measurement noise is estimated from
        % the data
        mParams = [0.15, measurementA, measurementB];
        
        measuredSpeciesIdx = get(handles.SpeciesMenu, 'Value');
        OutputWeightMatrix = zeros(NumSpecies, 1);
        OutputWeightMatrix(measuredSpeciesIdx) = 1;
        
        MeasurementTime = cells{1}.MeasurementTime;
        
        % Initialize mixed-effect state-space model
        model = InitializeMEGSM(kineticModel, MeasurementTime, ...
            OutputWeightMatrix, mParams, mDens, targetReactionIdx, ...
            hetReactionIdx, aPrior, bPrior, ...
            MuHyper, SigmaHyper, [], MuMorph, SigmaMorph, ...
            @GetPieceWiseConstantInput, inputParams);
        
        % Create plot options (path statistics)
        qPlotOptions = CreatePlotOptions([]);
        qPlotOptions.cellIdx = 1;
        
        % Monitor mRNA and Protein
        mRNAIdx = find(strcmp(kineticModel.SpeciesNames, 'mRNA'));
        proteinIdx = find(strcmp(kineticModel.SpeciesNames, 'Protein'));
        qPlotOptions.stateIdx = [mRNAIdx, proteinIdx];
        
        options = CreateDPPOptions();
        options.BlockSize = 200;
        

        % Store all cells that where selected in the textbox
        options.StorePathCellIdx = 1:NumCells;%cellIdx;
        options.StorePaths = 1 + get(handles.SwapCheckbox, 'Value');
        
        
        % Get number of particles
        numParticles = str2num(get(handles.NumberOfParticlesEdit, 'String'));
        options.M = numParticles;
        options.burnIn = ...
            round(numParticles * str2num(get(handles.BurnInEdit, 'String')) / 100);
        
        % Create plot options (parameter statistics). Pass an empty handle
        % to prefent a new figure being created.
        histPlotOptions = CreateHistogramOptions([]);
        
        % Initialize particle distribution at time zero.
        pDist = InitializeParticleDistributionME(options.M, model, NumCells);
        
        % % Get handle to monitor GUI and set properties
        
        % Set to the original cellIdx selected in the textbox.
        MonitorData.CellsIdx = cellIdx;%options.StorePathCellIdx;
        MonitorData.Cells = cells;
        MonitorData.qPlotOptions = qPlotOptions;
        MonitorData.Parameters = 'c1';
        MonitorData.Running = 1;
        MonitorData.Posterior = pDist;
        MonitorData.TimeIndex = 0;
        MonitorData.UpdateParameterDistribution = 1;
        MonitorData.UpdateStateDistribution = 0;
        MonitorData.Model = model;
        H = Monitor('UserData', MonitorData);
        
        Data.MonitorHandle = H;
        
        % Run DPP algorithm
        options.CallBackTime = ...
            @(Posterior, timeIndex)callbackTime(handles, Posterior, timeIndex);
        options.CallBackIteration = @()callbackIteration(handles);
        
        Data.Running = 1;
        set(handles.DPP, 'UserData', Data);
        
        try
            pDist = RunDPP(pDist, model, cells, options, histPlotOptions, ...
                qPlotOptions);
            
            Data.Posterior = pDist;
            Data.Running = 0;
            Data.Model = model;
            set(handles.DPP, 'UserData', Data);
            
            fprintf('Inference completed successfully. Press the ''Export'' button to save you results.\n');
            
        catch err
            if strcmp(err.identifier, 'DPP:STOP')
                fprintf('Algorithm stopped by user...\n');
            else
                error(err.identifier, err.message);
            end
        end
        set(hObject, 'String', 'Run!');
    catch err
        fprintf('ERROR occured while running DPP!\n')
        fprintf('ID: %s\n', err.identifier);
        fprintf('Msg: %s\n', err.message);
        
        Data.Running = 0;
        set(handles.DPP, 'UserData', Data);
        set(hObject, 'String', 'Run!');
    end
    
    if (H~=0)
        MonitorData = get(H, 'UserData');
        MonitorData.Running = 0;
        set(H, 'UserData', MonitorData);
    end
    
else % "Off-Botton"
    Data = get(handles.DPP, 'UserData');
    Data.Running = 0;
    set(handles.DPP, 'UserData', Data);
    set(hObject, 'String', 'Stopping...');
end


% --- Executes on button press in ExportButton.
function ExportButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = get(handles.DPP, 'UserData');

if (isfield(Data, 'Posterior'))
    try
        [FileName,PathName,FilterIndex] = uiputfile('*.mat', 'Export inference results');
        
        if (FilterIndex == 0)
            return;
        end
        
        Stat = GetParameterStatistics(Data.Posterior, Data.Model);
        
        save([PathName FileName], 'Stat');
        
        fprintf('Successfully exported results to %s!\n', [PathName FileName]);
    catch err
        fprintf('An error occurred while exporting inference results.\n');
    end
else
    fprintf('No inference performed!\n');
end
