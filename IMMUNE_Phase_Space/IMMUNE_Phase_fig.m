function varargout = IMMUNE_Phase_fig(varargin)
% IMMUNE_Phase_fig MATLAB code for IMMUNE_Phase_fig.fig
%      IMMUNE_Phase_fig, by itself, creates a new IMMUNE_Phase_fig or raises the existing
%      singleton*.
%
%      H = IMMUNE_Phase_fig returns the handle to a new IMMUNE_Phase_fig or the handle to
%      the existing singleton*.
%
%      IMMUNE_Phase_fig('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMMUNE_Phase_fig.M with the given input arguments.
%
%      IMMUNE_Phase_fig('Property','Value',...) creates a new IMMUNE_Phase_fig or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IMMUNE_Phase_fig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IMMUNE_Phase_fig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IMMUNE_Phase_fig

% Last Modified by GUIDE v2.5 12-Apr-2018 20:17:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IMMUNE_Phase_fig_OpeningFcn, ...
                   'gui_OutputFcn',  @IMMUNE_Phase_fig_OutputFcn, ...
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


% --- Executes just before IMMUNE_Phase_fig is made visible.
function IMMUNE_Phase_fig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IMMUNE_Phase_fig (see VARARGIN)

% Choose default command line output for IMMUNE_Phase_fig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IMMUNE_Phase_fig wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = IMMUNE_Phase_fig_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%----------------------------------------------------------------------------------------------
% Custom functions
%----------------------------------------------------------------------------------------------
function slider_Update(handles)
C0     = get(handles.C0,    'Value');
Islope = get(handles.Imult, 'Value');
ICmult = get(handles.ICmult,'Value');
Smult  = get(handles.Smult, 'Value');

% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
IMMUNE_Phase(handles)

function save_button_Callback(hObject, eventdata, handles)
IMMUNE_Phase(handles)

function delete_button_Callback(hObject, eventdata, handles)
h  = handles.initialValues;
hS = h.String;
hS = hS(setdiff(1:size(hS,1),h.Value(1)),:);
handles.initialValues.String = hS;

function clear_button_Callback(hObject, eventdata, handles)
handles.reset_button.Value = true;
handles.clear_button.Value = true;
handles.initialValues.String = '';
IMMUNE_Phase(handles)

function initialValues_Callback(hObject, eventdata, handles)
handles.initialValues.UserData = true;
IMMUNE_Phase(handles)
handles.initialValues.UserData = false;
% Hints: contents = cellstr(get(hObject,'String')) returns initialValues contents as cell array
%        contents{get(hObject,'Value')} returns selected item from initialValues

function read_button_Callback(hObject, eventdata, handles)
handles.initialValues.String = char(...
    ' 1 0.00 0.00 0.21 0.00 1.84 0.91 1.00 ',... % no SSc C->0
    ' 2 0.00 0.00 0.21 0.00 0.70 0.91 1.00 ',... %           equil
    ' 3 0.00 0.00 0.21 0.00 0.40 0.91 1.00 ',... %           Inf
    ' 4 0.00 0.00 0.21 1.18 0.60 0.39 1.00 ',... %    SSc C->0: decrease ICmult or Imult, C->Inf, SSc increase
    ' 5 0.00 0.00 0.21 0.86 0.54 0.75 1.00 ' ... %    similar, but with C->equil
    );  

function unusedBox_Callback(hObject, eventdata, handles)

% --- Executes on slider movement.
function C0_Callback    (hObject, eventdata, handles)
slider_Update(handles)
IMMUNE_Phase(handles)

function Imult_Callback(hObject, eventdata, handles)
slider_Update(handles)
IMMUNE_Phase(handles)

function ICmult_Callback(hObject, eventdata, handles)
slider_Update(handles)
IMMUNE_Phase(handles)

function Smult_Callback (hObject, eventdata, handles)
slider_Update(handles)
IMMUNE_Phase(handles)

%----------------------------------------------------------------------------------------------
% Other functions
%----------------------------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C0_CreateFcn     (hObject, eventdata, handles)
function Imult_CreateFcn  (hObject, eventdata, handles)
function C0_txt_CreateFcn (hObject, eventdata, handles)
function figure1_CreateFcn(hObject, eventdata, handles)
function ICmult_CreateFcn (hObject, eventdata, handles)
function Smult_CreateFcn  (hObject, eventdata, handles)
function initialValues_CreateFcn(hObject, eventdata, handles)
