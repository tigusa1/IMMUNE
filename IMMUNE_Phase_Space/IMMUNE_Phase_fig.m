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

% Last Modified by GUIDE v2.5 23-Nov-2018 13:48:09

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
% function slider_Update(handles)
% C0     = get(handles.C0,    'Value');
% Islope = get(handles.Ilim0, 'Value');
% ICmult = get(handles.Clim0, 'Value');
% Smult  = get(handles.Slim0, 'Value');

% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
IMMUNE_Phase(handles)
handles.reset_button.Value = false;

function save_button_Callback(hObject, eventdata, handles)
IMMUNE_Phase(handles)
handles.save_button.Value = false;

function delete_button_Callback(hObject, eventdata, handles)
h  = handles.initialValues;
hS = h.String;
hS = hS(setdiff(1:size(hS,1),h.Value(1)),:);
handles.initialValues.String = hS;
handles.delete_button.Value = false;

function clear_button_Callback(hObject, eventdata, handles)
handles.reset_button.Value = true;
handles.initialValues.String = '';
IMMUNE_Phase(handles)
handles.reset_button.Value = false;
handles.clear_button.Value = false;

function initialValues_Callback(hObject, eventdata, handles)
handles.initialValues.UserData = true;
IMMUNE_Phase(handles)
handles.initialValues.UserData = false;
% Hints: contents = cellstr(get(hObject,'String')) returns initialValues contents as cell array
%        contents{get(hObject,'Value')} returns selected item from initialValues

function read_button_Callback(hObject, eventdata, handles)
handles.initialValues.String = char(...
    );  
handles.read_button.Value = false;

function unusedBox_Callback(hObject, eventdata, handles)

% --- Executes on slider movement.
function C0_Callback    (hObject, eventdata, handles)
% slider_Update(handles)
IMMUNE_Phase(handles)

function Ilim0_Callback (hObject, eventdata, handles)
% slider_Update(handles)
IMMUNE_Phase(handles)

function Clim0_Callback (hObject, eventdata, handles)
% slider_Update(handles)
IMMUNE_Phase(handles)

function Slim0_Callback (hObject, eventdata, handles)
% slider_Update(handles)
IMMUNE_Phase(handles)

%----------------------------------------------------------------------------------------------
% Other functions
%----------------------------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C0_CreateFcn     (hObject, eventdata, handles)
function Ilim0_CreateFcn  (hObject, eventdata, handles)
function C0_txt_CreateFcn (hObject, eventdata, handles)
function figure1_CreateFcn(hObject, eventdata, handles)
function Clim0_CreateFcn  (hObject, eventdata, handles)
function Slim0_CreateFcn  (hObject, eventdata, handles)
function initialValues_CreateFcn(hObject, eventdata, handles)
