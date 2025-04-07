function varargout = Project(varargin)
% PROJECT MATLAB code for Project.fig
%      PROJECT, by itself, creates a new PROJECT or raises the existing
%      singleton*.
%
%      H = PROJECT returns the handle to a new PROJECT or the handle to
%      the existing singleton*.
%
%      PROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT.M with the given input arguments.
%
%      PROJECT('Property','Value',...) creates a new PROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Project_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Project_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Project

% Last Modified by GUIDE v2.5 21-Jun-2020 22:02:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Project_OpeningFcn, ...
                   'gui_OutputFcn',  @Project_OutputFcn, ...
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


% --- Executes just before Project is made visible.
function Project_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Project (see VARARGIN)

% Choose default command line output for Project
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Project wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Project_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global X
global fs
global Ns
global m
global n

%   Reading Testing file
[X,fs] = audioread('FSA1.wav');
X      = X - mean(X);                   % Remove the mean value to normalize
Ns     = length(X);
[m,n]  = size(X);


% --- Executes on button press in Process.
function Process_Callback(hObject, eventdata, handles)
% hObject    handle to Process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global X
global fs
global Ns
global m
global n
global loge_test
global zcr_test
global frametime
global VUV

%   Loading Training Data
load Dataset

% Initializing vectors
loge_test = [];
zcr_test  = [];
framedur  = 30;                       % Frame duration
shift     = 0.5;                      % Shift rate

framelength = framedur * fs/1000;       % Frame length
shiftlength = shift*framedur * fs/1000; % Shift length
frameNb     = floor(Ns - framelength)/shiftlength; % Number of short-time frames

% Loop to get audio signal from one ear 
if (n==2)
    X1 = X(:,1);
    clear X
    X = X1;
end

%   Extracting Feature vectors
for i=1:16
    Feature(i,:) = [loge_Mean(i,:), zcr_Mean(i,:), loge_Std(i,:), zcr_Std(i,:)];
    Label (i,:) = Threshold(i,:);
end

% Taking Transpose of Feature and Label
Feature = Feature.';
Label = Label.';

% Defining Neural Network
net = feedforwardnet(1);

% Training Neural Network
[nett,tr] = train(net,Feature,Label);

% Processing Testing signal
for k=1:frameNb
    frame       = X((k-1)*shiftlength+1:(k-1)*shiftlength+framelength); % frame
    frametime(k)=((k-1)*shiftlength+1)*1000/fs; % Time index for each frame
    loge_test   = [loge_test 10*log10(sum(frame.^2))/framelength];
    zcr_test    = [zcr_test sum(abs(diff(sign(frame))))/(2*framelength)];
end
     
% Extracting Feature vectors from testing signal 
% Taking Mean and standard deviation  of log of energy
loge_test_Mean = sum(loge_test)/length(loge_test);
loge_test_Std = sqrt(sum(loge_test.^2)/length(loge_test)-loge_test_Mean^2);

% Taking Mean and standard deviation  of zero crossings rate
zcr_test_Mean = sum(zcr_test)/length(zcr_test);
zcr_test_Std = sqrt(sum(zcr_test.^2)/length(zcr_test)-zcr_test_Mean^2);    

% Combining all features
test = [loge_test_Mean, zcr_test_Mean, loge_test_Std, zcr_test_Std];
test = test.';

% predicting output
threshold = floor(nett(test));

% Getting VUV using threshold
for k=1:frameNb
    VUV(k)= zcr_test(k)<(threshold/framelength) & loge_test(k)> -0.02; % if short-time ZCR is lower
end


% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global X
global loge_test
global zcr_test
global nframe
global frametime
global fs
global Ns
global VUV

% Plotting Logarithm of Energy
axes(handles.axes1);
plot(frametime/1000,loge_test,'r','LineWidth',1.5);
title('Logarithm of Energy')

% Plotting Zero crossing rate
axes(handles.axes2);
plot(frametime/1000,zcr_test,'g','LineWidth',1.5);
title('Zero Crossing Rate')

% Plotting speech signal and VUV
axes(handles.axes3);
t = (0 : Ns - 1)/fs;
plot(t,X,'r','LineWidth',1.5);
hold on 
plot(frametime/1000,0.75*VUV)
xlabel('Time(sec)');
legend('speech signal','VUV')
title('Voicing activity in speech signal')