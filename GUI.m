function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 29-Mar-2018 15:36:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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

% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GUI.
%if strcmp(get(hObject,'Visible'),'off')
%    plot(rand(5));
%end
set(handles.axes1, 'Visible', 'off');


% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in gen_tr_set_btn.
function gen_tr_set_btn_Callback(hObject, eventdata, handles)
% hObject    handle to gen_tr_set_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = str2double(get(handles.h, 'String')); % шаг дискретизации области ответственности 
N = str2double(get(handles.N, 'String')); % количество приемных пунктов

Xmin = str2double(get(handles.Xmin, 'String'));
Xmax = str2double(get(handles.Xmax, 'String'));

Ymin = str2double(get(handles.Ymin, 'String'));
Ymax = str2double(get(handles.Ymax, 'String'));

Zmin = str2double(get(handles.Zmin, 'String'));
Zmax = str2double(get(handles.Zmax, 'String'));

R = str2double(get(handles.mrs_radius, 'String'));   % радиус базы

sigma_tau = 25 * 10^(-9);   % СКО измерения времени задержки
c = 3 * 10^5;               % скорость света
Rerror = c * sigma_tau / sqrt(2);

% границы области ответственности МРС
D = [Xmin Xmax;
    Ymin Ymax;
    Zmin Zmax];

% размеры области ответственности
Dx = D(1,2) - D(1,1);
Dy = D(2,2) - D(2,1);
Dz = D(3,2) - D(2,1);

% количество интервалов разбиения
Nx = round(Dx / h);
Ny = round(Dy / h);
Nz = round(Dz / h);

method = get(handles.method, 'Value');

% центр области ответственности
Xc = (Xmax - Xmin) / 2;     
Yc = (Ymax - Ymin) / 2;

if (method == 1)
    % Разностно - дальномерный метод
    
    % задаем координаты приемных пунктов на окружности
    PP_Coord(:, 1) = [Xc Yc 0];
    
    for i = 1 : N - 1
        X = Xc + R * cos(2 * pi * i / (N - 1));
        Y = Yc + R * sin(2 * pi * i / (N - 1));
        Z = 0.04 * (Zmax - Zmin) * rand(1);
        PP_Coord(:, i + 1) = [X Y Z];
    end
    
    % записываем нужную информацию в файлы
    dlmwrite('D_daln.dat', D, ' ');
    dlmwrite('PP_Coord_daln.dat', PP_Coord, ' ');
    
    % генерируем обучающую выборку по равномерному закону
    t = 1;
    for i = 1 : Nx
        X = i * h;
        for j = 1 : Ny
            Y = j * h;
            for k = 1 : Nz
                Z = k * h;
                % расстояние до ведущего ПП
                R1 = Rasst(X, Y, Z, PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1)) + Rerror * rand(1);
                m = 1;
                for n = 2 : N
                    R = Rasst(X, Y, Z, PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n)) + Rerror * rand(1);
                    r = Rasst(PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1), PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n));
                    Batch_X(m) = R1 + r - R;
                    m = m + 1;
                end
                TrSet_X(t, :) = Batch_X;
                TrSet_Y(t, :) = [X, Y, Z];
                t = t + 1;
            end
        end
    end

    dlmwrite('TrSet_X_daln.dat', TrSet_X', ' ');
    dlmwrite('TrSet_Y_daln.dat', TrSet_Y', ' ');
else
    if (method == 2)
        % Угломерный метод

        % задаем координаты приемных пунктов
        for i = 1 : N
            X = (i - 1) * (Xmax - Xmin) / (N - 1);
            Y = 0;
            Z = 0.03 * Zmax * rand(1);
            PP_Coord(:, i) = [X Y Z];
        end

        % записываем нужную информацию в файлы
        dlmwrite('D_ugl.dat', D, ' ');
        dlmwrite('PP_Coord_ugl.dat', PP_Coord, ' ');

        skoip = 1; % СКО измерения пеленга, град

        skoias = skoip; % СКО измерения азимута, град
        skoium = skoip; % СКО измерения угла места, град

        % перевод СКО измерения азимута из градусов в радианы
        skoias_r = skoias * (pi / 180); 
        % перевод СКО измерения угла места из градусов в радианы
        skoium_r = skoium*(pi / 180);

        % генерируем обучающую выборку по равномерному закону
        t = 1;
        for i = 1 : Nx
            X = i * h;
            for j = 1 : Ny
                Y = j * h;
                for k = 1 : Nz
                    Z = k * h;
                    m = 1;
                    for n = 1 : N
                        % генерация текущих измеренных значений азимута и угла места
                        as = Asimut(PP_Coord(1, n), PP_Coord(2, n), X, Y) + skoias_r * randn(1);
                        um = Ugmest(PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n), X, Y, Z) + skoium_r * randn(1);
                        Batch_X(m, :) = [sin(as) cos(as) sin(um) cos(um)];
                        m = m + 1;
                    end
                    TrSet_X(t, :) = reshape(Batch_X, 4 * length(Batch_X), 1);
                    TrSet_Y(t, :) = [X, Y, Z];
                    t = t + 1;
                end
            end
        end

        dlmwrite('TrSet_X_ugl.dat', TrSet_X', ' ');
        dlmwrite('TrSet_Y_ugl.dat', TrSet_Y', ' ');
    else
        % угломерно-разностно-дальномерный метод
        
        % задаем координаты приемных пунктов
        for i = 1 : N
            X = (i - 1) * (Xmax - Xmin) / (N - 1);
            Y = 0;
            Z = 0.03 * Zmax * rand(1);
            PP_Coord(:, i) = [X Y Z];
        end

        % записываем нужную информацию в файлы
        dlmwrite('D_comb.dat', D, ' ');
        dlmwrite('PP_Coord_comb.dat', PP_Coord, ' ');

        skoip = 1; % СКО измерения пеленга, град

        skoias = skoip; % СКО измерения азимута, град
        skoium = skoip; % СКО измерения угла места, град

        % перевод СКО измерения азимута из градусов в радианы
        skoias_r = skoias * (pi / 180); 
        % перевод СКО измерения угла места из градусов в радианы
        skoium_r = skoium*(pi / 180);

        % генерируем обучающую выборку по равномерному закону
        t = 1;
        for i = 1 : Nx
            X = i * h;
            for j = 1 : Ny
                Y = j * h;
                for k = 1 : Nz
                    Z = k * h;
                    R1 = Rasst(X, Y, Z, PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1)) + Rerror * rand(1);
                    m = 1;
                    for n = 2 : N
                        R = Rasst(X, Y, Z, PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n)) + Rerror * rand(1);
                        r = Rasst(PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1), PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n));
                        Batch_X1(m) = R1 + r - R;
                        m = m + 1;
                    end
                    m = 1;
                    for n = 1 : N
                        % генерация текущих измеренных значений азимута и угла места
                        as = Asimut(PP_Coord(1, n), PP_Coord(2, n), X, Y) + skoias_r * randn(1);
                        um = Ugmest(PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n), X, Y, Z) + skoium_r * randn(1);
                        Batch_X2(m, :) = [sin(as) cos(as) sin(um) cos(um)];
                        m = m + 1;
                    end
                    Batch = [reshape(Batch_X2, 1, 4 * length(Batch_X2)), Batch_X1];
                    TrSet_X(t, :) = Batch;
                    TrSet_Y(t, :) = [X, Y, Z];
                    t = t + 1;
                end
            end
        end

        dlmwrite('TrSet_X_comb.dat', TrSet_X', ' ');
        dlmwrite('TrSet_Y_comb.dat', TrSet_Y', ' ');
        
    end
end

h = msgbox('Выборка сгенерирована!');

function h_Callback(hObject, eventdata, handles)
% hObject    handle to h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h as text
%        str2double(get(hObject,'String')) returns contents of h as a double


% --- Executes during object creation, after setting all properties.
function h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ymin_Callback(hObject, eventdata, handles)
% hObject    handle to Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ymin as text
%        str2double(get(hObject,'String')) returns contents of Ymin as a double


% --- Executes during object creation, after setting all properties.
function Ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ymax_Callback(hObject, eventdata, handles)
% hObject    handle to Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ymax as text
%        str2double(get(hObject,'String')) returns contents of Ymax as a double


% --- Executes during object creation, after setting all properties.
function Ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xmax_Callback(hObject, eventdata, handles)
% hObject    handle to Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xmax as text
%        str2double(get(hObject,'String')) returns contents of Xmax as a double


% --- Executes during object creation, after setting all properties.
function Xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xmin_Callback(hObject, eventdata, handles)
% hObject    handle to Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xmin as text
%        str2double(get(hObject,'String')) returns contents of Xmin as a double


% --- Executes during object creation, after setting all properties.
function Xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zmin_Callback(hObject, eventdata, handles)
% hObject    handle to Zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zmin as text
%        str2double(get(hObject,'String')) returns contents of Zmin as a double


% --- Executes during object creation, after setting all properties.
function Zmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zmax_Callback(hObject, eventdata, handles)
% hObject    handle to Zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zmax as text
%        str2double(get(hObject,'String')) returns contents of Zmax as a double


% --- Executes during object creation, after setting all properties.
function Zmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over gen_tr_set_btn.
function gen_tr_set_btn_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to gen_tr_set_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Test_N_Callback(hObject, eventdata, handles)
% hObject    handle to Test_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Test_N as text
%        str2double(get(hObject,'String')) returns contents of Test_N as a double


% --- Executes during object creation, after setting all properties.
function Test_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Test_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gen_test_set_btn.
function gen_test_set_btn_Callback(hObject, eventdata, handles)
% hObject    handle to gen_test_set_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N = str2double(get(handles.Test_N, 'String'));     % количество точек тестовой выборки
sigmaMin = 1 * 10^(-6);     % 1 мкс
sigmaMax = 5.5 * 10^(-6);   % 5.5 мкс
h_tau = 0.5 * 10^(-6);    % шаг 1 мкс
c = 3 * 10^5;           % скорость света


one_point_test = get(handles.one_point, 'Value');
method = get(handles.method, 'Value');

if (method == 1)
    % Разностно - дальномерный метод
    
    D = dlmread('D_daln.dat');   % границы области ответственности
    PP_Coord = dlmread('PP_Coord_daln.dat'); % координаты ПП
    
    if (one_point_test)
        % центр области ответственности
        X = D(1, 2) / 2;
        Y = D(2, 2) / 2;
        Z = D(3, 2) / 2;
        M = round((sigmaMax - sigmaMin) / h_tau);
        K = 100;     % количество точек на каждую ошибку

        t = 1;
        for i = 1 : M + 1
            for k = 1 : K
                Rerror = c * i * h_tau / sqrt(2);
                % расстояние до ведущего ПП
                R1 = Rasst(X, Y, Z, PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1)) + Rerror * rand(1);
                m = 1;
                for n = 2 : length(PP_Coord)
                    R = Rasst(X, Y, Z, PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n)) + Rerror * rand(1);
                    r = Rasst(PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1), PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n));
                    Batch_X(m) = R1 + r - R;
                    m = m + 1;
                end
                TestSet_X(t, :) = Batch_X;
                TestSet_Y(t, :) = [X, Y, Z];
                t = t + 1;
            end
            sigma(i) = i * h_tau;
        end
        dlmwrite('TestSet_X_daln.dat', TestSet_X', ' ');
        dlmwrite('TestSet_Y_daln.dat', TestSet_Y', ' ');
        dlmwrite('sigma_daln.dat', sigma, ' ');
    else
        Rerror = c * 25 * 10^(-9) / sqrt(2);
        t = 1;
        for i = 1 : N
            X = D(1, 2) * rand(1);
            Y = D(2, 2) * rand(1);
            Z = 0.04 * D(3, 2) + D(3, 2) * rand(1);
            % расстояние до ведущего ПП
            R1 = Rasst(X, Y, Z, PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1)) + Rerror * rand(1);
            m = 1;
            for n = 2 : length(PP_Coord)
                R = Rasst(X, Y, Z, PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n)) + Rerror * rand(1);
                r = Rasst(PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1), PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n));
                Batch_X(m) = R1 + r - R;
                m = m + 1;
            end
            TestSet_X(t, :) = Batch_X;
            TestSet_Y(t, :) = [X, Y, Z];
            t = t + 1;
        end

        dlmwrite('TestSet_X_daln.dat', TestSet_X', ' ');
        dlmwrite('TestSet_Y_daln.dat', TestSet_Y', ' ');
    end
else
    if (method == 2)
        
        % угломерный метод
        
        D = dlmread('D_ugl.dat');   % границы области ответственности
        PP_Coord = dlmread('PP_Coord_ugl.dat'); % координаты ПП
        
        if (one_point_test)
            % центр области ответственности
            X = D(1, 2) / 2;
            Y = D(2, 2) / 2;
            Z = D(3, 2) / 2;
            K = 100;     % количество точек на каждую ошибку

            skoip = 1; % СКО измерения пеленга, град

            t = 1;
            for i = skoip : 10

                skoias = i; % СКО измерения азимута, град
                skoium = i; % СКО измерения угла места, град

                % перевод СКО измерения азимута из градусов в радианы
                skoias_r = skoias * (pi / 180); 
                % перевод СКО измерения угла места из градусов в радианы
                skoium_r = skoium * (pi / 180);

                for k = 1 : K
                    m = 1;
                    for n = 1 : length(PP_Coord)
                        % генерация текущих измеренных значений азимута и угла места
                        as = Asimut(PP_Coord(1, n), PP_Coord(2, n), X, Y) + skoias_r * randn(1);
                        um = Ugmest(PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n), X, Y, Z) + skoium_r * randn(1);
                        Batch_X(m, :) = [sin(as) cos(as) sin(um) cos(um)];
                        m = m + 1;
                    end
                    TestSet_X(t, :) = reshape(Batch_X, 4 * length(Batch_X), 1);
                    TestSet_Y(t, :) = [X, Y, Z];
                    t = t + 1;
                end
                sigma(i) = i;
            end

            dlmwrite('TestSet_X_ugl.dat', TestSet_X', ' ');
            dlmwrite('TestSet_Y_ugl.dat', TestSet_Y', ' ');
            dlmwrite('sigma_ugl.dat', sigma, ' ');
        else
            skoip = 1; % СКО измерения пеленга, град

            skoias = skoip; % СКО измерения азимута, град
            skoium = skoip; % СКО измерения угла места, град

            % перевод СКО измерения азимута из градусов в радианы
            skoias_r = skoias * (pi / 180); 
            % перевод СКО измерения угла места из градусов в радианы
            skoium_r = skoium*(pi / 180);
            t = 1;
            for i = 1 : N
                X = D(1, 2) * rand(1);
                Y = D(2, 2) * rand(1);
                Z = D(3, 2) * rand(1);
                m = 1;
                for n = 1 : length(PP_Coord)
                    % генерация текущих измеренных значений азимута и угла места
                    as = Asimut(PP_Coord(1, n), PP_Coord(2, n), X, Y) + skoias_r * randn(1);
                    um = Ugmest(PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n), X, Y, Z) + skoium_r * randn(1);
                    Batch_X(m, :) = [sin(as) cos(as) sin(um) cos(um)];
                    m = m + 1;
                end

                TestSet_X(t, :) = reshape(Batch_X, 4 * length(Batch_X), 1);
                TestSet_Y(t, :) = [X, Y, Z];
                t = t + 1;
            end

            dlmwrite('TestSet_X_ugl.dat', TestSet_X', ' ');
            dlmwrite('TestSet_Y_ugl.dat', TestSet_Y', ' ');
        end
        
    else
        % угломерно-разностно-дальномерный метод
        
        D = dlmread('D_comb.dat');   % границы области ответственности
        PP_Coord = dlmread('PP_Coord_comb.dat'); % координаты ПП
        Rerror = c * 25 * 10^(-9) / sqrt(2);
        
        if(one_point_test)
            % центр области ответственности
            X = D(1, 2) / 2;
            Y = D(2, 2) / 2;
            Z = D(3, 2) / 2;
            K = 100;     % количество точек на каждую ошибку

            skoip = 1; % СКО измерения пеленга, град

            t = 1;
            for i = skoip : 10

                skoias = i; % СКО измерения азимута, град
                skoium = i; % СКО измерения угла места, град

                % перевод СКО измерения азимута из градусов в радианы
                skoias_r = skoias * (pi / 180); 
                % перевод СКО измерения угла места из градусов в радианы
                skoium_r = skoium * (pi / 180);

                for k = 1 : K
                    R1 = Rasst(X, Y, Z, PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1)) + Rerror * rand(1);
                    m = 1;
                    for n = 2 : length(PP_Coord)
                        R = Rasst(X, Y, Z, PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n)) + Rerror * rand(1);
                        r = Rasst(PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1), PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n));
                        Batch_X1(m) = R1 + r - R;
                        m = m + 1;
                    end
                    m = 1;
                    for n = 1 : length(PP_Coord)
                        % генерация текущих измеренных значений азимута и угла места
                        as = Asimut(PP_Coord(1, n), PP_Coord(2, n), X, Y) + skoias_r * randn(1);
                        um = Ugmest(PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n), X, Y, Z) + skoium_r * randn(1);
                        Batch_X2(m, :) = [sin(as) cos(as) sin(um) cos(um)];
                        m = m + 1;
                    end
                    Batch = [reshape(Batch_X2, 1, 4 * length(Batch_X2)), Batch_X1];
                    TestSet_X(t, :) = Batch;
                    TestSet_Y(t, :) = [X, Y, Z];
                    t = t + 1;
                end
                sigma(i) = i;
            end

            dlmwrite('TestSet_X_comb.dat', TestSet_X', ' ');
            dlmwrite('TestSet_Y_comb.dat', TestSet_Y', ' ');
            dlmwrite('sigma_comb.dat', sigma, ' ');
            
        else
            
            
            skoip = 1; % СКО измерения пеленга, град

            skoias = skoip; % СКО измерения азимута, град
            skoium = skoip; % СКО измерения угла места, град

            % перевод СКО измерения азимута из градусов в радианы
            skoias_r = skoias * (pi / 180); 
            % перевод СКО измерения угла места из градусов в радианы
            skoium_r = skoium * (pi / 180);
            
            t = 1;
            for i = 1 : N
                X = D(1, 2) * rand(1);
                Y = D(2, 2) * rand(1);
                Z = D(3, 2) * rand(1);
                
                R1 = Rasst(X, Y, Z, PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1)) + Rerror * rand(1);
                m = 1;
                for n = 2 : length(PP_Coord)
                    R = Rasst(X, Y, Z, PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n)) + Rerror * rand(1);
                    r = Rasst(PP_Coord(1, 1), PP_Coord(2, 1), PP_Coord(3, 1), PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n));
                    Batch_X1(m) = R1 + r - R;
                    m = m + 1;
                end
                m = 1;
                for n = 1 : length(PP_Coord)
                    % генерация текущих измеренных значений азимута и угла места
                    as = Asimut(PP_Coord(1, n), PP_Coord(2, n), X, Y) + skoias_r * randn(1);
                    um = Ugmest(PP_Coord(1, n), PP_Coord(2, n), PP_Coord(3, n), X, Y, Z) + skoium_r * randn(1);
                    Batch_X2(m, :) = [sin(as) cos(as) sin(um) cos(um)];
                    m = m + 1;
                end
                Batch = [reshape(Batch_X2, 1, 4 * length(Batch_X2)), Batch_X1];
                TestSet_X(t, :) = Batch;
                TestSet_Y(t, :) = [X, Y, Z];
                t = t + 1;
            end

            dlmwrite('TestSet_X_comb.dat', TestSet_X', ' ');
            dlmwrite('TestSet_Y_comb.dat', TestSet_Y', ' ');
            
            
        end
    end
end

m = msgbox('Тестовая выборка сгенерирована!');


% --- Executes on button press in train_network_btn.
function train_network_btn_Callback(hObject, eventdata, handles)
% hObject    handle to train_network_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
epochs = 500;    % задаем количество эпох обучения

Num_of_Layers = str2double(get(handles.hidden_layers_num, 'String'));
Num_of_Neurons = str2double(get(handles.num_of_neurons, 'String'));
fromFile = get(handles.trainFromFile, 'Value');

%neural_numbers(1) = 3;

for i = 1 : Num_of_Layers
    neural_numbers(i) = Num_of_Neurons;
end

%neural_numbers(Num_of_Layers + 2) = 10;
for i = 1 : length(neural_numbers)
    activation_funcs(i) = {'tansig'};
end

%neural_numbers(Num_of_Layers + 3) = 10;

%activation_funcs(length(activation_funcs) + 1) = {'purelin'};
%activation_funcs(length(activation_funcs)) = {'purelin'};

method = get(handles.method, 'Value');

if (method == 1)
    % разностно-дальномерный метод
    X = dlmread('TrSet_X_daln.dat');
    Y = dlmread('TrSet_Y_daln.dat');
    filename = 'model_daln';
else
    if (method == 2)
        %угломерный метод
        X = dlmread('TrSet_X_ugl.dat');
        Y = dlmread('TrSet_Y_ugl.dat');
        filename = 'model_ugl';
    else
        %угломерно-разностно-дальномерный метод
        X = dlmread('TrSet_X_comb.dat');
        Y = dlmread('TrSet_Y_comb.dat');
        filename = 'model_comb';
    end
end


if (fromFile)
    if(method == 1)
        load model_daln;
        net = checkpoint.net;
    else
        if(method == 2)
            load model_ugl;
            net = checkpoint.net;
        else
            load model_comb;
            net = checkpoint.net;
        end
    end
else
    net = newff(X, Y, neural_numbers, activation_funcs);
end
net.trainParam.epochs = epochs;
% обучаем нейросеть
%[net, tr] = train(net, X, Y, 'CheckpointFile','MyCheckpoint');
net = train(net, X, Y, 'CheckpointFile', filename);


function hidden_layers_num_Callback(hObject, eventdata, handles)
% hObject    handle to hidden_layers_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hidden_layers_num as text
%        str2double(get(hObject,'String')) returns contents of hidden_layers_num as a double


% --- Executes during object creation, after setting all properties.
function hidden_layers_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hidden_layers_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_of_neurons_Callback(hObject, eventdata, handles)
% hObject    handle to num_of_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_of_neurons as text
%        str2double(get(hObject,'String')) returns contents of num_of_neurons as a double


% --- Executes during object creation, after setting all properties.
function num_of_neurons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_of_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in test_model_btn.
function test_model_btn_Callback(hObject, eventdata, handles)
% hObject    handle to test_model_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

method = get(handles.method, 'Value');

if (method == 1)
    % разностно-дальномерный метод
    testX = dlmread('TestSet_X_daln.dat');
    testY = dlmread('TestSet_Y_daln.dat');
    load model_daln;
    net = checkpoint.net;
    PP_Coord = dlmread('PP_Coord_daln.dat');
    D = dlmread('D_daln.dat');
else
    if (method == 2)
        %угломерный метод
        testX = dlmread('TestSet_X_ugl.dat');
        testY = dlmread('TestSet_Y_ugl.dat');
        load model_ugl;
        net = checkpoint.net;
        PP_Coord = dlmread('PP_Coord_ugl.dat');
        D = dlmread('D_ugl.dat');
    else
        %угломерно-разностно-дальномерный метод
        testX = dlmread('TestSet_X_comb.dat');
        testY = dlmread('TestSet_Y_comb.dat');
        load model_comb;
        net = checkpoint.net;
        PP_Coord = dlmread('PP_Coord_comb.dat');
        D = dlmread('D_comb.dat');
    end
end

C1 = clock;
outputs = net(testX);
C2 = clock;

deltaC = C2(6) - C1(6);

set(handles.all_points_time, 'String', deltaC);
set(handles.one_point_time, 'String', deltaC / length(testX));
set(handles.mean_location_error, 'String', immse(outputs, testY));

one_point_test = get(handles.one_point, 'Value');

if (one_point_test)
    for i = 1 : 10
        m = 1;
        for j = 1 : 100
            outVec = outputs(:, i*j);
            trueVec = testY(:, i*j);
            errors(m) = immse(outVec, trueVec);
            m = m + 1;
        end
        Err(i) = sum(errors) / length(errors);
    end
    dlmwrite('NN_Errors.dat', Err);
end

% отображение результатов
plot3(outputs(1,:),outputs(2,:),outputs(3,:),'b+', testY(1,:),testY(2,:),testY(3,:),'r*', PP_Coord(1,:),PP_Coord(2,:),PP_Coord(3,:),'go'); grid on;
xlabel('X'); xlim([D(1, 1), D(1, 2)]);
ylabel('Y'); ylim([D(2, 1), D(2, 2)]); 
zlabel('Z'); zlim([D(3, 1), D(3,2)]);




function test_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to test_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of test_sigma as text
%        str2double(get(hObject,'String')) returns contents of test_sigma as a double


% --- Executes during object creation, after setting all properties.
function test_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to test_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trainFromFile.
function trainFromFile_Callback(hObject, eventdata, handles)
% hObject    handle to trainFromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trainFromFile


% --- Executes on selection change in method.
function method_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method


% --- Executes during object creation, after setting all properties.
function method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in one_point.
function one_point_Callback(hObject, eventdata, handles)
% hObject    handle to one_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of one_point



function mrs_radius_Callback(hObject, eventdata, handles)
% hObject    handle to mrs_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mrs_radius as text
%        str2double(get(hObject,'String')) returns contents of mrs_radius as a double


% --- Executes during object creation, after setting all properties.
function mrs_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mrs_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
