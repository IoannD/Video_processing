clear
clc
close all
%% ИНИЦИАЛИЗАЦИЯ ПРОГРАММЫ 
VAR.FileName = 'default name'; 
VAR.frame_period = 1; % длительность по умолчанию
VAR.pixel_size = 1; 
VAR.size_for_average = 3; % размер квадратного усредняющего фильтра 

% бинаризация 
VAR.Sens = 0.9; %sensitivity v binarizacii
VAR.gamma = 5; %parameter in imadjust
% VAR.k = 1; % счетчик обработки
image_binor =zeros(10,10); 
VAR.resize_parameter = 50; % величина расширения промоугольного окна
% в котором помещается капля - необходимо в функцие "get_resized_image"
X = 0; Y = 0; 

VAR.frame_period = 1/1000; % длительность кадра 
VAR.pixel_size = 8/791; %pixels per mm по умолчанию 

%% ПОТОКОВАЯ ОБРАБОТКА ВСЕГО ВИДЕО 
% файл для разработки и отладки функции получения кадров из avi-файлаclc


% ТУТ ПИСАТЬ АДРЕСС ВИДЕО ДЛЯ ОБРАБОТКИ
VAR.FileName = 'C:\Users\idobr\YandexDisk\3 НИР\7 семестр\12 Работа с видео\Эксперимент\4.avi';

% ------------------------------------------------------------------------
Video = VideoReader( VAR.FileName); 
k=1; 
X = 0; Y = 0; 

VAR = get_gamma_corr(VAR);

while  hasFrame(Video) 

    image = readFrame(Video);
    image_gray = rgb2gray(image);
    
    [resized_image,X,Y] = get_resized_image(image_gray, image_binor, X,Y, VAR); 
    image_gray = resized_image;
    
    % aлгоритм 1 - два раза применяется усредняющий фильтр
    after_filtr = average_filter(image_gray,VAR);
    after_filtr = average_filter(after_filtr,VAR);
    % алгорит 2 - без усреднения изображения
%     after_filtr = image_gray;

    image_binor = binarization(after_filtr, VAR);

    image_binor = clear_image(image_binor);

    boundary = plot_doundary(image_binor, image_gray, 1);
    % ЕСЛИ НУЖЕН ВЫВОД ГРАНИЦЫ НА КАЖДОМ КАДРЕ, ТО ЗАМЕНИТЬ 0 НА 1

    [half_axis_big, half_axis_small] = get_axis(image_binor);

        % условие на неразрывность деформации
        if half_axis_small>0
                hor_ax_small(k) = half_axis_small;
                vert_ax_big(k) = half_axis_big;
        else
                hor_ax_small(k) = hor_ax_small(k-1);
                vert_ax_big(k) = vert_ax_big(k-1);
        end
                VAR.volume(k) = 4/3*pi*half_axis_small^2*half_axis_big*VAR.pixel_size^3;

    k= k+1; 
end 

%%
VAR = get_deformation(hor_ax_small, vert_ax_big, VAR); 

VAR = get_radius(VAR); 

plot_deformation(VAR)
plot_volume(VAR)

100*VAR.Error_Radius/VAR.Radius
figure
plot(VAR.time, VAR.R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ФУНКЦИИ
%% АВТОМАТИЧЕСКОЕ ОПРЕДЕЛЕНИЕ ПАРАМЕТРОВ ГАММА-КОРРЕКЦИИ
function VAR = get_gamma_corr(VAR)
image_first = get_frame(VAR, 'first');
VAR.gamma_correction = [0 0];
n = 2;  
Idouble = im2double(image_first); 
avg = mean2(Idouble); 
sigma = std2(Idouble);
VAR.gamma_correction = [avg-n*sigma avg+n*sigma];
end
%% ФУНЦКЦИЯ ДЛЯ ПОЛУЧЕНИЕ ПЕРВОГО ИЛИ ПОСЛЕДНЕГО КАДРА ИЗ ВИДЕО
% ПРИМЕР ВЫЗОВА: 
% image_last = get_frame(VAR, 'last')
% image_last = get_frame(VAR, 'last')
function image = get_frame(VAR, string) 
    Video = VideoReader( VAR.FileName ); % получение с указателя на видео
    % Video - содержит также и информацю о самом видео: длительность,
    % частота кадров, длительность и тд
    if strcmp(string , 'first') 
        image = readFrame(Video); % получение первого кадра
        
    elseif strcmp(string, 'last') 
        time = Video.Duration - 1/Video.FrameRate; 
        Video = VideoReader(VAR.FileName, 'CurrentTime', time); % создается указатель 
        % на видео с длительностью в один, последний кадр
        image = readFrame(Video);% получение последнего кадра
        
    else % блок вывода ошибок 
        disp('Error: function get_frame have wrong input argument')
        disp('Argument must be "first" or "last"')
    end
end
%% ПРИМЕР ВЫЗОВА: VAR.frame_period = get_frame_period(VAR.FileName)
function [frame_period] = get_frame_period(FileName)
    % функция принимается на вход название файла с видео  
    % пример 'VideoName.avi'
    % strcat - Concatenate strings horizontally
    ending = ('.info');
    
    file_for_read =strcat( FileName,ending); % образуется имя файла .info
   
    FileID = fopen( file_for_read, 'r'); % открываем для чтения файл
    formatSpec = '%s %f %s'; % подготовка атрибута чтения из файла 
    %   %s - Read as a cell array of character vectors.
    %   %f - Floating-point number	
    cell_out = textscan(FileID, formatSpec); % чтения из файла
    
    [M, N] = size(cell_out{1,1}); 
    i = 1;
    % strcmp (string1, string2) - compare 2 strings
    while strcmp(cell2mat(cell_out{1,1}(i)) , 'Период:') % поиск номера ячейки
        % в которой указан период 
        i = i + 1; 
    end
    
    time = cell_out{1,2}(i);    % запись значения из cell величины периоды
    prefix = cell_out{1,3}(i);  % запись приставки СИ для времени в переменную
    
    % поиск совпадения приставки с базовыми приставками СИ 
    % и определение множетеля - factor 
    if strcmp(prefix , 'мкс')
        factor = 1e-6; 
    elseif strcmp(prefix , 'с') 
        factor = 1; 
    elseif strcmp(prefix , 'мс')
        factor = 1e-3; 
    elseif strcmp(prefix , 'нс')
        factor = 1e-9; 
    else 
        disp('Error - did not found the time prefics in video info-file')
        disp('Check the prefics or function "get_frame_period"')
    end 
    
    % вычисление длительность одного кадра 
    frame_period = time * factor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_gray_image - получает кадр с необходимым номером, при условии
% что задано текстовое начало файла и директория
function image_gray = get_gray_image(k, VAR) 
    num = num2str(k-1,'%i');    %определение имени файла
    full_str = [VAR.str num '.bmp'];
    image = imread(full_str);   %открытие изображения
    image_gray = rgb2gray(image);
end

% прямоугольный усредняющий фильтр
function image = average_filter(image,VAR)
    w = fspecial('average', [VAR.size_for_average]);
    % 'average' fspecial('average', [r c]). Прямоугольный усредняющий 
    % фильтр размера rxc. По умолчанию 3x3. Одно число на месте [r c] означает квадратный фильтр
    image = imfilter(image,w,'replicate');
    % replicate - Размер изображения увеличивается повторением величин на его боковых границах.
end

% binarization - преобразование с адаптивным порогом с применением 
% преобразования яркости с гамма коррекцией
function image_binar = binarization(image, VAR) 
    b=im2double(image);   
    
    %преобразование интенсивностей
%         b = imadjust(b,VAR.gamma_correction,[]); % автоматический выбор настройки
     b=imadjust(b,[0 0.5],[],VAR.gamma); % ручная настройка
    
    image_binar = imbinarize(b,'adaptive', 'Sensitivity', VAR.Sens, 'ForegroundPolarity','bright');
    
end 

function image_binor = clear_image(image) 

    BW1 = edge(image,'Canny',0.15, 4); %выделение границ
    
    BW1 = bwmorph(BW1, 'bridge'); %связывание границ 
    
    %узнать подробнее
    bw = imfill(BW1,'holes'); % заполнение области

    se = strel('disk',35);% стирание
    bw2 = imopen(bw,se); % лишнего шума
    image_binor = bw2;  
end

function [half_axis_big, half_axis_small] = get_axis(image_binor)
    stats = regionprops('table',image_binor,'Centroid',...
    'MajorAxisLength','MinorAxisLength'); %определение свойств капли
    half_axis_big = stats.MajorAxisLength/2;
    half_axis_small = stats.MinorAxisLength/2;
end  

function boundary = plot_doundary(image_binor, image_gray, flag)
% flag == 0 - do not plot 
% flag == 1 - plot 
    boundary = bwperim(image_binor,8);
    if flag == 1     
    imshowpair(boundary, image_gray, 'falsecolor');
    % 'falsecolor'	Creates a composite RGB image showing A and B 
    ... overlaid in different color bands. Gray regions in the 
    ...composite image show where the two images have the same 
    ...intensities. Magenta and green regions show where the 
    ...intensities are different. This is the default method.
    end  
end

function [VAR] = get_deformation(hor_ax_small, vert_ax_big, VAR)
    l = vert_ax_big*VAR.pixel_size; 
    d = hor_ax_small*VAR.pixel_size;
%     VAR.R = (l+d)/2; %расчет радиуса капли
%     Volume = d^2*l;
    VAR.R = (d.^2.*l).^(1/3);
    
    VAR.D = (l-d)./(l+d)*100; %расчет деформации
    VAR.time = linspace(VAR.frame_period,length(VAR.D)*VAR.frame_period,length(VAR.D)); 
    %время деформации
end

function [] = plot_deformation(VAR)
    % Построение графика деформации от времени
    figure1 = figure('Color',[1 1 1]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    ylabel('D, %');
    xlabel('Time, s');
    grid(axes1,'on');
    set(axes1,'FontSize',14);
    hold on
    plot(VAR.time,VAR.D,'.-')
end

function [] = plot_volume(VAR)
    % Построение графика выполнения закона сохранения
    % массы от времени
    figure1 = figure('Color',[1 1 1]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    ylabel('Volume, %');
    xlabel('Time, s');
    grid(axes1,'on');
    set(axes1,'FontSize',14);
    hold on
    
    referenceVolume = (4/3)*pi*VAR.Radius^3; 
    relativeVolume = (referenceVolume - VAR.volume)/referenceVolume;
    plot(VAR.time,relativeVolume,'.-')
end

function VAR = get_radius(VAR)

    k = 1;
    meanRadius = 0;
    while (VAR.D(k) < 0.4) && (k < length(VAR.R))
        meanRadius = meanRadius + VAR.R(k); 
        k = k + 1; 
    end
    VAR.Radius = meanRadius / k; 
    VAR.Error_Radius = std(VAR.R)/sqrt(k)  
end
%% ФУНКЦИЯ МАСШТАБИРОВАНИЯ ПОЛУТОНОВОГО ИЗОБРАЖЕНИЯ 
function [resized_image,X,Y] = get_resized_image(image_gray, image_binor, X,Y, VAR) 
%     k = VAR.k;
    [N_row ,N_column] = size(image_gray);
    % Column - столбец - X-coordinate 
    % Row - строка - Y-coordinate
    
    if (X~=0)     
        s = regionprops(image_binor, 'BoundingBox');
        % s.BoundingBox(1)- абцисса левого верхнего
        ... конца ограничивающего прямоугольника
        % s.BoundingBox(2) - ордината левого конца ограничивающего прямоугольника
        % s.BoundingBox(3) - width ограничивающего прямоугольника
        % s.BoundingBox(4) - height ограничивающего прямоугольника
        % fix(a) - округляет до ближайшего целого числа в сторону нуля
        X_left = fix(s.BoundingBox(1)) + X - VAR.resize_parameter; 
        Y_left = fix(s.BoundingBox(2)) + Y - VAR.resize_parameter;
            
        wight = s.BoundingBox(3);
        height = s.BoundingBox(4);
            
        X_r = fix(s.BoundingBox(1)) + X + wight + VAR.resize_parameter;
        Y_r = fix(s.BoundingBox(2)) + Y + height + VAR.resize_parameter; 

        if (X_left < 0) || (X_left == 1) 
            X_left = 1; 
        end
        
        if (Y_left < 0) || (Y_left == 1)
            Y_left = 1;     
        end
        
        if (X_r > N_column) || (X_r == N_column)
            X_r = N_column;    
        end
        
        
        if (Y_r > N_row) || (Y_r == N_row)
            Y_r = N_row;    
        end
        
        resized_image = image_gray(Y_left:Y_r, X_left:X_r);
        X = X_left;
        Y = Y_left;    
%             imshow(resized_image)   
    end
        
    if (X==0)&&(Y==0) % выполняется при первом вызове функции get_resized_image
        X = 1;
        Y = 1;
        resized_image = image_gray;
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% БЛОК ФУНКЦИЙ ДЛЯ АППРОКСИМАЦИИ ЭЛЛИПСОМ

function [X,Y] = get_boundary_coordinates(boundary)
    
    [M, N] = size(boundary); 
    t = 1;
    for i=1:M
        for j=1:N
%             if boudnary(i,j) ~= 0 
            if boundary(i , j) ~=0
                X(t) = j; 
                Y(t) = i;
                t = t+1; 
            end
        end
    end
    X = X';
    Y = Y'; 
 end
 
function ell = detect_ellipse(x, y)
    % получает на вход x и y - векторы строки строки с координатами точек, 
    % которые аппроксимируются
    x2 = x.*x; y2 = y.*y; xy = x.*y; one = ones(length(x),1);
    data = [x2 y2 xy x y one];
    X = slau_mnk(data);
    a11 = X(1); a22 = X(2); a12 = X(3)/2; a13 = X(4)/2; a23 = X(5)/2; a33 = 1;
    M = [a11 a12 a13; a12 a22 a23; a13 a23 a33];
    d = det(M);
    D = a11*a22 - a12*a12;
    I = a11 + a22;

    Mo.X = (a12*a23 - a13*a22)/D;
    Mo.Y = (a13*a12 - a11*a23)/D;
    tg = 2*a12/(a11-a22);
    k = (sqrt(tg*tg+1)-1)/tg;
    alpha = atan(k);

    [issol lamda1 lamda2] = quadratic_equation(1, -I, D);

    a = sqrt(-d/(lamda1*lamda2*lamda2));
    b = sqrt(-d/(lamda1*lamda1*lamda2));
    c = sqrt(abs(a*a-b*b));

    if a < b 
        k = -1/k; 
    end
    p = c/sqrt(1+k*k);
    F1.X = Mo.X - p; F1.Y = Mo.Y - k*p;
    F2.X = Mo.X + p; F2.Y = Mo.Y + k*p;

    ell.Mo = Mo; ell.F1 = F1; ell.F2 = F2;
    ell.alpha = alpha; ell.alphagrad = 180/pi*alpha;
    ell.a = a; ell.b = b; ell.c = c; ell.k = k; ell.M = M;
end

function x = slau_mnk(M)
    [rows cols] = size(M); Z = [];
    if cols-1 < rows
        for col = 1:cols-1
            N = [];
            for i = 1:rows
                N = [N; M(i,col)*M(i, : )];
            end
            Z = [Z; sum(N)];
        end
    end
    b = -Z(:,end);
    A = Z(1:end,1:end-1);
    x = A\b; % решение приведенной системы средствами матлаба
end

function [issol x1 x2] = quadratic_equation(A, B, C)
    D = B*B - 4*A*C; 
    if D < 0 issol = 0; x1 = 0; x2 = 0; 
        return; 
    else issol = 1; 
    end
    x1 = (-B + sqrt(D))/(2*A); x2 = (-B - sqrt(D))/(2*A);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
