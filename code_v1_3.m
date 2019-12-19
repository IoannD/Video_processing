clear
clc
close all
%% ������������� ��������� 
VAR.FileName = 'default name'; 
VAR.frame_period = 1; % ������������ �� ���������
VAR.pixel_size = 1; 
VAR.size_for_average = 3; % ������ ����������� ������������ ������� 

% ����������� 
VAR.Sens = 0.9; %sensitivity v binarizacii
VAR.gamma = 5; %parameter in imadjust
% VAR.k = 1; % ������� ���������
image_binor =zeros(10,10); 
VAR.resize_parameter = 50; % �������� ���������� �������������� ����
% � ������� ���������� ����� - ���������� � ������� "get_resized_image"
X = 0; Y = 0; 

VAR.frame_period = 1/1000; % ������������ ����� 
VAR.pixel_size = 8/791; %pixels per mm �� ��������� 

%% ��������� ��������� ����� ����� 
% ���� ��� ���������� � ������� ������� ��������� ������ �� avi-�����clc


% ��� ������ ������ ����� ��� ���������
VAR.FileName = 'C:\Users\idobr\YandexDisk\3 ���\7 �������\12 ������ � �����\�����������\4.avi';

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
    
    % a������� 1 - ��� ���� ����������� ����������� ������
    after_filtr = average_filter(image_gray,VAR);
    after_filtr = average_filter(after_filtr,VAR);
    % ������� 2 - ��� ���������� �����������
%     after_filtr = image_gray;

    image_binor = binarization(after_filtr, VAR);

    image_binor = clear_image(image_binor);

    boundary = plot_doundary(image_binor, image_gray, 1);
    % ���� ����� ����� ������� �� ������ �����, �� �������� 0 �� 1

    [half_axis_big, half_axis_small] = get_axis(image_binor);

        % ������� �� ������������� ����������
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
%% �������
%% �������������� ����������� ���������� �����-���������
function VAR = get_gamma_corr(VAR)
image_first = get_frame(VAR, 'first');
VAR.gamma_correction = [0 0];
n = 2;  
Idouble = im2double(image_first); 
avg = mean2(Idouble); 
sigma = std2(Idouble);
VAR.gamma_correction = [avg-n*sigma avg+n*sigma];
end
%% �������� ��� ��������� ������� ��� ���������� ����� �� �����
% ������ ������: 
% image_last = get_frame(VAR, 'last')
% image_last = get_frame(VAR, 'last')
function image = get_frame(VAR, string) 
    Video = VideoReader( VAR.FileName ); % ��������� � ��������� �� �����
    % Video - �������� ����� � ��������� � ����� �����: ������������,
    % ������� ������, ������������ � ��
    if strcmp(string , 'first') 
        image = readFrame(Video); % ��������� ������� �����
        
    elseif strcmp(string, 'last') 
        time = Video.Duration - 1/Video.FrameRate; 
        Video = VideoReader(VAR.FileName, 'CurrentTime', time); % ��������� ��������� 
        % �� ����� � ������������� � ����, ��������� ����
        image = readFrame(Video);% ��������� ���������� �����
        
    else % ���� ������ ������ 
        disp('Error: function get_frame have wrong input argument')
        disp('Argument must be "first" or "last"')
    end
end
%% ������ ������: VAR.frame_period = get_frame_period(VAR.FileName)
function [frame_period] = get_frame_period(FileName)
    % ������� ����������� �� ���� �������� ����� � �����  
    % ������ 'VideoName.avi'
    % strcat - Concatenate strings horizontally
    ending = ('.info');
    
    file_for_read =strcat( FileName,ending); % ���������� ��� ����� .info
   
    FileID = fopen( file_for_read, 'r'); % ��������� ��� ������ ����
    formatSpec = '%s %f %s'; % ���������� �������� ������ �� ����� 
    %   %s - Read as a cell array of character vectors.
    %   %f - Floating-point number	
    cell_out = textscan(FileID, formatSpec); % ������ �� �����
    
    [M, N] = size(cell_out{1,1}); 
    i = 1;
    % strcmp (string1, string2) - compare 2 strings
    while strcmp(cell2mat(cell_out{1,1}(i)) , '������:') % ����� ������ ������
        % � ������� ������ ������ 
        i = i + 1; 
    end
    
    time = cell_out{1,2}(i);    % ������ �������� �� cell �������� �������
    prefix = cell_out{1,3}(i);  % ������ ��������� �� ��� ������� � ����������
    
    % ����� ���������� ��������� � �������� ����������� �� 
    % � ����������� ��������� - factor 
    if strcmp(prefix , '���')
        factor = 1e-6; 
    elseif strcmp(prefix , '�') 
        factor = 1; 
    elseif strcmp(prefix , '��')
        factor = 1e-3; 
    elseif strcmp(prefix , '��')
        factor = 1e-9; 
    else 
        disp('Error - did not found the time prefics in video info-file')
        disp('Check the prefics or function "get_frame_period"')
    end 
    
    % ���������� ������������ ������ ����� 
    frame_period = time * factor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_gray_image - �������� ���� � ����������� �������, ��� �������
% ��� ������ ��������� ������ ����� � ����������
function image_gray = get_gray_image(k, VAR) 
    num = num2str(k-1,'%i');    %����������� ����� �����
    full_str = [VAR.str num '.bmp'];
    image = imread(full_str);   %�������� �����������
    image_gray = rgb2gray(image);
end

% ������������� ����������� ������
function image = average_filter(image,VAR)
    w = fspecial('average', [VAR.size_for_average]);
    % 'average' fspecial('average', [r c]). ������������� ����������� 
    % ������ ������� rxc. �� ��������� 3x3. ���� ����� �� ����� [r c] �������� ���������� ������
    image = imfilter(image,w,'replicate');
    % replicate - ������ ����������� ������������� ����������� ������� �� ��� ������� ��������.
end

% binarization - �������������� � ���������� ������� � ����������� 
% �������������� ������� � ����� ����������
function image_binar = binarization(image, VAR) 
    b=im2double(image);   
    
    %�������������� ��������������
%         b = imadjust(b,VAR.gamma_correction,[]); % �������������� ����� ���������
     b=imadjust(b,[0 0.5],[],VAR.gamma); % ������ ���������
    
    image_binar = imbinarize(b,'adaptive', 'Sensitivity', VAR.Sens, 'ForegroundPolarity','bright');
    
end 

function image_binor = clear_image(image) 

    BW1 = edge(image,'Canny',0.15, 4); %��������� ������
    
    BW1 = bwmorph(BW1, 'bridge'); %���������� ������ 
    
    %������ ���������
    bw = imfill(BW1,'holes'); % ���������� �������

    se = strel('disk',35);% ��������
    bw2 = imopen(bw,se); % ������� ����
    image_binor = bw2;  
end

function [half_axis_big, half_axis_small] = get_axis(image_binor)
    stats = regionprops('table',image_binor,'Centroid',...
    'MajorAxisLength','MinorAxisLength'); %����������� ������� �����
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
%     VAR.R = (l+d)/2; %������ ������� �����
%     Volume = d^2*l;
    VAR.R = (d.^2.*l).^(1/3);
    
    VAR.D = (l-d)./(l+d)*100; %������ ����������
    VAR.time = linspace(VAR.frame_period,length(VAR.D)*VAR.frame_period,length(VAR.D)); 
    %����� ����������
end

function [] = plot_deformation(VAR)
    % ���������� ������� ���������� �� �������
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
    % ���������� ������� ���������� ������ ����������
    % ����� �� �������
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
%% ������� ��������������� ������������ ����������� 
function [resized_image,X,Y] = get_resized_image(image_gray, image_binor, X,Y, VAR) 
%     k = VAR.k;
    [N_row ,N_column] = size(image_gray);
    % Column - ������� - X-coordinate 
    % Row - ������ - Y-coordinate
    
    if (X~=0)     
        s = regionprops(image_binor, 'BoundingBox');
        % s.BoundingBox(1)- ������� ������ ��������
        ... ����� ��������������� ��������������
        % s.BoundingBox(2) - �������� ������ ����� ��������������� ��������������
        % s.BoundingBox(3) - width ��������������� ��������������
        % s.BoundingBox(4) - height ��������������� ��������������
        % fix(a) - ��������� �� ���������� ������ ����� � ������� ����
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
        
    if (X==0)&&(Y==0) % ����������� ��� ������ ������ ������� get_resized_image
        X = 1;
        Y = 1;
        resized_image = image_gray;
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���� ������� ��� ������������� ��������

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
    % �������� �� ���� x � y - ������� ������ ������ � ������������ �����, 
    % ������� ����������������
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
    x = A\b; % ������� ����������� ������� ���������� �������
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
