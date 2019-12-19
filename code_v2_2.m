clear
clc
close all

%% ������������� ��������� 
VAR.FileName = 'default name'; 
VAR.frame_period = 1; % the duration of frame
VAR.pixel_size = 1; 
VAR.size_for_average = 3; % ������ ����������� ������������ ������� 
VAR.method = ('automatic');
% VAR.method = ('constant');
VAR.threshold = [0.15, 0.19];

% ����������� 
VAR.Sens = 0.9; % ���������������� � �����������
% VAR.Sens = 0.99; % ���������������� � �����������
VAR.gamma = 5; % �������� � �������������� ������� - imadjust
% VAR.k = 1; % ����o��� ���������
image_binor =zeros(10,10); 
% VAR.resize_parameter = 50; % �������� ���������� �������������� ����
VAR.resize_parameter = 50; % �������� ���������� �������������� ����
% � ������� ���������� ����� - ���������� � "get_resized_image"
X = 0; Y = 0; 


VAR.frame_period = 1/1000; % ������������ ����� 
VAR.pixel_size = 8/791; %pixels per mm �� ��������� 

% ��������� ��������� ����� ����� 
    
% ��� ������ ������ ����� ��� ���������
VAR.FileName = 'C:\Users\idobr\YandexDisk\3 ���\7 �������\12 ������ � �����\�����������\5.avi';
% VAR.FileName = 'C:\Users\idobr\Desktop\��������� �����������\1_1.avi';
% VAR.FileName = 'C:\Users\idobr\Desktop\����� �� ������ ����������\2017_04_04\T_500_def.avi';
% ------------------------------------------------------------------------
Video = VideoReader(VAR.FileName); 
% Video = VideoReader(VAR.FileName, 'CurrentTime', 23);
k=1; 
X = 0; Y = 0; 

VAR = get_gamma_corr(VAR);

while  hasFrame(Video) 


        image = readFrame(Video);

        image_gray = rgb2gray(image);

        [resized_image,X,Y] = get_resized_image(image_gray, image_binor, X,Y, VAR); 
        image_gray = resized_image;

        % a������� 1 - ��� ���� ����������� ����������� ������
    %     after_filtr = average_filter(image_gray,VAR);
    %     after_filtr = average_filter(after_filtr,VAR);
        % ������� 2 - ��� ���������� �����������
        after_filtr = image_gray;

        image_binor = binarization(after_filtr, VAR);

%          image_binor = clear_image(image_binor);
        image_binor = clear_image_2(image_binor,VAR);


            boundary = plot_doundary(image_binor, after_filtr, 1);
        % ���� ����� ����� ������� �� ������ �� ���, �� �������� 0 �� 1


            % �������� 2 - ������������ ��������
        [X_boundary,Y_boundary] = get_boundary_coordinates(boundary);

        ell = detect_ellipse(X_boundary, Y_boundary);
        Deformtion_from_ell_aprox(k) = 100*(ell.a - ell.b)/(ell.a + ell.b);
        Volume_from_ell_aprox(k) = 4/3*pi*ell.b^2*ell.a*VAR.pixel_size^3;

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
hold on
plot(VAR.time, Deformtion_from_ell_aprox)

plot_volume(VAR)

 referenceVolume = (4/3)*pi*VAR.Radius^3; 
hold on 
plot( VAR.time, 100*(-1)*(Volume_from_ell_aprox - referenceVolume)/referenceVolume); 
100*VAR.Error_Radius/VAR.Radius
figure
plot(VAR.time, VAR.R)
figure
plot(VAR.time, hor_ax_small,'b') 
figure
plot(VAR.time,vert_ax_big,'r')

figure
plot(VAR.time, 100*(max(VAR.volume)- VAR.volume)/max(VAR.volume))

