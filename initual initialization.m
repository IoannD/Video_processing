%% БЛОК ИНИАЦЛИЗАЦИИ ОКОН ДЛЯ НАСТРОЙКИ ПАРАМЕТРОВ ОБРАБОТКИ

% загрузка двух изображений из файла с видео
image_first = get_frame(VAR, 'first');
image_last = get_frame( VAR, 'last');

% первичная обработка - получения обрезанного изображения 
image_first = take_the_intial_processing(image_first, VAR); % вызов функции
image_last = take_the_intial_processing(image_last, VAR);

%% ЧАСТЬ ИСПОЛНЯЕМАЯ ПОСЛЕ ИЗМЕНЕНИЯ НАСТРОЕК БИНАРИЗАЦИИ
    % скорее всего нужно переписывать функцию вывода в окно в приложении
    % первое изображение
    image_binor_first = binarization(image_first, VAR);
    image_binor_first = clear_image(image_binor_first);
    boundary = plot_doundary(image_binor_first, image_gray_first, 1);
    % второе изображение
    image_binor_last = binarization(image_last, VAR);
    image_binor_last = clear_image(image_binor_last);
    boundary = plot_doundary(image_binor_last, image_gray_last, 1);   
