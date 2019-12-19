%% РАБОТА С ВИДЕО РАЗБИТЫМИ НА КАДРЫ
% для вызова нужна переменная с префиксом имени файла : 
% VAR.str = 'prefix';
% get_gray_image - получает кадр с необходимым номером, при условии
% что задано текстовое начало файла и директория
function image_gray = get_gray_image(k, VAR) 
    num = num2str(k-1,'%i');    %определение имени файла
    full_str = [VAR.str num '.bmp'];
    image = imread(full_str);   %открытие изображения
    image_gray = rgb2gray(image);
end