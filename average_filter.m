% прямоугольный усредняющий фильтр
function image = average_filter(image,VAR)
    w = fspecial('average', [VAR.size_for_average]);
    % 'average' fspecial('average', [r c]). Прямоугольный усредняющий 
    % фильтр размера rxc. По умолчанию 3x3. 
    % Одно число на месте [r c] означает квадратный фильтр
    image = imfilter(image,w,'replicate');
    % replicate - Размер изображения увеличивается повторением величин 
    % на его боковых границах.
end