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