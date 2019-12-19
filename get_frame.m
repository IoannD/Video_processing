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