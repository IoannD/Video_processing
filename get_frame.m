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