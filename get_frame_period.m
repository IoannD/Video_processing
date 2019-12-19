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