% ������������� ����������� ������
function image = average_filter(image,VAR)
    w = fspecial('average', [VAR.size_for_average]);
    % 'average' fspecial('average', [r c]). ������������� ����������� 
    % ������ ������� rxc. �� ��������� 3x3. 
    % ���� ����� �� ����� [r c] �������� ���������� ������
    image = imfilter(image,w,'replicate');
    % replicate - ������ ����������� ������������� ����������� ������� 
    % �� ��� ������� ��������.
end