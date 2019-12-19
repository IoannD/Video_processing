%% ������ � ����� ��������� �� �����
% ��� ������ ����� ���������� � ��������� ����� ����� : 
% VAR.str = 'prefix';
% get_gray_image - �������� ���� � ����������� �������, ��� �������
% ��� ������ ��������� ������ ����� � ����������
function image_gray = get_gray_image(k, VAR) 
    num = num2str(k-1,'%i');    %����������� ����� �����
    full_str = [VAR.str num '.bmp'];
    image = imread(full_str);   %�������� �����������
    image_gray = rgb2gray(image);
end