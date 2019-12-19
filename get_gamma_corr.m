%% �������������� ����������� ���������� �����-���������
function VAR = get_gamma_corr(VAR)
% ���� - �����, �������� ������������� ������� ����� � ��������� 
% (average - n*standard_deviation, average + n*standard_deviation)
% ����������� � ���� �������� �������������
image_first = get_frame(VAR, 'first');
VAR.gamma_correction = [0 0];
n = 2;  
Idouble = im2double(image_first); 
avg = mean2(Idouble); % Average of matrix elements
sigma = std2(Idouble); % Standard deviation of matrix elements
VAR.gamma_correction = [avg-n*sigma avg+n*sigma];
end