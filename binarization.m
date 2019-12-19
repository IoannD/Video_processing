%  binarization - �������������� � ���������� ������� � ����������� 
% �������������� ������� � ����� ����������
function image_binar = binarization(image, VAR) 
    b=im2double(image);   
    
    
%     n = 2;  
%     avg = mean2(b); % Average of matrix elements
%     sigma = std2(b); % Standard deviation of matrix elements
%     VAR.gamma_correction = [avg-n*sigma avg+n*sigma];

    %�������������� ��������������
%     b = imadjust(b,VAR.gamma_correction,[]); % �������������� ����� ���������
%      b=imadjust(b,[0 0.5],[],VAR.gamma); % ������ ���������
     
%      imshow(b)
%      pause(1)
    
    image_binar = imbinarize(b,'adaptive', 'Sensitivity', VAR.Sens, 'ForegroundPolarity','bright');
%     imshow(b)
    %     imshow( image_binar)
%     pause(1)
end 