%% ������� ��������������� ������������ ����������� 
function [resized_image,X,Y] = get_resized_image(image_gray, image_binor, X,Y, VAR) 
%     k = VAR.k;
    [N_row ,N_column] = size(image_gray);
    % Column - ������� - X-coordinate 
    % Row - ������ - Y-coordinate
    
    if (X~=0)     
        s = regionprops(image_binor, 'BoundingBox');
        % s.BoundingBox(1)- ������� ������ ��������
        ... ����� ��������������� ��������������
        % s.BoundingBox(2) - �������� ������ ����� ��������������� ��������������
        % s.BoundingBox(3) - width ��������������� ��������������
        % s.BoundingBox(4) - height ��������������� ��������������
        % fix(a) - ��������� �� ���������� ������ ����� � ������� ����
        X_left = fix(s.BoundingBox(1)) + X - VAR.resize_parameter; 
        Y_left = fix(s.BoundingBox(2)) + Y - VAR.resize_parameter;
            
        wight = s.BoundingBox(3);
        height = s.BoundingBox(4);
            
        X_r = fix(s.BoundingBox(1)) + X + wight + VAR.resize_parameter;
        Y_r = fix(s.BoundingBox(2)) + Y + height + VAR.resize_parameter; 

        if (X_left < 0) || (X_left == 1) 
            X_left = 1; 
        end
        
        if (Y_left < 0) || (Y_left == 1)
            Y_left = 1;     
        end
        
        if (X_r > N_column) || (X_r == N_column)
            X_r = N_column;    
        end
        
        
        if (Y_r > N_row) || (Y_r == N_row)
            Y_r = N_row;    
        end
        
        resized_image = image_gray(Y_left:Y_r, X_left:X_r);
        X = X_left;
        Y = Y_left;    
%             imshow(resized_image)   
    end
        
    if (X==0)&&(Y==0) % ����������� ��� ������ ������ ������� get_resized_image
        X = 1;
        Y = 1;
        resized_image = image_gray;
    end  
end