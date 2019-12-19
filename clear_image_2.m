function image_binor = clear_image_2(image,VAR) 
    
% image = edge(image,'Canny',0.15, 4); %��������� ������
image = edge(image,'Canny',VAR.threshold, 4);
% image = edge(image,'Canny',0.15, sqrt(2)); %��������� ������
%     image = edge(image,'log',2.9573e-05) ;
    
    BW1 = bwmorph(image, 'bridge'); %���������� ������ 

    %������ ���������
    image = imfill(image,'holes'); % ���������� �������
  
    max = 0;
    number = 0; 
    cc = bwconncomp(image);
    stats = regionprops(cc,'Area','PixelList');

    for i = 1:length([stats(:).Area])
        if max < stats(i).Area
            max = stats(i).Area;
            number = i; 
        end
    end

%     cc = bwconncomp(image); 

    selected = ismember(labelmatrix(cc), number);
    
    image_binor = selected;  
end