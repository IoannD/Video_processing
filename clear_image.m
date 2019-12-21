function image_binor = clear_image(image) 

    
    BW1 = edge(image,'Canny',0.15, 4); %��������� ������
%     BW1 = edge(image,'Canny',0.1, 6); %��������� ������
%     imshow(BW1)
%     pause(1)
    
    BW1 = bwmorph(BW1, 'bridge'); %���������� ������ 
%         imshow(BW1)
%     pause(1)
    %������ ���������
    bw = imfill(BW1,'holes'); % ���������� �������
%     imshow(bw)
%     pause(1) 
    
    max = 0;
    number = 0; 
    
%     imshow(bw) 


%     se = strel('disk',15);% ��������
%   bw = imopen(bw,se); % ������� ����
% %     image_binor = bw2;  

    cc = bwconncomp(bw);
    stats = regionprops(cc,'Area','PixelList');

    for i = 1:length([stats(:).Area])
        if max < stats(i).Area
            max = stats(i).Area;
            number = i; 
        end
    end

    cc = bwconncomp(bw); 

    BW2 = ismember(labelmatrix(cc), number);
    
    

    
    image_binor = BW2;
    
%     imshow(BW2)
    
    
end

% function image_binor = clear_image(image) 
% 
% %     BW1 = edge(image,'Canny',0.15, 4); %��������� ������
%         BW1 = edge(image,'Canny',0.10, 4); %��������� ������
% %     imshow(BW1)
% %     pause(1)
%     
%     BW1 = bwmorph(BW1, 'bridge'); %���������� ������ 
% %         imshow(BW1)
% %     pause(1)
%     %������ ���������
%     bw = imfill(BW1,'holes'); % ���������� �������
% %         imshow(bw)
% %     pause(1)
% %     se = strel('disk',35);% ��������
%     
%         se = strel('disk',10);% ��������
%     
%     bw2 = imopen(bw,se); % ������� ����
%     image_binor = bw2; 
%     
% %     imshow(bw2)
% %     pause(1)
% end
