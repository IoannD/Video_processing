function image_binor = clear_image(image) 

    
    BW1 = edge(image,'Canny',0.15, 4); %выделение границ
%     BW1 = edge(image,'Canny',0.1, 6); %выделение границ
%     imshow(BW1)
%     pause(1)
    
    BW1 = bwmorph(BW1, 'bridge'); %связывание границ 
%         imshow(BW1)
%     pause(1)
    %узнать подробнее
    bw = imfill(BW1,'holes'); % заполнение области
%     imshow(bw)
%     pause(1) 
    
    max = 0;
    number = 0; 
    
%     imshow(bw) 


%     se = strel('disk',15);% стирание
%   bw = imopen(bw,se); % лишнего шума
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
% %     BW1 = edge(image,'Canny',0.15, 4); %выделение границ
%         BW1 = edge(image,'Canny',0.10, 4); %выделение границ
% %     imshow(BW1)
% %     pause(1)
%     
%     BW1 = bwmorph(BW1, 'bridge'); %связывание границ 
% %         imshow(BW1)
% %     pause(1)
%     %узнать подробнее
%     bw = imfill(BW1,'holes'); % заполнение области
% %         imshow(bw)
% %     pause(1)
% %     se = strel('disk',35);% стирание
%     
%         se = strel('disk',10);% стирание
%     
%     bw2 = imopen(bw,se); % лишнего шума
%     image_binor = bw2; 
%     
% %     imshow(bw2)
% %     pause(1)
% end
