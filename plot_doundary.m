function boundary = plot_doundary(image_binor, image_gray, flag)
% flag == 0 - do not plot 
% flag == 1 - plot 
    boundary = bwperim(image_binor,8);
    if flag == 1     
    imshowpair(boundary, image_gray, 'falsecolor');
%     imshowpair (boundary,image_gray,'diff')
    % 'falsecolor'	Creates a composite RGB image showing A and B 
    ... overlaid in different color bands. Gray regions in the 
    ...composite image show where the two images have the same 
    ...intensities. Magenta and green regions show where the 
    ...intensities are different. This is the default method.
    end  
end
