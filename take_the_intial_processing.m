function [image_output] =   take_the_intial_processing(image_input, VAR)
    X = 0; Y = 0;
    image_gray = rgb2gray(image_input);
    [image_gray, X, Y] = get_resized_image(image_gray, image_binor, X, Y, VAR);
    after_filtr = average_filter(image_gray, VAR);
    image_binor = binarization(after_filtr, VAR);
    image_binor = clear_image(image_binor);
    boundary = plot_doundary(image_binor, image_gray, 0);
    % второй проход по изображению 
    [image_gray,X,Y] = get_resized_image(image_gray, image_binor, X, Y, VAR);
    after_filtr = average_filter(image_gray, VAR);
    image_binor = binarization(after_filtr, VAR);
    image_binor = clear_image(image_binor);
    boundary = plot_doundary(image_binor, image_gray, 1);
    
    image_output = image_binor;
    end