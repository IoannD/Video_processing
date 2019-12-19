function [VAR] = get_deformation(hor_ax_small, vert_ax_big, VAR)
    l = vert_ax_big*VAR.pixel_size; % переход от длины в пикселях к мм
    d = hor_ax_small*VAR.pixel_size;
    
%     VAR.R = (l+d)/2; % метод определение радиуса на основе усреднения 
%     Volume = d^2*l; % и большой полуосей

    VAR.R = (d.^2.*l).^(1/3); % метод определения радиуса капли - как радиуса 
    % шара с объемом равным объему эллипсу. Работает как с деформированной 
    % изначально, так и в случае отсутствия первичной деформации
    
    VAR.D = (l-d)./(l+d)*100; %расчет деформации
    VAR.time = linspace(VAR.frame_period,length(VAR.D)*VAR.frame_period,length(VAR.D)); 
    %время деформации
end
