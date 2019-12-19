%% АВТОМАТИЧЕСКОЕ ОПРЕДЕЛЕНИЕ ПАРАМЕТРОВ ГАММА-КОРРЕКЦИИ
function VAR = get_gamma_corr(VAR)
% ИДЕЯ - набор, пикселей интенсивность которых лежит в диапазоне 
% (average - n*standard_deviation, average + n*standard_deviation)
% отображется в весь диапазон интесивностей
image_first = get_frame(VAR, 'first');
VAR.gamma_correction = [0 0];
n = 2;  
Idouble = im2double(image_first); 
avg = mean2(Idouble); % Average of matrix elements
sigma = std2(Idouble); % Standard deviation of matrix elements
VAR.gamma_correction = [avg-n*sigma avg+n*sigma];
end