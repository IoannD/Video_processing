function [issol x1 x2] = quadratic_equation(A, B, C)
    D = B*B - 4*A*C; 
    if D < 0 issol = 0; x1 = 0; x2 = 0; 
        return; 
    else issol = 1; 
    end
    x1 = (-B + sqrt(D))/(2*A); x2 = (-B - sqrt(D))/(2*A);
end