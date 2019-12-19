function ell = detect_ellipse(x, y)
    % получает на вход x и y - векторы строки строки с координатами точек, 
    % которые аппроксимируются
    x2 = x.*x; y2 = y.*y; xy = x.*y; one = ones(length(x),1);
    data = [x2 y2 xy x y one];
    X = slau_mnk(data);
    a11 = X(1); a22 = X(2); a12 = X(3)/2; a13 = X(4)/2; a23 = X(5)/2; a33 = 1;
    M = [a11 a12 a13; a12 a22 a23; a13 a23 a33];
    d = det(M);
    D = a11*a22 - a12*a12;
    I = a11 + a22;

    Mo.X = (a12*a23 - a13*a22)/D;
    Mo.Y = (a13*a12 - a11*a23)/D;
    tg = 2*a12/(a11-a22);
    k = (sqrt(tg*tg+1)-1)/tg;
    alpha = atan(k);

    [issol lamda1 lamda2] = quadratic_equation(1, -I, D);

    a = sqrt(-d/(lamda1*lamda2*lamda2));
    b = sqrt(-d/(lamda1*lamda1*lamda2));
    c = sqrt(abs(a*a-b*b));

    if a < b 
        k = -1/k; 
    end
    p = c/sqrt(1+k*k);
    F1.X = Mo.X - p; F1.Y = Mo.Y - k*p;
    F2.X = Mo.X + p; F2.Y = Mo.Y + k*p;

    ell.Mo = Mo; ell.F1 = F1; ell.F2 = F2;
    ell.alpha = alpha; ell.alphagrad = 180/pi*alpha;
    ell.a = a; ell.b = b; ell.c = c; ell.k = k; ell.M = M;
end
