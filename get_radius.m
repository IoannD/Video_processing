function VAR = get_radius(VAR)

    k = 1;
    meanRadius = 0;
    while (VAR.D(k) < 0.4) && (k < length(VAR.R))
        meanRadius = meanRadius + VAR.R(k); 
        k = k + 1; 
    end
    VAR.Radius = meanRadius / k; 
    VAR.Error_Radius = std(VAR.R)/sqrt(k)  
end