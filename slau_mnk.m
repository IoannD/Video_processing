function x = slau_mnk(M)
    [rows cols] = size(M); Z = [];
    if cols-1 < rows
        for col = 1:cols-1
            N = [];
            for i = 1:rows
                N = [N; M(i,col)*M(i, : )];
            end
            Z = [Z; sum(N)];
        end
    end
    b = -Z(:,end);
    A = Z(1:end,1:end-1);
    x = A\b; % решение приведенной системы средствами матлаба
end