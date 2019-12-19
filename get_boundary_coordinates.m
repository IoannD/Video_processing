%% акнй тсмйжхи дкъ юоопнйяхлюжхх щккхоянл

function [X,Y] = get_boundary_coordinates(boundary)
    
    [M, N] = size(boundary); 
    t = 1;
    for i=1:M
        for j=1:N
%             if boudnary(i,j) ~= 0 
            if boundary(i , j) ~=0
                X(t) = j; 
                Y(t) = i;
                t = t+1; 
            end
        end
    end
    X = X';
    Y = Y'; 
 end