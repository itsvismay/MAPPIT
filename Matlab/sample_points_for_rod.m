function [P,J] = sample_points_for_rod(X, E, varargin)
    s = 10; %samples per edge
    %for each edge
    P = zeros(s*(size(X,1)-1)+1, 3);
    J = zeros(3*size(P,1), 3*size(X,1));
    for i = 1:size(E,1)
        tl = X(E(i,1),3);
        tr = X(E(i,2),3);
        xl = X(E(i,1),1);
        xr = X(E(i,2),1);
        yl = X(E(i,1),2);
        yr = X(E(i,2), 2);
        for j =1:s
            P(s*(i-1)+j,3) = tl + (j-1)*(tr - tl)/s;
            P(s*(i-1)+j,1) = xl + (j-1)*(xr - xl)/s;
            P(s*(i-1)+j,2) = yl + (j-1)*(yr - yl)/s;
            J1 = 1 -(j-1)/s;
            J2 = (j-1)/s;
            i_p = 3*(s*(i-1)+j - 1) + 1;
            i_x = 3*(E(i,1)-1) + 1;
            J(i_p + 0, i_x  + 0) = J1;
            J(i_p + 0, i_x + 3 + 0) = J2;
            
            J(i_p + 1, i_x  + 1) = J1;
            J(i_p + 1, i_x + 3 + 1) = J2;
            
            J(i_p + 2, i_x + 2) = J1;
            J(i_p + 2, i_x + 3 + 2) = J2;
        end
    end
    P(end,:) = X(end,:);
    %final point jacobian coordinates is set manually
    J(end-2,end-2) = 1;
    J(end-1,end-1) = 1;
    J(end,end) = 1;
    
  
end