function [P,Edges,J] = sample_points_for_rod(X, E, varargin)
    s = 10; %5 samples per edge
    %for each edge
    P = zeros(s*(size(X,1)-1)+1, 3);
    J = zeros(size(P,1), size(X,1));
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
            J(s*(i-1)+j,E(i,1)) = 1 -(j-1)/s;
            J(s*(i-1)+j,E(i,2)) = (j-1)/s;
        end
    end
    P(end,:) = X(end,:);
    J(end,end-1) = 0;
    J(end,end) = 1;
    Edges = [(1:size(P,1)-1)' (2:size(P,1))'];
  
end