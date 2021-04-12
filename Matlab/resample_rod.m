function [P,Edges,J] = resample_rod(X, E, s)
    %INPUTS
    % s - vector of ints. Points to sample points between each edge
    %for each edge
    P = zeros(sum(s), 3);
    J = zeros(size(P,1), size(X,1));
    seg_count =0;
    ind = 0;
    for i = 1:size(E,1)
        tl = X(E(i,1),3);
        tr = X(E(i,2),3);
        xl = X(E(i,1),1);
        xr = X(E(i,2),1);
        yl = X(E(i,1),2);
        yr = X(E(i,2),2);
        for j =1:s(i)
            ind = ind + 1;
            seg_count = seg_count +1;
            P(ind,3) = tl + (j-1)*(tr - tl)/s(i);
            P(ind,1) = xl + (j-1)*(xr - xl)/s(i);
            P(ind,2) = yl + (j-1)*(yr - yl)/s(i);
            J(ind,E(i,1)) = 1 -(j-1)/s(i);
            J(ind,E(i,2)) = (j-1)/s(i);
        end
    end
    P(end,:) = X(end,:);
    J(end,end-1) = 0;
    J(end,end) = 1;
    Edges = [(1:size(P,1)-1)' (2:size(P,1))'];
  
end
