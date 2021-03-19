function [Dist, Path, newA, newA_visited] = mydijk3d(VV, EE, newA, newA_visited, s, t, BV, Bind)
    % set the edge weights
    newA = set_edge_weights_directed(VV, newA, BV, newA_visited);
    
    % 3d dijkstra
    n = size(newA, 1);
    Qlabel = linspace(1, size(VV,1), size(VV,1));
    
    D = Inf*ones(n,1); 
    D(s) = 0;
    P = -1*ones(n,1);

	while length(Qlabel)>0
        [d, ind]= min(D(Qlabel));
        u = Qlabel(ind);
        Qlabel(ind) = [];
    
        [neighbors,kA,Aj] = find(newA(:,u));
        for vi = neighbors'
            edge_weight = 1.0/newA(vi, u);%min_edge_to_boundary_dist(Q(u,:), Q(vi,:), BV);
            
            alt_dist = D(u) + edge_weight*norm(VV(u,:) - VV(vi,:));
            if alt_dist < D(vi)
                D(vi) = alt_dist;
                P(vi) = u;
            end
        end
    end

    
    Dist = D(t);
    b = t;
    Path = [];
    
    %Vertex Based
    [xx,vv] = find(newA(:,b));
    for num=1:length(xx)
        newA_visited(xx(num),b) = 1e-7;
        newA_visited(b,xx(num)) = 1e-7;
    end

    %unwind
    while P(b)>0 && b ~= s
        Path = [b Path];
        %% vertices based
        [xx,vv] = find(newA(:,P(b)));
        for num=1:length(xx)
            newA_visited(xx(num),P(b)) = 1e-7;
            newA_visited(P(b), xx(num)) = 1e-7;
        end
        b = P(b);
        %% edges based
%         newA_visited(b,P(b)) = 1e-3;
%         newA_visited(P(b),b) = 1e-3;
%         b = P(b);
    end
    Path = [s Path];
    
end

function [E] = set_edge_weights_directed(Q, A, BV, A_visited)
    [ii,jj,ss] = find(A);
    for k=1:length(ii)
       %// A nonzero element of A: ss(k) = S(ii(k),jj(k))
       d = min_edge_to_boundary_dist(Q(ii(k),:), Q(jj(k),:), BV);
       ss(k) = d;
    end
    ii = [ii; size(A,1)];
    jj = [jj; size(A,2)];
    ss = [ss; 0];
    
    spA = sparse(ii,jj,ss);
    E = A_visited.*spA;
    
end

function [d] = min_edge_to_boundary_dist(P1, P2, M)
    v = P2 - P1;
    %quadrature sample points along edgee (P1, s1, s2,..., P2)
    
    s1 = P1 + (1/4)*v;
    s2 = P1 + (2/4)*v;
    s3 = P1 + (3/4)*v;
    
    
    D = [sqrt(sum((M - P1).^2, 2)) ... 
        sqrt(sum((M - s1).^2, 2)) ...
        sqrt(sum((M - s2).^2, 2)) ...
        sqrt(sum((M - s3).^2, 2)) ...
        sqrt(sum((M - P2).^2, 2))];
    %d = rand*min(min(D)) + 1e-7;
    d = min(min(D)) + 1e-7;
end