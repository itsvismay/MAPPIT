function [Dist, Path, VV, EE] = mydijk3d(Q, A, A_visited, s, t, BV, Bind)
    % adjacency matrix of 2d graph
    A = set_edge_weights(Q, A, BV, A_visited);
    
    % build up 3d graph
    % find all the edges in 2d graph
    A_lt = tril(A);
    [edge_s,edge_t] = find(A_lt);
    edge = zeros(length(edge_s),2);
    edge(:,1) = edge_s;
    edge(:,2) = edge_t;
    
    % set the num of layer and get the spacetime graph 
    % (vertices: VV and edges: EE)
    nLayer = 2;
    time = linspace(0,10,nLayer)';
    [ver,EE] = spacetime_graph(Q,edge,time);
    VV = zeros(length(ver),3);
    VV(:,1) = ver(:,1);
    VV(:,2) = ver(:,2);
    VV(:,3) = ver(:,4);
    tsurf(EE,VV);
    
    % adjacency matrix of 3d graph
    newA = adjacency_matrix(EE);
    vPerLayer = length(ver) / nLayer;
    
    newA_visited = newA;
    
    % set the edge weights
    newA = set_edge_weights(VV, newA, BV, newA_visited);
    
     % make the graph directed along the time dimension
    [ii,jj,ss] = find(newA);
    for k=1:length(ii)
       s_temp = ii(k);
       t_temp = jj(k);
       if VV(s_temp,3) > VV(t_temp,3)
           %newA(s_temp, t_temp) = 0;
       end
    end
    
    % replace the orginal end point with the one on the top layer
    t = t + (nLayer-1)*vPerLayer;
   
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
    %unwind
    while P(b)>0 && b ~= s
        Path = [b Path];
        newA_visited(P(b), b) = 1e-3;
        newA_visited(b, P(b)) = 1e-3;
        b = P(b);
    end
    Path = [s Path];
    
end

function [E] = set_edge_weights(Q, A, BV, A_visited)
    A_lt = tril(A);
    [ii,jj,ss] = find(A_lt);
    for k=1:length(ii)
       %// A nonzero element of A: ss(k) = S(ii(k),jj(k))
       d = min_edge_to_boundary_dist(Q(ii(k),:), Q(jj(k),:), BV);
       ss(k) = d;
    end
    ii = [ii; size(A,1)];
    jj = [jj; size(A,2)];
    ss = [ss; 0];
    
    spA = sparse(ii,jj,ss); %triangular matrix. needs to be made symmetric
    E = A_visited.*(spA + spA');
    
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