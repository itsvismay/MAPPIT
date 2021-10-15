% VV - vertices
% EE - edges
% newA - Adjacency Matrix used to compute current path
% newA_visited - Adj Mat used to keep track of previously visited vertices
% s, t - start and end
% BV, Bind - Boundary verts and boundary indices into VV
function [Dist, Path, newA, newA_visited] = mydijk3d(VV, EE, newA, newA_visited, s, t, BV, Bind, agentRadius)
    global space_time_diags
    agent_radius = 2*agentRadius;
    
    % set the edge weights
    newA = set_edge_weights_directed(VV, newA, BV,Bind, newA_visited, agent_radius);
    
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
            edge_weight = edge_weight + 100*norm(VV(vi,:) - VV(t,:));
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
        newA_visited(xx(num),b) = 0;
        newA_visited(b,xx(num)) = 0;
    end

    %unwind
    hold on;
    while P(b)>0 && b ~= s
        Path = [b Path];
        
        %% Cross out edges within agent_radius from Adj Mat
        Idx_neighbors_in_radius = cell2mat(rangesearch(VV, VV(b,:), agent_radius));
        %plot3(VV(Idx_neighbors_in_radius,1), VV(Idx_neighbors_in_radius,2), VV(Idx_neighbors_in_radius,3), 'ko')
        Idx_neighbors_in_radius = Idx_neighbors_in_radius(2:end);%ignores node b
        %for each nearby neighbor, wipe out edges
        for idx_n = 1:length(Idx_neighbors_in_radius)
            neighbor_idx = Idx_neighbors_in_radius(idx_n);
            [xx,vv] = find(newA(:,neighbor_idx));
            for num=1:length(xx)
                newA_visited(xx(num),neighbor_idx) = 0;
                newA_visited(neighbor_idx, xx(num)) = 0;
            end
        end
        
 
            
        %% vertices based, cross out edges connected to vertices already visited from Adj mat
        [xx,vv] = find(newA(:,P(b)));
        for num=1:length(xx)
            newA_visited(xx(num),P(b)) = 0;
            newA_visited(P(b), xx(num)) = 0;
        end
        b = P(b);
        %% edges based
%         newA_visited(b,P(b)) = 1e-3;
%         newA_visited(P(b),b) = 1e-3;
%         b = P(b);
    end
    Path = [s Path];
    
    
end

function [E] = set_edge_weights_directed(Q, A, BV, Bind, A_visited, ar)   
    
     %% Cross out edges within agent_radius from Boundary
%     for i=1:length(Bind)
%         b = Bind(i);
%         Idx_neighbors_in_radius = cell2mat(rangesearch(Q, Q(b,:), ar));
%         %plot3(VV(Idx_neighbors_in_radius,1), VV(Idx_neighbors_in_radius,2), VV(Idx_neighbors_in_radius,3), 'ko')
%         Idx_neighbors_in_radius = Idx_neighbors_in_radius(2:end);%ignores node b
%         %for each nearby neighbor, wipe out edges
%         for idx_n = 1:length(Idx_neighbors_in_radius)
%             neighbor_idx = Idx_neighbors_in_radius(idx_n);
%             [xx,vv] = find(A(:,neighbor_idx));
%             for num=1:length(xx)
%                 A_visited(xx(num),neighbor_idx) = 0;
%                 A_visited(neighbor_idx, xx(num)) = 0;
%             end
%         end
%     end
%     E = A_visited;
    
    [ii,jj,ss] = find(A);
    for k=1:length(ii)
        ss(k) = 1;
    end
    
    for k=1:length(ii)
       %// A nonzero element of A: ss(k) = S(ii(k),jj(k))
       d = min_edge_to_boundary_dist(Q(ii(k),:), Q(jj(k),:), BV);
       
       % if d < ar, then set weight to 0
       if d<ar
           ss(k) = 0; 
       else
           ss(k) = d-ar;
       end
       
       
    end
    
    
    ii = [ii; size(A,1)];
    jj = [jj; size(A,2)];
    ss = [ss; 0];
    
    spA = sparse(ii,jj,ss);
    E = A_visited.*spA;
    
end

%Find edges in adjacency matrix A of mesh Q
%near node n within radius r
function [d] = find_edges_within_radius_of_node(n, Q, A, r)
    Q_idx = rangesearch(Q, n, r); %indexes of Q 
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