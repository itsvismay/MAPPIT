function [Dist,Path, A, A_visited] = mydijk(Q, A, A_visited, s, t, BV, Bind)
    A = set_edge_weights(Q, A, BV, A_visited);
    
    n = size(A, 1);
    Qlabel = linspace(1, size(Q,1), size(Q,1));
    %CDATA = full(sum(A))';
    D = Inf*ones(n,1); 
    D(s) = 0;
    P = -1*ones(n,1);

	while length(Qlabel)>0
        [d, ind]= min(D(Qlabel));
        u = Qlabel(ind);
        Qlabel(ind) = [];
    
        [neighbors,kA,Aj] = find(A(:,u));
        for vi = neighbors'
            edge_weight = 1.0/A(vi, u);%min_edge_to_boundary_dist(Q(u,:), Q(vi,:), BV);
            
            alt_dist = D(u) + edge_weight*norm(Q(u,:) - Q(vi,:));
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
        A_visited(P(b), b) = 1e-3;
        A_visited(b, P(b)) = 1e-3;
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