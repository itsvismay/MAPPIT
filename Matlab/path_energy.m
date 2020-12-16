function [f,g] = path_energy(q, num_agents, e, surf_anim)
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    G = zeros(size(Q));
    
    T = 0;
    V = 0;
    V2 = 0;
    
    g = zeros(size(q));

    for i=1:num_agents
        q_i = Q(:, i); %3*nodes
        
        %various fun path energies, I'll use principle of least action becase I like it
        %kinetic energy of curve integrated along piecewise linear segment is 
        % 0.5*dx'*dx./dt 
        m = 1;

        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        T = T + sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        
        dEdq_left = zeros(numel(q_i)/3, 3);
        dEdq_right = zeros(numel(q_i)/3, 3);
        
        dEdq_left(1:end-1,1) = -m*dx(:,1)./dx(:,3);
        dEdq_left(1:end-1,2) = -m*dx(:,2)./dx(:,3);
        dEdq_left(1:end-1,3) = 0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));
        
        dEdq_right(2:end, 1) = m*dx(:,1)./dx(:,3);
        dEdq_right(2:end, 2) = m*dx(:,2)./dx(:,3);
        dEdq_right(2:end, 3) = -0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));
        
        dEdq = dEdq_left + dEdq_right;
        
        G(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
    end
    g = reshape(G, size(G,1)*size(G,2),1);
    
    %find minimum energy curve
    f = T;
    
%     T = abs(T)
%     PV = reshape(q, 3, numel(q)/3)';
%     PE = e;
%     [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
%     surf_anim.Vertices = CV;
%     drawnow;
    

end