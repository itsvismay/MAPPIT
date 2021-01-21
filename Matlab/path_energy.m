function [f,g] = path_energy(q_i, UserTols, num_agents, scene, e, surf_anim)
    q = q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = q_i(end-num_agents+1:end);
    f = 0;
    T = 0;
    B = 0;
    W = 0;
    K = 0;
    
    g = zeros(size(q_i));
    GT = zeros(size(Q));
    gW = zeros(num_agents,1);
    GB = zeros(size(Q));
    GK = zeros(size(Q));

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
        
        GT(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1, J1] = sample_points_for_rod(A1, 50);
        for j =1+1:num_agents
            if j==i
                continue;
            end
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A2, J2] = sample_points_for_rod(A2, 50);
            dist_is_good = 0;
            alpha_count = 10;
            alpha_val = 50;
            while dist_is_good==0
                [D,G1] = soft_distance(alpha_val,A2, A1);
                [D,G2] = soft_distance(alpha_val,A1, A2);
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
            JG1 = J1'*G1;
            JG2 = J2'*G2;
            
            tol = Tols(i) + Tols(j);
            A = -1*log(-tol + D);

            B = B + -1*log(-tol + D);
            GB(:,i) = GB(:,i)+ (-1/(-tol + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
            GB(:,j) = GB(:,j)+ (-1/(-tol + D))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
            gW(i) = gW(i) + (1/(-tol+D));
            gW(j) = gW(j) + (1/(-tol+D));
            
%             B = B + 1/D;
%             GB(:,i) = GB(:,i)+ (-1/(D^2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
%             GB(:,j) = GB(:,j)+ (-1/(D^2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);

        end
        A1(:, 3) = zeros(size(A1,1),1);
        [P, JP] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
        [D,GP] = soft_distance(50,P, A1);
        JG1 = J1'*GP;
        GB(:,i) = GB(:,i)+ (-1/(-tol + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
        
        B = B + -1*log(-UserTols(i) + D);
        
        %V_i = reshape(q_i, 3, numel(q_i)/3)';
        [ak, gk] = acceleration_energy(q_i);
        K = K + 1000*ak;
        GK(:,i) = GK(:,i) + 1000*gk;
%         
        
        %f = f + V_i(end,3);
    end
    
    W = 1000*0.5*(Tols - UserTols')'*(Tols - UserTols');
    gW = gW + 1000*(Tols - UserTols');
    
    %find minimum energy curve
    f = f + T+ B + W + K;
    T
    B
    W
    K
    g(1:end-num_agents) =reshape(GB, size(GB,1)*size(GB,2),1);
    gkk = reshape(GK, size(GK,1)*size(GK,2),1);
    gtt = reshape(GT, size(GT,1)*size(GT,2),1);
    g(1:end-num_agents) = g(1:end-num_agents) + gkk;
    g(1:end-num_agents) = g(1:end-num_agents) + gtt;
    g(end-num_agents+1:end) = gW;
%     
%     QQ = reshape(q, 3, numel(q)/3)';
%     GG = reshape(gkk, 3, (numel(gkk))/3)';
%     Xquiver = QQ(:,1);
%     Yquiver = QQ(:,2);
%     Zquiver = QQ(:,3);
%     Uquiver = GG(:,1);
%     Vquiver = GG(:,2);
%     Wquiver = GG(:,3);
%     quiver3(Xquiver, Yquiver, Zquiver, Uquiver, Vquiver, Wquiver);
    
   
%     PV = reshape(q, 3, numel(q)/3)';
%     PE = e;
%     plot3(PV(:,1), PV(:,2), PV(:,3), '-ok', 'LineWidth',5);
%     %[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 4);
%     %surf_anim.Vertices = CV;
%     drawnow;
end


