function [f, g] = path_energy_hessian(q_i, UserTols, num_agents, scene, e, surf_anim, constraint_set)
    
    q = q_i;% q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    g_full = zeros(size(q_i));

    H_full = zeros(size(q_i,1),size(q_i,1));
    
    Tols = q_i(end-num_agents+1:end);
    
    %Weights
    K_agent = 0*scene.coeff_matrix(1,:);
    K_tol =   0*scene.coeff_matrix(2,:); %don't touch
    K_accel = 0*scene.coeff_matrix(3,:);
    K_map =   0*scene.coeff_matrix(4,:);
    K_ke =    1*scene.coeff_matrix(5,:);
    K_pv =    0*scene.coeff_matrix(6,:);
    K_reg =   1;
    

%      [e_agent, g_agent] = agent_agent_energy(Q, Tols, scene, K_agent);
%      [e_tol, g_tol] = tolerance_energy(Tols, UserTols, K_tol);
%      [e_accel, g_accel] = acceleration_energy(Q, scene, K_accel);
%      [e_map, g_map] = agent_map_energy( Q,Tols, UserTols, scene, K_map);
     [e_ke, g_ke] = kinetic_energy(Q, scene, K_ke, constraint_set);
     %[e_ke, H_ke] = kinetic_energy_with_hessian(Q, scene, K_ke, constraint_set);
%      [e_pv, g_pv] = preferred_time_energy(Q, scene, K_pv);
%      [e_rg, g_rg] = regularizer_energy(Q, scene, K_reg);
    
    f = e_ke;%e_rg + e_agent + e_map + e_tol + e_ke + e_accel + e_pv;

    g = g_full + g_ke;% + g_agent;
    %g(1:end-num_agents) = g(1:end-num_agents) + g_ke + g_accel + g_pv + g_map + g_rg;
    %g(end-num_agents+1:end) = g(end-num_agents+1:end) + g_tol;
    
%     for ii=1:num_agents
%         H_full(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H_ke(:,:,ii);
%         %H_full(end-ii+1,end-ii+1) = 1; %assign the last block to identity
%     end
    
%     H = H_full;
    
    plottings(surf_anim, q, e, g_full(1:end-3));
end
function [e, g] = regularizer_energy(Q, scene, K)
    if sum(K)==0
        e=0;
        g = zeros(numel(Q),1);
        return;
    end
    GT = zeros(size(Q));
    e=0;
    
    
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        dt = dx(:,3)+1e-6;%add epsilon to make sure there is never a divide by 0
        endtime = q_i(end);
        segments = numel(q_i)/3-1;
        kt = (endtime/segments)*ones(size(dx,1),1);%regular time intervals over the rod;
        e = e + K*(sum(kt./dt));
        dEdq_left = zeros(numel(q_i)/3, 3);
        dEdq_right = zeros(numel(q_i)/3, 3);
        dEdq_left(1:end-1,3) = kt./(dt.^2);
        dEdq_right(2:end, 3) = -kt./(dt.^2);
        
        %e = e + K*0.5*(sum((dt - kt).^2));
        %dEdq_left = zeros(numel(q_i)/3, 3);
        %dEdq_right = zeros(numel(q_i)/3, 3);
        %dEdq_left(1:end-1,3) = -(dt - kt);
        %dEdq_right(2:end, 3) = (dt - kt);
        
        dEdq = dEdq_left + dEdq_right;
        GT(:,i) = K*reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
    end
    g = reshape(GT, size(GT,1)*size(GT,2),1);
end
function [e, g] = preferred_time_energy(Q, scene, K)
    if sum(K)==0
        e=0;
        g = zeros(numel(Q),1);
        return;
    end
    GT = zeros(size(Q));
    e=0;
    pv = 5;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        
        e = e + K(i)*0.5*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv) ).^2;

%         e = e + K*0.5*(q_i(end).^2 - (sum((dx(:,1).^2  + dx(:,2).^2))/(pv*pv)) ).^2;
        
        gradE = ones(size(q_i));
        gradE = gradE*K(i)*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv));
        
        dedx2 = -(dx(:,1).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;
        dedx1 = (dx(:,1).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;
        dedy2 = -(dx(:,2).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;
        dedy1 = (dx(:,2).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;

        dEdq_left = zeros(numel(q_i)/3, 3);
        dEdq_right = zeros(numel(q_i)/3, 3);
        dEdq_left(1:end-1,1) = dedx1;
        dEdq_left(1:end-1,2) = dedy1;
        dEdq_right(2:end, 1) = dedx2;
        dEdq_right(2:end, 2) = dedy2;

        dEdq = dEdq_left + dEdq_right;
        flatdEdq = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        gradE(1:end-1) = gradE(1:end-1).*flatdEdq(1:end-1);
        
        GT(:,i) = gradE;
        
    end
    g = reshape(GT, size(GT,1)*size(GT,2),1);
end
function [e, g] = agent_agent_energy(Q, Tols, scene, K)
    num_agents = numel(scene.agents);
    if sum(K)==0
        e=0;
        g = zeros(numel(Q) + num_agents, 1);
        return;
    end
    GB = zeros(size(Q));
    gW = zeros(num_agents,1);
    e=0;
    g = zeros(numel(Q) + num_agents, 1);
   
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A11, E11, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        for j =i+1:num_agents
            if j==i
                continue;
            end
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A22, E22, J2] = sample_points_for_rod(A2, scene.agents(j).e);
            dist_is_good = 0;
            alpha_count =2;
            alpha_val = 5;
            while dist_is_good==0
                %[~,G1] = soft_distance(alpha_val,A2, A1);
                %[D,G2] = soft_distance(alpha_val,A1, A2);
                %[~, G2] = smooth_min_distance(A1,[],alpha_val,A2,[],alpha_val);
                %[D, G1] = smooth_min_distance(A2,[],alpha_val,A1,[],alpha_val);
                
                [B1,I1] = build_distance_bvh(A1,scene.agents(i).e);
                [B2,I2] = build_distance_bvh(A2,scene.agents(j).e);
                [~, G1] = smooth_min_distance(A2,scene.agents(j).e,B2,I2,alpha_val,A11,[],alpha_val);
                [D, G2] = smooth_min_distance(A1,scene.agents(i).e,B1,I1,alpha_val,A22,[],alpha_val);
                
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
            JG1 = J1'*G1;
            JG2 = J2'*G2;
            
            tol = Tols(i) + Tols(j);
            
            e = e + -K(i)*log(-tol + D);
            if(ismember(j,scene.agents(i).friends))
                e = e + -K(i)*log(-D + 2);
                GB(:,i) = GB(:,i)+ K(i)*(-1/(-D + 2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
                GB(:,j) = GB(:,j)+ K(i)*(-1/(-D + 2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
            end
            
            GB(:,i) = GB(:,i)+ K(i)*(-1/(-tol + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
            GB(:,j) = GB(:,j)+ K(j)*(-1/(-tol + D))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
            gW(i) = gW(i) + K(i)*(1/(-tol+D));
            gW(j) = gW(j) + K(j)*(1/(-tol+D));
            
%             B = B + 1/D;
%             GB(:,i) = GB(:,i)+ (-1/(D^2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
%             GB(:,j) = GB(:,j)+ (-1/(D^2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);

        end        
    end
    g(1:end-num_agents) = reshape(GB, size(GB,1)*size(GB,2),1);
    g(end-num_agents+1:end) = gW;
end
function [e, g] = tolerance_energy(Tols, UserTols, K)
    e = 0.5*(K'.*(Tols - UserTols'))'*(Tols - UserTols');
    g = K'.*(Tols - UserTols');
end
function [e, g] = acceleration_energy(Q, scene, K)
    if sum(K)==0
        e=0;
        g = zeros(numel(Q),1);
        return;
    end
    GK = zeros(size(Q));
    e=0;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        V1 = dx(1:end-1,:);
        V2= dx(2:end,:);
        V1xV2 = cross(V1, V2, 2);
        V1dV2 = dot(V1, V2, 2);
        V1xV2norm = sqrt(V1xV2(:,1).^2+V1xV2(:,2).^2+V1xV2(:,3).^2);
        V1xV2norm(V1xV2norm < eps)= eps; %making sure norms are bigger than epsilon
        Z = V1xV2 ./ V1xV2norm; 
        V1norm = sqrt(sum(V1.^2,2));
        V2norm = sqrt(sum(V2.^2,2));
        Y = dot(V1xV2, Z, 2);
        X = V1norm.*V2norm + V1dV2;
        angle = 2*atan2(Y,X);
        
        %bending energy
        e = e + 0.5*K(i)*sum((angle - 0).^2);

        %vectorized gradient computation
        %end points have zero gradient (no bending energy applied)
        gr = [0 0 0; ...
            K(i).*angle.*((cross(V2,Z)./(V2norm.*V2norm)) + (cross(V1,Z)./(V1norm.*V1norm))); ...
            0 0 0]; 
        GK(:,i) = reshape(gr', numel(gr),1);
    end
    g = reshape(GK, size(GK,1)*size(GK,2),1); 
end        
function [e, g] = agent_map_energy( Q, Tols, UserTols, scene, K)
    if sum(K)==0
        e=0;
        g = zeros(numel(Q),1);
        return;
    end
    GB = zeros(size(Q));
    e=0;
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1,E1, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        
        A1(:, 3) = zeros(size(A1,1),1);
        %[P, ~] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
        P = scene.terrain.BV;
        [D,GP] = soft_distance(100,P, A1);
        
        JG1 = J1'*GP;
        GB(:,i) = K(i)*(-1/(-UserTols(i) + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
        e = e + -K(i)*log(-UserTols(i) + D);
    end
    g = reshape(GB, size(GB,1)*size(GB,2),1);
    
end
function [e, g] = kinetic_energy(Q, scene, K, constraint_set)
    if sum(K)==0
        e=0;
        g = zeros(numel(Q),1);
        return;
    end
    %various fun path energies, I'll use principle of least action becase I like it
    %kinetic energy of curve integrated along piecewise linear segment is 
    % 0.5*dx'*dx./dt
    GT = zeros(size(Q));
    e=0;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        m = scene.agents(i).mass;
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        
         dEdq_left = zeros(numel(q_i)/3, 3);
         dEdq_right = zeros(numel(q_i)/3, 3);
   
        if constraint_set == 1
            %spring energy
            e = e + K(i)*sum(0.5*m*sum(dx(:, 1:3).*dx(:,1:3),2)); %kinetic energy
            
            dEdq_left(1:end-1,1) = -m*dx(:,1);
            dEdq_left(1:end-1,2) = -m*dx(:,2);
            dEdq_left(1:end-1,3) = -m*dx(:,3);
            
            dEdq_right(2:end, 1) = m*dx(:,1);
            dEdq_right(2:end, 2) = m*dx(:,2);
            dEdq_right(2:end, 3) = m*dx(:,3);

        
        else
            e = e + K(i)*sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
           
            dEdq_left(1:end-1,1) = -m*dx(:,1)./dx(:,3);
            dEdq_left(1:end-1,2) = -m*dx(:,2)./dx(:,3);
            dEdq_left(1:end-1,3) = 0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));

            dEdq_right(2:end, 1) = m*dx(:,1)./dx(:,3);
            dEdq_right(2:end, 2) = m*dx(:,2)./dx(:,3);
            dEdq_right(2:end, 3) = -0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));

        end
        
        
        dEdq = dEdq_left + dEdq_right;
        
        GT(:,i) = K(i)*reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
    end
    g = reshape(GT, size(GT,1)*size(GT,2),1);
end
function [s] = plottings(surf_anim, q, e, g)
%     QQ = reshape(q, 3, numel(q)/3)';
%     GG = reshape(g, 3, (numel(g))/3)';
%     Xquiver = QQ(:,1);
%     Yquiver = QQ(:,2);
%     Zquiver = QQ(:,3);
%     Uquiver = GG(:,1);
%     Vquiver = GG(:,2);
%     Wquiver = GG(:,3);
%     quiver3(Xquiver, Yquiver, Zquiver, Uquiver, Vquiver, Wquiver);
%     drawnow;

    PV = reshape(q, 3, numel(q)/3)';
    PE = e;
%     plot3(PV(:,1), PV(:,2), PV(:,3), '-ok', 'LineWidth',5);
    [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.75, 'PolySize', 4);
    surf_anim.Vertices = CV;
    drawnow;
end



