function [f,g] = path_energy(q_i, UserTols, num_agents, scene, e, surf_anim)
    
    q = q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = q_i(end-num_agents+1:end);
    
    %Weights
    K_agent = 1;
    K_tol = 1000; %don't touch
    K_accel = 0;
    K_map = 0;
    K_ke = 1;
    K_pv = 0;
    
    [e_agent, g_full] = agent_agent_energy(Q, Tols, scene, K_agent);
    [e_tol, g_tol] = tolerance_energy(Tols, UserTols, K_tol);
    [e_accel, g_accel] = acceleration_energy(Q, scene, K_accel);
    [e_map, g_map] = agent_map_energy( Q,Tols, UserTols, scene, K_map);
    [e_ke, g_ke] = kinetic_energy(Q, scene, K_ke);
    [e_pv, g_pv] = preferred_time_energy(Q, scene, K_pv);
    
    f = e_agent + e_tol + e_map + e_ke + e_accel + e_pv;
    
    g = g_full;
    g(1:end-num_agents) = g(1:end-num_agents) + g_map + g_ke + g_accel + g_pv;
    g(end-num_agents+1:end) = g(end-num_agents+1:end) + g_tol;
    
    %plottings(surf_anim, q, e, g);
end

function [e, g] = preferred_time_energy(Q, scene, K)
    GT = zeros(size(Q));
    e=0;
    pv = 0.001;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        v = dx(:,1:2)./dx(:,3);
        % 0.5* (s  - m/pv)^2
        
        e = e + K*0.5*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv) ).^2;
        
        dEdq = zeros(size(q_i));
        
        for j=1:size(dx,1)
            r_inds = [3*(j-1) + 1 : 3*(j-1)+6];
            q_i1 = q_i(r_inds(1));
            q_i2 = q_i(r_inds(2));
            q_i3 = q_i(r_inds(3));
            q_i4 = q_i(r_inds(4));
            q_i5 = q_i(r_inds(5));
            q_i6 = q_i(r_inds(6));
            
            dEdq(r_inds) = dEdq(r_inds) + ...
            [ -(K*(2*q_i1 - 2*q_i4)*(q_i6 - ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/pv))/(2*pv*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)),...
            -(K*(2*q_i2 - 2*q_i5)*(q_i6 - ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/pv))/(2*pv*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)),...
            0,...
            (K*(2*q_i1 - 2*q_i4)*(q_i6 - ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/pv))/(2*pv*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)),...
            (K*(2*q_i2 - 2*q_i5)*(q_i6 - ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/pv))/(2*pv*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)),...
            (K*(2*q_i6 - (2*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2))/pv))/2]';

        end
        
        GT(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
    end
    g = K*reshape(GT, size(GT,1)*size(GT,2),1);
end

function [e, g] = agent_agent_energy(Q, Tols, scene, K)
    num_agents = numel(scene.agents);
    GB = zeros(size(Q));
    gW = zeros(num_agents,1);
    e=0;
    g = zeros(numel(Q) + num_agents, 1);
   
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        for j =1+1:num_agents
            if j==i
                continue;
            end
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A2, J2] = sample_points_for_rod(A2, scene.agents(j).e);
            dist_is_good = 0;
            alpha_count = 10;
            alpha_val = 50;
            while dist_is_good==0
                [~,G1] = soft_distance(alpha_val,A2, A1);
                [D,G2] = soft_distance(alpha_val,A1, A2);
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
            JG1 = J1'*G1;
            JG2 = J2'*G2;
            
            tol = Tols(i) + Tols(j);

            e = e + -K*log(-tol + D);
            
            GB(:,i) = GB(:,i)+ (-1/(-tol + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
            GB(:,j) = GB(:,j)+ (-1/(-tol + D))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
            gW(i) = gW(i) + (1/(-tol+D));
            gW(j) = gW(j) + (1/(-tol+D));
            
%             B = B + 1/D;
%             GB(:,i) = GB(:,i)+ (-1/(D^2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
%             GB(:,j) = GB(:,j)+ (-1/(D^2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);

        end        
    end
    g(1:end-num_agents) = K*reshape(GB, size(GB,1)*size(GB,2),1);
    g(end-num_agents+1:end) = K*gW;
end
function [e, g] = tolerance_energy(Tols, UserTols, K)
    e = K*0.5*(Tols - UserTols')'*(Tols - UserTols');
    g = K*(Tols - UserTols');
end
function [e, g] = acceleration_energy(Q, scene, K)
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
        e = e + 0.5*K*sum((angle - 0).^2);

        %vectorized gradient computation
        %end points have zero gradient (no bending energy applied)
        gr = [0 0 0; ...
            K.*angle.*((cross(V2,Z)./(V2norm.*V2norm)) + (cross(V1,Z)./(V1norm.*V1norm))); ...
            0 0 0]; 
        GK(:,i) = reshape(gr', numel(gr),1);
    end
    g = reshape(GK, size(GK,1)*size(GK,2),1); 
end        
function [e, g] = agent_map_energy( Q, Tols, UserTols, scene, K)
    GB = zeros(size(Q));
    e=0;
    
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        
        A1(:, 3) = zeros(size(A1,1),1);
        %[P, ~] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
        P = scene.terrain.BV;
        [D,GP] = soft_distance(100,P, A1);
        
        JG1 = J1'*GP;
        GB(:,i) = (-1/(-UserTols(i) + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
        e = e + -K*log(-UserTols(i) + D);
    end
    g = K*reshape(GB, size(GB,1)*size(GB,2),1);
    
end
function [e, g] = kinetic_energy(Q, scene, K)
    %various fun path energies, I'll use principle of least action becase I like it
    %kinetic energy of curve integrated along piecewise linear segment is 
    % 0.5*dx'*dx./dt
    GT = zeros(size(Q));
    e=0;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        m = 1;
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';

        e = e + K*sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
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
    end
    g = K*reshape(GT, size(GT,1)*size(GT,2),1);
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


    PV = reshape(q, 3, numel(q)/3)';
    PE = e;
%     plot3(PV(:,1), PV(:,2), PV(:,3), '-ok', 'LineWidth',5);
    [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 4);
    surf_anim.Vertices = CV;
    drawnow;
end



