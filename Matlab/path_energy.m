function [f,g] = path_energy(q_i, UserTols, num_agents, scene, e, surf_anim)
    global simple_sd mu_barrier;
    q = q_i;%q_i(1:end-num_agents); %TODO
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = [];%q_i(end-num_agents+1:end); %TODO
    g = zeros(size(q_i));
    
    %Weights
    K_agent = 1*scene.coeff_matrix(1,:);
    %K_tol =   0*scene.coeff_matrix(2,:); %don't touch
    K_accel = 1*scene.coeff_matrix(3,:);
    K_map =   1*scene.coeff_matrix(4,:);
    K_ke =    1*scene.coeff_matrix(5,:);
    K_pv =    0*scene.coeff_matrix(6,:);
    K_reg =   1*scene.coeff_matrix(7,:);
    
 
    [e_agent, g_full] = agent_agent_energy(Q, Tols, scene, K_agent);
    %[e_tol, g_tol] = tolerance_energy(Tols, UserTols, K_tol);
    [e_map, g_map] = agent_map_energy( Q,Tols, UserTols, scene, K_map);
    
 
    %[e_accel, g_accel] = acceleration_energy(Q, scene, K_accel);
    %[e_ke, g_ke] = kinetic_energy(Q, scene, K_ke);
    [e_pv, g_pv] = preferred_time_energy(Q, scene, K_pv);
    [e_rg, g_rg] = regularizer_energy(Q, scene, K_reg);
    
    e_accel = mex_accel_energy(Q(:), num_agents, size(Q, 1)/3, K_accel);
    e_ke = mex_kinetic_energy(Q(:), num_agents, size(Q, 1)/3, K_ke);
    e_tol = 0;% mex_tol_energy(Tols, UserTols, K_tol); TODO
    g_accel = mex_accel_gradient(Q(:), num_agents, size(Q,1)/3, K_accel);
    g_ke = mex_kinetic_gradient(Q(:), num_agents, size(Q,1)/3, K_ke);
    %g_tol = mex_tol_gradient(Tols, UserTols, K_tol);
    
    f = 0;
    f = f + mu_barrier*e_agent + e_tol + mu_barrier*e_map;
    
    f = f + e_ke + e_accel + e_pv + e_rg;
    
    g = g + mu_barrier*g_full;
    g = g + g_map + g_ke + g_accel + g_pv + g_rg + mu_barrier*g_map;
    
    %g(end-num_agents+1:end) = g(end-num_agents+1:end) + g_tol; TODO
    %g(1:end-num_agents) = g(1:end-num_agents) + g_ke + g_accel + g_pv + g_rg;
    
    plottings(surf_anim, q, e, g_accel);
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
        e = e + K(i)*(sum(kt./dt));
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
        GT(:,i) = K(i)*reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
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
    global simple_sd;
    num_agents = numel(scene.agents);
    if sum(K)==0
        e=0;
        g = zeros(numel(Q), 1); %zeros(numel(Q) + num_agents, 1); %TODO: tols
        return;
    end
    GB = zeros(size(Q));
    gW = zeros(num_agents,1);
    e=0;
    g = zeros(numel(Q), 1); %zeros(numel(Q) + num_agents, 1); %TODO: tols
   
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A11, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        for j =i+1:num_agents
            if j==i
                continue;
            end
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A22, J2] = sample_points_for_rod(A2, scene.agents(j).e);
            dist_is_good = 0;
            if simple_sd
                alpha_count = 10;
                alpha_val = 10;
                beta_val = 0.0;
            else
                alpha_count = 1;
                alpha_val = 1;
                beta_val = 0.0;
            end
            while dist_is_good==0
                if simple_sd
                    [~,G1] = soft_distance(alpha_val,A22, A11);
                    [D,G2] = soft_distance(alpha_val,A11, A22);
                else
                    [B1,I1] = build_distance_bvh(A11,[]);
                    [B2,I2] = build_distance_bvh(A22,[]);
                    [~, G1] = smooth_min_distance(A22,[],B2,I2,alpha_val,A11,[],alpha_val, beta_val);
                    [D, G2] = smooth_min_distance(A11,[],B1,I1,alpha_val,A22,[],alpha_val, beta_val);
                end
                
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
            
         
            JG1 = J1'*reshape(G1', size(G1,1)*size(G1,2), 1);
            JG2 = J2'*reshape(G2', size(G2,1)*size(G2,2), 1);     
            
            tol = scene.agents(i).radius;%Tols(i) + Tols(j);
            
            % Old energy
            %e = e + -K(i)*log((-tol + D).^2);
            %GB(:,i) = GB(:,i)+ K(i)*(-2/((-tol + D)))*JG1;
            %GB(:,j) = GB(:,j)+ K(j)*(-2/((-tol + D)))*JG2;
            
            %New energy
             e = e + -K(i)*log((-tol + D));
             GB(:,i) = GB(:,i)+ -(K(i)*JG1)/(-tol + D);
             GB(:,j) = GB(:,j)+ -(K(j)*JG2)/(-tol + D);
            
%             if(ismember(j,scene.agents(i).friends))
%                 e = e + -K(i)*log(-D + 2);
% %                 GB(:,i) = GB(:,i)+ K(i)*(-1/(-D + 2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
% %                 GB(:,j) = GB(:,j)+ K(i)*(-1/(-D + 2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
%                 GB(:,i) = GB(:,i)+ K(i)*(1/(-D + 2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
%                 GB(:,j) = GB(:,j)+ K(j)*(1/(-D + 2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
%             end
            
            
            
%             gW(i) = gW(i) + K(i)*(2/((-tol + D))); TODO
%             gW(j) = gW(j) + K(j)*(2/((-tol + D)));
            
        end        
    end
    g = reshape(GB, size(GB,1)*size(GB,2),1); %g(1:end-num_agents) = reshape(GB, size(GB,1)*size(GB,2),1); %TODO
    %g(end-num_agents+1:end) = gW; TODO
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
        [A11, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        
        A11(:, 3) = zeros(size(A11,1),1);
        %[P, ~] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
        P = scene.terrain.BV;
        
        dist_is_good = 0;
        alpha_count = 10;
        alpha_val = 10;

        while dist_is_good==0
            [D,GP] = soft_distance(alpha_val,P, A11);

            if(D>-1e-8)
                dist_is_good =1;
            end
            alpha_val = alpha_val+alpha_count;
        end
        
        %check if the values for D make any sense
        %make sure E, G, H are good
        % make sure softdistance alphas are in agreement
        
        tol = scene.agents(i).radius;
        JG1 = J1'*reshape(GP', size(GP,1)*size(GP,2), 1);
        
        
        ei = -K(i)*log((-tol + D));
%         if(~isreal(ei))
%             ei = inf;
%         end
        e = e+ ei;
        GB(:,i) = GB(:,i)+ -(K(i)*JG1)/(-tol + D);
        %e = e + -K(i)*log((-tol + D).^2);
        %GB(:,i) = GB(:,i)+ K(i)*(-2/((-tol + D)))*JG1;
        
    end
    g = reshape(GB, size(GB,1)*size(GB,2),1);
    
end
function [e, g] = kinetic_energy(Q, scene, K)
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

        e = e + K(i)*sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        dEdq_left = zeros(numel(q_i)/3, 3);
        dEdq_right = zeros(numel(q_i)/3, 3);

        dEdq_left(1:end-1,1) = -m*dx(:,1)./dx(:,3);
        dEdq_left(1:end-1,2) = -m*dx(:,2)./dx(:,3);
        dEdq_left(1:end-1,3) = 0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));

        dEdq_right(2:end, 1) = m*dx(:,1)./dx(:,3);
        dEdq_right(2:end, 2) = m*dx(:,2)./dx(:,3);
        dEdq_right(2:end, 3) = -0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));

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
%     quiver3(Xquiver, Yquiver, Zquiver, Uquiver, Vquiver, Wquiver, 5);
%     drawnow;

    PV = reshape(q, 3, numel(q)/3)';
    PE = e;
%     plot3(PV(:,1), PV(:,2), PV(:,3), '-ok', 'LineWidth',5);
    [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1.5, 'PolySize', 4);
    surf_anim.Vertices = CV;
    drawnow;
end



