function [f,g] = path_energy(q_i, UserTols, num_agents, e, surf_anim)
    global scene simple_sd mu_barrier;
    egTic = tic;
    
    q = q_i;%q_i(1:end-num_agents); %TODO
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = [];%q_i(end-num_agents+1:end); %TODO
    g = zeros(size(q_i));
    
    %% Weights and parameters
    K_agent = 1*scene.coeff_matrix(1,:);
    K_tol =   1*scene.coeff_matrix(2,:); %friendship
    K_accel = 1*scene.coeff_matrix(3,:);
    K_map =   1*scene.coeff_matrix(4,:);
    K_ke =    1*scene.coeff_matrix(5,:);
    K_pv =    1*scene.coeff_matrix(6,:);
    K_reg =   1*scene.coeff_matrix(7,:);
    A_mass = zeros(numel(scene.agents), 1)';
    A_pv = zeros(numel(scene.agents), 1)';
    KDTMdls = cell(1,numel(scene.agents));
    %KDTEnv = KDTreeSearcher(scene.terrain.BV);
    for i=1:numel(scene.agents)
        A_mass(i) = scene.agents(i).mass;
        A_pv(i) = scene.agents(i).preferred_end_time;
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A11, ~] = sample_points_for_rod(A1, scene.agents(i).e);
        KDTMdls{i} = KDTreeSearcher(A11);
    end
    
    
    %% fxns
    oneTic = tic;
    [e_agent, g_full] = agent_agent_energy(Q, Tols, scene, K_agent, K_tol, KDTMdls);
    scene.timings.iterations(end).egAgent = scene.timings.iterations(end).egAgent + toc(oneTic);
    
    oneTic = tic;
    [e_map, g_map] = agent_map_energy( Q,Tols, UserTols, scene, K_map);
    scene.timings.iterations(end).egMap = scene.timings.iterations(end).egMap + toc(oneTic);
    
    oneTic = tic;
    e_pv = mex_pv_energy(Q(:), K_pv', A_pv', num_agents, size(Q, 1)/3);
    scene.timings.iterations(end).ePv = scene.timings.iterations(end).ePv + toc(oneTic);
    
    oneTic = tic;
    e_accel = mex_accel_energy(Q(:), K_accel', num_agents, size(Q, 1)/3);
    scene.timings.iterations(end).eAcc = scene.timings.iterations(end).eAcc + toc(oneTic);
    
    oneTic = tic;
    e_ke = mex_kinetic_energy(Q(:), K_ke', A_mass', num_agents, size(Q, 1)/3);
    %[e_ke, g_ke] = kinetic_energy(Q, scene, K_ke);
    scene.timings.iterations(end).eKE = scene.timings.iterations(end).eKE + toc(oneTic);
    
    oneTic = tic;
    e_rg = mex_reg_energy(Q(:), K_reg', num_agents, size(Q, 1)/3);
    %[e_rg, g_rg] = regularizer_energy(Q, scene, K_reg);
    scene.timings.iterations(end).eReg = scene.timings.iterations(end).eReg + toc(oneTic);
    
    oneTic = tic;
    g_pv = mex_pv_gradient(Q(:),K_pv', A_pv', num_agents, size(Q,1)/3);
    scene.timings.iterations(end).ePv = scene.timings.iterations(end).ePv + toc(oneTic);
    
    oneTic = tic;
    g_accel = mex_accel_gradient(Q(:),K_accel', num_agents, size(Q,1)/3);
    scene.timings.iterations(end).gAcc = scene.timings.iterations(end).gAcc + toc(oneTic);
    
    oneTic = tic;
    g_ke = mex_kinetic_gradient(Q(:), K_ke', A_mass', num_agents, size(Q,1)/3);
    scene.timings.iterations(end).gKE = scene.timings.iterations(end).gKE + toc(oneTic);
    
    oneTic = tic;
    g_rg = mex_reg_gradient(Q(:), K_reg', num_agents, size(Q,1)/3);
    scene.timings.iterations(end).gReg = scene.timings.iterations(end).gReg + toc(oneTic);
    
    f = 0;
    f = f + mu_barrier*e_agent  + mu_barrier*e_map;
    
    f = f + e_ke + e_accel + e_rg + e_pv;
    
    g = g + mu_barrier*g_full;
    g = g + g_map + g_ke + g_accel + g_pv + g_rg + mu_barrier*g_map;
    scene.timings.iterations(end).egTotal = scene.timings.iterations(end).egTotal + toc(egTic);
    
    plottings(scene, q, [])
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

function [e, g] = agent_agent_energy(Q, Tols, scene, K, Ktol, KDTMdls)
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
            %if j is not included in the i's collision_interactions list,
            %continue
            if(~ismember(j,scene.agents(i).collision_interactions))
                continue;
            end

            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A22, J2] = sample_points_for_rod(A2, scene.agents(j).e);
            tol = scene.agents(i).radius+scene.agents(j).radius;%Tols(i) + Tols(j);
            dist_is_good = 0;
            if simple_sd
                alpha_count = scene.smoothDistAlpha/10;
                alpha_val = scene.smoothDistAlpha;
                beta_val = 0.0;
            else
                alpha_count = scene.smoothDistAlpha/10;
                alpha_val = scene.smoothDistAlpha;
                beta_val = 0.0;
            end
            while dist_is_good==0
                if simple_sd
                    [D,G1] = soft_distance(alpha_val,A22, A11, KDTMdls{i}, 10*tol);
                    [D,G2] = soft_distance(alpha_val,A11, A22,  KDTMdls{j}, 10*tol);
                else
                   return;
                end
                
                if(D>1e-8)
                    dist_is_good =1;
                end
                if(isinf(D))
                    sprintf("Error: path_energy, agent energy, D is inf")
                    alpha_val = alpha_val-alpha_count;
                else
                    alpha_val = alpha_val+alpha_count;
                end
       
            end
            
         
            JG1 = J1'*reshape(G1', size(G1,1)*size(G1,2), 1);
            JG2 = J2'*reshape(G2', size(G2,1)*size(G2,2), 1);     
            
            % Old energy
            %e = e + -K(i)*log((-tol + D).^2);
            %GB(:,i) = GB(:,i)+ K(i)*(-2/((-tol + D)))*JG1;
            %GB(:,j) = GB(:,j)+ K(j)*(-2/((-tol + D)))*JG2;
            
            %New energy
            
            e = e + -K(i)*log((-tol + D));
        
            GB(:,i) = GB(:,i)+ -(K(i)*JG1)/(-tol + D);
            GB(:,j) = GB(:,j)+ -(K(j)*JG2)/(-tol + D);
            
            if(ismember(j,scene.agents(i).friends))
                %friendship_radius = 4*max(scene.agents(i).radius, scene.agents(j).radius);
                e = e + Ktol(i)*K(i)*(D - tol)^2;
                GB(:,i) = GB(:,i)+ 2*Ktol(i)*K(i)*(D-tol)*JG1;% -(K(i)*JG1)/(-D + friendship_radius);
                GB(:,j) = GB(:,j)+ 2*Ktol(i)*K(j)*(D-tol)*JG2;%-(K(j)*JG2)/(-D + friendship_radius);
            end
                   
        end        
    end
    g = reshape(GB, size(GB,1)*size(GB,2),1); %g(1:end-num_agents) = reshape(GB, size(GB,1)*size(GB,2),1); %TODO
    %g(end-num_agents+1:end) = gW; TODO
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
        KDTMdli = KDTreeSearcher(A11);
        %[P, ~] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
        P = scene.terrain.BV;
        
        tol = scene.agents(i).radius;
        dist_is_good = 0;
        alpha_count =  scene.smoothDistAlpha/10;
        alpha_val = scene.smoothDistAlpha;
        
        while dist_is_good==0
            [D,GP] = soft_distance(alpha_val,P, A11, KDTMdli, 10*tol);
            %[D,GP] = soft_distance_original(alpha_val,P, A11);
            if(D>-1e-8)
                dist_is_good =1;
            end
            
            if(isinf(D))
                sprintf("Error: Path_energy, map energy, D is inf");
            end
            alpha_val = alpha_val+alpha_count;
        end
        
        %check if the values for D make any sense
        %make sure E, G, H are good
        % make sure softdistance alphas are in agreement
        
        
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


