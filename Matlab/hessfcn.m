function [H] = hessfcn(q_i,lambda)
    global num_agents scene simple_sd mu_barrier;
    hessTic = tic;
    q = q_i;%TODO: (1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = [];%TODO%Tols = q_i(end-num_agents+1:end);
    
    %Weights
    K_agent = 1*scene.coeff_matrix(1,:);
    K_tol =   1*scene.coeff_matrix(2,:); %friendship
    K_accel = 1*scene.coeff_matrix(3,:);
    K_map =   1*scene.coeff_matrix(4,:);
    K_ke =    1*scene.coeff_matrix(5,:);
    K_reg =   1*scene.coeff_matrix(7,:);
    A_mass = zeros(numel(scene.agents), 1)';
    A_pv = zeros(numel(scene.agents), 1)';
    for i=1:numel(scene.agents)
        A_mass(i) = scene.agents(i).mass;
        A_pv(i) = scene.agents(i).preferred_end_time;
    end
    
    H = sparse(size(q_i,1),size(q_i,1));

    oneHessTic = tic;
    if simple_sd
        [e_agent, HB_agent, hw_agent] = agent_agent_energy_with_hessian(Q, Tols, scene, K_agent, K_tol);
    else
        [e_agent, HB_agent, hw_agent] = agent_agent_energy_with_hessian_abhisheks(Q, Tols, scene, K_agent); 
    end
    scene.timings.iterations(end).hAgent = scene.timings.iterations(end).hAgent + toc(oneHessTic);
    
    oneHessTic = tic;
    [e_map, H_map] = map_energy_with_hessian(Q, Tols, scene, K_map);
    scene.timings.iterations(end).hMap = scene.timings.iterations(end).hMap + toc(oneHessTic);
    
    for ii=1:num_agents
%         for jj=1:num_agents
%             H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(jj-1)*size(Q,1):jj*size(Q,1)) = H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(jj-1)*size(Q,1):jj*size(Q,1))...
%                                                                                 + mu_barrier*HB_agent{ii,jj};
%         end
        H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1))...
                                                                            + mu_barrier*H_map{ii}...
                                                                            + mu_barrier*HB_agent{ii,ii};
    end
    
    oneHessTic = tic;
    [H_acc] = mex_accel_hessian(Q(:),K_accel', num_agents, size(Q,1)/3);
    scene.timings.iterations(end).hAcc = scene.timings.iterations(end).hAcc + toc(oneHessTic);
    H = H + H_acc;
     
     
    oneHessTic = tic;
    [H_ke] = mex_kinetic_hessian(Q(:),K_ke', A_mass', num_agents, size(Q, 1)/3);
    %[e_ke1, H_ke1] = kinetic_energy_with_hessian(Q, scene, K_ke);
    scene.timings.iterations(end).hKE = scene.timings.iterations(end).hKE + toc(oneHessTic);
    H = H + H_ke;
%     for ii=1:num_agents
%         H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1))...
%                                                                             + H_ke1(:,:,ii);
%     end
    
    oneHessTic = tic;
    [H_reg] = mex_reg_hessian(Q(:),K_reg', num_agents, size(Q, 1)/3);
    %[e_reg1, H_reg1] = regularizer_energy_with_hessian(Q, scene, K_reg);
    scene.timings.iterations(end).hReg = scene.timings.iterations(end).hReg + toc(oneHessTic);
    H = H + H_reg;

    scene.timings.iterations(end).hTotal = scene.timings.iterations(end).hTotal + toc(hessTic);
end

function [e, HT] = regularizer_energy_with_hessian(Q, scene, K)
    if sum(K)==0
        e=0;
        HT = zeros(size(Q,1), size(Q,1), size(Q,2));
        return;
    end
    HT = zeros(size(Q,1), size(Q,1), size(Q,2));
    e=0;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        dt = dx(:,3)+1e-6;%add epsilon to make sure there is never a divide by 0
        endtime = q_i(end);
        segments = numel(q_i)/3-1;
        kt = (endtime/segments)*ones(size(dx,1),1);%regular time intervals over the rod;
        e = e + K(i)*(sum(kt./dt));
       
        d2Edq2 = zeros(numel(q_i),numel(q_i));
        
        for kk=1:numel(q_i)/3-1
            d2Edq2((kk-1)*3+3,(kk-1)*3+3) = d2Edq2((kk-1)*3+3,(kk-1)*3+3) + 2*kt(kk,1)./(dt(kk,1).^3);

            d2Edq2((kk-1)*3+3,(kk-1)*3+6) = d2Edq2((kk-1)*3+3,(kk-1)*3+6) - 2*kt(kk,1)./(dt(kk,1).^3);

            d2Edq2((kk-1)*3+6,(kk-1)*3+3) = d2Edq2((kk-1)*3+6,(kk-1)*3+3) - 2*kt(kk,1)./(dt(kk,1).^3);

            d2Edq2((kk-1)*3+6,(kk-1)*3+6) = d2Edq2((kk-1)*3+6,(kk-1)*3+6) + 2*kt(kk,1)./(dt(kk,1).^3);
        end
               
        HT(:,:,i) = K(i)*d2Edq2;
            
    end
end

function [e, HB, hW] = agent_agent_energy_with_hessian(Q, Tols, scene, K, Ktol)
    HB = cell(size(Q,2));
    for i=1:numel(scene.agents)
        HB{i,i} = sparse(size(Q,1), size(Q,1));
        for j =i+1:numel(scene.agents)
            HB{i,j} = sparse(size(Q,1), size(Q,1));
            HB{j,i} = sparse(size(Q,1), size(Q,1));
        end
    end
    if sum(K)==0
        e=0;
        hW = zeros(size(Q,2), size(Q,2));
        return;
    end
   
    
    hW = zeros(size(Q,2), size(Q,2));
    e=0;
   
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A11, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        for j =i+1:numel(scene.agents)
            if j==i
                continue;
            end
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A22, J2] = sample_points_for_rod(A2, scene.agents(j).e);
            
            tol = scene.agents(i).radius+scene.agents(j).radius;%Tols(i) + Tols(j);
            dist_is_good = 0;
            alpha_count = scene.smoothDistAlpha/10;
            alpha_val = scene.smoothDistAlpha;
            while dist_is_good==0
                [~,G1] = soft_distance(alpha_val,A22, A11);
                [D,G2] = soft_distance(alpha_val,A11, A22);
                      
                if(D>-1-8)
                    dist_is_good =1;
                end
                
                if(isinf(D))
                    sprintf("Error: hessfcn, agent energy, D is inf");
                end
                
                alpha_val = alpha_val+alpha_count;
            end
            JG1 = J1'*reshape(G1', size(G1,1)*size(G1,2), 1);
            JG2 = J2'*reshape(G2', size(G2,1)*size(G2,2), 1);  
            
            %Old Energy
            %e = e + -K(i)*log((-tol + D).^2);
            %Hi = (K(i)*(2/((-tol + D).^2))) * (JG1 * JG1');
            %Hj = (K(j)*(2/((-tol + D).^2))) * (JG2 * JG2');
            
            %New Energy
            e = e + -K(i)*log((-tol + D));
            Hi = (K(i)*(1/((-tol + D).^2))) * (JG1 * JG1');
            Hj = (K(j)*(1/((-tol + D).^2))) * (JG2 * JG2');
            Hij = (K(j)*(1/((-tol + D).^2))) * (JG2 * JG1');
            if(ismember(j,scene.agents(i).friends))
                %friendship_radius = 3*max(scene.agents(i).radius, scene.agents(j).radius);
                e = e + Ktol(i)*K(i)*(D - tol)^2;
                %GB(:,i) = GB(:,i)+ 2*K(i)*(D-tol)*JG1;
                %GB(:,j) = GB(:,j)+ 2*K(j)*(D-tol)*JG2;
                Hi = Hi + 2*Ktol(i) * K(i) * (JG1 * JG1');
                Hj = Hj + 2*Ktol(i) * K(j) * (JG2 * JG2');
                %Hi = Hi + (K(i)*(1/((-D + friendship_radius).^2))) * (JG1 * JG1');
                %Hj = Hj + (K(j)*(1/((-D + friendship_radius).^2))) * (JG2 * JG2');
            end
            
            [spHi, spHj] = sparsify_dense_hessian(A1, A2, Hi, Hj, 0.2);
            [~, spHij] = sparsify_dense_hessian(A1, A2, Hij, Hij, 0.2);
            
            HB{i,i} = HB{i,i} + spHi;
            HB{j,j} = HB{j,j} + spHj;
            
            HB{i,j} = HB{i,j} + spHij;
            HB{j,i} = HB{j,i} - spHij;
            
            %hW(i) = hW(i) + K(i)*(2/((-tol+D).^2));TODO
            %hW(j) = hW(j) + K(j)*(2/((-tol+D).^2));
            
        end        
    end
end

function [e, HB, hW] = map_energy_with_hessian(Q, Tols, scene, K)
    if sum(K)==0
        e=0;
        HB(1:size(Q,2)) = {sparse(size(Q,1), size(Q,1))};
        hW = zeros(size(Q,2), size(Q,2));
        return;
    end
   
    HB(1:size(Q,2)) = {sparse(size(Q,1), size(Q,1))};
    hW = zeros(size(Q,2), size(Q,2));
    e=0;
   
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A11, J1] = sample_points_for_rod(A1, scene.agents(i).e);
       
        A11(:, 3) = zeros(size(A11,1),1);
        %[P, ~] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
        P = scene.terrain.BV;
        
        tol = scene.agents(i).radius;%Tols(i) + Tols(j);
        dist_is_good = 0;
        alpha_count = scene.smoothDistAlpha/10;
        alpha_val = scene.smoothDistAlpha;

        while dist_is_good==0
            [D,GP] = soft_distance(alpha_val,P, A11);

            if(D > -1e-8)
                dist_is_good =1;
            end
            if(isinf(D))
                sprintf("Error: hessfcn, map energy, D is inf");
            end
            alpha_val = alpha_val+alpha_count;
        end

        JG1 = J1'*reshape(GP', size(GP,1)*size(GP,2), 1);       
        

        e = e + -K(i)*log((-tol + D));
        Hi = (K(i)*(1/((-tol + D).^2))) * (JG1 * JG1');
        %e = e + -K(i)*log((-tol + D).^2);
        %Hi = (K(i)*(2/((-tol + D).^2))) * (JG1 * JG1');


        %double up on Hi because I don't want to write a new fxn
        [~, spHi] = sparsify_dense_hessian(A1, A1, Hi, Hi, 0.2); 

        HB{i} = HB{i} + spHi;
            
             
    end
end

function [spHi, spHj] = sparsify_dense_hessian(X, V, Hi, Hj, cutoff_percent)
    %plot soft distance to get a sense of what is going on
     dx = (X(:,1) - V(:,1)');
     dy = (X(:,2) - V(:,2)');
     dz = (X(:,3) - V(:,3)');
    
     d = sqrt(dx.^2 + dy.^2 + dz.^2);
     Dl = sparse(d < max(max(d))*cutoff_percent);
     spHi = sparse(Hi.* kron(Dl, ones(3)));
     spHj = sparse(Hj.* kron(Dl, ones(3)));
     
end
    

