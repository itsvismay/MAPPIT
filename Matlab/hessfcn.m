function [H] = hessfcn(q_i,lambda)
    global num_agents scene simple_sd mu_barrier;
    q = q_i;%TODO: (1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = [];%TODO%Tols = q_i(end-num_agents+1:end);
    
    %Weights
    K_agent = 1*scene.coeff_matrix(1,:);
    K_tol =   0*scene.coeff_matrix(2,:);
    K_accel = 0*scene.coeff_matrix(3,:);
    K_map =   0*scene.coeff_matrix(4,:);
    K_ke =    10*scene.coeff_matrix(5,:);
    K_pv =    0*scene.coeff_matrix(6,:);
    K_reg =   1;
    
    if simple_sd
        [e_agent, HB_agent, hw_agent] = agent_agent_energy_with_hessian(Q, Tols, scene, K_agent);
    else
        [e_agent, HB_agent, hw_agent] = agent_agent_energy_with_hessian_abhisheks(Q, Tols, scene, K_agent); 
    end
    
    %[e_map, H_map] = map_energy_with_hessian(Q, Tols, scene, K_map);
    [e_reg, H_reg] = regularizer_energy_with_hessian(Q, scene, K_reg);
    [e_ke1, H_ke1] = kinetic_energy_with_hessian(Q, scene, K_ke);
    [H_acc] = mex_accel_hessian(Q(:), num_agents, size(Q,1)/3, K_accel);
    
    H = sparse(size(q_i,1),size(q_i,1));
    for ii=1:num_agents
        H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H_ke1(:,:,ii)...
                                                                            + H_reg(:,:,ii)...
                                                                            + mu_barrier*HB_agent{ii};
                                                                            %+ mu_barrier*H_map{ii};

        % tolerance_energy_with_hessian
       %H(end-ii+1,end-ii+1) = K_tol(end-ii+1) + hw_agent(end-ii+1, end-ii+1); 
    end
   H = H + sparse(H_acc);
   
%     [H_ke] = mex_kinetic_hessian(Q(:), num_agents, size(Q, 1)/3, K_ke);
%     unused = [0;0];
%     [H_tol] = mex_tol_hessian(Tols, unused, K_tol);
    
%     spH = 1e-5*eye(size(q_i,1),size(q_i,1));
%     spH(1:end-num_agents, 1:end-num_agents) = spH(1:end-num_agents, 1:end-num_agents)+  H_ke;
%     spH(end-num_agents+1:end, end-num_agents+1:end) = spH(end-num_agents+1:end, end-num_agents+1:end) + H_tol;
%     
%     for ii=1:num_agents
%         spH(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = spH(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) + H_reg(:,:,ii);
%         %H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) +  HB_agent(:,:,ii);
%         % tolerance_energy_with_hessian
%         %H(end-ii+1,end-ii+1) = K_tol(end-ii+1) + hw_agent(end-ii+1, end-ii+1); 
%     end
end

function [e, HT] = acceleration_energy_with_hessian(Q, scene, K)
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
        
        HT(:,:,i) = K(i)*d2Edq2;
    end
end


function [e, HT] = kinetic_energy_with_hessian(Q, scene, K)
    if sum(K)==0
        e=0;
        HT = zeros(size(Q,1), size(Q,1), size(Q,2));
        return;
    end
    HT = zeros(size(Q,1), size(Q,1), size(Q,2));
    e=0;
    for i=1:numel(scene.agents)

        q_i = Q(:, i); %3*nodes
        m = scene.agents(i).mass;
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        d2Edq2 = zeros(numel(q_i),numel(q_i));
        e = e + K*sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        for kk=1:numel(q_i)/3-1
            d2Edq2((kk-1)*3+1,(kk-1)*3+1) = d2Edq2((kk-1)*3+1,(kk-1)*3+1) + m./dx(kk,3);
            d2Edq2((kk-1)*3+2,(kk-1)*3+2) = d2Edq2((kk-1)*3+2,(kk-1)*3+2) + m./dx(kk,3);
            d2Edq2((kk-1)*3+3,(kk-1)*3+3) = d2Edq2((kk-1)*3+3,(kk-1)*3+3) + m*(dx(kk,1).*dx(kk,1) + dx(kk,2).*dx(kk,2))./(dx(kk,3).*dx(kk,3).*dx(kk,3));

            d2Edq2((kk-1)*3+1,(kk-1)*3+4) = d2Edq2((kk-1)*3+1,(kk-1)*3+4) - m./dx(kk,3);
            d2Edq2((kk-1)*3+2,(kk-1)*3+5) = d2Edq2((kk-1)*3+2,(kk-1)*3+5) - m./dx(kk,3);
            d2Edq2((kk-1)*3+3,(kk-1)*3+6) = d2Edq2((kk-1)*3+3,(kk-1)*3+6) - m*(dx(kk,1).*dx(kk,1) + dx(kk,2).*dx(kk,2))./(dx(kk,3).*dx(kk,3).*dx(kk,3));

            d2Edq2((kk-1)*3+4,(kk-1)*3+1) = d2Edq2((kk-1)*3+4,(kk-1)*3+1) - m./dx(kk,3);
            d2Edq2((kk-1)*3+5,(kk-1)*3+2) = d2Edq2((kk-1)*3+5,(kk-1)*3+2) - m./dx(kk,3);
            d2Edq2((kk-1)*3+6,(kk-1)*3+3) = d2Edq2((kk-1)*3+6,(kk-1)*3+3) - m*(dx(kk,1).*dx(kk,1) + dx(kk,2).*dx(kk,2))./(dx(kk,3).*dx(kk,3).*dx(kk,3));

            d2Edq2((kk-1)*3+4,(kk-1)*3+4) = d2Edq2((kk-1)*3+4,(kk-1)*3+4) + m./dx(kk,3);
            d2Edq2((kk-1)*3+5,(kk-1)*3+5) = d2Edq2((kk-1)*3+5,(kk-1)*3+5) + m./dx(kk,3);
            d2Edq2((kk-1)*3+6,(kk-1)*3+6) = d2Edq2((kk-1)*3+6,(kk-1)*3+6) + m*(dx(kk,1).*dx(kk,1) + dx(kk,2).*dx(kk,2))./(dx(kk,3).*dx(kk,3).*dx(kk,3));
        end
               
        HT(:,:,i) = K(i)*d2Edq2;
            
    end
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
        e = e + K*(sum(kt./dt));
       
        d2Edq2 = zeros(numel(q_i),numel(q_i));
        
        for kk=1:numel(q_i)/3-1
            d2Edq2((kk-1)*3+3,(kk-1)*3+3) = d2Edq2((kk-1)*3+3,(kk-1)*3+3) + 2*kt(kk,1)./(dt(kk,1).^3);

            d2Edq2((kk-1)*3+3,(kk-1)*3+6) = d2Edq2((kk-1)*3+3,(kk-1)*3+6) - 2*kt(kk,1)./(dt(kk,1).^3);

            d2Edq2((kk-1)*3+6,(kk-1)*3+3) = d2Edq2((kk-1)*3+6,(kk-1)*3+3) - 2*kt(kk,1)./(dt(kk,1).^3);

            d2Edq2((kk-1)*3+6,(kk-1)*3+6) = d2Edq2((kk-1)*3+6,(kk-1)*3+6) + 2*kt(kk,1)./(dt(kk,1).^3);
        end
               
        HT(:,:,i) = K*d2Edq2;
            
    end
end

function [e, HB, hW] = agent_agent_energy_with_hessian_abhisheks(Q, Tols, scene, K)
    if sum(K)==0
        e=0;
        HB = zeros(size(Q,1), size(Q,1), size(Q,2));
        hW = zeros(size(Q,2), size(Q,2));
        return;
    end
   
    HB = zeros(size(Q,1), size(Q,1), size(Q,2));
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
            dist_is_good = 0;
            
            alpha_count = 1;
            alpha_val = 1;
            beta_val = 0.0;
            while dist_is_good==0     
                [B1,I1] = build_distance_bvh(A11,[]);
                [B2,I2] = build_distance_bvh(A22,[]);
                [~, G1] = smooth_min_distance(A22,[],B2,I2,alpha_val,A11,[],alpha_val, beta_val);
                [D, G2] = smooth_min_distance(A11,[],B1,I1,alpha_val,A22,[],alpha_val, beta_val);
                
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end

%             newH1 = reshape(H1', 3, size(H1,1)*3)';
%             newH2 = reshape(H2', 3, size(H2,1)*3)';
%             H1blocks = mat2cell(newH1, repmat(3, size(H1,1),1), [3]);
%             H2blocks = mat2cell(newH2, repmat(3, size(H2,1),1), [3]);
%             spH1 = sparse(blkdiag(H1blocks{:}));
%             spH2 = sparse(blkdiag(H2blocks{:}));
%             JH1J = J1'*spH1*J1;
%             JH2J = J2'*spH2*J2;
            JG1 = J1'*reshape(G1', size(G1,1)*size(G1,2), 1);
            JG2 = J2'*reshape(G2', size(G2,1)*size(G2,2), 1);
            
            tol = Tols(i) + Tols(j);
            
            e = e + -K(i)*log((-tol + D)^2);
            
            %HB(:,:,i) = HB(:,:,i) + (K(i)*(2/((-tol + D).^2))) * JH1J;
            %HB(:,:,j) = HB(:,:,j) + (K(j)*(2/((-tol + D).^2))) * JH2J;
            HB(:,:,i) = HB(:,:,i)+ JG1' * (K(i)*(2/((-tol + D).^2))) * JG1;
            HB(:,:,j) = HB(:,:,j)+ JG2' * (K(j)*(2/((-tol + D).^2))) * JG2;
            
            
            hW(i) = hW(i) + K(i)*(-1/((-tol+D).^2));
            hW(j) = hW(j) + K(j)*(-1/((-tol+D).^2));
            
        end        
    end
end
   
function [e, HB, hW] = agent_agent_energy_with_hessian(Q, Tols, scene, K)
    if sum(K)==0
        e=0;
        HB = zeros(size(Q,1), size(Q,1), size(Q,2));
        hW = zeros(size(Q,2), size(Q,2));
        return;
    end
   
    HB(1:size(Q,2)) = {sparse(size(Q,1), size(Q,1))};
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
            dist_is_good = 0;
            alpha_count = 10;
            alpha_val = 10;
            while dist_is_good==0
                [~,G1] = soft_distance(alpha_val,A22, A11);
                [D,G2] = soft_distance(alpha_val,A11, A22);
                      
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
            JG1 = J1'*reshape(G1', size(G1,1)*size(G1,2), 1);
            JG2 = J2'*reshape(G2', size(G2,1)*size(G2,2), 1);  
            
            tol = 0.75;%Tols(i) + Tols(j);
            
            %Old Energy
            %e = e + -K(i)*log((-tol + D).^2);
            %Hi = (K(i)*(2/((-tol + D).^2))) * (JG1 * JG1');
            %Hj = (K(j)*(2/((-tol + D).^2))) * (JG2 * JG2');
            
            %New Energy
            e = e + -K(i)*log((-tol + D));
            Hi = (K(i)*(1/((-tol + D).^2))) * (JG1 * JG1');
            Hj = (K(j)*(1/((-tol + D).^2))) * (JG2 * JG2');
            
            
            
            [spHi, spHj] = sparsify_dense_hessian(A1, A2, Hi, Hj, 0.2);
            
            HB{i} = HB{i} + spHi;
            HB{j} = HB{j} + spHj;
            
            %hW(i) = hW(i) + K(i)*(2/((-tol+D).^2));TODO
            %hW(j) = hW(j) + K(j)*(2/((-tol+D).^2));
            
        end        
    end
end

function [e, HB, hW] = map_energy_with_hessian(Q, Tols, scene, K)
    if sum(K)==0
        e=0;
        HB = zeros(size(Q,1), size(Q,1), size(Q,2));
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

        JG1 = J1'*reshape(GP', size(GP,1)*size(GP,2), 1);       
        tol = 0.75;%Tols(i) + Tols(j);

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
    

