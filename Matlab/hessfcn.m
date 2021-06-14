function [H] = hessfcn(q_i,lambda)
    global num_agents scene simple_sd;
      q = q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = q_i(end-num_agents+1:end);
    
    %Weights
    K_agent = 1*scene.coeff_matrix(1,:);
    K_tol =   1*scene.coeff_matrix(2,:);
    K_accel = 0*scene.coeff_matrix(3,:);
    K_map =   0*scene.coeff_matrix(4,:);
    K_ke =    1*scene.coeff_matrix(5,:);
    K_pv =    0*scene.coeff_matrix(6,:);
    K_reg =   1;
    
   
    [e_ke, H_ke] = kinetic_energy_with_hessian(Q, scene, K_ke);
    [e_reg, H_reg] = regularizer_energy_with_hessian(Q, scene, K_reg);
   
    if simple_sd
        [e_agent, HB_agent, hw_agent] = agent_agent_energy_with_hessian(Q, Tols, scene, K_agent);
    else
        [e_agent, HB_agent, hw_agent] = agent_agent_energy_with_hessian_abhisheks(Q, Tols, scene, K_agent); 
    end
   
    H = zeros(size(q_i,1),size(q_i,1));
    for ii=1:num_agents
        H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H_ke(:,:,ii) + H_reg(:,:,ii);
        H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) +  HB_agent(:,:,ii);

            % tolerance_energy_with_hessian
            H(end-ii+1,end-ii+1) = K_tol(end-ii+1) + hw_agent(end-ii+1, end-ii+1); 
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
            
            alpha_count = 10;
            alpha_val = 10;
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
            
            
            
            tol = Tols(i) + Tols(j);
            
            e = e + -K(i)*log((-tol + D).^2);
%             if(ismember(j,scene.agents(i).friends))
%                 e = e + -K(i)*log(-D + 2);
%                 HB(:,:,i) = HB(:,:,i) + reshape_JG1' * (K(i)*(1/((-D + 2).^2))) * reshape_JG1;
%                 HB(:,:,j) = HB(:,:,j) + reshape_JG2' * (K(j)*(1/((-D + 2).^2))) * reshape_JG2;
%             end
            
            HB(:,:,i) = HB(:,:,i)+ JG1' * (K(i)*(2/((-tol + D).^2))) * JG1;
            HB(:,:,j) = HB(:,:,j)+ JG2' * (K(j)*(2/((-tol + D).^2))) * JG2;
            
            hW(i) = hW(i) + K(i)*(2/((-tol+D).^2));
            hW(j) = hW(j) + K(j)*(2/((-tol+D).^2));
            
        end        
    end
end
    
    

