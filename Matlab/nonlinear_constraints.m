function [c, ceq, gc, gceq] = nonlinear_constraints(q_i, UserTols, scene)
    ceq = [];
    gceq = [];
    
    num_agents = numel(scene.agents);
    q = q_i(1:end-numel(scene.agents));
    Q = reshape(q, numel(q)/num_agents, num_agents);
    %hard coded tolerance, just like in Abhisheks version
    Tols = q_i(end-num_agents+1:end); 
    
    %doing multiple constraints since I *know* this works in Abhisheks
    %implementation
    c = [0];
    gc = [zeros(numel(Q) + num_agents, 1)];
      
    interaction_count = 1;
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
            alpha_count = 1;
            alpha_val = 10;
            while dist_is_good==0
                % vismay code
                [~,G1] = soft_distance(alpha_val,A22, A11);
                [D,G2] = soft_distance(alpha_val,A11, A22);
                %use sampled rods
%                 [B1,I1] = build_distance_bvh(A11,E11);
%                 [B2,I2] = build_distance_bvh(A22,E22);
%                 [~, G1] = smooth_min_distance(A22,[],B2,I2,alpha_val,A11,[],alpha_val);
%                 [D, G2] = smooth_min_distance(A11,[],B1,I1,alpha_val,A22,[],alpha_val);
                % abhisheks code
%                 [B1,I1] = build_distance_bvh(A1,scene.agents(i).e);
%                 [B2,I2] = build_distance_bvh(A2,scene.agents(j).e);
%                 [~, G1] = smooth_min_distance(A2,scene.agents(j).e,B2,I2,alpha_val,A11,[],alpha_val);
%                 [D, G2] = smooth_min_distance(A1,scene.agents(i).e,B1,I1,alpha_val,A22,[],alpha_val);
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
                        
            tol = 1 + 1;%UserTols(i) + UserTols(j); %hard coded tolerances to 1
            c(interaction_count) = tol - D;
            %tol-D
            if isnan(sum(sum(G1))) | isnan(sum(sum(G2)))
                D
                alpha_val
            end
            
            JG1 = J1'*G1;
            JG2 = J2'*G2;
            
            gW = zeros(num_agents,1);%tolerance gradient
            GB = zeros(size(Q)); %agent-agent interaction gradient
            GB(:,i) = GB(:,i)+ -1*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
            GB(:,j) = GB(:,j)+ -1*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
            gW(i) = 1;%hard coded tolerances to 1
            gW(j) = 1;%hard coded tolerances to 1
            
            g_col = zeros(numel(Q) + num_agents, 1);
            g_col(1:end-num_agents) = reshape(GB, size(GB,1)*size(GB,2),1);
            g_col(end-num_agents+1:end) = gW;
            if isnan(sum(g_col))
                D
            end
            gc(:, interaction_count) = g_col;
            interaction_count = interaction_count + 1;
        end        
    end
    % Do a smooth maximum of tol-D
    % TODO: What's a good way to modify alpha here?
    alpha = 5;
    exp_weights = exp(alpha*c);
    total = sum(exp_weights);
    gc = (gc.*exp_weights)./total;
    gc = sum(gc,2);
    if isnan(sum(gc))
        exp_weights;
    end
    c = 1/alpha*log(total);
   
end

