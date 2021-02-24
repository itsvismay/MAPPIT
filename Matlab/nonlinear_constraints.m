function [c,ceq, GC,GCeq] = nonlinear_constraints(q_i,scene)
    c = 0;
    w = 0;
    num_agents = numel(scene.agents);
    q = q_i(1:end-numel(scene.agents));
   
    Q = reshape(q, numel(q)/num_agents, num_agents);  
    GQ = zeros(size(Q));
    Tols = q_i(end-num_agents+1:end);
    
    %intersection constraints
    for i=1:num_agents
        for j=i+1:num_agents
            A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
                
            [A1, J1] = sample_points_for_rod(A1, scene.agents(i).e);
            [A2, J2] = sample_points_for_rod(A2, scene.agents(j).e);

            %true min D
            [D, G2]= soft_distance( 200, A1, A2);
            [D, G1]= soft_distance( 200, A2, A1);
  
            tol = Tols(i) + Tols(j);
            
            c = max(c, tol-D);
            
            JG1 = -J1'*G1;
            JG2 = -J2'*G2;
            GQ(:, i) = GQ(:, i) + reshape(JG1', numel(JG1),1);
            GQ(:, j) = GQ(:, j) + reshape(JG2', numel(JG2),1);
        end
       
        P = scene.terrain.BV;
        [P, J1] = sample_points_for_rod(P, scene.agents(i).e);
        %add radius (-x,  to make signed_
        A1(:, 3) = zeros(size(A1,1),1);
        P = scene.terrain.BV;
        [D,GP] = soft_distance(100,P, A1);
%       JG1 = J1'*GP;

        tol = Tols(j);
        if(tol-D>0)
            w = max(w, tol - D);
        end
        
    end
    
    %GC = [reshape(GQ, size(GQ,1)*size(GQ,2),1); 0];
    c= max(c,w);
    
    ceq = [];
    GCeq = [];
end

