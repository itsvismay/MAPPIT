function [c, ceq] = nonlinear_constraints(q, scene)
    c = 0;
    w = 0;
    num_agents = numel(scene.agents);
    Q = reshape(q, numel(q)/num_agents, num_agents);
    GQ = zeros(size(Q));
    
    
    %intersection constraints
    for i=1:num_agents
        for j=i+1:num_agents
            A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
            B1 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            GA = zeros(size(A1));
            GB = zeros(size(B1));
            
            %if time is monotically increasing
            %sometimes, its not, for whatever reason
            if(all(diff(A1(:,3))>0) && all(diff(B1(:,3))>0))
                [A1, E1, lindA1] = sample_points_for_rod(A1, 30);
                [A2, E2, lindA2] =  sample_points_for_rod(B1, 30);
%                 A1 = A1;
%                 A2 = B1;
%                 lindA1 = linspace(1,size(A1,1), size(A1,1))';
%                 lindA2 = linspace(1,size(A2,1), size(A2,1))';
%                 E1 = [(1:size(A1,1)-1)' (2:size(A1,1))'];
%                 E2 = [(1:size(A2,1)-1)' (2:size(A2,1))'];
                % min dist >= tol ----> 0 >= tol- mindist
                tol = scene.agents(i).radius + scene.agents(j).radius;

                %[D, G]= soft_distance(1, A1, A2);
                %D
                %c = c + tol-D;
                
%                 GQ(:, i) = G;
%                 GQ(:, j) = -G;
                
                d = (A1(:,1) - A2(:,1)').^2 + (A1(:,2) - A2(:,2)').^2 + (A1(:,3) - A2(:,3)').^2; 
                min_d = min(d(:));
                c = c +  tol^2 - min_d; % tol<=D -> tol - D <=
%                 GA(lindA1,:) = -G;
%                 GB(lindA2,:) = G;
%                 GQ(:,i) = GQ(:,i)+ reshape(GA', size(GA,1)*size(GA,2), 1);
%                 GQ(:,j) = GQ(:,j)+ reshape(GB', size(GB,1)*size(GB,2), 1);
                
                
            else
                c = c + 100;
            end

           
        end
    end
    GC = reshape(GQ, size(GQ,1)*size(GQ,2),1);
    c = c+ w;
    c
    ceq = [];
    GCeq = [];
    
end

function [c, G] = terrain_constraints()
%     %terrain constraints
    % for each rod point, check if within the terrain
    % project into the  terrain if outside
%      for j=1:num_agents
%         P = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
%         if(all(diff(P(:,3))>0))
%             P = sample_points_for_rod(P, 10, 'cylinders', scene.agents(j).radius);
%             %add radius (-x,  to make signed_
%             P(:, 3) = zeros(size(P,1),1);
%             W = signed_distance(P, scene.terrain.V, scene.terrain.F, 'SignedDistanceType', 'winding_number');
%             w = w + sum(W);
%         else
%             w = w + 100;
%         end
%     end
end

