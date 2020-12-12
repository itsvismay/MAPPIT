%read in from the scene file
fname = "scene_2/";
setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.agents = [];

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];

a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 5;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    %edges
    r1e = [(1:agent.segments)' (2:(agent.segments+1))'];
    agent.e = r1e;
    r1e = r1e + size(v,1);
    e = [e; r1e];
    
    % %vertices
    r1v = [linspace(agent.xse(1,1),agent.xse(end,1), agent.segments+1)', ... %x
            linspace(agent.xse(1,2),agent.xse(end,2),agent.segments+1)', ... %y
            linspace(agent.xse(1,3),agent.xse(end,3),agent.segments+1)'];    %t
    
    %wiggles the rod start so that they aren't intersecting
    endtime = r1v(end,3);
    r1v(:,3) = sort(rand(1,size(r1v,1))*(endtime));
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];
    
    %fix end points
    num_constraints = size(agent.xse,1)*2 +1;
    Aeq_rows = [1:num_constraints];
    Aeq_cols = [1:3 reshape([3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 1; 
                3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 2], 1,[])];
    Aeq_vals = ones(1,num_constraints);
    
    A1eq = sparse(Aeq_rows, Aeq_cols, Aeq_vals, num_constraints, numel(r1v));
    b1eq = [agent.xse(1,:)'; reshape(agent.xse(2:end,1:2)', 1, [])'];
    
    Aeq = [Aeq zeros(size(Aeq,1), size(A1eq,2)); 
            zeros(size(A1eq,1), size(Aeq,2)) A1eq];
        
    beq = [beq; b1eq];

    %time is monotonic constraints 
    %t_i+1 - t_i > 0 --> in this format --> A1*q <= 0
    A1 = sparse(reshape([1:agent.segments; 1:agent.segments], 1, 2*agent.segments),...
                reshape([3:3:(numel(r1v)-3); 6:3:(numel(r1v))], 1, 2*agent.segments),...
                reshape([ones(1, agent.segments); -ones(1,agent.segments)], 1, 2*agent.segments),... 
                agent.segments, numel(r1v));
    b1 = zeros(agent.segments, 1);

    %max time via inequality constraint (comment out to turn off, but make sure
    %you uncomment the end time potential energy in cost)
    A1 = [A1; zeros(1, 3*agent.segments+2) 1];
    b1 = [b1; agent.max_time];
    
    A = [A zeros(size(A,1), size(A1,2));
        zeros(size(A1,1), size(A,2)) A1];
    b = [b; b1];
    
    scene.agents = [scene.agents agent];   
end

PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
axis equal;
drawnow;

%minimize here

options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true);

options.MaxFunctionEvaluations = 1e4;

[q, fval, exitflag, output] = fmincon(@(x) path_energy(x, numel(scene.agents), e, surf_anim), reshape(v', numel(v),1), A,b,Aeq,beq,[],[], @(x)constraints(x, scene), options);

PV = reshape(q, 3, numel(q)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',2, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;


% %PRINT AGENTS 
print_agents(fname+"agents_test.json",scene, q);

function p = print_agents(fname, scene, q)
    Q = reshape(q, numel(q)/numel(scene.agents), numel(scene.agents));

    for i=1:size(Q,2)
        scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
    end
    jsontext = jsonencode(scene);
    jsontext = strrep(jsontext, ',', sprintf(',\n'));
    jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
    jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
    fileID = fopen(fname, 'w');
    fprintf(fileID, jsontext);
end

function [f,g] = path_energy(q, num_agents, e, surf_anim)
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    G = zeros(size(Q));
    
    T = 0;
    V = 0;
    V2 = 0;
    
    g = zeros(size(q));

    for i=1:num_agents
        q_i = Q(:, i); %3*nodes
        
        %various fun path energies, I'll use principle of least action becase I like it
        %kinetic energy of curve integrated along piecewise linear segment is 
        % 0.5*dx'*dx./dt 
        m = 1;

        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        T = T + sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        
        dEdq_left = zeros(numel(q_i)/3, 3);
        dEdq_right = zeros(numel(q_i)/3, 3);
        
        dEdq_left(1:end-1,1) = -m*dx(:,1)./dx(:,3);
        dEdq_left(1:end-1,2) = -m*dx(:,2)./dx(:,3);
        dEdq_left(1:end-1,3) = 0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));
        
        dEdq_right(2:end, 1) = m*dx(:,1)./dx(:,3);
        dEdq_right(2:end, 2) = m*dx(:,2)./dx(:,3);
        dEdq_right(2:end, 3) = -0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));
        
        dEdq = dEdq_left + dEdq_right;
        
        G(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
    end
    g = reshape(G, size(G,1)*size(G,2),1);
    
    %find minimum energy curve
    f = T;
    
%     T = abs(T)
%     PV = reshape(q, 3, numel(q)/3)';
%     PE = e;
%     [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
%     surf_anim.Vertices = CV;
%     drawnow;
    

end

function [c, ceq] = constraints(q, scene)
    c = 0;
    DC = zeros(size(q));
    w = 0;
    num_agents = numel(scene.agents);
    Q = reshape(q, numel(q)/num_agents, num_agents);
    GQ = zeros(size(Q));
    
    
%     %terrain constraints
    % for each rod point, check if within the terrain
    % project into the  terrain if outside
%      for j=1:num_agents
%         P = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
%         if(all(diff(P(:,3))>0))
%             P = SamplePointsForRod(P, 10, 'cylinders', scene.agents(j).radius);
%             %add radius (-x,  to make signed_
%             P(:, 3) = zeros(size(P,1),1);
%             W = signed_distance(P, scene.terrain.V, scene.terrain.F, 'SignedDistanceType', 'winding_number');
%             w = w + sum(W);
%         else
%             w = w + 100;
%         end
%     end
    
    
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
                [A1 lindA1] = SamplePointsForRod(A1, 30);
                [B1 lindB1] =  SamplePointsForRod(B1, 30);
                %lindA1 = linspace(1,size(A1,1), size(A1,1))';
                %lindB1 = linspace(1,size(B1,1), size(B1,1))';
                % min dist >= tol ----> 0 >= tol- mindist
                tol = scene.agents(i).radius + scene.agents(j).radius;

                %[D, G]= smooth_min_distance(A1, B1);
                d = (A1(:,1) - B1(:,1)').^2 + (A1(:,2) - B1(:,2)').^2 + (A1(:,3) - B1(:,3)').^2;
                d = sqrt(d);
                D = min(d(:));
                
                c = c +  tol - D; % tol<=D -> tol - D <=0
                    
                

%                 GA(lindA1,:) = -G;
%                 GB(lindB1,:) = G;
% 
%                 GQ(:,i) = GQ(:,i)+ reshape(GA', size(GA,1)*size(GA,2), 1);
%                 GQ(:,j) = GQ(:,j)+ reshape(GB', size(GB,1)*size(GB,2), 1);
                
                
            else
                c = c + 100;
            end

           
        end
    end
    DC = reshape(GQ, size(GQ,1)*size(GQ,2),1);
    c = c+ w;
    
    ceq = [];
    DCeq = [];
    
end

function [P, lind] = SamplePointsForRod(X, samples, varargin)
    ts = linspace(X(1,3)+0.001, X(end,3)-0.001, samples);
    PV = [interp1(X(:,3), X(:, 1:2), ts, 'linear') ts'];

    if numel(varargin)>0 && strcmp('cylinders',varargin{1})
        rad = varargin{2}; 
        PE = [(1:size(PV,1)-1)' (2:size(PV,1))'];
        [P,~,~,~] = edge_cylinders(PV,PE, 'Thickness',2*rad, 'PolySize', 4);
    else
        % samples along centerline
        ts = linspace(X(1,3)+0.001, X(end,3)-0.001, samples);
        P = PV;
    end
    [trash, trash, lind] = histcounts(ts, X(:,3)');
   
end

