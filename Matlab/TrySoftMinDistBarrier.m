%read in from the scene file
addpath("../external/smooth-distances/build/");
fname = "../Scenes/output_results/eight_agents/agent_circle/";
%fname = "../Scenes/output_results/three_agents/test/";

setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BF = boundary_faces(tF);
scene.terrain.BVind = unique(scene.terrain.BF);
scene.terrain.BV = tV(scene.terrain.BVind,:);
scene.agents = [];

global Kw Kt Ka
Kw = 10; %tolerancee
Kt = 1; % kinetic
Ka = 1; %agent-agent

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];
UserTols = [];
AdjM = adjacency_matrix(scene.terrain.F);
AdjM_visited = AdjM;

a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.mass = getfield(setup_params.agents, a{i}).mass;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 10;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    
    [r1e, r1v, AdjM, AdjM_visited] = set_path(AdjM, AdjM_visited, agent, scene);
    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    
    e = [e; r1e];
    
    %wiggles the rod start so that they aren't intersecting
    endtime = r1v(end,3);
    r1v(:,3) = r1v(:,3)/i;%sort(rand(1,size(r1v,1))*(endtime));
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];
    agent.rest_edge_lengths = r1el;
    agent.rest_region_lengths = [0; r1el(1:size(r1el)-1) + r1el(2:size(r1el))];
    
    %sets up the agent bvh
    [B,I] = build_distance_bvh(agent.v,[]);
    agent.bvh.B = B;
    agent.bvh.I = I;
    
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
    UserTols = [UserTols agent.radius];
end
Aeq = [Aeq zeros(size(Aeq,1),numel(scene.agents))];
A = [A zeros(size(A,1),numel(scene.agents))];


PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
axis equal;
drawnow;

%minimize here
options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient',true, 'Display', 'iter', 'UseParallel', false, 'HessianApproximation', 'lbfgs');
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = 1e3;
q_i  = [reshape(v', numel(v),1);  -1e-8*ones(numel(scene.agents),1)];
[q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x, UserTols, numel(scene.agents), scene, e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            [], options);
q = q_i(1:end-numel(scene.agents));
                        
PV = reshape(q, 3, numel(q)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',1, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;

function [f,g] = path_energy(q_i, UserTols, num_agents, scene, e, surf_anim)
    global Kw Kt Ka;

    q = q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    Tols = q_i(end-num_agents+1:end);   

    [e_agent, g_full] = agent_agent_energy(Q, Tols, scene, Ka);
    [e_tol, g_tol] = tolerance_energy(Tols, UserTols, Kw);
    [e_ke, g_ke] = kinetic_energy(Q, scene, Kt);
    
    f = e_agent + e_tol + e_ke;
    
    g = g_full;
    g(1:end-num_agents) = g(1:end-num_agents) + g_ke;
    g(end-num_agents+1:end) = g(end-num_agents+1:end) + g_tol;
    
end

function [e, g] = agent_agent_energy(Q, Tols, scene, K)
    num_agents = numel(scene.agents);
    if sum(K)==0
        e=0;
        g = zeros(numel(Q) + num_agents, 1);
        return;
    end
    GB = zeros(size(Q));
    gW = zeros(num_agents,1);
    e=0;
    g = zeros(numel(Q) + num_agents, 1);
   
    for i=1:numel(scene.agents)
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1, E1, J1] = sample_points_for_rod(A1, scene.agents(i).e);
        for j =i+1:num_agents
            if j==i
                continue;
            end
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A2,E2, J2] = sample_points_for_rod(A2, scene.agents(j).e);
            dist_is_good = 0;
            alpha_count = 5;
            alpha_val = 50;
            while dist_is_good==0
                [~,G1] = soft_distance(alpha_val,A2, A1);
                [D,G2] = soft_distance(alpha_val,A1, A2);
                %[~, G2] = smooth_min_distance(A1,[],alpha_val,A2,[],alpha_val);
                %[D, G1] = smooth_min_distance(A2,[],alpha_val,A1,[],alpha_val);
                %[~,G2] = smooth_min_distance(A1,[],scene.agents(i).bvh.B, scene.agents(i).bvh.I,alpha_val,A2,[],alpha_val);
                %[D,G1] = smooth_min_distance(A2,[],scene.agents(j).bvh.B, scene.agents(j).bvh.I,alpha_val,A1,[],alpha_val);
                if(D>-1e-8)
                    dist_is_good =1;
                end
                alpha_val = alpha_val+alpha_count;
            end
            JG1 = J1'*G1;
            JG2 = J2'*G2;
            
            tol = Tols(i) + Tols(j);
            
            e = e + -K*log(-tol + D);
        
            GB(:,i) = GB(:,i)+ K*(-1/(-tol + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
            GB(:,j) = GB(:,j)+ K*(-1/(-tol + D))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);
            gW(i) = gW(i) + K*(1/(-tol+D));
            gW(j) = gW(j) + K*(1/(-tol+D));
            
%             B = B + 1/D;
%             GB(:,i) = GB(:,i)+ (-1/(D^2))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
%             GB(:,j) = GB(:,j)+ (-1/(D^2))*reshape(JG2', size(JG2,1)*size(JG2,2), 1);

        end        
    end
    g(1:end-num_agents) = reshape(GB, size(GB,1)*size(GB,2),1);
    g(end-num_agents+1:end) = gW;
end
function [e, g] = tolerance_energy(Tols, UserTols, K)
    e = 0.5*(K*(Tols - UserTols'))'*(Tols - UserTols');
    g = K*(Tols - UserTols');
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
        
        GT(:,i) = K*reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
    end
    g = reshape(GT, size(GT,1)*size(GT,2),1);
end
