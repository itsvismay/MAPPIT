%read in from the scene file
addpath("../external/smooth-distances/build/");
%fname = "../Scenes/output_results/eight_agents/agent_circle/";
%fname = "../Scenes/output_results/three_agents/test/";
fname = "../Scenes/output_results/scaling_tests/2_agents/";
%fname = "../Scenes/output_results/scaling_tests/test/";

setup_params = jsondecode(fileread(fname+"setup.json"));
global scene;
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BF = boundary_faces(tF);
scene.terrain.BVind = unique(scene.terrain.BF);
scene.terrain.BV = tV(scene.terrain.BVind,:);
scene.agents = [];


v = [];
e = [];
el = [];
Aeq1 = [];
beq1 = [];
Aeq2 = [];
beq2 = [];
A = [];
b = [];
UserTols = [];
coefficients_matrix = zeros(6, numel(fieldnames(setup_params.agents)));

%% SETUP DIJIKSTRAS
AdjM = adjacency_matrix(scene.terrain.F);
AdjM_visited = AdjM;
nLayer = 15;
nTotalVer = nLayer * length(scene.terrain.V(:,1));
% build up 3d graph
% find all the edges in 2d graph
A_lt = tril(AdjM);
[edge_s,edge_t] = find(A_lt);
edge = zeros(length(edge_s),2);
edge(:,1) = edge_s;
edge(:,2) = edge_t;
    
% set the num of layer and get the spacetime graph 
% (vertices: VV and edges: EE)
time = linspace(0,1,nLayer)';
[VV,EE, newA] = spacetime_graph(scene.terrain.V,edge,time);


% make the graph directed along the time dimension
[ii,jj,ss] = find(newA);
for k=1:length(ii)
   s_temp = ii(k);
   t_temp = jj(k);
   if VV(s_temp,3) > VV(t_temp,3)
       %newA(s_temp, t_temp) = 0;
   end
end

newA_visited = newA;
%%

a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.mass = getfield(setup_params.agents, a{i}).mass;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    i
    agent.seg_per_waypoint = 50;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    agent.mass = getfield(setup_params.agents, a{i}).mass;
    agent.friends = getfield(setup_params.agents, a{i}).friends;
    agent.mesh = getfield(setup_params.agents, a{i}).mesh;
    agent.animation_cycles = getfield(setup_params.agents, a{i}).animation_cycles;
    coefficients_matrix(1, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_agent;
    coefficients_matrix(2, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_tol;
    coefficients_matrix(3, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_accel;
    coefficients_matrix(4, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_map;
    coefficients_matrix(5, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_ke;
    coefficients_matrix(6, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_pv;
    
    
    %[r1e, r1v, newA, newA_visited] = set_path3d(newA, newA_visited, agent, scene, VV, EE, nLayer, nTotalVer);
    [r1e, r1v, newA, newA_visited] = set_path3d(newA, newA_visited, agent, scene, VV, EE, nLayer, nTotalVer, agent.segments);
    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    
    e = [e; r1e];
    
    %wiggles the rod start so that they aren't intersecting
    starttime = r1v(1,3);
    endtime = r1v(end,3);
    %r1v(:,3) = ones();%sort(rand(1,size(r1v,1))*(endtime));%r1v(:,3)/i;%
    % smoothing_eps = k*[1,2,3....]
    %adds time to make sure no div by 0 (flat paths)
    smoothing_eps = 1e-1*linspace(1,size(r1v,1), size(r1v,1))'; 
    r1v(:,3) = r1v(:,3) + smoothing_eps;
    r1v(1,3) = starttime;
    %r1v(end,3) = endtime;
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
    
    %fix end points constraint set 1
        num_constraints = size(agent.xse,1)*3;
        Aeq_rows = [1:num_constraints];
        Aeq_cols = [1:3 reshape([3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 1; 
                                3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 2 ;
                                3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 3], 1,[])];
        Aeq_vals = ones(1,num_constraints);
        A1eq = sparse(Aeq_rows, Aeq_cols, Aeq_vals, num_constraints, numel(r1v));
        b1eq = [r1v(1,:)'; reshape(r1v(end,:)', 1, [])'];
        Aeq1 = [Aeq1 zeros(size(Aeq1,1), size(A1eq,2)); 
                zeros(size(A1eq,1), size(Aeq1,2)) A1eq];
        beq1 = [beq1; b1eq];
    
    %fix end points constraint set 2
        num_constraints = size(agent.xse,1)*2 +1;
        Aeq_rows = [1:num_constraints];
        Aeq_cols = [1:3 reshape([3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 1; 
                                3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 2], 1,[])];
        Aeq_vals = ones(1,num_constraints);
        A2eq = sparse(Aeq_rows, Aeq_cols, Aeq_vals, num_constraints, numel(r1v));
        b2eq = [r1v(1,:)'; reshape(r1v(end,1:2)', 1, [])'];
        
        Aeq2 = [Aeq2 zeros(size(Aeq2,1), size(A2eq,2)); 
                zeros(size(A2eq,1), size(Aeq2,2)) A2eq];
        beq2 = [beq2; b2eq];

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
%Aeq1 = [Aeq1 zeros(size(Aeq1,1),numel(scene.agents))];
%Aeq2 = [Aeq2 zeros(size(Aeq2,1),numel(scene.agents))];
%A = [A zeros(size(A,1),numel(scene.agents))];
scene.coeff_matrix = coefficients_matrix;
global num_agents
num_agents = numel(scene.agents);

%%
PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
hold on;
axis equal;
drawnow;
qn = reshape(v', numel(v),1);

options = optimoptions('fmincon','Algorithm','interior-point',...
    "SpecifyConstraintGradient",true,...
    "SpecifyObjectiveGradient",true,...
    'HessianFcn',@hessfcn);

%minimize with KE here
%user constraint set 2
[q_i, fval,exitflag,output,lambda, grad,hessian] = fmincon(@(x) path_energy_hessian(x,UserTols, numel(scene.agents),scene, e, surf_anim, 2),... 
                            qn, ...
                            A,b,Aeq2,beq2,[],[], ...
                            [], options);
qn = q_i;%q_i(1:end-numel(scene.agents));
Q = reshape(qn, numel(qn)/numel(scene.agents), numel(scene.agents));
scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
PV = reshape(qn, 3, numel(qn)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',1, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;



print_agents(fname+"agents.json", scene, qn)

%path_energy(q_i,UserTols, numel(scene.agents),scene, e, surf_anim)
