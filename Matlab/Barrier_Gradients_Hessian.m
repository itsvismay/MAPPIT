%read in from the scene file
addpath("../external/smooth-distances/build/");
addpath("../CrowdSolverCpp/matlab/");

%Scaling Tests
fname = "../Scenes/output_results/scaling_tests/2_agents/"; nLayer = 3; 
    num_segments = 30; max_iters = 5; num_inside_iters = 30;
%fname = "../Scenes/output_results/scaling_tests/10_agents/"; nLayer = 7; 
    %num_segments = 50; max_iters = 5; num_inside_iters = 30;
%fname = "../Scenes/output_results/scaling_tests/20_agents/"; nLayer = 8;
    %num_segments = 50; max_iters = 5; num_inside_iters = 30;
%fname = "../Scenes/output_results/scaling_tests/30_agents/"; nLayer = 7;
    %num_segments = 50; max_iters = 5; num_inside_iters = 30;
%fname = "../Scenes/output_results/scaling_tests/60_agents/"; nLayer = 7;
    %num_segments = 50; max_iters = 5; num_inside_iters = 30;
%fname = "../Scenes/output_results/3x_3_agents/test/"; nLayer = 3; 
    %num_segments = 30; max_iters = 10; num_inside_iters = 50;
%fname = "../Scenes/output_results/complex_maze/square_maze/one_agent/";nLayer = 5; 
    %num_segments = 250; max_iters = 10; num_inside_iters = 20;
%fname = "../Scenes/output_results/complex_maze/square_maze/three_agents/";nLayer = 5; 
    %num_segments = 250; max_iters = 10; num_inside_iters = 20;


setup_params = jsondecode(fileread(fname+"setup.json"));
global scene num_agents simple_sd mu_barrier;
simple_sd = 1;
mu_barrier= 1;
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BF = boundary_faces(tF);
scene.terrain.BVind = unique(scene.terrain.BF);
scene.terrain.BV = tV(scene.terrain.BVind,:);
scene.agents = [];
tV(scene.terrain.BVind,3) = 1;
tsurf(tF, tV);hold on;

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];
UserTols = [];
coefficients_matrix = zeros(6, numel(fieldnames(setup_params.agents)));

%% SETUP DIJIKSTRAS
AdjM = adjacency_matrix(scene.terrain.F);
AdjM_visited = AdjM;
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
a = fieldnames(setup_params.agents);
max_time = 0;
for i = 1:numel(a)
    if getfield(setup_params.agents, a{i}).xse(end, end)> max_time
        max_time = getfield(setup_params.agents, a{i}).xse(end, end)
    end
end
if max_time == 0
    "Something is broken, max time for scene shouldn't be 0"
    return
end
time = linspace(0,max_time,nLayer)';
[VV,EE, AdjM, BVind] = spacetime_graph(scene.terrain.V,edge,time, scene.terrain.BVind);
BV = VV(BVind,:);
plot3(BV(:,1), BV(:,2), BV(:,3), 'ko');hold on;

% make the graph directed along the time dimension
[ii,jj,ss] = find(AdjM);
for k=1:length(ii)
   s_temp = ii(k);
   t_temp = jj(k);
   if VV(s_temp,3) > VV(t_temp,3)
       %newA(s_temp, t_temp) = 0;
   end
end

AdjM_visited = AdjM;

for i = 1:numel(a)
    agent.id = i
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = num_segments;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    agent.mass = getfield(setup_params.agents, a{i}).mass;
    agent.friends = getfield(setup_params.agents, a{i}).friends;
    agent.mesh = getfield(setup_params.agents, a{i}).mesh;
    agent.animation_cycles = getfield(setup_params.agents, a{i}).animation_cycles;
    %agent.mesh_direction_forward_up = getfield(setup_params.agents, a{i}).mesh_direction_forward_up;
    
    coefficients_matrix(1, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_agent;
    coefficients_matrix(2, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_tol;
    coefficients_matrix(3, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_accel;
    coefficients_matrix(4, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_map;
    coefficients_matrix(5, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_ke;
    coefficients_matrix(6, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_pv;
        
    %set path new
    %[r1e, r1v, AdjM, AdjM_visited] = set_path(AdjM, AdjM_visited, agent, scene);
    [r1e, r1v, AdjM, AdjM_visited] = set_path3d(AdjM, AdjM_visited, agent, scene, VV, EE, BV, BVind, nLayer, nTotalVer, agent.segments);

    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    e = [e; r1e]; 
    
    %wiggles the rod start so that they aren't intersecting
    starttime = r1v(1,3);
    smoothing_eps = 1e-3*linspace(1,size(r1v,1), size(r1v,1))'; 
    r1v(:,3) = r1v(:,3) + smoothing_eps;
    r1v(1,3) = starttime;
    endtime = r1v(end,end);
    
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
    b1 = [b1; endtime];
    
    A = [A zeros(size(A,1), size(A1,2));
        zeros(size(A1,1), size(A,2)) A1];
    b = [b; b1];
    
    scene.agents = [scene.agents agent]; 
    UserTols = [UserTols agent.radius];
    
end
scene.coeff_matrix = coefficients_matrix;
%Aeq = [Aeq zeros(size(Aeq,1),numel(scene.agents))];
%A = [A zeros(size(A,1),numel(scene.agents))];
num_agents = numel(scene.agents);

PV = v;
PE = e;
hold on;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1.5, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
hold on;
axis equal;

drawnow;

%minimize here
options = optimoptions('fmincon', ...
                        'SpecifyObjectiveGradient', true, ...
                        'Display', 'iter',...
                        'UseParallel', true,...
                        'HessianFcn',@hessfcn);
                    
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = num_inside_iters;
q_i  = [reshape(v', numel(v),1)];
qn = reshape(v', numel(v),1);
print_agents(fname+"initial.json", scene, qn);
for  iter=1:max_iters
    [q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x,UserTols, numel(scene.agents),scene, e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            [], options);
    qn = q_i;

    %print_agents(fname+"agents.json", scene, qn) %TODO

    Q = reshape(qn, numel(qn)/numel(scene.agents), numel(scene.agents));
    scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
    PV = reshape(qn, 3, numel(qn)/3)';
    [CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness', 1.5, 'PolySize', 4);
    surf_anim.Vertices = CV;
    drawnow;
    mu_barrier = mu_barrier*0.5;
end

%path_energy(q_i,UserTols, numel(scene.agents),scene, e, surf_anim)
