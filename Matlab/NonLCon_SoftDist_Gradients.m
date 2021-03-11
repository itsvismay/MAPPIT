%read in from the scene file
addpath("../external/smooth-distances/build/");
%fname = "../Scenes/output_results/eight_agents/agent_circle/";
% fname = "../Scenes/output_results/three_agents/friends_asymmetric_grouping/";
% fname = "../Scenes/output_results/three_agents/mass_asymmetric_collisions/";
% fname = "../Scenes/output_results/three_agents/no_collisions/";
% 
% fname = "../Scenes/output_results/three_agents/size_asymmetric_collisions/";
% fname = "../Scenes/output_results/three_agents/size_mass_asymmetric_collisions/";
% fname = "../Scenes/output_results/three_agents/symmetric_collisions/";
%fname = "../Scenes/roomba_maze/scene_3/";
fname = "../Scenes/output_results/scaling_tests/10_agents/";

setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BF = boundary_faces(tF);
scene.terrain.BVind = unique(scene.terrain.BF);
scene.terrain.BV = tV(scene.terrain.BVind,:);
scene.agents = [];

% surf_anim = tsurf(scene.terrain.F, scene.terrain.V); 
% hold on;
% axis equal;

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];
UserTols = [];
coefficients_matrix = zeros(6, numel(fieldnames(setup_params.agents)));
AdjM = adjacency_matrix(scene.terrain.F);
AdjM_visited = AdjM;
a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 6;
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
        
    
    [r1e, r1v, AdjM, AdjM_visited] = set_path(AdjM, AdjM_visited, agent, scene);
    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    
    e = [e; r1e];
    
    %wiggles the rod start so that they aren't intersecting
    endtime = r1v(end,3);
   
    r1v(:,3) = r1v(:,3)/(i);% sort(rand(1,size(r1v,1))*(endtime));%%
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];
    agent.rest_edge_lengths = r1el;
    agent.rest_region_lengths = [0; r1el(1:size(r1el)-1) + r1el(2:size(r1el))];
    
    %sets up the agent bvh
    [B,I] = build_distance_bvh(agent.v, []);
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
scene.coeff_matrix = coefficients_matrix;
Aeq = [Aeq zeros(size(Aeq,1),numel(scene.agents))];
A = [A zeros(size(A,1),numel(scene.agents))];


PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
hold on;
axis equal;
drawnow;

%minimize here
options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient',true, 'Display', 'iter', 'UseParallel', false);
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = 1e3;
otions.StepTolerance = 1e-4;
q_i  = [reshape(v', numel(v),1); -1e-8*ones(numel(scene.agents),1)];
[q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x,UserTols, numel(scene.agents),scene, e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            @(x)nonlinear_constraints(x, UserTols, scene), options);
qn = q_i(1:end-numel(scene.agents));

print_agents(fname+"agents.json", scene, qn)

Q = reshape(qn, numel(qn)/numel(scene.agents), numel(scene.agents));
 scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
PV = reshape(qn, 3, numel(qn)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',0.75, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;

%path_energy(q_i,UserTols, numel(scene.agents),scene, e, surf_anim)
