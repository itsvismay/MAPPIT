%read in from the scene file
addpath("../external/smooth-distances/build/");
addpath("../CrowdSolverCpp/matlab/");

global scene num_agents simple_sd mu_barrier  surf_anim;

mu_barrier= 1;
mu_barrier_decrease_factor = 0.5;

%looks at the output folder and gets the latest test run number

% 
% fname = "../Scenes/1_input_scenes/bottleneck/dense_bottleneck/";nLayer = 5; 
%     num_segments = 200; max_iters = 2; num_inside_iters = 100; mu_barrier= 1; smoothing_eps_coeff = 1e-2;

%run = 6;
% fname = "../Scenes/1_input_scenes/circle_maze/seven_agents/";nLayer = 10; 
%     num_segments = 300; max_iters = 5; num_inside_iters = 30;mu_barrier= 1; smoothing_eps_coeff = 1e-2;

runId = 0;
fname = "../Scenes/1_input_scenes/battlefield/denser/";nLayer = 3; 
    num_segments = 100; max_iters = 5; num_inside_iters = 30;mu_barrier= 1; smoothing_eps_coeff = 1e-2;

printedScene = jsondecode(fileread(strrep(fname,"1_input_scenes","2_output_results")+"run"+runId+"/initial.json"));

setup_params = jsondecode(fileread(fname+"setup.json"));
simple_sd = 1;
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);

scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BF = boundary_faces(tF);
scene.terrain.BVind = unique(scene.terrain.BF);
scene.terrain.BV = tV(scene.terrain.BVind,:);
scene.terrain.bvhBV = 0;%makePointBVH(scene.terrain.BV(:, 1:2), [1:size(scene.terrain.BV,1)]);
scene.agents = [];
scene.smoothDistAlpha = setup_params.smoothDistAlpha;
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
scene.coeff_matrix = printedScene.energyMultiplierMatrix;

a = fieldnames(setup_params.agents);
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
    agent.preferred_end_time = agent.xse(end, end);
    if (isfield(getfield(setup_params.agents, a{i}), "preferred_end_time"))
        agent.preferred_end_time = getfield(setup_params.agents, a{i}).preferred_end_time;
    end
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
    coefficients_matrix(7, i) = getfield(setup_params.agents, a{i}).energy_coefficients.K_rg;

    %set path new
    r1e = printedScene.agents(i).e;
    r1v = printedScene.agents(i).v;
    
    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    e = [e; r1e]; 
    
    
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

%Aeq = [Aeq zeros(size(Aeq,1),numel(scene.agents))];
%A = [A zeros(size(A,1),numel(scene.agents))];
num_agents = numel(scene.agents);

scene.timings.preprocess = toc(pre_tic);

PV = v;
PE = e;
hold on;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.25, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
hold on;
axis equal;
axis tight manual;
drawnow;



%minimize here
options = optimoptions('fmincon', ...
                        'OutputFcn', @stoppingCriteria,...
                        'SpecifyObjectiveGradient', true, ...
                        'Display', 'iter',...
                        'UseParallel', true,...
                        'HessianFcn',@hessfcn);
                    
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = num_inside_iters;
q_i  = [reshape(v', numel(v),1)];
qn = reshape(v', numel(v),1);
scene.timings.iterations = [];
plottings(scene, qn, []);

for  iter=1:max_iters
    iterTic = tic;
    iterOutput = struct;
    iterOutput.egAgent = 0;
    iterOutput.egMap = 0;
    iterOutput.eAcc = 0;
    iterOutput.eKE = 0;
    iterOutput.eReg = 0;
    iterOutput.ePv = 0;
    iterOutput.gAcc = 0;
    iterOutput.gKE = 0;
    iterOutput.gReg = 0;
    iterOutput.gPv = 0;
    iterOutput.egTotal = 0;
    iterOutput.hAgent = 0;
    iterOutput.hMap = 0;
    iterOutput.hAcc = 0;
    iterOutput.hKE = 0;
    iterOutput.hReg = 0;
    iterOutput.hPv = 0;
    iterOutput.hTotal = 0;
    iterOutput.output = struct;
    iterOutput.totalTime =0;
    scene.timings.iterations = [scene.timings.iterations iterOutput];
    [q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x,UserTols, numel(scene.agents), e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            [], options);
    
    qn = q_i;
    output.message = strrep(output.message, newline, "");
    scene.timings.iterations(end).output = output;
    scene.timings.iterations(end).totalTime = toc(iterTic);
    print_agents(strrep(fname,"1_input_scenes","2_output_results") +"run"+runId+"/", "agents", qn) %TODO
    plottings(scene, qn, []);
   
    mu_barrier = mu_barrier*mu_barrier_decrease_factor;
end

%path_energy(q_i,UserTols, numel(scene.agents),scene, e, surf_anim)
