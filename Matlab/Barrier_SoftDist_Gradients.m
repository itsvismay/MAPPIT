%read in from the scene file
fname = "../tunnel_maze/scene_1/";
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

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];
UserTols = [];

a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 10;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    
    [r1e, r1v] = set_path(agent, scene);
    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    
    e = [e; r1e];
    
    %wiggles the rod start so that they aren't intersecting
    endtime = r1v(end,3);
    r1v(:,3) = r1v(:,3)/2;%sort(rand(1,size(r1v,1))*(endtime));
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];
    agent.rest_edge_lengths = r1el;
    agent.rest_region_lengths = [0; r1el(1:size(r1el)-1) + r1el(2:size(r1el))];
    
    
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


% PV = v;
% PE = e;
% [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 4);
% surf_anim = tsurf(CF, CV); 
% hold on;
% axis equal;
% drawnow;

A1 = scene.agents(1).v;
[A1, J1] = sample_points_for_rod(A1, scene.agents(1).e);
plot3(A1(:,1), A1(:,2), A1(:,3), '-ok')
hold on;
A1(:, 3) = zeros(size(A1,1),1);
[P, JP] = sample_points_for_rod(scene.terrain.V, scene.terrain.BF);
plot3(P(:,1), P(:,2), P(:,3), '*k')

% [D,GP] = soft_distance(50,P, A1);
% JG1 = J1'*GP;
% gkk = (-1/(-scene.agents(1).radius + D))*reshape(JG1', size(JG1,1)*size(JG1,2), 1);
% QQ = A1;
% GG = reshape(gkk, 3, (numel(gkk))/3)';
% Xquiver = scene.agents(1).v(:,1);
% Yquiver = scene.agents(1).v(:,2);
% Zquiver = scene.agents(1).v(:,3);
% Uquiver = GG(:,1);
% Vquiver = GG(:,2);
% Wquiver = GG(:,3);
% quiver3(Xquiver, Yquiver, Zquiver, Uquiver, Vquiver, Wquiver);
% plot3(scene.agents(1).v(:,1), scene.agents(1).v(:,2), scene.agents(1).v(:,3), '-ok', 'LineWidth',5);

% min d from every point on surface to rod
minX = min(scene.terrain.V(:,1));
maxX = max(scene.terrain.V(:,1));
minY = min(scene.terrain.V(:,2));
maxY = max(scene.terrain.V(:,2));
[X,Y]  = meshgrid(minX:0.1:maxX, minY:0.1:maxY);
C = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        [D,GP] = soft_distance(50,P, [X(i,j), Y(i,j),0]);
        C(i,j) = D;
    end
end
surf(X,Y,zeros(size(X))-0.5,C);

%minimize here
options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'iter', 'UseParallel', true, 'HessianApproximation', 'lbfgs');
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = 1e4;
q_i  = [reshape(v', numel(v),1); -1e-8*ones(numel(scene.agents),1)];
[q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x,UserTols, numel(scene.agents),scene, e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            [], options);
[q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x,UserTols, numel(scene.agents),scene, e, surf_anim),... 
    q_i, ...
    A,b,Aeq,beq,[],[], ...
    [], options);
qn = q_i(1:end-numel(scene.agents));

print_agents(fname+"agents.json", scene, qn)

Q = reshape(qn, numel(qn)/numel(scene.agents), numel(scene.agents));
 scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
PV = reshape(qn, 3, numel(qn)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',1, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;

%path_energy(q_i,UserTols, numel(scene.agents),scene, e, surf_anim)
