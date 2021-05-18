%read in from the scene file
addpath("../external/smooth-distances/build/");
%fname = "../Scenes/output_results/eight_agents/agent_circle/";
%fname = "../Scenes/output_results/three_agents/test/";
fname = "../Scenes/output_results/scaling_tests/test/";
%fname = "../Scenes/output_results/scaling_tests/10_agents/";

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
Kw = 100;
Kt = 1;
Ka = 1;

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

nLayer = 20;
nTotalVer = nLayer * length(scene.terrain.V(:,1));
% build up 3d graph
% find all the edges [[i,j],...] in 2d graph
A_lt = tril(AdjM);
[edge_s,edge_t] = find(A_lt);
edge = zeros(length(edge_s),2);
edge(:,1) = edge_s;
edge(:,2) = edge_t;

    
% set the num of layer and get the spacetime graph 
% (vertices: VV and edges: EE)
time = linspace(0, 1, nLayer)';
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

a = fieldnames(setup_params.agents);
PATHS = [];
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.mass = getfield(setup_params.agents, a{i}).mass;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 50;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    % set_path
    [r1e, r1v, newA, newA_visited] = set_path3d(newA, newA_visited, agent, scene, VV, EE, nLayer, nTotalVer);
    
    % edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    e = [e; r1e];
    
    % vertices
    endtime = r1v(end,3);
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
end

PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.75, 'PolySize', 10);
figure;
surf_anim = tsurf(CF, CV); 
axis equal;
drawnow;
hold on;
[ii,jj,ss] = find(newA_visited);
EEE = []
for k=1:length(ii)
    EEE = [EEE; ii(k), jj(k)];
end
tsurf(EEE, VV);

