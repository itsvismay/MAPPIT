%%read in from the scene file
addpath("../external/smooth-distances/build/");
%fname = "../Scenes/output_results/eight_agents/agent_circle/";
%fname = "../Scenes/output_results/three_agents/test/";
fname = "../Scenes/output_results/scaling_tests/1_agents/";
%fname = "../Scenes/output_results/scaling_tests/test/";
setup_params = jsondecode(fileread(fname+"setup.json"));

%% initialize global variables into `scene`

%% initialize the agent paths


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