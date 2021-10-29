fname = "../../../1_input_scenes/antelopes/200/";
inputVerts = readOBJ(fname + "input_verts.obj");
outputVerts = readOBJ(fname + "output_verts.obj");

xy_starts = [inputVerts(:,1), inputVerts(:,3), zeros(size(inputVerts,1),1)];
xy_ends = [outputVerts(:,1), outputVerts(:,3), 100*ones(size(inputVerts,1),1)];
scene = struct;
scene.name = "antelopes";
agents = struct;
for i=1:size(xy_starts,1)
    a = struct;
    a.id = i;
    a.xse = [xy_starts(i,:)' xy_ends(i,:)']';
    a.mesh = "deer";
    a.animation_cycles = "pass";
    a.friends = [];
    a.radius = 1;
    a.mass = 1;
    ec = struct;
    ec.K_agent = 1; % agent
    ec.K_tol = 0; % friends
    ec.K_accel = 0;
    ec.K_map = 1;
    ec.K_ke = 1;
    ec.K_pv = 0;
    ec.K_rg = 10;
    a.energy_coefficients = ec;
    agents = setfield(agents, "agent"+string(i), a);
end
t = struct;
t.mesh = "terrain.obj";
scene.agents = agents;
scene.terrain = t;
scene.obstacles= struct;
scene.smoothDistAlpha = 30;
scene.useSpaceTimeGraphDiags = 0;

jsontext = jsonencode(scene);
jsontext = strrep(jsontext, ',', sprintf(',\n'));
jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
fileID = fopen(fname+"setup.json", 'w');
fprintf(fileID, jsontext);