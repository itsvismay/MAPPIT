fname = "../Scenes/1_input_scenes/bottleneck/30Agents/";
% inputVerts = readOBJ(fname + "start_verts.obj");
% outputVerts = readOBJ(fname + "end_verts.obj");

% xy_starts = [inputVerts(:,1), inputVerts(:,3), outputVerts(:,2)];
% xy_ends = [outputVerts(:,1), outputVerts(:,3), abs(outputVerts(:,2))]; %100*ones(size(inputVerts,1),1)
% 
% %sort so that its in order of smallest to largest by axis
% %lets me control which agents start/end and order
% [B, I] = sort(xy_starts(:,1), 'ascend');
% xy_starts = xy_starts(I,:);
% xy_starts(:,3) =  zeros(size(inputVerts,1),1);
% 
% [B, I] = sort(xy_ends(:,1), 'ascend');
% xy_ends = xy_ends(I,:);
% %xy_ends(:,3) = 30*ones(size(inputVerts,1),1);



scene = struct;
scene.name = "bottleneck";
agents = struct;
for i=1:15
    a = struct;
    a.id = i;
    a.xse = [[6*(i+1), (2*mod(i,2)-1)*-37.5, 0]' [6*(i+1), (2*mod(i,2)-1)*37.5, 100]']';
    a.mesh = "elephant";
    a.preferred_end_time = 60;
    a.animation_cycles = "pass";
    a.friends = [];
    a.collision_interactions = [];
%     if(i==1)
%         a.collision_interactions(end+1) = [2];
%     elseif(i==102)
%         a.collision_interactions(end+1) = [101];
%     else
%         a.collision_interactions(end+1) = i-1;
%         a.collision_interactions(end+1) = i+1;
%     end
    
    a.radius = 0.35;
    a.mass = 100;
    ec = struct;
    ec.K_agent = 10; % agent
    ec.K_tol = 0; % friends
    ec.K_accel = 0;
    ec.K_map = 1;
    ec.K_ke = 1;
    ec.K_pv = 0;
    ec.K_rg = 100;
    a.energy_coefficients = ec;
    agents = setfield(agents, "agent"+string(i), a);
end
t = struct;
t.mesh = "terrain.obj";
scene.agents = agents;
scene.terrain = t;
scene.obstacles= struct;
scene.smoothDistAlpha = 40;
scene.useSpaceTimeGraphDiags = 0;

jsontext = jsonencode(scene);
jsontext = strrep(jsontext, ',', sprintf(',\n'));
jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
fileID = fopen(fname+"setup.json", 'w');
fprintf(fileID, jsontext);