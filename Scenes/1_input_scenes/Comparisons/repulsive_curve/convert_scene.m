folder_name = "repulsive_curve/repulsive_curve/";

[V1,F1] =readOBJ(folder_name+'agent1.obj');
[V2, F2] = readOBJ(folder_name+'agent2.obj');


maxTime = max([A(1), B(1), C(1), D(1), E(1), F(1), G(1), H(1),]);
A = [A; repmat(A(end-2:end), [maxTime-A(1),1])];
B = [B; repmat(B(end-2:end), [maxTime-B(1),1])];
C = [C; repmat(C(end-2:end), [maxTime-C(1),1])];
D = [D; repmat(D(end-2:end), [maxTime-D(1),1])];
E = [E; repmat(E(end-2:end), [maxTime-E(1),1])];
F = [F; repmat(F(end-2:end), [maxTime-F(1),1])];
G = [G; repmat(G(end-2:end), [maxTime-G(1),1])];
H = [H; repmat(H(end-2:end), [maxTime-H(1),1])];
Q = [A(2:end) B(2:end) C(2:end) D(2:end) E(2:end) F(2:end) G(2:end) H(2:end)];

scene = struct;
num_agents = 3;
scene.agents = [];
fname = "../../2_output_results/Comparisons/"+folder_name + "agents";
%Function needs to print the:
% 1. Scene and agents (without the environment, or BVH)
% 2. The figure (matlab figure of trajectory on map)
% 3. Timings

   

%%Print JSON
for i=1:size(Q,2)
    a = struct;
    a.v = reshape(Q(:,i), 3, size(Q,1)/3)';
    a.radius = 0.5;
    a.id = i;
    a.mesh = "agent";
    segments = size(a.v)-1;
    a.e = [(1:segments)' (2:(segments+1))'];
    scene.agents = [scene.agents a];
end

jsontext = jsonencode(scene);
jsontext = strrep(jsontext, ',', sprintf(',\n'));
jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
fileID = fopen(fname+".json", 'w');
fprintf(fileID, jsontext);
    
%% Figure
colors = get_colors();
fig = figure;
% hold the map
for i=1:numel(scene.agents)
    PV = scene.agents(i).v;
    PE = scene.agents(i).e;
    rad = scene.agents(i).radius;
    [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness', 2.0*rad, 'PolySize', 10);
    tsurf(CF, CV,'FaceColor', colors(:, i)'); hold on;
end


function colors = get_colors()
colors = [[0.19154,      0.78665,      0.56299]', ...
          [1.0,      0.31997,      0.20648]',...
          [0.24929,       0.4175,      0.87439]',...
          [0.18907,      0.61792,      0.73418]',...
          [0.19126,      0.50687,      0.84305]',...
          [0.20127,      0.55259,      0.78732]',...
          [0.35756,      0.34505,      0.83857]',...
          [0.22557,      0.28247,       1.0]',...
          [0.13424,       1.0,      0.31652]',...
          [0.66463,      0.19804,      0.67851]',...
          [0.19239,      0.44728,       0.9015]',...
          [0.95037,      0.20389,      0.38692]',...
          [0.77345,      0.61533,       0.1524]',...
          [0.30151,      0.34759,      0.89208]',...
          [0.99527,      0.17835,      0.36756]'];
end