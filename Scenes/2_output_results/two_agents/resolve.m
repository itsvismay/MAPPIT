%resolve_collision
fname = "agents.json";
val = jsondecode(fileread(fname));

% PLOT AGENTS 
figure('Name', 'Simple collision'); hold on;
for i = 1:size(val.agents,1)
    a = val.agents(i);
    PV = a.v;
    PE = a.e;
    [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.001, 'PolySize', 10);
    tsurf(CF, CV); 
    axis equal;
end
hold off;