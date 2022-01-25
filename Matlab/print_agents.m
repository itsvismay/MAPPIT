function p = print_agents(path, fname, q)
    global scene
    %Function needs to print the:
    % 1. Scene and agents (without the environment, or BVH)
    % 2. The figure (matlab figure of trajectory on map)
    % 3. Timings
    
    if ~exist(path, 'dir')
        mkdir(path)
    end
    Q = reshape(q, numel(q)/numel(scene.agents), numel(scene.agents));
    
    %%Print JSON
    for i=1:size(Q,2)
        scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
    end
    
    for i=1:numel(scene.timings.iterations)
        scene.timings.iterations(i).output.message = strrep(scene.timings.iterations(i).output.message, " ", "");
        scene.timings.iterations(i).output.message = strrep(scene.timings.iterations(i).output.message, ",", "");
    end
    printableScene = struct;
    printableScene.agents = scene.agents;
    printableScene.energyMultiplierMatrix = scene.coeff_matrix;
    printableScene.timings = scene.timings;
    
    printableScene.agents = rmfield(printableScene.agents,"bvh");
    
    jsontext = jsonencode(printableScene);
    jsontext = strrep(jsontext, ',', sprintf(',\n'));
    jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
    jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
    fileID = fopen(path+fname+".json", 'w');
    fprintf(fileID, jsontext);
    
    %% Figure
    colors = get_colors();
    fig = figure;
    % hold the map
    tsurf(scene.terrain.F, scene.terrain.V); 
    hold on;
    axis equal;
    for i=1:numel(scene.agents)
        PV = scene.agents(i).v;
        PE = scene.agents(i).e;
        rad = scene.agents(i).radius;
        [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness', 2.0*rad, 'PolySize', 2);
        tsurf(CF, CV,'FaceColor', colors(:, i)'); hold on;
    end
    saveas(fig, path+"fig-"+fname);
    close(fig);
end
