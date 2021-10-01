function p = print_agents(path, fname, scene, q)
    %Function needs to print the:
    % 1. Scene and agents (without the environment, or BVH)
    % 2. The figure (matlab figure of trajectory on map)
    % 3. Timings
    
    if ~exist(path, 'dir')
        mkdir(path)
    end
    Q = reshape(q, numel(q)/numel(scene.agents), numel(scene.agents));
    

    for i=1:size(Q,2)
        scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
    end
    
    printableScene = struct;
    printableScene.agents = scene.agents;
    printableScene.energyMultiplierMatrix = scene.coeff_matrix;
    
    printableScene.agents = rmfield(printableScene.agents,"bvh");
    
    jsontext = jsonencode(printableScene);
    jsontext = strrep(jsontext, ',', sprintf(',\n'));
    jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
    jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
    fileID = fopen(path+fname, 'w');
    fprintf(fileID, jsontext);
end