function p = print_agents(fname, scene, q)
    Q = reshape(q, numel(q)/numel(scene.agents), numel(scene.agents));

    for i=1:size(Q,2)
        scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
    end
    jsontext = jsonencode(scene);
    jsontext = strrep(jsontext, ',', sprintf(',\n'));
    jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
    jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
    fileID = fopen(fname, 'w');
    fprintf(fileID, jsontext);
end