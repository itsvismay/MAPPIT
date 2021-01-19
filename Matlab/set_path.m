function [re, rv] = set_path(agent, scene)
    s = [agent.xse(1,1),agent.xse(1,2),0];
    t = [agent.xse(end,1),agent.xse(end,2),0];
    [I1, minD, VI] = snap_points([s; t], scene.terrain.V);
    agent.xse(:, 1:2) = VI(:, 1:2);
    
    AdjM = adjacency_matrix(scene.terrain.F);
    
    [D1, P1] = mydijk(scene.terrain.V, AdjM, I1(1), I1(2), scene.terrain.BV, scene.terrain.BVind);
    
    v = scene.terrain.V(P1,:);
    v(:,3) = linspace(0,10, size(v,1));
    
    re = [(1:agent.segments)' (2:(agent.segments+1))'];
    rv = interp1(v(:,3), v(:,1:2), linspace(agent.xse(1,3),agent.xse(end,3),agent.segments+1));
    rv = [rv linspace(agent.xse(1,3),agent.xse(end,3),agent.segments+1)'];
    
end