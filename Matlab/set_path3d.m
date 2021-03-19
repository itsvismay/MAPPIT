function [re, rv, AdjM, A_visited] = set_path3d(AdjM, A_visited, agent, scene, VV, EE, nLayer, nTotalVer)
%SET_PATH3D Summary of this function goes here
%   Detailed explanation goes here
% set_path
    s = [agent.xse(1,1),agent.xse(1,2),0];
    t = [agent.xse(end,1),agent.xse(end,2),0];
    [I1, minD, VI] = snap_points([s; t], scene.terrain.V);
    agent.xse(:, 1:2) = VI(:, 1:2);
    % 3d dijkstra
    % P1 contains all vertices on the path 
    % (note: the index of P1 corresponds to the one in 3d graph VV)
    I1(2) = I1(2) + (nLayer-1)*(nTotalVer/nLayer);
    I1
    [D1, P1, AdjM, A_visited] = mydijk3d(VV, EE, AdjM, A_visited, I1(1), I1(2), scene.terrain.BV, scene.terrain.BVind);
    
    % P1 are the indexes from djikstra's path
    % read out the vertex values into vv
    vv = VV(P1,:);
    vv(:,3) = linspace(0,agent.xse(end,end), size(vv,1));
    re = [(1:agent.segments)' (2:(agent.segments+1))'];
    r1v = interp1(vv(:,3), vv(:,1:2), linspace(agent.xse(1,3),agent.xse(end,3),agent.segments+1));
    rv = [r1v linspace(agent.xse(1,3),agent.xse(end,3),agent.segments+1)'];
    
end

