function [re, rv, AdjM, A_visited] = set_path3d(AdjM, A_visited, agent, scene, VV, EE, BV, BVind, nLayer, nTotalVer, segments)
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
    
    [D1, P1, AdjM, A_visited] = mydijk3d(VV, EE, AdjM, A_visited, I1(1), I1(2), BV, BVind, agent.radius);
    
    % P1 are the indexes from djikstra's path
    % read out the vertex values into vv
    % P1 are the indexes from djikstra's path
    % read out the vertex values into vv
    if(segments < length(P1))
        segments
        length(P1)
        sprintf("Error: Set_path3d: add more segments");
        return
    end
    
    vv = VV(P1,:);
    vv(:,3) = vv(:,3);
    
    re = [(1:segments)' (2:(segments+1))'];
    if(size(P1)==1)
        PV = scene.agents(end).v;
        PE = scene.agents(end).e;
        rad = scene.agents(end).radius;
        [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness', 2.0*rad, 'PolySize', 4);
        tsurf(CF, CV, 'FaceColor', [0.19154,      0.78665,      0.56299]'); hold on;
        sprintf("Error: Something bad happened. Djikstras wasn't unwound properly")
        %if this happens:
        % check agent start and end points. Make sure they're valid
    end
    rv = interpolate_rod_segments(vv, segments);
    %re = [(1:size(vv,1)-1)' (2:size(vv,1))'];
    %rv = vv;
end

