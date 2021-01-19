[I1, minD, VI] = snap_points([0 -8.75 0; 8.75 0 0], scene.terrain.V);
[I2, minD, VI] = snap_points([-8.75 0 0; 0 8.75 0], scene.terrain.V);

AdjM = adjacency_matrix(scene.terrain.F);
[D1, P1, CDATA1] = mydijk(scene.terrain.V, AdjM, I1(1), I1(2), scene.terrain.BVind);
D1
[D2, P2, CDATA2] = mydijk(scene.terrain.V, AdjM, I2(1), I2(2), scene.terrain.BVind);
D2
% NODES = [linspace(1,size(scene.terrain.V,1), size(scene.terrain.V,1))' scene.terrain.V];
% SEGMENTS = [linspace(1,size(scene.terrain.F,1), size(scene.terrain.F,1))' scene.terrain.F];
% SID = I1(1);
% FID = I1(2);
% [D,P] = dijkstra(NODES, SEGMENTS, FID, SID);

% v = scene.terrain.V(P1,:);
% v(:,3) = linspace(0,10, size(v,1));
% rv = interp1(v(:,3), v(:,1:2), linspace(v(1,3),v(end,3),100));
% rv = [rv linspace(v(1,3),v(end,3),100)'];
% e = [(1:size(rv, 1)-1)' (2:size(rv, 1))'];
% PV = rv;
% PE = e;
% [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness', 0.5, 'PolySize', 10);
% surf_anim = tsurf(CF, CV); 
% hold on

% v = scene.terrain.V(P2,:);
% v(:,3) = linspace(0,10, size(v,1));
% rv = interp1(v(:,3), v(:,1:2), linspace(v(1,3),v(end,3),100));
% rv = [rv linspace(v(1,3),v(end,3),100)'];
% e = [(1:size(rv, 1)-1)' (2:size(rv, 1))'];
% PV = rv;
% PE = e;
% [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness', 0.5, 'PolySize', 10);
% surf_anim = tsurf(CF, CV); 
% hold on

tV = scene.terrain.V;
tV(scene.terrain.BVind,3) = 0.5;
surf_anim = tsurf(scene.terrain.F, tV); 
axis equal;
drawnow;