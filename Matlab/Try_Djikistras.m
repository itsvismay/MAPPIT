[I1, minD, VI] = snap_points([0 -8.75 0; 8.75 0 0], scene.terrain.V);
[I2, minD, VI] = snap_points([-8.75 0 0; 0 8.75 0], scene.terrain.V);

AdjM = adjacency_matrix(scene.terrain.F);
[D1, P1, CDATA1] = mydijk(scene.terrain.V, AdjM,AdjM, I1(1), I1(2), scene.terrain.BVind);
D1
[D2, P2, CDATA2] = mydijk(scene.terrain.V, AdjM, AdjM, I2(1), I2(2), scene.terrain.BVind);
D2

tV = scene.terrain.V;
tV(scene.terrain.BVind,3) = 1;
surf_anim = tsurf(scene.terrain.F, tV); 
axis equal;
drawnow;