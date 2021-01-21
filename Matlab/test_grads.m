% V = [ -1.0278       5.5347            0;
%      -0.82776       4.0588          0.5;
%      -0.77776       2.5829            1;
%      -0.77776       1.1069          1.5;
%      -0.77776     -0.36898            2;
%      -0.77776      -1.8449          2.5;
%      -0.97776      -1.8449            3;
%      -0.87776      -2.9518          3.5;
%      -0.87776      -3.6898            4;
%       -1.0278      -4.0588          4.5;
%       -1.0278      -5.5347           10]; %l shape
% E = [1:size(V,1)-1; 2:size(V,1)]';
% 
% 
% V1 = V(1:(end-2),:) - V(2:(end-1),:);
% V2 = V(3:end,:) - V(2:(end-1),:);
% %V2 = V(E(:,2),:);
% V12 = cross(V1,V2);
% V12norm = sqrt(V12(:,1).^2+V12(:,2).^2+V12(:,3).^2);
% V1norm = sqrt(V1(:,1).^2+V1(:,2).^2+V1(:,3).^2);
% V2norm = sqrt(V2(:,1).^2+V2(:,2).^2+V2(:,3).^2);
% 
% if(V12norm < eps)
%     V12norm = eps;
% end
% 
% zhat = V12./V12norm;
% 
% 
% theta = 2*atan2(dot(V12, zhat, 2), V1norm.*V2norm + dot(V1,V2,2));
% theta0 = pi;
% %bending energy 
% k = 100;
% e = 0.5*100*(theta-theta0).^2;
% 
% %bending energy gradient
% grad = [k*theta.*(cross(V1,zhat)./(V1norm.*V1norm)); ...
%         -k*(theta-0).*((cross(V2,zhat)./(V2norm.*V2norm)) - (cross(V1,zhat)./(V1norm.*V1norm))); ...
%        k*theta.*(cross(V2,zhat)./(V2norm.*V2norm))]; 
%     
% 
% hold on;
% plot(V(:,1), V(:,2), 'r');
% plot(V(:,1), V(:,2), 'r*');
% line([V(:,1)'; V(:,1)' + 0.01*grad(:,1)'], [V(:,2)'; V(:,2)' + 0.01*grad(:,2)'], 'Color', [0 0 1]);
% hold off;
% axis equal;

fname = "../tunnel_maze/scene_1/";
setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BVind = unique(boundary_faces(tF));
scene.terrain.BV = tV(scene.terrain.BVind,:);
scene.agents = [];
surf_anim = tsurf(scene.terrain.F, scene.terrain.V); 
hold on;

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];
UserTols = [];

a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 10;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    
    [r1e, r1v] = set_path(agent, scene);
    %edges
    agent.e = r1e;
    r1e = r1e + size(v,1);
    
    e = [e; r1e];
    
    %wiggles the rod start so that they aren't intersecting
    endtime = r1v(end,3);
    r1v(:,3) = r1v(:,3)/2;%sort(rand(1,size(r1v,1))*(endtime));
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];
    agent.rest_edge_lengths = r1el;
    agent.rest_region_lengths = [0; r1el(1:size(r1el)-1) + r1el(2:size(r1el))];
    
    
    %fix end points
    num_constraints = size(agent.xse,1)*2 +1;
    Aeq_rows = [1:num_constraints];
    Aeq_cols = [1:3 reshape([3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 1; 
                3*agent.seg_per_waypoint*linspace(1,agent.waypoints,agent.waypoints) + 2], 1,[])];
    Aeq_vals = ones(1,num_constraints);
    
    A1eq = sparse(Aeq_rows, Aeq_cols, Aeq_vals, num_constraints, numel(r1v));
    b1eq = [agent.xse(1,:)'; reshape(agent.xse(2:end,1:2)', 1, [])'];
    
    Aeq = [Aeq zeros(size(Aeq,1), size(A1eq,2)); 
            zeros(size(A1eq,1), size(Aeq,2)) A1eq];
        
    beq = [beq; b1eq];

    %time is monotonic constraints 
    %t_i+1 - t_i > 0 --> in this format --> A1*q <= 0
    A1 = sparse(reshape([1:agent.segments; 1:agent.segments], 1, 2*agent.segments),...
                reshape([3:3:(numel(r1v)-3); 6:3:(numel(r1v))], 1, 2*agent.segments),...
                reshape([ones(1, agent.segments); -ones(1,agent.segments)], 1, 2*agent.segments),... 
                agent.segments, numel(r1v));
    b1 = zeros(agent.segments, 1);

    %max time via inequality constraint (comment out to turn off, but make sure
    %you uncomment the end time potential energy in cost)
    A1 = [A1; zeros(1, 3*agent.segments+2) 1];
    b1 = [b1; agent.max_time];
    
    A = [A zeros(size(A,1), size(A1,2));
        zeros(size(A1,1), size(A,2)) A1];
    b = [b; b1];
    
    scene.agents = [scene.agents agent]; 
    UserTols = [UserTols agent.radius];
end
Aeq = [Aeq zeros(size(Aeq,1),numel(scene.agents))];
A = [A zeros(size(A,1),numel(scene.agents))];
