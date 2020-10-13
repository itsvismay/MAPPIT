%% THE SETUP
xse = [0 0 0;
       10 0 5;
       0 7 8;
       10 10 10];

max_time = 10;
waypoints = size(xse,1)-1;
seg_per_waypoint = 200;
segments = seg_per_waypoint*waypoints; %total_segments = segments per waypoint * num waypoints

%vertices
r1v = [linspace(xse(1,1),xse(end,1), segments+1)', linspace(xse(1,2),xse(end,2),segments+1)', linspace(xse(1,3),xse(end,3),segments+1)'];
%edges
r1e = [(1:segments)' (2:(segments+1))'];
    
%rest edge lenths
r1el = sqrt(sum((r1v(r1e(:,2),:) - r1v(r1e(:,1))).^2,2));

%compute kb (disable, not in use)
%r1kb = ComputeKB(r1v, r1e, r1el);

%rest region lenghts (disable, don't know what this is doing)
%r1l = [0; r1el(1:size(r1el)-1) + r1el(2:size(r1el))];

%make indexing less irritating in the optimization functions
flat_index = @(I, j) 3*I + (j-3); 

%fix end points
num_constraints = size(xse,1)*2 +1;
Aeq_rows = [1:num_constraints];
Aeq_cols = [1:3 reshape([3*seg_per_waypoint*linspace(1,waypoints,waypoints) + 1; 3*seg_per_waypoint*linspace(1,waypoints,waypoints) + 2], 1,[])];
Aeq_vals = ones(1,num_constraints);
Aeq = sparse(Aeq_rows, Aeq_cols, Aeq_vals, num_constraints, numel(r1v));
beq = [xse(1,:)'; reshape(xse(2:end,1:2)', 1, [])'];

%time is monotonic constraints 
A = sparse(reshape([1:segments; 1:segments], 1, 2*segments), reshape([3:3:(numel(r1v)-3); 6:3:(numel(r1v))], 1, 2*segments), reshape([ones(1, segments); -ones(1,segments)], 1, 2*segments), segments, numel(r1v));
b = zeros(segments, 1);

%max time via inequality constraint (comment out to turn off, but make sure
%you uncomment the end time potential energy in cost)
A = [A; zeros(1, 3*segments+2) 1];
b = [b; max_time];

%minimize here
options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 1e6;
q = fmincon(@(x) path_energy(x,max_time), reshape(r1v', numel(r1v),1), A,b,Aeq,beq,[],[],[], options);

r1v = reshape(q, 3, size(r1v,1))';
PV = r1v;
PE = r1e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.5, 'PolySize', 10);
% %PLOT AGENTS
figure('Name', 'Simple collision');
tsurf(CF, CV);
axis equal


function c = path_energy(q, max_time) 
    
    %various fun path energies, I'll use principle of least action becase
    %I like it
    %kinetic energy of curve integrated along piecewise linear segment is 
    % 0.5*dx'*dx./dt 
    dx = reshape(q(4:end) - q(1:end -3), 3, numel(q)/3-1)';
    T = sum(0.5*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
    
    %arrive at a particular time
    V = 0.5*(q(end) - max_time).^2;
    
    V2 = sum(sqrt(sum(dx.*dx,2))); %sum of squared distances
    
    %find minimum energy curve
    c = T+V+V2;
    
end

