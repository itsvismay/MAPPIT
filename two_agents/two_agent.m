%% THE SETUP
agent1.id = 1;
agent1.xse =[10 0 0;
             0 10 10];
agent1.max_time = 10;
agent1.waypoints = size(agent1.xse,1)-1;
agent1.seg_per_waypoint = 10;
agent1.segments = agent1.seg_per_waypoint*agent1.waypoints; %total_segments = segments per waypoint * num waypoints
agent1.v = 0;
agent1.e = 0;

agent2.id = 2;
agent2.xse =[0 0 0;
            10 10 10];
agent2.max_time = 10;
agent2.waypoints = size(agent2.xse,1)-1;
agent2.seg_per_waypoint = 10;
agent2.segments = agent2.seg_per_waypoint*agent2.waypoints;
agent2.v = 0;
agent2.e = 0;

scene.agents = [agent1, agent2];

i=1;
v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];
while i <= size(scene.agents,2)
    a = scene.agents(i);
    
    %edges
    r1e = [(1:a.segments)' (2:(a.segments+1))'];
    scene.agents(i).e = r1e;
    r1e = r1e + size(v,1);
    e = [e; r1e];
    
    % %vertices
    r1v = [linspace(a.xse(1,1),a.xse(end,1), a.segments+1)', ...
            linspace(a.xse(1,2),a.xse(end,2),a.segments+1)', ...
            linspace(a.xse(1,3),a.xse(end,3),a.segments+1)'];
    scene.agents(i).v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];

    %fix end points
    num_constraints = size(a.xse,1)*2 +1;
    Aeq_rows = [1:num_constraints];
    Aeq_cols = [1:3 reshape([3*a.seg_per_waypoint*linspace(1,a.waypoints,a.waypoints) + 1; 
                3*a.seg_per_waypoint*linspace(1,a.waypoints,a.waypoints) + 2], 1,[])];
    Aeq_vals = ones(1,num_constraints);
    
    A1eq = sparse(Aeq_rows, Aeq_cols, Aeq_vals, num_constraints, numel(r1v));
    b1eq = [a.xse(1,:)'; reshape(a.xse(2:end,1:2)', 1, [])'];
    
    Aeq = [Aeq zeros(size(Aeq,1), size(A1eq,2)); 
            zeros(size(A1eq,1), size(Aeq,2)) A1eq];
        
    beq = [beq; b1eq];

    %time is monotonic constraints 
    A1 = sparse(reshape([1:a.segments; 1:a.segments], 1, 2*a.segments), reshape([3:3:(numel(r1v)-3); 6:3:(numel(r1v))], 1, 2*a.segments), reshape([ones(1, a.segments); -ones(1,a.segments)], 1, 2*a.segments), a.segments, numel(r1v));
    b1 = zeros(a.segments, 1);

    %max time via inequality constraint (comment out to turn off, but make sure
    %you uncomment the end time potential energy in cost)
    A1 = [A1; zeros(1, 3*a.segments+2) 1];
    b1 = [b1; a.max_time];
    
    A = [A zeros(size(A,1), size(A1,2));
        zeros(size(A1,1), size(A,2)) A1];
    b = [b; b1];
    
    i= i + 1;
end

%minimize here
options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 1e6;
q = fmincon(@(x) path_energy(x), reshape(v', numel(v),1), A,b,Aeq,beq,[],[], @(x)constraints(x), options);


v = reshape(q, 3, size(v,1))';
PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.1, 'PolySize', 10);
tsurf(CF, CV); 
axis equal;

scene.agents(1).v = v( 1:scene.agents(1).segments+1, :);
scene.agents(2).v = v( scene.agents(1).segments+2 : 2+scene.agents(1).segments+scene.agents(2).segments, :);
    
% %PRINT AGENTS 
print_agents("agents.json",scene);

function p = print_agents(fname, scene)
    jsontext = jsonencode(scene);
    jsontext = strrep(jsontext, ',', sprintf(',\n'));
    jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
    jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
    fileID = fopen(fname, 'w');
    fprintf(fileID, jsontext);
end

function c = path_energy(q) 

    %various fun path energies, I'll use principle of least action becase I like it
    %kinetic energy of curve integrated along piecewise linear segment is 
    % 0.5*dx'*dx./dt 
    dx = reshape(q(4:end) - q(1:end -3), 3, numel(q)/3-1)';
    T = sum(0.5*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy

    %arrive at a particular time
    %V = 0.5*(q(end) - max_time).^2;

    V2 = sum(sqrt(sum(dx.*dx,2))); %sum of squared distances

    %find minimum energy curve
    c = T;%+V+V2;

end

function [c, ceq] = constraints(q)
    c = 0;
    
    A = q(1:numel(q)/2);
    B = q(numel(q)/2+1:numel(q));
    
    A1 = reshape(A, 3, numel(A)/3)';
    B1 = reshape(B, 3, numel(B)/3)';
    
    D = (A1(:,1) - B1(:,1)').^2 + (A1(:,2) - B1(:,2)').^2 + (A1(:,3) - B1(:,3)').^2;
    D = sqrt(D) + eye(size(D));
    % min dist >= tol ----> 0 >= tol- mindist
    c = 2.95-min(D(:));
    
    ceq = [];
    
end

function p = WhereAmIAt(t, x)
    [d, ind] = min(abs(x(:,3) - t));

    %if t rounds up, then use [ind-1 , ind]
    %if t rounds down, use [ind, ind+1]
    if ind>1 &&(t<=round(t)|| ind==size(x,1))
        lind = ind-1;
        uind = ind;
    else
        lind = ind;
        uind = ind+1;
    end    
    
    s = (t - x(lind,3))/(x(uind,3) - x(lind,3));
    qp0 = x(lind, 1:2);
    qp1 = x(uind, 1:2);
    
    p = qp0 + s*(qp1 - qp0);
    
end

