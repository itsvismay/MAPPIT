%read in from the scene file
fname = "../roomba_maze/scene_2/";
setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.agents = [];

v = [];
e = [];
el = [];
Aeq = [];
beq = [];
A = [];
b = [];

a = fieldnames(setup_params.agents);
for i = 1:numel(a)
    agent.id = i;
    agent.xse = getfield(setup_params.agents, a{i}).xse;
    agent.max_time = agent.xse(end, end);
    agent.waypoints = size(agent.xse,1)-1;
    agent.seg_per_waypoint = 15;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    agent.radius = getfield(setup_params.agents, a{i}).radius;
    
    %edges
    r1e = [(1:agent.segments)' (2:(agent.segments+1))'];
    agent.e = r1e;
    r1e = r1e + size(v,1);
    e = [e; r1e];
    
    % %vertices
    r1v = [linspace(agent.xse(1,1),agent.xse(end,1), agent.segments+1)', ... %x
            linspace(agent.xse(1,2),agent.xse(end,2),agent.segments+1)', ... %y
            linspace(agent.xse(1,3),agent.xse(end,3),agent.segments+1)'];    %t
    
    %wiggles the rod start so that they aren't intersecting
    endtime = r1v(end,3);
    r1v(:,3) = sort(rand(1,size(r1v,1))*(endtime));
    r1v(end,3) = endtime;
    agent.v = r1v;            
    v = [v;r1v];
    
    
    %rest edge lenths
    r1el = sqrt(sum((v(r1e(:,2),:) - v(r1e(:,1))).^2,2));
    el = [el; r1el];
    
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
end

PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
axis equal;
drawnow;

UserTols = 2*ones(size(scene.agents(1).v,1),size(scene.agents(2).v,1));
Aeq = [Aeq zeros(size(Aeq,1),numel(UserTols))];
A = [A zeros(size(A,1),numel(UserTols))];
%minimize here

options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 1e3;
q = [reshape(v', numel(v), 1); 0*reshape(UserTols',numel(UserTols),1)];
mu = 1;
for k = 1:10
    [q, fval, exitflag, output] = fmincon(@(x) path_energy(x, numel(scene.agents), e, surf_anim, numel(v), UserTols, mu), q, A,b,Aeq,beq,[],[], [],options);
    mu = mu/2.0;
    PV = reshape(q(1:numel(v)), 3, numel(v)/3)';
    [CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',2, 'PolySize', 4);
    surf_anim.Vertices = CV;
    drawnow;
end



function [f] = path_energy(q, num_agents, e, surf_anim, num_v, UserTols, mu)
    Q = reshape(q(1:num_v), num_v/num_agents, num_agents); %3*nodes x agents
    Tols = reshape(q(num_v+1:end), size(UserTols))';
    
    T = 0;
    for i=1:num_agents
        q_i = Q(:, i); %3*nodes
        
        m = 1;

        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        T = T + sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
            
    end
    
    B = 0;
    U = 0;
    for i=1:num_agents
        for j=i+1:num_agents
            A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            d = (A1(:,1) - A2(:,1)').^2 + (A1(:,2) - A2(:,2)').^2 + (A1(:,3) - A2(:,3)').^2;
            d = sqrt(d);
            B = sum(sum(-mu*log(-Tols + d)));
            
            U = sum(sum((Tols -UserTols).^2));
        end
    end
    
    
    %find minimum energy curve
    T
    B
    U
    f = T+B+U;
    

end

