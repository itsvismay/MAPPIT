%read in from the scene file
fname = "../../three_agents/scene_1/";
setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;
[tV, tF] = readOBJ(fname+setup_params.terrain.mesh);
scene.terrain.V = tV;
scene.terrain.F = tF;
scene.terrain.BV = tV(unique(boundary_faces(tF)),:);
scene.agents = [];

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
    agent.seg_per_waypoint = 6;
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
    r1v(:,3) = r1v(:,3)-0.5*i;%sort(rand(1,size(r1v,1))*(endtime));
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
    UserTols = [UserTols agent.radius];
end
Aeq = [Aeq zeros(size(Aeq,1),numel(scene.agents))];
A = [A zeros(size(A,1),numel(scene.agents))];


PV = v;
PE = e;
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
surf_anim = tsurf(CF, CV); 
axis equal;
drawnow;

%minimize here
options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 1e5;
q_i  = [reshape(v', numel(v),1); zeros(numel(scene.agents),1)];
[q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x,UserTols, numel(scene.agents),scene, e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            @(x)nonlinear_constraints(x, scene), options);
q = q_i(1:end-numel(scene.agents));
                        
PV = reshape(q, 3, numel(q)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',2, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;

function [c, ceq] = nonlinear_constraints(q_i, scene)
    c = zeros(numel(scene.agents),1);
    w = zeros(numel(scene.agents),1);
    r = zeros(numel(scene.agents),1)
    ;
    num_agents = numel(scene.agents);
    q = q_i(1:end-numel(scene.agents));
   
    Q = reshape(q, numel(q)/num_agents, num_agents);  
    GQ = zeros(size(Q));
    Tols = q_i(end-num_agents+1:end);
    
    %intersection constraints
    for i=1:num_agents
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1, J1] = sample_points_for_rod(A1, 50);
        A1(:, 3) = zeros(size(A1,1),1);
        P = scene.terrain.BV;
        [D,G] = soft_distance(200,P, A1);
        
        tol = Tols(i);
        if(1-D>0)
            w(i) = max(w(i), 1 - D);
        end
        c(i) = max(w(i), r(i));
    end
    c
    ceq = [];
   
end


function [f,g] = path_energy(q_i, UserTols, num_agents, scene, e, surf_anim)
    q = q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    G = zeros(size(Q));
    Tols = q_i(end-num_agents+1:end);

    T = 0;
    B = 0;
    W = 0;
    
    g = zeros(size(q));

    for i=1:num_agents
        q_i = Q(:, i); %3*nodes
        
        %various fun path energies, I'll use principle of least action becase I like it
        %kinetic energy of curve integrated along piecewise linear segment is 
        % 0.5*dx'*dx./dt 
        m = 1;

        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        T = T + sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        
        A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
        [A1, J1] = sample_points_for_rod(A1, 50);
        for j =i+1:num_agents
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            [A2, J2] = sample_points_for_rod(A2, 50);
            [D,G] = soft_distance(20,A1, A2);

            tol = Tols(i) + Tols(j);
            
            B = B + -1*log(-tol + D);
        end
%         A1(:, 3) = zeros(size(A1,1),1);
%         P = scene.terrain.BV;
%         [D,G] = soft_distance(200,scene.terrain.BV, A1);
%         tol = Tols(i);
%         B = B + -1*log(-UserTols(i) + D);
     
        
        W = 0.5*(Tols - UserTols')'*(Tols - UserTols');
        W = 10*W;
        
    end
    
    %find minimum energy curve
    f = T+W+B;
    f
    
%     PV = reshape(q, 3, numel(q)/3)';
%     PE = e;
%     [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
%     surf_anim.Vertices = CV;
%     drawnow;
    

end
