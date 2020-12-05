%read in from the scene file
fname = "scene_1/";
setup_params = jsondecode(fileread(fname+"setup.json"));
scene = struct;

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
    agent.seg_per_waypoint = 10;
    agent.segments = agent.seg_per_waypoint*agent.waypoints;
    agent.v = 0;
    
    
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
    r1v(:,3) = r1v(:,3)/(2.0*i);
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
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.5, 'PolySize', 10);
surf_anim = tsurf(CF, CV); 
axis equal;
drawnow;

%minimize here
options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 1e5;
q = fmincon(@(x) path_energy(x, numel(scene.agents), e, surf_anim), reshape(v', numel(v),1), A,b,Aeq,beq,[],[], @(x)constraints(x, numel(scene.agents), A), options);

PV = reshape(q, 3, numel(q)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.5, 'PolySize', 10);
surf_anim.Vertices = CV;
drawnow;


% %PRINT AGENTS 
print_agents(fname+"agents.json",scene, q);

function p = print_agents(fname, scene, q)
    Q = reshape(q, numel(q)/numel(scene.agents), numel(scene.agents));

    for i=1:size(Q,2)
        scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
    end
    jsontext = jsonencode(scene);
    jsontext = strrep(jsontext, ',', sprintf(',\n'));
    jsontext = strrep(jsontext, '[{', sprintf('[\r{\r'));
    jsontext = strrep(jsontext, '}]', sprintf('\r}\r]'));
    fileID = fopen(fname, 'w');
    fprintf(fileID, jsontext);
end

function c = path_energy(q, num_agents, e, surf_anim) 
    
    Q = reshape(q, numel(q)/num_agents, num_agents);
    T = 0;
    V = 0;
    V2 = 0;
    %for 
    for i=1:num_agents
        q_i = Q(:, i);
        
        %various fun path energies, I'll use principle of least action becase I like it
        %kinetic energy of curve integrated along piecewise linear segment is 
        % 0.5*dx'*dx./dt 

        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        T = T + sum(0.5*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        
        %arrive at a particular time
        V = 0.5*(q_i(end) - 10).^2;
        V2 = sum(sqrt(sum(dx.*dx,2))); %sum of squared distances
    end

    %find minimum energy curve
    c = T;%+V2;%+V+V2;
    
    if T<0
        T
        PV = reshape(q, 3, numel(q)/3)';
        PE = e;
        [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',0.5, 'PolySize', 10);
        surf_anim.Vertices = CV;
        drawnow;
    end
    T
    

end

function [c, ceq] = constraints(q, num_agents, A)
    c = 0;
    Q = reshape(q, numel(q)/num_agents, num_agents);
    
    for i=1:num_agents
        for j=i+1:num_agents
            A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
            B1 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
            %if time is monotically increasing
            %sometimes, its not, for whatever reason
            if(all(diff(A1(:,3))>0) && all(diff(B1(:,3))>0))
                A1 = SamplePointsForRod(A1);
                B1 = SamplePointsForRod(B1);
            else
                c = c + 100;
            end
            D = (A1(:,1) - B1(:,1)').^2 + (A1(:,2) - B1(:,2)').^2 + (A1(:,3) - B1(:,3)').^2;
            D = sqrt(D);
            % min dist >= tol ----> 0 >= tol- mindist
            if  min(D(:)) < 1
                c = c+(1 - min(D(:)));
            end
        end
    end
    
    ceq = [];
    
end

function p = SamplePointsForRod(x)
    ts = linspace(x(1,3)+0.001, x(end,3)-0.001, 100); %time slots to sample
    
    %get lower index for time and upper index
    
    [trash, trash, lind] = histcounts(ts, x(:,3)');
    uind = lind +1;
    s = (ts'- x(lind,3))./(x(uind,3) - x(lind,3));
    qp0 = x(lind, :);
    qp1 = x(uind, :);

    p = qp0 + s.*(qp1 - qp0);
end

