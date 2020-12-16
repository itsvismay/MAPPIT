%read in from the scene file
fname = "../roomba_maze/scene_2/";
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
[q_i, fval, exitflag, output] = fmincon(@(x) path_energy(x, numel(scene.agents), e, surf_anim),... 
                            q_i, ...
                            A,b,Aeq,beq,[],[], ...
                            @(x)nonlinear_constraints(x, scene), options);
q = q_i(1:end-numel(scene.agents));
                        
PV = reshape(q, 3, numel(q)/3)';
[CV,CF,CJ,CI] = edge_cylinders(PV,e, 'Thickness',2, 'PolySize', 4);
surf_anim.Vertices = CV;
drawnow;

function [c, ceq] = nonlinear_constraints(q_i, scene)
    c = 0;
    w = 0;
    num_agents = numel(scene.agents);
    q = q_i(1:end-numel(scene.agents));
   
    Q = reshape(q, numel(q)/num_agents, num_agents);  
    GQ = zeros(size(Q));
    Tols = q_i(end-num_agents+1:end);
    
    %intersection constraints
    for i=1:num_agents
        for j=i+1:num_agents
            A1 = reshape(Q(:,i), 3, numel(Q(:,i))/3)';
            A2 = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
                
            [A1, J1] = sample_points_for_rod(A1, 50);
            [A2, J2] = sample_points_for_rod(A2, 50);

            %true min D
            d = sqrt((A1(:,1) - A2(:,1)').^2 + (A1(:,2) - A2(:,2)').^2 + (A1(:,3) - A2(:,3)').^2); 
            true_min_d = min(d(:));

%             [D, G2]= soft_distance( 200, A1, A2);
%             [D, G1]= soft_distance( 200, A2, A1);

            D = true_min_d;
            D
            
            tol = Tols(i) + Tols(j);
            
            c = max(c, tol-D);

%             JG1 = -J1'*G1;
%             JG2 = -J2'*G2;
%             
%             GQ(:, i) = GQ(:, i) + reshape(JG1', numel(JG1),1);
%             GQ(:, j) = GQ(:, j) + reshape(JG2', numel(JG2),1);
        end
        
        P = reshape(Q(:,j), 3, numel(Q(:,j))/3)';
        
        [P, J1] = sample_points_for_rod(P, 50);
        %add radius (-x,  to make signed_
        P(:, 3) = zeros(size(P,1),1);

%         [D,G] = soft_distance(10,scene.terrain.BV, P);
        A1 = scene.terrain.BV;
        A2 = P;
        d = sqrt((A1(:,1) - A2(:,1)').^2 + (A1(:,2) - A2(:,2)').^2 + (A1(:,3) - A2(:,3)').^2); 
        true_min_d = min(d(:));
%         JG = -J1'*G;
        D = true_min_d;
        tol = Tols(j);
        if(tol-D>0)
            w = max(w, tol - D);
%             GQ(:, j) = GQ(:, j) + reshape(JG, numel(JG),1);
        end
    end
%     GC = [reshape(GQ, size(GQ,1)*size(GQ,2),1); 0];
    c
    w
    c= max(c,w);
    c
    ceq = [];
    GCeq = [];

    
end

function [P,J] = sample_points_for_rod(X, samples, varargin)
%     ts = linspace(X(1,3)+0.001, X(end,3)-0.001, samples);
%     PV = [interp1(X(:,3), X(:, 1:2), ts, 'linear') ts'];
%     PE = [(1:size(PV,1)-1)' (2:size(PV,1))'];
%     if numel(varargin)>0 && strcmp('cylinders',varargin{1})
%         rad = varargin{2}; 
%         
%         [P,J,~,~] = edge_cylinders(PV,PE, 'Thickness',2*rad, 'PolySize', 10);
%     else
%         % samples along centerline
%         P = PV;
%         J = PE;
%     end

    s = 10; %5 samples per edge
    %for each edge
    P = zeros(s*(size(X,1)-1)+1, 3);
    J = zeros(size(P,1), size(X,1));
    for i = 1:size(X,1)-1
        tl = X(i,3);
        tr = X(i+1,3);
        xl = X(i,1);
        xr = X(i+1,1);
        yl = X(i,2);
        yr = X(i+1, 2);
        for j =1:s
            P(s*(i-1)+j,3) = tl + (j-1)*(tr - tl)/s;
            P(s*(i-1)+j,1) = xl + (j-1)*(xr - xl)/s;
            P(s*(i-1)+j,2) = yl + (j-1)*(yr - yl)/s;
            J(s*(i-1)+j,i) = 1 -(j-1)/s;
            J(s*(i-1)+j,i+1) = (j-1)/s;
        end
    end
    P(end,:) = X(end,:);
    J(end,end-1) = 0;
    J(end,end) = 1;
    
  
end

% soft distance function
function [D, G] = soft_distance(alpha, X, V)

     %plot soft distance to get a sense of what is going on
     dx = (X(:,1) - V(:,1)');
     dy = (X(:,2) - V(:,2)');
     dz = (X(:,3) - V(:,3)');
    
     d = sqrt(dx.^2 + dy.^2 + dz.^2); 
     diff_all_d = exp(-alpha*d);
     diff_d = sum(sum(diff_all_d));
     
     D = -1./alpha.*log(diff_d);
     
     G = zeros(size(V));
     %loop over all mesh vertices
     for i=1:size(V,1)
         
        G(i, 1) = -(1./diff_d).*sum(diff_all_d(:,i).*(dx(:,i)./d(:,i)));
        G(i, 2) = -(1./diff_d).*sum(diff_all_d(:,i).*(dy(:,i)./d(:,i)));
        G(i, 3) = -(1./diff_d).*sum(diff_all_d(:,i).*(dz(:,i)./d(:,i)));
     end
     
                
end

function [f,g] = path_energy(q_i, num_agents, e, surf_anim)
    q = q_i(1:end-num_agents);
    Q = reshape(q, numel(q)/num_agents, num_agents); %3*nodes x agents
    G = zeros(size(Q));
    Tols = q_i(end-num_agents+1:end);

    T = 0;
    V = 0;
    V2 = 0;
    
    g = zeros(size(q));

    for i=1:num_agents
        q_i = Q(:, i); %3*nodes
        
        %various fun path energies, I'll use principle of least action becase I like it
        %kinetic energy of curve integrated along piecewise linear segment is 
        % 0.5*dx'*dx./dt 
        m = 1;

        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        T = T + sum(0.5*m*sum(dx(:, 1:2).*dx(:,1:2),2)./dx(:,3)); %kinetic energy
        
%         dEdq_left = zeros(numel(q_i)/3, 3);
%         dEdq_right = zeros(numel(q_i)/3, 3);
%         
%         dEdq_left(1:end-1,1) = -m*dx(:,1)./dx(:,3);
%         dEdq_left(1:end-1,2) = -m*dx(:,2)./dx(:,3);
%         dEdq_left(1:end-1,3) = 0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));
%         
%         dEdq_right(2:end, 1) = m*dx(:,1)./dx(:,3);
%         dEdq_right(2:end, 2) = m*dx(:,2)./dx(:,3);
%         dEdq_right(2:end, 3) = -0.5*m*(dx(:,1).*dx(:,1) + dx(:,2).*dx(:,2))./(dx(:,3).*dx(:,3));
        
%         dEdq = dEdq_left + dEdq_right;
        
        W = 0.5*(Tols - [2;2])'*(Tols - [2;2]);
        W = 100*W;
%         G(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
    end
%     g = [reshape(G, size(G,1)*size(G,2),1); q_i(end) - 2];
    
    %find minimum energy curve
    f = T+W;
    
    PV = reshape(q, 3, numel(q)/3)';
    PE = e;
    [CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',2, 'PolySize', 4);
    surf_anim.Vertices = CV;
    drawnow;
    

end



