%single kinetic energy minimizing rod
global nverts;
global nsegments;

%tend
tn = 19.8;
%start point and end point
p0 = [0 0 0];
pn = [1 1 tn]; %x, y, sum(1/dt)*segments

%curve discretization is staggered
%store 2d positions at nodes and reciprocal time on each section
%kinetic energy is linear in reciprocal delta time, makes the sensitivity of the
% gradient better. t

%make weird sigmoid curve for testing

%flat line, vertical line, flat line 
nsegments = 99;
nverts = nsegments + 1;

%useful matrices
%ONE = ones(


%inverse time array
it = [repmat(1e3, 1, nsegments/3) repmat(3.0, 1, nsegments/3) repmat(1e3, 1, nsegments/3)]';

%position away
x = [linspace(p0(1), pn(1)/2, 33)' linspace(p0(2), pn(2)/2, 33)'; linspace(pn(1)/2, pn(1)/2, 33)' linspace(pn(2)/2, pn(2)/2, 33)';  linspace(pn(1)/2, pn(1), 34)' linspace(pn(2)/2, pn(2), 34)' ]

fig = figure;

hold on;
h = plot_curve(x, it);
hold off;

options = optimoptions(@fmincon,'Algorithm','interior-point',...
                       'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false,...
                       'PlotFcn',@optplotcurve,...
                       'MaxIterations', 1e9,...
                       'MaxFunctionEvaluations', 1e9, ...
                       'Display', 'iter', ...
                       'HessianApproximation', 'finite-difference',...
                       'SubproblemAlgorithm', 'cg');


%bounds on times
lb = -inf(2*nverts+nsegments,1);
lb((2*nverts+1):end) = 0;

% inequality constraint
A = zeros(1, 2*nverts+nsegments);
A((2*nverts+1):end) = 1;

b = sum(it);
%end point equality constraints
Aeq = zeros(5, 2*nverts + nsegments); 
Aeq(1,1) = 1;
Aeq(2,nverts+1) = 1;
Aeq(3,nverts) = 1;
Aeq(4,2*nverts) = 1;
Aeq(5, (2*nverts+1):end) = 1;
beq = [p0(1:2) pn(1:2) pn(3)]';
q0 = [x(:); it];
q_out = fmincon(@energy, q0, A, b, Aeq, beq, lb, [], [], options); 

x = [q_out(1:nverts) q_out((nverts+1):2*nverts)];
it = q_out((2*nverts+1):end);
h = plot_curve(x, it);
    
%kinetic energy integrated over time (assume constant kinetic energy over
%each segment)
%q stacked vector of xy positions then per segment times
function [T, grad] = energy(q)

    global nverts;
    global nsegments;

    r = 10;
    u = [q(1:nverts) q((nverts+1):2*nverts)];
    it = q((2*nverts+1):end);
    
    du = u(2:end,:) - u(1:(end-1),:);
    %T = 0.5*sum(sum((du.*du).*(it*r)));

    
    %l = du + 
    %need some inverse length weighting to stabilize
    T = 0.5*sum(sum((du.*du)*r)) + 0.5*sum(it);
    
    e = ones(nsegments,1);
    D = spdiags([-e e],0:1,nsegments,nverts);
    if nargout > 1
        grad = r*(D'*D)*u;
        grad = [grad(:);sign(it)]; %real gradient
    end

end

%plot curve (don't forget to unwind and integrate the time 
function h = plot_curve(pos, itime) 

 time = cumsum([0; 1./itime]);
 
 h = plot3(pos(:,1), pos(:,2), time);
 
 hold on; 
 plot3(pos(:,1), pos(:,2), time, 'r*');
 hold off;
 view(0, 30);
 xlabel('x');
 ylabel('y');
 zlabel('t');
 
end

function stop = optplotcurve(q, optimValues, state)

    global nverts;
    global nsegments;
    
    u = [q(1:nverts) q((nverts+1):2*nverts)];
    it = q((2*nverts+1):end);
    plot_curve(u, it);
    
    stop = false;
end


