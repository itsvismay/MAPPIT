% models motion using spring forces%
clear;
[x,y] = meshgrid(0:1:10);
z = zeros(size(x,1));
goals = [10,9,0;
        10,1,0];
agents = [0,6,0;
          0,7,0];
goals_q = reshape(goals',[],1);
agents_q = reshape(agents',[],1);
agents_v = agents_q*0;
agents_k = kron(eye(size(agents,1)), eye(3));
agents_m = kron(eye(size(agents,1)), eye(3));
%F = ma
%F = kx
%a = kx/m
dt = 0.1;
damping = 10;
options=optimoptions('fmincon','Algorithm','interior-point');

for t=1:150
    clf;
    surf(x,y,z, x.*y*0);
    hold on;
    F = agents_k*(goals_q - agents_q) - damping*agents_v;
    q0 = agents_q;
    v0 = agents_v;
    M = agents_m;
    fun = @(q)objective(q,q0,v0,F,M, dt);
    %A = [1,1,1,-1,-1,-1;
         %-1,-1,-1,1,1,1];
    A = kron([1,-1;-1,1], eye(3));
    b = [1,1];
    
    nonlncon = @(q)constraints(q, A, b, agents);
    q = fmincon(fun, q0, [],[],[],[],[],[], nonlncon,options);
    agents_q = q;
    for i=1:size(agents,1)
        ax = agents_q(3*(i-1)+1);
        ay = agents_q(3*(i-1)+2);
        if i==1
            plot3(ax,ay,0,'bo', 'LineWidth',3,'MarkerSize', 15);
        else
            plot3(ax,ay,0,'go', 'LineWidth',3,'MarkerSize', 15);
        end
        plot3(goals(i,1), goals(i,2),0,'ro', 'LineWidth',3,'MarkerSize', 15)
    end
    drawnow;
    %plotting, imaging
    [imind, cm] = rgb2ind(frame2im(getframe(gcf)), 256);
    if t == 1 
        imwrite(imind, cm,'./crappycrowds.gif','gif', 'DelayTime',0.1, 'Loopcount',inf); 
    else
        imwrite(imind, cm,'./crappycrowds.gif','gif','DelayTime',0.1, 'WriteMode','append');
    end
end

function [e] = objective(q,q0, v0, F, M,dt)
% T - V
    e = norm((1/(dt*dt))*(M*(q - q0 - v0)) - F)^2;

end
function [c,ceq] = constraints(q, A, b,agents)
    %|Aq|
    v = A*q;
    c = 0;
    for i=1:size(agents,1)
       d = norm(v(3*(i-1)+1 : 3*(i-1)+3));
       if(d)<1
           c = c+ 1/d;
       end
    end
    
    ceq = [];
end