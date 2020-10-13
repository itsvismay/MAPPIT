x_field = 10;
y_field = 10;
max_time = 10;
segments = 100;

%vertices
r1v = [linspace(0,x_field,segments+1)', linspace(0,y_field,segments+1)', linspace(0,max_time,segments+1)'];
%edges
r1e = [linspace(1,segments,segments)', linspace(2,segments+1,segments)'];
%rest edge lenths
r1el = sqrt(sum((r1v(r1e(:,2),:) - r1v(r1e(:,1))).^2,2));
%compute kb
r1kb = ComputeKB(r1v, r1e, r1el);
%rest region lenghts
r1l = [0; r1el(1:size(r1el)-1) + r1el(2:size(r1el))];

r2v = [linspace(x_field,0,segments+1)', linspace(0,y_field,segments+1)', linspace(0,max_time,segments+1)'];
r2e = [ linspace(1,segments,segments)', linspace(2,segments+1,segments)'];
r2el = sqrt(sum((r2v(r2e(:,2),:) - r2v(r2e(:,1))).^2,2));
r2kb = ComputeKB(r2v, r2e, r2el);
r2l = [0; r2el(1:size(r2el)-1) + r2el(2:size(r2el))];



PV = [r1v;r2v];
PE = [r1e;r2e + [size(r1v, 1),  size(r1v, 1)]];
[CV,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness',1, 'PolySize', 10);

%SOLVE RODS FOR COLLISIONS

%%%%%

%PLOT AGENTS
figure('Name', 'Simple collision');
for ti=1:segments
    clf;
    tsurf(CF, CV);
    hold on;
    
    %agent 1
    ax = r1v(ti,1);
    ay = r1v(ti,2);
    plot3(ax,ay,0,'go', 'LineWidth',3,'MarkerSize', 15);
    hold on;
    
    %agent 2
    ax = r2v(ti,1);
    ay = r2v(ti,2);
    plot3(ax,ay,0,'bo', 'LineWidth',3,'MarkerSize', 15);
    drawnow;
    
    %plotting, imaging
    [imind, cm] = rgb2ind(frame2im(getframe(gcf)), 256);
    if ti == 1 
        imwrite(imind, cm,'./simple_collision.gif','gif', 'DelayTime',0.1, 'Loopcount',inf); 
    else
        imwrite(imind, cm,'./simple_collision.gif','gif','DelayTime',0.1, 'WriteMode','append');
    end
    
end

function [be] = BendEnergy(restRegionLengths, kb)
    %bending modulus 
    alpha = 1;
    be = 0;
    for i = 2:size(kb)
        be = be + (alpha*(kb(i)*kb(i))/restRegionLengths(i))
    end
    
end

function [se] = StretchEnergy(verts, edges, restEL)
    ym = 1e4;
    rod_rad = 0.5;
    se = 0;
    
    for i = 1:size(edges,1)
        l = sqrt(sum((verts(edges(:,2),:) - verts(edges(:,1))).^2,2));
        eps = (l - restEL(i))/restEL(i);
        se = se + 0.5*ym*rod_rad*rod_rad*pi*eps*eps;
    end
end

function [e] = TotalEnergy(verts, edges, restEL, restRL, kb)
    e = 0;
    e = e + BendEnergy(restRL, kb);
    e = e + StretchEnergy(verts, edges, restEL);
end

function [kb] = ComputeKB(verts, edges, restEdgeLengths)
    E = (verts(edges(:,2),:) - verts(edges(:,1)));
    e_i_1 = E(1:size(edges,1)-1,:);
    e_i = E(2:size(edges,1),:);
    l_i_1 = restEdgeLengths(1:size(edges,1)-1);
    l_i = restEdgeLengths(2:size(edges,1));
    kb = [[0 0 0]; cross(e_i_1, e_i)./(l_i_1.*l_i + dot(e_i, e_i_1,2))];
end

