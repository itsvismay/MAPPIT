function [s] = plottings(scene, q, g)
    global surf_anim
%     QQ = reshape(q, 3, numel(q)/3)';
%     GG = reshape(g, 3, (numel(g))/3)';
%     Xquiver = QQ(:,1);
%     Yquiver = QQ(:,2);
%     Zquiver = QQ(:,3);
%     Uquiver = GG(:,1);
%     Vquiver = GG(:,2);
%     Wquiver = GG(:,3);
%     quiver3(Xquiver, Yquiver, Zquiver, Uquiver, Vquiver, Wquiver, 5);
%     drawnow;
    
    Q = reshape(q, numel(q)/numel(scene.agents), numel(scene.agents));
    
    CV = [];
    for i=1:numel(scene.agents)
        scene.agents(i).v = reshape(Q(:,i), 3, size(Q,1)/3)';
        PV = scene.agents(i).v;
        PE = scene.agents(i).e;
        rad = scene.agents(i).radius;
        [CVi,CF,CJ,CI] = edge_cylinders(PV,PE, 'Thickness', 2.0*rad, 'PolySize', 2);
        CV = [CV; CVi];
    end
    surf_anim.Vertices = CV;
    drawnow;
end