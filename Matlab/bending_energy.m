function [be] = bending_energy(agent, verts)
    edges = agent.e;
    
    %compute kb
    E = (verts(edges(:,2),:) - verts(edges(:,1)));
    e_i_1 = E(1:size(edges,1)-1,:);
    e_i = E(2:size(edges,1),:);
    l_i_1 = agent.rest_edge_lengths(1:size(edges,1)-1);
    l_i = agent.rest_edge_lengths(2:size(edges,1));
    kb = [[0 0 0]; cross(e_i_1, e_i)./(l_i_1.*l_i + dot(e_i, e_i_1,2))];
    
    %bending modulus 
    alpha = 10000;
    be = 0;
    for i = 2:size(kb)
        be = be + (alpha*(kb(i)*kb(i))/agent.rest_region_lengths(i));
    end
    
end