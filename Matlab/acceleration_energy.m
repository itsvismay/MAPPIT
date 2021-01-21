function [a, g] = acceleration_energy(q_i)

    %vertices Vx3
    %V = reshape(q_i, 3, numel(q_i)/3)'
    dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
    V1 = dx(1:end-1,:);
    V2= dx(2:end,:);
    V1xV2 = cross(V1, V2, 2);
    V1dV2 = dot(V1, V2, 2);
    V1xV2norm = sqrt(V1xV2(:,1).^2+V1xV2(:,2).^2+V1xV2(:,3).^2);
    V1xV2norm(V1xV2norm < eps)= eps; %making sure norms are bigger than epsilon
    Z = V1xV2 ./ V1xV2norm; 
    V1norm = sqrt(sum(V1.^2,2));
    V2norm = sqrt(sum(V2.^2,2));
    Y = dot(V1xV2, Z, 2);
    X = V1norm.*V2norm + V1dV2;
    angle = 2*atan2(Y,X);
    %bending energy
    k =1e-2;
    a = 0.5*k*sum((angle - 0).^2);
    
    %vectorized gradient computation
    %end points have zero gradient (no bending energy applied)
    gr = [0 0 0; ...
        k.*angle.*((cross(V2,Z)./(V2norm.*V2norm)) + (cross(V1,Z)./(V1norm.*V1norm))); ...
        0 0 0]; 
    g = reshape(gr', numel(gr),1);
end
