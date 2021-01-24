function [a, g] = acceleration_energy(q_i)

%     %vertices Vx3
%     %V = reshape(q_i, 3, numel(q_i)/3)'
%     dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
%     V1 = dx(1:end-1,:);
%     V2= dx(2:end,:);
%     V1xV2 = cross(V1, V2, 2);
%     V1dV2 = dot(V1, V2, 2);
%     V1xV2norm = sqrt(V1xV2(:,1).^2+V1xV2(:,2).^2+V1xV2(:,3).^2);
%     V1xV2norm(V1xV2norm < eps)= eps; %making sure norms are bigger than epsilon
%     Z = V1xV2 ./ V1xV2norm; 
%     V1norm = sqrt(sum(V1.^2,2));
%     V2norm = sqrt(sum(V2.^2,2));
%     Y = dot(V1xV2, Z, 2);
%     X = V1norm.*V2norm + V1dV2;
%     angle = 2*atan2(Y,X);
%     %bending energy
%     k =1e-2;
%     a = 0.5*k*sum((angle - 0).^2);
%     
%     %vectorized gradient computation
%     %end points have zero gradient (no bending energy applied)
%     gr = [0 0 0; ...
%         k.*angle.*((cross(V2,Z)./(V2norm.*V2norm)) + (cross(V1,Z)./(V1norm.*V1norm))); ...
%         0 0 0]; 
    g = reshape(gr', numel(gr),1);
end

function [e, g] = preferred_time_energy(Q, scene, K)
    GT = zeros(size(Q));
    e=0;
    pv= 10;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        v = dx(:,1:2)./dx(:,3);
        % 0.5* (s  - m/pv)^2
        
        e = e + K*0.5*(sum(sqrt(dx(:,1).^2 + dx(:,2).^2))./dx(end,3) - pv).^2;
        
        jac_left = zeros(numel(q_i)/3, 3);
        jac_right = zeros(numel(q_i)/3, 3);

%         jac_left(1:end-1,1) = ; %
%         jac_left(1:end-1,2) = ;
%         jac_left(1:end-1,3) = ;
% 
%         jac_right(2:end, 1) = ;
%         jac_right(2:end, 2) = ;
%         jac_right(2:end, 3) = ;
        
        dEdq = (jac_left + jac_right);
        
        for j=1:size(dx,1)
            r_inds = [3*(j-1) + 1 : 3*(j-1)+6];
            q_i1 = q_i(r_inds(1));
            q_i2 = q_i(r_inds(2));
            q_i3 = q_i(r_inds(3));
            q_i4 = q_i(r_inds(4));
            q_i5 = q_i(r_inds(5));
            q_i6 = q_i(r_inds(6));
            
            dEdq(r_inds) = dEdq(r_inds) + ...
            [ ((2*q_i1 - 2*q_i4)*(pv + ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/(q_i3 - q_i6)))/(2*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)*(q_i3 - q_i6)), ((2*q_i2 - 2*q_i5)*(pv + ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/(q_i3 - q_i6)))/(2*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)*(q_i3 - q_i6)), -(((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)*(pv + ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/(q_i3 - q_i6)))/(q_i3 - q_i6)^2, -((2*q_i1 - 2*q_i4)*(pv + ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/(q_i3 - q_i6)))/(2*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)*(q_i3 - q_i6)), -((2*q_i2 - 2*q_i5)*(pv + ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/(q_i3 - q_i6)))/(2*((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)*(q_i3 - q_i6)), (((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)*(pv + ((q_i1 - q_i4)^2 + (q_i2 - q_i5)^2)^(1/2)/(q_i3 - q_i6)))/(q_i3 - q_i6)^2];

        end
        
        GT(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
    end
    g = K*reshape(GT, size(GT,1)*size(GT,2),1);
end

function [e, g] = preferred_velocity_energy(Q, scene, K)
    GT = zeros(size(Q));
    e=0;
    pv= 10;
    for i=1:numel(scene.agents)
        q_i = Q(:, i); %3*nodes
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        v = dx(:,1:2)./dx(:,3);
        e = e + K*0.5*sum((sum(v.^2,2) - pv^2).^2);
        
        jac_left = zeros(numel(q_i)/3, 3);
        jac_right = zeros(numel(q_i)/3, 3);

        jac_left(1:end-1,1) = -2*dx(:,1)./(dx(:,3).^2); %
        jac_left(1:end-1,2) = -2*dx(:,2)./(dx(:,3).^2);
        jac_left(1:end-1,3) = (2*(dx(:,1).^2 + dx(:,2).^2))./(dx(:,3).^3);
        jac_left(1:end-1,:) = jac_left(1:end-1,:).*(sum(v.^2,2) - pv^2);

        jac_right(2:end, 1) = 2*dx(:,1)./(dx(:,3).^2);
        jac_right(2:end, 2) = 2*dx(:,2)./(dx(:,3).^2);
        jac_right(2:end, 3) = -(2*(dx(:,1).^2 + dx(:,2).^2))./(dx(:,3).^3);
        jac_right(2:end,:) = jac_right(2:end,:).*(sum(v.^2,2) - pv^2);
        
        dEdq = (jac_left + jac_right);
        
        GT(:,i) = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        
    end
    g = K*reshape(GT, size(GT,1)*size(GT,2),1);
end

