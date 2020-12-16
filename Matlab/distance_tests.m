%Plots soft distance function and optionally checks gradients
%check_gradients_on - boolean specifying whether or not to check gradients,
%
% gradient check is really slow so don't do it by default
% function distance_tests(check_gradients_on)

% 
% %                 J_a(3*i-2) = -(1./diff_d).*diff_all_d(:,i).*(dx(:,i)./d(:,i));
% %                 J_a(3*i-1) = -(1./diff_d).*diff_all_d(:,i).*(dy(:,i)./d(:,i));
% %                 J_a(3*i-0) = -(1./diff_d).*diff_all_d(:,i).*(dz(:,i)./d(:,i));
% 


V1 = [linspace(0,10, 3+1)', ... %x
        linspace(0,10, 3+1)', ... %y
        linspace(0,10, 3+1)'];    %t

V2 = [linspace(10,0, 3+1)', ... %x
    linspace(0,10, 3+1)', ... %y
    linspace(0,10, 3+1)'];    %t

Gfd = soft_distance_fd_gradient(50, V1, V2);
[D, G] = soft_distance_gradient(50, V1, V2);

function [Gfd] = soft_distance_fd_gradient(alpha, X, V)
    fd_tol = 1e-6;
    Gfd = zeros(size(V));
    for i=1:size(V,1)
                
        %central finite differences
        for j = 1:3
            Vp1 = V;
            Vm1 = V;

            Vp1(i,j) = Vp1(i,j) + fd_tol;
            Vm1(i,j) = Vm1(i,j) - fd_tol;

            Gfd(i,j) = (soft_distance(alpha, X, Vp1) - soft_distance(alpha, X, Vm1))./(2*fd_tol);
        end
    end
end

function [D, G] = soft_distance_gradient(alpha, X, V)
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
     
     
     
%      G = zeros(size(X,1), numel(V));
%      %loop over all mesh vertices
%      for i=1:size(V,1)
%         G(:, 3*i-2) = -(1./diff_d).*diff_all_d(:,i).*(dx(:,i)./d(:,i));
%         G(:, 3*i-1) = -(1./diff_d).*diff_all_d(:,i).*(dy(:,i)./d(:,i));
%         G(:, 3*i-0) = -(1./diff_d).*diff_all_d(:,i).*(dz(:,i)./d(:,i));
%      end
   
end

function [D] = soft_distance(alpha, X, V)

     %plot soft distance to get a sense of what is going on
     dx = (X(:,1) - V(:,1)');
     dy = (X(:,2) - V(:,2)');
     dz = (X(:,3) - V(:,3)');
    
     d = sqrt(dx.^2 + dy.^2 + dz.^2); 
     diff_all_d = exp(-alpha*d);
     diff_d = sum(sum(diff_all_d));
     
     D = -1./alpha.*log(diff_d);
                
end
