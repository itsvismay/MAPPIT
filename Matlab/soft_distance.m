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
     if(isinf(D) || isnan(D))
         D = 1e1;
         return;
     end
     %loop over all mesh vertices
     for i=1:size(V,1)
        G(i, 1) = -(1./diff_d).*sum(diff_all_d(:,i).*(dx(:,i)./d(:,i)));
        G(i, 2) = -(1./diff_d).*sum(diff_all_d(:,i).*(dy(:,i)./d(:,i)));
        G(i, 3) = -(1./diff_d).*sum(diff_all_d(:,i).*(dz(:,i)./d(:,i)));
        if(sum(isnan(G(i,:)))>0)
            G(i,:) = 1e-3*rand(3,1);
        end
     end
     
                
end

