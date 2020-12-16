% soft distance function
function [D, G] = soft_distance(alpha, X, V)

     %plot soft distance to get a sense of what is going on
     dx = (X(:,1)' - V(:,1)');
     dy = (X(:,2)' - V(:,2)');
     dz = (X(:,3)' - V(:,3)');
    
     d = sqrt(dx.^2 + dy.^2 + dz.^2); 
     diff_all_d = exp(-alpha*d);
     diff_d = sum(diff_all_d,2);
     
     D = -1./alpha.*log(diff_d);
     
     G = zeros(numel(V),1);
     %loop over all mesh vertices
     for i=1:size(V,1)
        G((3*i-2)) = -(1./diff_d).*diff_all_d(:,i).*(dx(:,i)./d(:,i));
        G((3*i-1)) = -(1./diff_d).*diff_all_d(:,i).*(dy(:,i)./d(:,i));
        G((3*i-0)) = -(1./diff_d).*diff_all_d(:,i).*(dz(:,i)./d(:,i));
     end
                
end

