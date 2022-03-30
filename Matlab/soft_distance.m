% alpha is the smd alpha value
% X is current agent's curve
% HV is the hashed verts for 
function [D, G] = soft_distance(alpha, AllX, AllV, KDTreeV, r)
     %searches for all neighbors in KDTV within radius
     %of points in X. Returns indexes of X, AllV
%      [D, G] = soft_distance2(alpha, AllX, AllV);
%      return;
     %IDXV -> rows correspond to idx of X, 
     % each row has vector which corresponds to idx of V
     Indices = rangesearch(KDTreeV, AllX, r);
    
     %find non-empty cell functions
     IdxX = find(~cellfun(@isempty, Indices));
     T = Indices(IdxX);
     IdxV = unique(cast([T{:}], 'int16'));
     X = AllX(IdxX, :);   
     V = AllV(IdxV, :);
     %plot soft distance to get a sense of what is going on
     dx = (X(:,1) - V(:,1)');
     dy = (X(:,2) - V(:,2)');
     dz = (X(:,3) - V(:,3)');
    
     d = sqrt(dx.^2 + dy.^2 + dz.^2); 
     diff_all_d = exp(-alpha*d);
     diff_d = sum(sum(diff_all_d));
     if(diff_d<1e-100)
         diff_d = 1e-100;
     end
     GV = zeros(size(V));
     G = zeros(size(AllV));
%      if(diff_d<1e-6)
%          D = 1e1;
%          return;
%      end
     D = -1./alpha.*log(diff_d);
     if(isinf(D) || isnan(D))
         D = 1e1;
         return;
     end
     
    
     %loop over all mesh vertices
     for i=1:size(V,1)
        GV(i, 1) = -(1./diff_d).*sum(diff_all_d(:,i).*(dx(:,i)./d(:,i)));
        GV(i, 2) = -(1./diff_d).*sum(diff_all_d(:,i).*(dy(:,i)./d(:,i)));
        GV(i, 3) = -(1./diff_d).*sum(diff_all_d(:,i).*(dz(:,i)./d(:,i)));
        if(sum(isnan(GV(i,:)))>0)
            GV(i,:) = 1e-3*rand(3,1);
        end
     end
     G(IdxV,:) = GV;
                
end
