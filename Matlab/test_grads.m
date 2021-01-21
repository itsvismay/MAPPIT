V = [ -1.0278       5.5347            0;
     -0.82776       4.0588          0.5;
     -0.77776       2.5829            1;
     -0.77776       1.1069          1.5;
     -0.77776     -0.36898            2;
     -0.77776      -1.8449          2.5;
     -0.97776      -1.8449            3;
     -0.87776      -2.9518          3.5;
     -0.87776      -3.6898            4;
      -1.0278      -4.0588          4.5;
      -1.0278      -5.5347           10]; %l shape
E = [1:size(V,1)-1; 2:size(V,1)]';


V1 = V(1:(end-2),:) - V(2:(end-1),:);
V2 = V(3:end,:) - V(2:(end-1),:);
%V2 = V(E(:,2),:);
V12 = cross(V1,V2);
V12norm = sqrt(V12(:,1).^2+V12(:,2).^2+V12(:,3).^2);
V1norm = sqrt(V1(:,1).^2+V1(:,2).^2+V1(:,3).^2);
V2norm = sqrt(V2(:,1).^2+V2(:,2).^2+V2(:,3).^2);

if(V12norm < eps)
    V12norm = eps;
end

zhat = V12./V12norm;


theta = 2*atan2(dot(V12, zhat, 2), V1norm.*V2norm + dot(V1,V2,2));
theta0 = pi;
%bending energy 
k = 100;
e = 0.5*100*(theta-theta0).^2;

%bending energy gradient
grad = [k*theta.*(cross(V1,zhat)./(V1norm.*V1norm)); ...
        -k*(theta-0).*((cross(V2,zhat)./(V2norm.*V2norm)) - (cross(V1,zhat)./(V1norm.*V1norm))); ...
       k*theta.*(cross(V2,zhat)./(V2norm.*V2norm))]; 
    

hold on;
plot(V(:,1), V(:,2), 'r');
plot(V(:,1), V(:,2), 'r*');
line([V(:,1)'; V(:,1)' + 0.01*grad(:,1)'], [V(:,2)'; V(:,2)' + 0.01*grad(:,2)'], 'Color', [0 0 1]);
hold off;
axis equal;
