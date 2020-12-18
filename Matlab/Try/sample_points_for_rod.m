function [P,J] = sample_points_for_rod(X, samples, varargin)
%     ts = linspace(X(1,3)+0.001, X(end,3)-0.001, samples);
%     PV = [interp1(X(:,3), X(:, 1:2), ts, 'linear') ts'];
%     PE = [(1:size(PV,1)-1)' (2:size(PV,1))'];
%     if numel(varargin)>0 && strcmp('cylinders',varargin{1})
%         rad = varargin{2}; 
%         
%         [P,J,~,~] = edge_cylinders(PV,PE, 'Thickness',2*rad, 'PolySize', 10);
%     else
%         % samples along centerline
%         P = PV;
%         J = PE;
%     end

    s = 10; %5 samples per edge
    %for each edge
    P = zeros(s*(size(X,1)-1)+1, 3);
    J = zeros(size(P,1), size(X,1));
    for i = 1:size(X,1)-1
        tl = X(i,3);
        tr = X(i+1,3);
        xl = X(i,1);
        xr = X(i+1,1);
        yl = X(i,2);
        yr = X(i+1, 2);
        for j =1:s
            P(s*(i-1)+j,3) = tl + (j-1)*(tr - tl)/s;
            P(s*(i-1)+j,1) = xl + (j-1)*(xr - xl)/s;
            P(s*(i-1)+j,2) = yl + (j-1)*(yr - yl)/s;
            J(s*(i-1)+j,i) = 1 -(j-1)/s;
            J(s*(i-1)+j,i+1) = (j-1)/s;
        end
    end
    P(end,:) = X(end,:);
    J(end,end-1) = 0;
    J(end,end) = 1;
    
  
end