function [VV,EE, newA, BVind] = spacetime_graph(V,E,t, Bind)
  nt = size(t,1);
  nV = size(V,1);
  nE = size(E,1);
  ver = [repmat(V,nt,1) kron(t,ones(nV,1))];
  EE_space = repmat(E,nt,1) + kron(cumsum([0;nV*ones(nt-1,1)]),ones(nE,1));
  EE_time = [1:nV*(nt-1) ; nV+1:nV*nt]';
  EE = [ EE_space; EE_time ];
%   EE_cross = [ EE_space(1:nE*(nt-1),1) EE_space(nE+1:nE*nt,2) ];
%   EE = [ EE_space; EE_time; EE_cross ];

  BVind = repmat(Bind, nt, 1) + kron(linspace(0,nt-1,nt)', nV*ones(size(Bind,1),1));
  
  VV = zeros(length(ver),3);
  VV(:,1) = ver(:,1);
  VV(:,2) = ver(:,2);
  VV(:,3) = ver(:,4);
  %tsurf(EE,VV);
  
  % adjacency matrix of 3d graph
  A_space = adjacency_matrix(EE_space);
  A_time = adjacency_matrix(EE_time);
  %A_cross = adjacency_matrix(EE_cross);
  newA = A_space + tril(A_time); %+ tril(A_cross);
  %newA = adjacency_matrix(EE);
  
end
