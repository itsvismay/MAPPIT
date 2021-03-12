function [VV,EE] = spacetime_graph(V,E,t)
  nt = size(t,1);
  nV = size(V,1);
  nE = size(E,1);
  VV = [repmat(V,nt,1) kron(t,ones(nV,1))];
  EE_space = repmat(E,nt,1) + kron(cumsum([0;nV*ones(nt-1,1)]),ones(nE,1));
  EE_time = [1:nV*(nt-1) ; nV+1:nV*nt]';
  EE_cross = [ EE_space(1:nE*(nt-1),1) EE_space(nE+1:nE*nt,2) ];
  EE = [ EE_space ; EE_time ; EE_cross ];
end
