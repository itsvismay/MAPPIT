% payoffs for the Shapley game
R = [0 7; 
     2 6];
S = R';

% use celp to build constraints describing the region of correlated 
% equilibria
ntot = numel(R);
A = -1*[celp(R); celp(S); eye(ntot)];
b = zeros(size(A,1),1);
Aeq = ones(1,ntot);
beq = 1;


% find some random corners of the set of correlated equilibria and
% display them.
for i = 1:5
    %max f(x) = x'Hx/2 + c'x 
    %           s.t. Ax <= b, Aeq*x = beq
    %so..........
    %min -f(x)
    %           s.t. Ax >= b, Aeq*x = beq
    c = -R(:) - S(:);
    H = 1e-4*speye(length(c));
    [pr, val] = quadprog(H, c, A, b, Aeq, beq);
    P = reshape(pr, size(R,1), [])
    val
    fprintf('P1 value %g\n', P(:)'*R(:));
    fprintf('P2 value %g\n\n', P(:)'*S(:));
end


% Build the linear inequalities that represent the rationality constraints
% for a player with payoff matrix R.  Constraints are returned in the
% matrix A: if P is the probability distribution over joint actions (a
% matrix of the same size as R), and if X=P(:), the constraints are
% A * X >= 0.

function [A] = celp(R)

    [n, m] = size(R);
    A = zeros(n*(n-1), n*m);

    for i = 1:n
      for j = [1:i-1 i+1:n]
        constr = R(i,:) - R(j,:);
        A((i-1)*(n-1) + j - (j>i), i:m:end) = constr;
      end
    end
end