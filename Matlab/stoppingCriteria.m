function stop = outfun(x,optimValues,state)
%instructions here, function stops fmincon
%https://www.mathworks.com/help/optim/ug/output-function.html#f11454

stop = false;
optimValues.stepsize; %norm of step size

% Check whether directional derivative norm is less than .01.
% if it is that means the sum of the delta of agents trajectories is less
% than 1e-2, which means basically converged
if (optimValues.stepsize < 1e-2 || optimValues.firstorderopt<1e-2) && optimValues.iteration>0
    stop = true;
end

end 
