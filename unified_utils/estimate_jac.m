function [J] = estimate_jac(f,state)
%ESTIMATE_JAC Estimates jacobian of function f.
%   Estimate jacobian by perturbing state slightly and measuring changes.
% f: function handle in question
% state: vector of the parameters.

nin = numel(state);

step = 0.000001;

vi = f(state);
nout = size(vi, 1);

J = zeros(nout, nin);
for i=1:nin
    ds = zeros(size(state));
    ds(i) = step;
    
    argi = f(state + ds);
    J(:,i) = (argi - vi) / step;
end

end

