function [Psi_lin] = linearized_group_error(G, psi)
%LINEARIZED_GROUP_ERROR Linearized group error equation
%   Tangent space linearization of group error. Based on linear 
% approximation of exp(x) as I + hat(x)
%   psi - group error in tangent space

% Psi = exp(psi) approx I + psi
Psi_lin = eye(numel(psi)) + G.hat(psi);

end

