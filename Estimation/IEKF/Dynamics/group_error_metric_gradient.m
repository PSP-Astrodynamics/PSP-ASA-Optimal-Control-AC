function [dh, dh_approx] = group_error_metric_gradient(G, Psi, psi)
%GROUP_ERROR_METRIC_GRADIENT Gradient of the quadratic group error metric
%   Group error metric h is 1/2 * norm(psi, P) ^ 2 where P is some
%   weighting matrix.

dh = Psi * G.hat(psi);

Psi_approx = linearized_group_error(psi);
dh_approx = Psi_approx * G.hat(psi);

end

