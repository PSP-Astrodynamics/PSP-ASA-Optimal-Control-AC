function [dPsi] = group_error_dynamics(G, Psi, twist, twist_d)
%GROUP_ERROR_DYNAMICS Dynamical equation for group error
%   Differential equation describing the evolution of the group error
%   G - group
%   Psi - group error
%   twist - twist (velocity in body frame)
%   twist_d - desired twist

velocity_error = twist - G.Ad(Psi.inv) * twist_d;
dPsi = Psi.compose(G.hat(velocity_error));

end

