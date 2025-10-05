function [psidot] = linearized_group_error_dynamics(G, psi, twist, twist_d)
%LINEARIZED_GROUP_ERROR_DYNAMICS Dynamics of linearzed group error
%   Differential equation for tangent space linearization of group error. 
% Based on linear approximation of exp(x) as I + x

psidot = -G.ad(twist_d) * psi + twist - twist_d;

end

