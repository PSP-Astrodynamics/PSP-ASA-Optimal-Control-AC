function [F_PD] = GroupPD(G, twist, twist_d, psi, K_p, K_d)
%GROUPPD PD controller on Lie group
%   Proportional-derivative controller on Lie group using group error in
%   tangent space.

psidot = linearized_group_error_dynamics(G, psi, twist, twist_d);

F_PD = K_p * psi + K_d * psidot;

end

