function [F_PD] = GroupPID(G, twist, twist_d, psi, psi_integ, K_p, K_i, K_d)
%GROUPPID PID controller on Lie group
%   Proportional-integral-derivative controller on Lie group using group error in
%   tangent space.

psidot = linearized_group_error_dynamics(G, psi, twist, twist_d);

F_PD = K_p * psi + K_i * psi_integ + K_d * psidot;

end

