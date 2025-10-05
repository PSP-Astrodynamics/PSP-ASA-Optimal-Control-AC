function [X_hat_plus, Sigma_hat_l_plus] = LIEKF_disc_time_update(X_hat_minus, Sigma_hat_l_minus, Q, dt, u, g, X_true)
%LIEKF_DISC_TIME_UPDATE Summary of this function goes here
%   Detailed explanation goes here

STM_l = X_hat_minus.STM_l(u, dt);
Sigma_hat_l_plus = STM_l * Sigma_hat_l_minus * STM_l' + Q * dt;

% Constant input model 
X_hat_plus = X_hat_minus.approx_dynamics(u, g, dt, X_true);

end

