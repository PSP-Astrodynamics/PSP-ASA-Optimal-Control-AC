function [X_hat_plus, Sigma_hat_r_plus] = RIEKF_approxdisc_time_update(X_hat_minus, Sigma_hat_r_minus, f_u, A_r, Q, dt)
%RIEKF_DISC_TIME_UPDATE Summary of this function goes here
%   Detailed explanation goes here

STM_r = expm(A_r * dt);

X_hat_plus = X_hat_minus.rplus(f_u(X_hat_minus));
Sigma_hat_r_plus = STM_r * Sigma_hat_r_minus * STM_r' + X_hat.Ad * Q * X_hat.Ad' * dt ^ 2;

end

