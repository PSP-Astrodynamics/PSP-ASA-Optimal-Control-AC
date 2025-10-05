function [X_hat_plus, Sigma_hat_r_plus] = LIEKF_approxdisc_time_update(X_hat_minus, Sigma_hat_r_minus, A_r, Q, dt)
%LIEKF_DISC_TIME_UPDATE Summary of this function goes here
%   Detailed explanation goes here

% 1st order approximation of matrix exponential
STM_r = eye(size(A_r)) + A_r * dt;

X_hat_plus = X_hat_minus;
X_hat_plus.element = STM_r * X_hat_minus.element;
Sigma_hat_r_plus = STM_r * Sigma_hat_r_minus * STM_r' + Q * dt;

end

