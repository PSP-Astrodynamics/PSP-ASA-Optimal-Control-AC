function [x_hat_plus, Sigma_hat_plus] = EKF_measurement_update(x_hat_minus, Sigma_hat_minus, Y, b, H, N)
%LIEKF_MEASUREMENT_UPDATE Summary of this function goes here
%   Detailed explanation goes here

n = size(Sigma_hat_minus, 1);

S = H * Sigma_hat_minus * H' + N;
S_inv = inv(S);
L = Sigma_hat_minus * H' * S_inv;
delta_x = L * (Y - b);
x_hat_plus = x_hat_minus + delta_x;
% Joseph form (full form)
%Sigma_hat_plus = (eye(n) - L * H) * Sigma_hat_minus * (eye(n) - L * H)' + L * N * L';
% Simplified form (exact if L is optimal so no precision issues)
Sigma_hat_plus = (eye(n) - L * H) * Sigma_hat_minus;

end

