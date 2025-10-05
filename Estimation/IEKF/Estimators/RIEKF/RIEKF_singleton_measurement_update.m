function [X_hat_plus, Sigma_hat_r_plus] = RIEKF_singleton_measurement_update(X_hat_minus, Sigma_hat_r_minus, z, b, H, R, M, Pi)
%RIEKF_MEASUREMENT_UPDATE Summary of this function goes here
%   Detailed explanation goes here

n = size(Sigma_hat_r_minus, 1);

S = H * Sigma_hat_r_minus * H' + R * M * R';
K = Sigma_hat_r_minus * H' / S;
delta_X = K * (Pi * X_h * z - Pi * b);
X_hat_plus = exp(-delta_X) * X_hat_minus;
Sigma_hat_r_plus = (eye(n) - K * H) * Sigma_hat_r_minus;

end

