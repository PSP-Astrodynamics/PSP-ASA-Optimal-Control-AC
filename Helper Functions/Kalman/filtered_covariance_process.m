function [X_kp1] = filtered_covariance_process(A_k, X_k, B_k, S_k, L_kp1, Ptilde_minus_y_k)
%FILTERED_COVARIANCE_PROCESS Summary of this function goes here
%   Detailed explanation goes here
% S_k = sqrtm(P_u_k) should equal K_k X_k
% X_k = sqrtm(Phat_k)

X_kp1 = [A_k * X_k + B_k * S_k, L_kp1 * sqrtm(Ptilde_minus_y_k)];
end

