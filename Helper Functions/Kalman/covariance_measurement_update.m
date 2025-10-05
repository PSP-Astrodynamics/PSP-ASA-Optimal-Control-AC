function [Ptilde_k] = covariance_measurement_update(L_k, C_k, Ptilde_minus_k, D_k)
%COVARIANCE_MEASUREMENT_UPDATE Summary of this function goes here
%   Detailed explanation goes here

L_k_C_k = L_k * C_k;
I = eye(size(L_k_C_k));

Ptilde_k = (I - L_k_C_k) * Ptilde_minus_k * (I - L_k_C_k)' + L_k * D_k * D_k' * L_k';
end

