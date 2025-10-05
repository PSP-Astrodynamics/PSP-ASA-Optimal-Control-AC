function [L_k] = kalman_gain(Ptilde_minus_k, C_k, D_k)
%KALMAN_GAIN Summary of this function goes here
%   Detailed explanation goes here

L_k = Ptilde_minus_k * C_k' / (C_k * Ptilde_minus_k * C_k' + D_k * D_k');
end

