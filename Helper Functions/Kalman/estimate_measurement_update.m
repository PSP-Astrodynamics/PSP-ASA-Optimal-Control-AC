function [xhat_k] = estimate_measurement_update(xhat_minus_k, L_k, ytilde_minus_k)
%ESTIMATE_MEASUREMENT_UPDATE Summary of this function goes here
%   Detailed explanation goes here

xhat_k = xhat_minus_k + L_k * ytilde_minus_k;
end

