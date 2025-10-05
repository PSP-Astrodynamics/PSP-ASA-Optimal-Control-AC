function [Ptilde_minus_k] = covariance_time_update(A_km1, Ptilde_km1, G_km1)
%COVARIANCE_TIME_UPDATE Summary of this function goes here
%   Detailed explanation goes here

Ptilde_minus_k = A_km1 * Ptilde_km1 * A_km1' + G_km1 * G_km1';
end

