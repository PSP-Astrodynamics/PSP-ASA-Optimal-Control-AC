function [Pdot] = continuous_kalman_covariance_derivative(A, P, G)
%CONTINUOUS_KALMAN_COVARIANCE_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

Pdot = A * P + P * A' + G * G';
end

