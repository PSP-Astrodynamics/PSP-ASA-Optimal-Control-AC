function [P_ytilde_minus_k] = innovation_process_covariance(C_k, Ptilde_minus_k, D_k)
%INNOVATION_PROCESS_COVARIANCE Summary of this function goes here
%   Detailed explanation goes here

P_ytilde_minus_k = C_k * Ptilde_minus_k * C_k' + D_k * D_k';
end

