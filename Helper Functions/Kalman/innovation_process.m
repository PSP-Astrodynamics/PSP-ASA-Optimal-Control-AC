function [ytilde_minus_k] = innovation_process(y_k, C_k, xhat_minus_k)
%INNOVATION_PROCESS Summary of this function goes here
%   Detailed explanation goes here

ytilde_minus_k = y_k - C_k * xhat_minus_k;

% Should be equal to C_k xtilde_k + D_k v_k
end

