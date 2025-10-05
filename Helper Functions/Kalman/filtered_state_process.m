function [xhat_kp1] = filtered_state_process(xhat_k, u_k, A_k, B_k, c_k, L_kp1, ytilde_minus_kp1)
%FILTERED_STATE_PROCESS Summary of this function goes here
%   Detailed explanation goes here

xhat_kp1 = A_k * xhat_k + B_k * u_k + c_k + L_kp1 * ytilde_minus_kp1;
end

