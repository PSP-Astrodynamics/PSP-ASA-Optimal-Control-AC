function [t_sub, x_array, u_array] = propagate_cont_feedback_no_kalman_filter(x_0, x_ref, u_ref, K_k, f, G, t_k, N_sub, w_k_func, delta_t, tolerances)
%CONT_KALMAN_FILTER Continuous propagation of continuous-discrete Kalman
%filter
%   Note - control is constant between measurements

nx = numel(x_0);

t_sub = linspace(t_k(1), t_k(end), N_sub * (numel(t_k) - 1) + 1);

x_array = zeros([nx, numel(t_sub)]);
x_array(:, 1) = x_0;

u_array = zeros([size(u_ref, 1), numel(t_sub) - 1]);

for k = 1:(numel(t_k) - 1)
    k_0 = (k - 1) * N_sub;
    u_k = u_ref(:, k) + K_k(:, :, k) * (x_array(:, k_0 + 1) - x_ref(:, k));
    % Continuously propogate over the interval between measurements (with constant control)

    [t_m, x_m, u_m] = sode45(f, G, @(t, x) u_k, 0, w_k_func, t_sub(k_0 + (1:(N_sub + 1))), delta_t, x_array(:, k_0 + 1), tolerances, w_k_func = w_k_func);

    x_array(:, k_0 + (1:(N_sub + 1))) = x_m;
    u_array(:, min(k_0 + (1:(N_sub + 1)), numel(t_sub) - 1)) = u_m;
end
end