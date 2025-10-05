function [t_sub, x_array, xhat_array, Phat_array, u_array] = propagate_cont_feedback_kalman_filter(x_0, Ptilde_0, p, f, G, A_ref, B_ref, c_ref, x_ref, u_ref, K, L_k, C_k, D_k, f_0, g_0, t_k, tspan, N_sub, w_k_func, v_k, tolerances)
%CONT_KALMAN_FILTER Continuous propagation of continuous-discrete Kalman
%filter
%   Note - control is constant between measurements

nx = numel(x_0(:, 1));

t_sub = linspace(t_k(1), t_k(end), N_sub * (numel(t_k) - 1) + 1);
if t_k(end) ~= tspan(end)
    t_sub = [t_sub(1:(end - 1)), linspace(t_k(end), tspan(end), max(ceil((t_k(end) - tspan(end)) / (t_sub(2) - t_sub(1))), 3))];
end

x_array = zeros([nx, numel(t_sub)]);
x_array(:, 1) = x_0(:, 1);
xhat_array = zeros([nx, numel(t_sub)]);
xhat_array(:, 1) = x_0(:, 2);
Phat_array = zeros([nx, nx, numel(t_sub)]);
Phat_array(:, :, 1) = Ptilde_0;
u_array = zeros(numel(u_ref(:, 1)), numel(t_sub));

u_ref_func = @(t) interp1(t_k(1:size(u_ref, 2)), u_ref', t, "previous", "extrap")';

for k = 1:(numel(t_k) - 1)
    if k ~= (numel(t_k) - 1) || t_k(end) == tspan(end)
        k_0 = (k - 1) * N_sub;
        y_k = [reshape(x_ref(:, k), [], 1); reshape(x_array(:, k_0 + 1), [], 1); reshape(xhat_array(:, k_0 + 1), [], 1); reshape(Phat_array(:, :, k_0 + 1), [], 1)];
        u_array(:, k_0 + (1:(N_sub + 1))) = repmat(u_ref_func(t_k(k)) + K(:, :, k) * (xhat_array(:, k_0 + 1) - x_ref(:, k)), 1, N_sub + 1);
        % Continuously propogate over the interval between measurements (with constant control)
        [t_m, y_m] = ode45(@(t, y) cont_kalman_filter_DE(t, y, u_array(:, k_0 + 1), p, f, G, A_ref, B_ref, c_ref, u_ref_func, w_k_func, nx), t_sub(k_0 + (1:(N_sub + 1))), y_k, tolerances);
    
        x_m = y_m(:, nx + (1:nx))';
        xhat_minus = y_m(:, 2 * nx + (1:nx))';
        Phat_minus = permute(reshape(y_m(:, (3 * nx + 1):end), [N_sub + 1, nx, nx]), [2, 3, 1]);
    
        x_array(:, k_0 + (1:(N_sub + 1))) = x_m;
        xhat_array(:, k_0 + (1:(N_sub + 1))) = xhat_minus;
        Phat_array(:, :, k_0 + (1:(N_sub + 1))) = Phat_minus;
    
        % Perform discrete measurement update
        y_k = f_0(t_k(k), x_m(:, end), u_array(:, k_0 + 1)) + g_0(x_m(:, end), u_array(:, k_0 + 1)) * v_k(:, k);
        ytilde_minus_k = innovation_process(y_k, C_k(:, :, k), xhat_minus(:, end));
        xhat_array(:, k * N_sub + 1) = estimate_measurement_update(xhat_minus(:, end), L_k(:, :, k), ytilde_minus_k);
        Phat_array(:, :, k * N_sub + 1) = covariance_measurement_update(L_k(:, :, k), C_k(:, :, k), Phat_minus(:, :, end), D_k(:, :, k));
    else
        k_0 = (k - 1) * N_sub;
        [t_m, y_m] = ode45(@(t, y) cont_kalman_filter_DE(t, y, u_array(:, k_0 + 1), p, f, G, A_ref, B_ref, c_ref, u_ref_func, w_k_func, nx), t_sub(((k - 1) * N_sub + 1):end), y_0, tolerances);

        x_m = y_m(:, nx + (1:nx))';
        xhat_minus = y_m(:, 2 * nx + (1:nx))';
        Phat_minus = permute(reshape(y_m(:, (3 * nx + 1):end), [[], nx, nx]), [2, 3, 1]);
    
        x_array(:, ((k - 1) * N_sub + 1):end) = x_m;
        xhat_array(:, ((k - 1) * N_sub + 1):end) = xhat_minus;
        Phat_array(:, :, ((k - 1) * N_sub + 1):end) = Phat_minus; 
    end
end
end

function [ydot] = cont_kalman_filter_DE(t, y, u, p, f, G, A_ref, B_ref, c_ref, u_ref, w_func, nx)
    x_mean = y(1:nx);
    x = y(nx + (1:nx));
    xhat = y(2 * nx + (1:nx));
    Phat = reshape(y((3 * nx + 1):end), [nx, nx]);

    xdot_mean = f(t, x_mean, u_ref(t), p);
    xdot = f(t, x, u, p) + G(t, x, u, p) * w_func(t);
    xhatdot = A_ref(t, x_mean, u_ref(t), p) * xhat + B_ref(t, x_mean, u_ref(t), p) * u + c_ref(t, x_mean, u_ref(t), p);
    Phatdot = continuous_kalman_covariance_derivative(A_ref(t, x_mean, u_ref(t), p), Phat, G(t, x_mean, u_ref(t), p));

    ydot = [xdot_mean(:); xdot(:); xhatdot(:); Phatdot(:)];
end