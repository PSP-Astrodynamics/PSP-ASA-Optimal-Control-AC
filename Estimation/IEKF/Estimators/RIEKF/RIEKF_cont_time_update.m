function [X_hat_plus, Sigma_hat_r_plus] = RIEKF_cont_time_update(X_hat_minus, Sigma_hat_r_minus, f_u, A_r, Q, tspan, tolerances)
%RIEKF_CONT_TIME_UPDATE Right invariant extended Kalman filter time update
%   Time update step for error state Kalman filter on a Lie group with 
% right invariant error.
% TODO:
% - Use general group composition/action instead of assuming matrix group

% Continuously propogate over the interval between measurements
[t_m, y_m] = ode45(@(t, y) cont_RIEKF_DE(t, y, , tolerances);

x_m = y_m(:, nx + (1:nx))';
xhat_minus = y_m(:, 2 * nx + (1:nx))';
Phat_minus = permute(reshape(y_m(:, (3 * nx + 1):end), [N_sub + 1, nx, nx]), [2, 3, 1]);


end


function [ydot] = cont_RIEKF_DE(t, y, f_u, A_r, Q, G)
    x_hat = reshape(y(G.size), size(G.identity));
    Sigma_hat = reshape(y((G.size + 1):end), [G.nx, G.nx]);
    
    X_hat = G;
    X_hat.element = x_hat;

    X_hat_dot = f_u(X_hat);
    Sigma_hat_dot = A_r * Sigma_hat + Sigma_hat * A_r' + X_hat.Ad * Q * X_hat.Ad';

    ydot = [X_hat_dot(:); Sigma_hat_dot(:)];
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