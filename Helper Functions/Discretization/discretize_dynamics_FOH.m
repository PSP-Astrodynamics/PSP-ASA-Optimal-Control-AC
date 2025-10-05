function [A_k, B_k_plus, B_k_minus, E_k, c_k, Delta_vec] = discretize_dynamics_FOH(f, A, B, E, c, N, tspan, x_ref, u_ref, p_ref, tolerances)
    % Discretization of a dynamical system assuming FOH control
    % Make c optional? c = f - Ax - Bu

    nx = numel(x_ref(:, 1));
    nu = numel(u_ref(:, 1));
    np = numel(p_ref);

    t_k = linspace(tspan(1), tspan(2), N);
    A_k = zeros([nx, nx, N - 1]);
    B_k_plus = zeros([nx, nu, N - 1]);
    B_k_minus = zeros([nx, nu, N - 1]);
    E_k = zeros([nx, np, N - 1]);
    c_k = zeros([nx, 1, N - 1]);
    Delta_vec = zeros([nx, N - 1]);

    %u_ref = @(t) interp1(t_k(1:size(u_ref, 2)), u_ref', t, "linear", "extrap")';

    for k = 1:(N - 1)
        u_ref_k = [u_ref(:, k), u_ref(:, k + 1)];
        [A_k(:, :, k), B_k_plus(:, :, k), B_k_minus(:, :, k), E_k(:, :, k), c_k(:, :, k), x_kp1] = integrate_discrete_FOH(x_ref(:, k), A, B, E, c, f, u_ref_k, p_ref, [t_k(k), t_k(k + 1)], tolerances);
        Delta_vec(:, k) = x_kp1 - x_ref(:, k + 1);
    end
end