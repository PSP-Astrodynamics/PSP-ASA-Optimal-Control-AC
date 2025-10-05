function [A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, tspan, x_ref, u_ref, p_ref, tolerances)
    % Discretization of a dynamical system assuming ZOH control'''
    % Make c optional? c = f - Ax - Bu - Ep

    nx = numel(x_ref(:, 1));
    nu = numel(u_ref(:, 1));
    np = numel(p_ref);

    t_k = linspace(tspan(1), tspan(2), N);
    A_k = zeros([nx, nx, N - 1]);
    B_k = zeros([nx, nu, N - 1]);
    E_k = zeros([nx, np, N - 1]);
    c_k = zeros([nx, 1, N - 1]);
    Delta = zeros([nx, N - 1]);
    
    for k = 1:(N - 1)
        [A_k(:, :, k), B_k(:, :, k), E_k(:, :, k), c_k(:, :, k), x_kp1] = integrate_discrete_ZOH(x_ref(:, k), A, B, E, c, f, u_ref(:, k), p_ref, [t_k(k), t_k(k + 1)], tolerances);
        Delta(:, k) = x_kp1 - x_ref(:, k + 1);
    end
end