function [A_k, B_k, E_k, c_k, G_k, Delta] = discretize_stochastic_dynamics_ZOH(f, A, B, E, c, G, N, tspan, x_ref, u_ref, p_ref, tolerances)
    % Discretization of a dynamical system assuming ZOH control'''
    % Make c optional? c = f - Ax - Bu - Ep

    nx = numel(x_ref(:, 1));
    nu = numel(u_ref(:, 1));
    np = numel(p_ref);
    t_k = linspace(tspan(1), tspan(2), N);
    A_k = zeros([nx, nx, N]);
    B_k = zeros([nx, nu, N]);
    E_k = zeros([nx, np, N]);
    c_k = zeros([nx, 1, N]);
    G_k = zeros(size(A_k));
    Delta = zeros([1, N - 1]);
    
    for k = 1:(N - 1)
        [A_k(:, :, k), B_k(:, :, k), E_k(:, :, k), c_k(:, :, k), G_k(:, :, k), x_kp1] = integrate_stochastic_discrete_ZOH(x_ref(:, k), A, B, E, c, f, G, u_ref(:, k), p_ref, [t_k(k), t_k(k + 1)], tolerances);
        Delta(k) = norm(x_kp1 - x_ref(:, k + 1), 2);
    end
end