function [x_hat_plus, Sigma_hat_plus] = EKF_disc_time_update(x_hat_minus, Sigma_hat_minus, Q, dt, u, g, A, G, A_rot, rotate)
    theta_minus = x_hat_minus(7:9);
    v_minus = x_hat_minus(4:6);
    r_minus = x_hat_minus(1:3);

    %theta_plus = rk4_explicit(@(t, x, u, p) A_rot(x, u), 0, theta_minus, repmat(u(4:6), 1, 3), [], dt);
    theta_plus = theta_minus + A_rot(theta_minus, u(4:6)) * dt;
    v_plus = v_minus + rotate(theta_minus, u(1:3) * dt) + g * dt;
    r_plus = r_minus + v_minus * dt + 1 / 2 * rotate(theta_minus, u(1:3) * dt ^ 2) + 1 / 2 * g * dt ^ 2;

    x_hat_plus = [r_plus; v_plus; theta_plus];

    STM = expm(A(x_hat_minus, u) * dt);
    Phi_u = expm(G(x_hat_minus, u) * dt); % Correct?

    Sigma_hat_plus = STM * Sigma_hat_minus * STM' + Phi_u * Q * Phi_u';
end