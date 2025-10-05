function [x_hat_plus, Sigma_hat_plus] = ESKF_disc_time_update(x_hat_minus, Sigma_hat_minus, Q, dt, u, g, A, G, x_true)
    q_minus = x_hat_minus(7:10);
    v_minus = x_hat_minus(4:6);
    r_minus = x_hat_minus(1:3);

    q_plus = q_mul(q_minus, qexp(u(4:6) * dt));
    v_plus = v_minus - quat_rot(q_conj(q_minus), u(1:3) * dt) + g * dt;
    r_plus = r_minus + v_minus * dt - 1 / 2 * quat_rot(q_conj(q_minus), u(1:3) * dt ^ 2) + 1 / 2 * g * dt ^ 2;

    x_hat_plus = [r_plus; v_plus; q_plus];

    STM = expm(A(x_hat_minus, u) * dt);
    Phi_u = expm(G(x_hat_minus, u) * dt); % Correct?

    Sigma_hat_plus = STM * Sigma_hat_minus * STM' + Phi_u * Q * Phi_u';
end

function [q] = qexp(tau)
    theta = norm(tau);
    u = tau / theta;
    q = [u * sin(theta / 2); cos(theta / 2)];
end