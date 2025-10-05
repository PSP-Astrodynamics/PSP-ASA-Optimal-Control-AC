function [x_hat_plus, Sigma_hat_plus] = ESKF_disc_time_update_localrot(x_hat_minus, Sigma_hat_minus, Q, dt, u, g, A, G, x_true)
    q_minus = x_hat_minus(7:10);
    v_minus = x_hat_minus(4:6);
    r_minus = x_hat_minus(1:3);

    q_plus = q_mul(q_minus, qexp(u(4:6) * dt));
    v_plus = v_minus + quat_rot(q_minus, u(1:3) * dt) + g * dt;
    r_plus = r_minus + v_minus * dt + 1 / 2 * quat_rot(q_minus, u(1:3) * dt ^ 2) + 1 / 2 * g * dt ^ 2;

    R_plus = quat_rotmatrix(q_plus);

    S3_minus = SO3_quaternion(q_minus);
    R_minus2 = S3_minus.R;
    q_plus2 = S3_minus.rplus(u(4:6) * dt);
    R_minus = SO3_RotationMatrix(quat_rotmatrix(q_minus));
    q_minus2 = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, 0] * dcm2quat(R_minus.R)';
    R_kp1 = R_minus.rplus(u(4:6) * dt).element;



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

function [R] = quat_rotmatrix(q)
    w = q(4);
    v = q(1:3);
    x = v(1);
    y = v(2);
    z = v(3);
    R = [w ^ 2 + x ^ 2 - y ^ 2 - z ^ 2, 2 * (x * y - w * z), 2 * (x * z + w * y);
         2 * (x * y + w * z), w ^ 2 - x ^ 2 + y ^ 2 - z ^ 2, 2 * (x * z - w * x);
         2 * (x * z - w * y), 2 * (y * z + w * x), w ^ 2 - x ^ 2 - y ^ 2 + z ^ 2]';
end