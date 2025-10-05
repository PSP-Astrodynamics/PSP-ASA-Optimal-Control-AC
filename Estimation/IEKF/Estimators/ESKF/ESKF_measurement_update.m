function [x_hat_plus, Sigma_hat_plus] = ESKF_measurement_update(x_hat_minus, Sigma_hat_minus, Y, b, H, N, add, mult, x_true)
%LESKF_MEASUREMENT_UPDATE Summary of this function goes here
%   Detailed explanation goes here

n = size(Sigma_hat_minus, 1);

S = H * Sigma_hat_minus * H' + N;
S_inv = inv(S);
L = Sigma_hat_minus * H' * S_inv;

delta_x = L * (Y - b);% .* mult;
x_hat_plus = add(x_hat_minus, delta_x);
% Joseph form (full form)
%Sigma_hat_plus = (eye(n) - L * H) * Sigma_hat_minus * (eye(n) - L * H)' + L * N * L';
% Simplified form (exact if L is optimal so no precision issues)
Sigma_hat_plus = (eye(n) - L * H) * Sigma_hat_minus;

% Covariance of error needs to be updated now that the new orientation
% error is with respect to the new nominal orientation
G = blkdiag(eye(6), eye(3) - 1 / 2 * skew(delta_x(7:9)));
Sigma_hat_plus = G * Sigma_hat_plus * G';

norm(quat_left_error(x_true(7:10), x_hat_minus(7:10)));
norm(quat_left_error(x_true(7:10), x_hat_plus(7:10)));
end

function [tau] = qLog(q)
    N = size(q, 2);
    for k = 1 : N
        w = q(4, k);
        v = q(1:3, k);
        w = w * sign(w);
        v = v * sign(w);
        tau(:, k) = 2 * v * atan2(norm(v), w) / norm(v);
    end
end

function [xi_theta_l] = quat_left_error(q_true, q_est)
    % q is from inertial to body frame
    xi_theta_l = qLog(q_mul_array(q_true, q_conj(q_est)));
end

