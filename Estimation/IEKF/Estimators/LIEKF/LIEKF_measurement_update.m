function [X_hat_plus, Sigma_hat_l_plus] = LIEKF_measurement_update(X_hat_minus, Sigma_hat_l_minus, v, S_inv, H, N, X_true)
%LIEKF_MEASUREMENT_UPDATE Summary of this function goes here
%   Detailed explanation goes here

n = size(Sigma_hat_l_minus, 1);

L = Sigma_hat_l_minus * H' * S_inv;
delta_X = L * v;
X_hat_plus = X_hat_minus.rplus(-delta_X);
% Joseph form (full form)
%Sigma_hat_l_plus = (eye(n) - L * H) * Sigma_hat_l_minus * (eye(n) - L * H)' + L * N * L';
% Simplified form (exact if L is optimal so no precision issues)
Sigma_hat_l_plus = (eye(n) - L * H) * Sigma_hat_l_minus;

% Covariance of error needs to be updated now that the new orientation
% error is with respect to the new nominal orientation
G = blkdiag(eye(3) - 1 / 2 * skew(-delta_X(1:3)), eye(6));
Sigma_hat_l_plus = G * Sigma_hat_l_plus * G';

err_l_before = norm(left_invariant_error(X_true, X_hat_minus))
err_l_after = norm(left_invariant_error(X_true, X_hat_plus))

end

function [xi_l] = left_invariant_error(X_true, X_est)
    xi_l = X_true.inv.compose(X_est).Log();
end