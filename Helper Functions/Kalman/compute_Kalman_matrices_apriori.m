function [Ptilde_minus_k, Ptilde_k, L_k] = compute_Kalman_matrices_apriori(P0, A_k, G_k, C_k, D_k)
%COMPUTE_KALMAN_MATRICES_APRIORI Summary of this function goes here
%   Detailed explanation goes here

n_k = size(A_k, 3);

Ptilde_minus_k = zeros(size(A_k));
Ptilde_minus_k(:, :, 1) = P0;

L_k = zeros(size(D_k, 2), size(D_k, 1)); % L_k in M^(nx by ny)

Ptilde_k = zeros(size(A_k));
Ptilde_k(:, :, 1) = P0;

P_ytilde_minus_k = zeros(size(D_k, 1), size(D_k, 1)); % ny by ny
P_ytilde_minus_k(:, :, 1) = innovation_process_covariance(C_k(:, :, 1), Ptilde_minus_k(:, :, 1), D_k(:, :, 1));

for k = 1:(n_k - 1)
    Ptilde_minus_k(:, :, k + 1) = covariance_time_update(A_k(:, :, k), Ptilde_k(:, :, k), G_k(:, :, k));
    L_k(:, :, k + 1) = kalman_gain(Ptilde_minus_k(:, :, k + 1), C_k(:, :, k + 1), D_k(:, :, k + 1));
    Ptilde_k(:, :, k + 1) = covariance_measurement_update(L_k(:, :, k + 1), C_k(:, :, k + 1), Ptilde_minus_k(:, :, k + 1), D_k(:, :, k + 1));
    P_ytilde_minus_k(:, :, k + 1) = innovation_process_covariance(C_k(:, :, k + 1), Ptilde_minus_k(:, :, k + 1), D_k(:, :, k + 1));
    % Compute X_k, S_k, Phat, P_uk?
end
end

