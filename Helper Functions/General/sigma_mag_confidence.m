function [n_sigma] = sigma_mag_confidence(epsilon, nu)
%N_SIGMA How many sigmas to satisfy magnitude constraint with 1-epsilon
%confidence
%   Detailed explanation goes here

n_sigma = sqrt(chi2inv(1 - epsilon, nu));
end

