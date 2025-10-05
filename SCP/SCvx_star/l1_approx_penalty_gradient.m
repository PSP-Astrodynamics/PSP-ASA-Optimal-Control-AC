function [nabla_phi] = l1_approx_penalty_gradient(z, tau)
%L1_APPROX_PENALTY_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
arguments
    z
    tau = 1.1
end

nabla_phi = sign(z) .* abs(z) .^ (tau - 1) + z;

end