function [phi] = l1_approx_penalty(z, tau, variant)
%L1_APPROX_PENALTY Summary of this function goes here
%   Detailed explanation goes here
arguments
    z
    tau = 1.1
    variant = "standard"
end

if variant == "standard"
    phi = 1 / tau * pow_abs(z, tau) + 1 / 2 * z .^ 2;
elseif variant == "positive"
    phi = 1 / tau * pow_pos(z, tau) + 1 / 2 * z .^ 2;
end

end