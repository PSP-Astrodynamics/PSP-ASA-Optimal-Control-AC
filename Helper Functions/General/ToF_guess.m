function [t_go] = ToF_guess(r, v, g)
%TOF_GUESS Summary of this function goes here
%   Detailed explanation goes here
% https://arc.aiaa.org/doi/epdf/10.2514/6.1997-3709

Gamma = 0;

a = -2 * dot(v, v) / (Gamma + g ^ 2 / 2);
b = -12 * dot(v, r) / (Gamma + g ^ 2 / 2);
c = -18 * dot(r, r) / (Gamma + g ^ 2 / 2);

alpha = 1/3 * (3 * (a ^ 2 - 4 * c) - 4 * a ^ 2);
beta = 1/27 * (16 * a ^ 3 - 18 * a * (a ^ 2 - 4 * c) - 27 * b ^ 2);

delta = alpha ^ 3 / 26 + beta ^ 2 / 4;

Z = (-beta / 2 + sqrt(delta)) ^ (1 / 3) + (-beta / 2 - sqrt(delta)) ^ (1 / 3);

eta = Z - 2 * a / 3;

zeta = -b / (2 * eta);
xi = (a + eta) / 2;

branch_12 = eta - 4 * (xi - sqrt(eta) * zeta) > 0;
branch_34 = eta - 4 * (xi + sqrt(eta) * zeta) > 0;

t_go = [];
if branch_12
    t_go_12 = [sqrt(eta) + sqrt(eta - 4 * (xi - sqrt(eta) * zeta)), ...
               sqrt(eta) - sqrt(eta - 4 * (xi - sqrt(eta) * zeta))] / 2;
    t_go = t_go_12(t_go_12 > 0);
end
if branch_34
    t_go_34 = [-sqrt(eta) + sqrt(eta - 4 * (xi + sqrt(eta) * zeta)), ...
               -sqrt(eta) - sqrt(eta - 4 * (xi + sqrt(eta) * zeta))] / 2;
    t_go = [t_go, t_go_34(t_go_34 > 0)];
end
t_go = min(t_go);

check = (Gamma + g ^ 2 / 2) * t_go ^ 4 - 2 * dot(v, v) * t_go ^ 2 - 12 * dot(v, r) - t_go - 18 * dot(r, r);

p = [(Gamma + g ^ 2 / 2), 0, -2 * dot(v, v), -1, -12 * dot(v, r) - 18 * dot(r, r)];
r = roots(p);

end

