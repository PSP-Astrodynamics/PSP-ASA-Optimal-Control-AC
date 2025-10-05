function [J_p] = penalty_function(lambda, xi, mu, zeta, phi, w)
    J_p = lambda' * xi + sum(1 / w * phi(w * xi, "standard")) ...
          + mu' * zeta + sum(1 / w * phi(w * zeta, "positive")); 
end