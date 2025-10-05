function [x_kp1, x_est] = rk4_explicit(f, t, x, u, p, dt)

    timestep = dt / 2;
   
    k0 = f(t, x, u(:, 1), p);
    est1 = timestep .* k0 + x;

    k1 = f(t + timestep, est1, u(:, 2), p);
    est2 = timestep .* k1 + x;
    
    k2 = f(t + timestep, est2, u(:, 2), p);
    est3 = timestep .* k2 * 2 + x;

    k3 = f(t + 2 * timestep, est3, u(:, 3), p);

    x_kp1 = x + (dt/6) * (k0 + 2 * k1 + 2 * k2 + k3);

    x_est = [est1, est2, est3];
end