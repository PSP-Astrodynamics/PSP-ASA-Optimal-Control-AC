function x_kp1 = rk4(f, t, x, u, p, dt)

    timestep = dt / 2;
   
    k0 = f(t, x, u(t), p);
    est1 = timestep .* k0 + x;

    k1 = f(t + timestep, est1, u(t + timestep), p);
    est2 = timestep .* k1 + x;
    
    k2 = f(t + timestep, est2, u(t + timestep), p);
    est3 = timestep .* k2 * 2 + x;

    k3 = f(t + 2 * timestep, est3, u(t + timestep * 2), p);

    x_kp1 = x + (dt/6) * (k0 + 2 * k1 + 2 * k2 + k3);
end