function [] = plot_3DoFc_time_histories(t, x, u)
%PLOT_3DOF_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

Nu = size(u, 2);

tiledlayout(3, 2)

nexttile
plot(t, x(1:2, :));
legend("r_" + string(1:2))
title("Position vs Time")
xlabel("Time [s]")
ylabel("Position [km]")
grid on

nexttile
plot(t, x(3:4, :) * 1e3);
legend("v_" + string(1:2))
title("Velocity vs Time")
xlabel("Time [s]")
ylabel("Velocity [m / s]")
grid on

nexttile
plot(t, wrapTo180(rad2deg(x(5, :))));
title("Orientation vs Time")
xlabel("Time [s]")
ylabel("Orientation [deg]")
grid on

nexttile
plot(t, rad2deg(x(6, :)));
title("Angular Velocity vs Time")
xlabel("Time [s]")
ylabel("Angular Velocity [deg / s]")
grid on

nexttile
plot(t(1:Nu), u(3, :) .* exp(x(7, 1:Nu)))
title("Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [kg km / s2]")
grid on

nexttile
gimbal_angle = atan2d(u(2, :), u(1, :));
plot(t(1:Nu), gimbal_angle)
title("Gimbal Angle vs Time")
xlabel("Time [s]")
ylabel("Gimbal Angle [deg]")
grid on

sgtitle("3DoF Rocket Landing Time Histories")

end

