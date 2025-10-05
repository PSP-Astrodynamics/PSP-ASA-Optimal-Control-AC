function [] = plot_6DoF_time_histories(t, x, u)
%PLOT_3DOF_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

Nu = size(u, 2);

tiledlayout(3, 2)

nexttile
plot(t, x(1:3, :));
legend("r_" + string(1:3))
title("Position vs Time")
xlabel("Time [s]")
ylabel("Position [km]")
grid on

nexttile
plot(t, x(4:6, :) * 1e3);
legend("v_" + string(1:3))
title("Velocity vs Time")
xlabel("Time [s]")
ylabel("Velocity [m / s]")
grid on

nexttile
plot(t(1:Nu), u(1:3, :)); hold on
plot(t(1:Nu), u(4, :), LineStyle = "--"); hold off
legend(["v_" + string(1:3), "||v||"])
title("Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [kg km / s2]")
grid on

nexttile
[theta1, theta2, theta3] = dcm2angle(angle2dcm(x(7,:), x(8,:), x(9,:),"XYX"), "ZYX"); % [rad]
plot(t, rad2deg(theta1), t, rad2deg(theta2), t, rad2deg(theta3));
legend("Yaw", "Pitch", "Roll")
title("Orientation vs Time")
xlabel("Time [s]")
ylabel("Orientation [deg]")
grid on

nexttile
plot(t, rad2deg(x(10:12, :)));
legend("\omega_" + string(1:3))
title("Angular Velocity vs Time")
xlabel("Time [s]")
ylabel("Angular Velocity [deg / s]")
grid on

nexttile
gimbal_angle = acosd(u(1, :) ./ u(4, :));
plot(t(1:Nu), gimbal_angle)
title("Gimbal Angle vs Time")
xlabel("Time [s]")
ylabel("Gimbal Angle [deg]")
grid on

sgtitle("3DoF Rocket Landing Time Histories")

end

