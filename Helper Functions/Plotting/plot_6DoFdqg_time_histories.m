function [] = plot_6DoFdqg_time_histories(t, x, u)
%PLOT_3DOF_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

Nu = size(u, 2);

q = x(2:5, :);
r_N = 2 * eye(3, 4) * q_mul_array(x(6:9, :), q_conj(q));

tiledlayout(4, 2)

nexttile
plot(t, r_N);
legend("r_" + string(1:3))
title("Position vs Time")
xlabel("Time [s]")
ylabel("Position [km]")
grid on

v_N = quat_rot_array(q_conj(q), x(13:15, :));

nexttile
plot(t, v_N);
legend("v_" + string(1:3))
title("Velocity vs Time")
xlabel("Time [s]")
ylabel("Velocity [m / s]")
grid on


T = u(1, :);
gimbal_angle = u(2, :);
clock_angle = u(3, :);

u_B = [T .* sin(gimbal_angle) .* cos(clock_angle); T .* sin(gimbal_angle) .* sin(clock_angle); T .* cos(gimbal_angle)];

u_N = quat_rot_array(q_conj(q), u_B);


nexttile
plot(t(1:Nu), u_N(1:3, :)); hold on
plot(t(1:Nu), T, LineStyle = "--"); hold off
legend(["T_" + string(1:3), "||T||"])
title("Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [kg m / s2]")
grid on


nexttile
q_mag = vecnorm(q, 2);
[theta1, theta2, theta3] = dcm2angle(quat2rot(q), "ZYX"); % [rad]
plot(t, rad2deg(theta1), t, rad2deg(theta2), t, rad2deg(theta3), t, q_mag * 100);
legend("Yaw", "Pitch", "Roll", "100*||q||")
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
plot(t(1:Nu), rad2deg(gimbal_angle)); hold on
plot(t(1:Nu), clock_angle); hold off
legend("Gimbal Angle", "Clock Angle")
title("Gimbal and Clock Angles vs Time")
xlabel("Time [s]")
ylabel("Gimbal Angle [deg], Clock Angle [rad]")
grid on

nexttile
plot(t(1:Nu), u(4:6, :)); hold on
legend(["\tau_" + string(1:3)])
title("Moment vs Time")
xlabel("Time [s]")
ylabel("Moment [kg m2 / s2]")
grid on

nexttile
plot(t, x(1, :)); hold off
title("Mass vs Time")
xlabel("Time [s]")
ylabel("Mass [kg]")
grid on

sgtitle("6DoF Rocket Landing Time Histories")

end

function R = quat2rot(Q)
q_0 = Q(4, :);
q_1 = Q(1, :);
q_2 = Q(2, :);
q_3 = Q(3, :);
R(1,1,:) = q_0.*q_0 + q_1.*q_1 - q_2.*q_2 - q_3.*q_3;
R(1,2,:) = 2.*(q_1.*q_2 - q_0.*q_3);
R(1,3,:) = 2.*(q_0.*q_2 + q_1.*q_3);
R(2,1,:) = 2.*(q_0.*q_3 + q_1.*q_2);
R(2,2,:) = q_0.*q_0 - q_1.*q_1 + q_2.*q_2 - q_3.*q_3;
R(2,3,:) = 2.*(q_2.*q_3 - q_0.*q_1);
R(3,1,:) = 2.*(q_1.*q_3 - q_0.*q_2);
R(3,2,:) = 2.*(q_0.*q_1 + q_2.*q_3);
R(3,3,:) = q_0.*q_0 - q_1.*q_1 - q_2.*q_2 + q_3.*q_3;
end