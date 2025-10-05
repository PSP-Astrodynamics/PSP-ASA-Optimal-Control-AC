function [] = comparison_plot_6DoFdqg_time_histories(t, x, u, names, options)
arguments
    t
    x
    u
    names
    options.title = "6DoF Rocket Landing Time Histories"
    options.colorpallete = "sail"
    options.linestyle = []
end
%PLOT_3DOF_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

Nc = numel(x);

tiledlayout(3, 2)

ax = nexttile;
newcolors = repmat([1, 0, 0; 0, 1, 0; 0, 0, 1], Nc, 1);
         
legend_labels = [];
for c = 1:Nc
    q = x{c}(2:5, :);
    r_N = 2 * eye(3, 4) * q_mul_array(x{c}(6:9, :), q_conj(q));

    plot(ax, t{c}, r_N(1:3, :), LineWidth = 2, LineStyle = options.linestyle(c)); hold on
    legend_labels = [legend_labels, "$r_{" + string(1:3) + "_{" + names(c) + "}}$"];
end
hold off
legend(ax, legend_labels, Interpreter="latex")
title("Position vs Time")
xlabel("Time [s]")
ylabel("Position [m]")
grid on

colororder(ax, newcolors)

ax = nexttile;
         
legend_labels = [];
for c = 1:Nc
    q = x{c}(2:5, :);
    v_N = quat_rot_array(q_conj(q), x{c}(13:15, :));

    plot(ax, t{c}, v_N(1:3, :), LineWidth = 2, LineStyle = options.linestyle(c)); hold on
    legend_labels = [legend_labels, "$v_{" + string(1:3) + "_{" + names(c) + "}}$"];
end
legend(ax, legend_labels, Interpreter="latex")
title("Velocity vs Time")
xlabel("Time [s]")
ylabel("Velocity [m / s]")
grid on

colororder(ax, newcolors)

ax = nexttile;

for c = 1:Nc
    plot(ax, t{c}(1:numel(u{c}(1, :))), u{c}(1, :), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
end
legend(ax, names)
title("Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [kg m / s2]")
grid on

colororder(ax, options.colorpallete)

ax = nexttile;
for c = 1:Nc
    q = x{c}(2:5, :);
    [theta1, theta2, theta3] = dcm2angle(quat2rot(q), "ZYX"); % [rad]

    plot(ax, t{c}, wrapTo180(rad2deg([theta1, theta2, theta3])), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
    legend_labels = [legend_labels, "$\theta_{" + string(1:3) + "_{" + names(c) + "}}$"];
end
hold off
legend(ax, legend_labels, Interpreter="latex")
title("Orientation vs Time")
xlabel("Time [s]")
ylabel("Orientation [deg]")
grid on

colororder(ax, newcolors)

ax = nexttile;
for c = 1:Nc
    plot(ax, t{c}, rad2deg(x{c}(10:12, :)), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
    legend_labels = [legend_labels, "$\omega_{" + string(1:3) + "_{" + names(c) + "}}$"];
end
hold off
legend(ax, legend_labels, Interpreter="latex")
title("Angular Velocity vs Time")
xlabel("Time [s]")
ylabel("Angular Velocity [deg / s]")
grid on

colororder(ax, newcolors)

ax = nexttile;
for c = 1:Nc
    gimbal_angle = rad2deg(u{c}(2, :));
    plot(ax, t{c}(1:numel(gimbal_angle)), gimbal_angle, LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
end
legend(ax, names)
title("Gimbal Angle vs Time")
xlabel("Time [s]")
ylabel("Gimbal Angle [deg]")
grid on

colororder(ax, options.colorpallete)

sgtitle(options.title)

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