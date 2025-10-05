function [] = comparison_plot_6DoFq_time_histories(t, x, u, names, options)
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
    plot(ax, t{c}, x{c}(1:3, :), LineWidth = 2, LineStyle = options.linestyle(c)); hold on
    legend_labels = [legend_labels, "$r_{" + string(1:3) + "_{" + names(c) + "}}$"];
end
hold off
legend(ax, legend_labels, Interpreter="latex")
title("Position vs Time")
xlabel("Time [s]")
ylabel("Position [km]")
grid on

colororder(ax, newcolors)

ax = nexttile;
         
legend_labels = [];
for c = 1:Nc
    plot(ax, t{c}, x{c}(4:6, :) * 1e3, LineWidth = 2, LineStyle = options.linestyle(c)); hold on
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
    plot(ax, t{c}(1:numel(u{c}(3, :))), vecnorm(u{c}(1:3, :)), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
end
legend(ax, names)
title("Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [kg km / s2]")
grid on

colororder(ax, options.colorpallete)

legend_labels = [];
ax = nexttile;
for c = 1:Nc
    theta_1 = zeros([numel(t{c}), 1]);
    theta_2 = zeros([numel(t{c}), 1]);
    theta_3 = zeros([numel(t{c}), 1]);
    for n = 1 : numel(t{c})
        R = quat_rotmatrix(x{c}(7:10, n));
        [theta_1(n), theta_2(n), theta_3(n)] = dcm2angle(R, "YXZ");
    end
    plot(ax, t{c}, rad2deg(theta_1), t{c}, rad2deg(theta_2), t{c}, rad2deg(theta_3), LineWidth = 2, LineStyle = options.linestyle(c)); hold on
    legend_labels = [legend_labels, "${" + ["Pitch", "Yaw", "Roll"] + "}_{" + names(c) + "}$"];
end
hold off
legend(ax, legend_labels, Interpreter="latex")
title("Orientation vs Time")
xlabel("Time [s]")
ylabel("Orientation [deg]")
grid on

colororder(ax, newcolors)

legend_labels = [];

ax = nexttile;
for c = 1:Nc
    plot(ax, t{c}, rad2deg(x{c}(11:13, :)), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
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
    gimbal_angle = acosd(u{c}(1, :) ./ vecnorm(u{c}(1:3, :)));
    plot(ax, t{c}(1:numel(gimbal_angle)), gimbal_angle, LineWidth = 2, LineStyle = options.linestyle(c), DisplayName="gimbal" + "_{" + names(c) + "}"); hold on
    plot(ax, t{c}(1:numel(gimbal_angle)), rad2deg(u{c}(4, :)), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName="vane" + "_{" + names(c) + "}"); 
end
legend(ax)
title("Gimbal Angle and Vane Angle vs Time")
xlabel("Time [s]")
ylabel("Gimbal Angle and Vane Angle [deg]")
grid on

colororder(ax, options.colorpallete)

sgtitle(options.title)

end


function [R] = quat_rotmatrix(q)
    w = q(4);
    v = q(1:3);

    R = (w ^ 2 - v' * v) * eye(3) + 2 * v * v' + 2 * w * skew(v);
end