function [] = comparison_plot_6DoF_time_histories(t, x, u, names, options)
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
    plot(ax, t{c}(1:numel(u{c}(4, :))), u{c}(4, :), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
end
legend(ax, names)
title("Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [kg km / s2]")
grid on

colororder(ax, options.colorpallete)

ax = nexttile;
for c = 1:Nc
    plot(ax, t{c}, wrapTo180(rad2deg(x{c}(7:9, :))), LineWidth = 2, LineStyle = options.linestyle(c), DisplayName=names(c)); hold on
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
    gimbal_angle = acosd(u{c}(1, :) ./ u{c}(4, :));
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

