function [] = comparison_plot_6DoF_trajectory(x, names, glideslope_angle, options)
arguments
    x 
    names
    glideslope_angle
    options.step = 4
    options.title = "State Trajectory"
    options.colorpallete = "sail"
    options.linestyle = []
end
%COMPARISON_PLOT_3DOF_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

Nc = numel(x);

if isempty(options.linestyle)
    options.linestyle = strings([1, Nc]);
    options.linestyle(:) = "-";
end

x_cat = [x{:}];

z_max = max(x_cat(3, :));
x_lim = [min(x_cat(1, :)) - z_max / 10, max(x_cat(1, :)) + z_max / 10];
y_lim = [min(x_cat(2, :)) - z_max / 10, max(x_cat(2, :)) + z_max / 10];
z_lim = [min(x_cat(3, :)) - z_max / 10, z_max * 1.1];

for c = 1:Nc
    plot3(x{c}(1, :), x{c}(2, :), x{c}(3, :), LineWidth = 2, LineStyle = options.linestyle(c)); hold on
end
hold off
%line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle), 'Color', 'k', 'LineStyle', '--'); hold off
title(options.title)
legend([names, "Glideslope"], Location="eastoutside", Orientation="vertical")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
axis equal
grid on
xlim(x_lim)
ylim(y_lim)
zlim(z_lim)

colororder(options.colorpallete)

end

