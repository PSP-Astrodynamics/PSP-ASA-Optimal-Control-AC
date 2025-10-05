function [] = comparison_plot_3DoF_trajectory(x, names, glideslope_angle, options)
arguments
    x 
    names
    glideslope_angle
    options.step = 4
    options.title = "State Trajectory"
    options.colorpallete = "sail"
    options.linestyle = []
    options.STC_glideslope_angle_height = []
end
%COMPARISON_PLOT_3DOF_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

Nc = numel(x);

if isempty(options.linestyle)
    options.linestyle = strings([1, Nc]);
    options.linestyle(:) = "-";
end

x_cat = [x{:}];

y_max = max(x_cat(2, :));
x_lim = [min(x_cat(1, :)) - y_max / 10, max(x_cat(1, :)) + y_max / 10];
y_lim = [min(x_cat(2, :)) - y_max / 10, y_max * 1.1];

for c = 1:Nc
    plot(x{c}(1, :), x{c}(2, :), LineWidth = 2, LineStyle = options.linestyle(c)); hold on
end
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle), 'Color', 'k', 'LineStyle', '--'); hold on
if ~isempty(options.STC_glideslope_angle_height)
    line([-options.STC_glideslope_angle_height(2), 0, options.STC_glideslope_angle_height(2)] * tan(options.STC_glideslope_angle_height(1)), [options.STC_glideslope_angle_height(2), 0, options.STC_glideslope_angle_height(2)], Color = 'k', LineStyle = '--');
    yline(options.STC_glideslope_angle_height(2), LineStyle="--")
end
hold off
title(options.title)
legend([names, "Glideslope", "STC Glideslope", "STC Glideslope Trigger"], Location="eastoutside", Orientation="vertical")
xlabel("X [km]")
ylabel("Y [km]")
axis equal
grid on
xlim(x_lim)
ylim(y_lim)

colororder(options.colorpallete)

end

