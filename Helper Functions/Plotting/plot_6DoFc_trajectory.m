function [] = plot_6DoFc_trajectory(t, x, u, glideslope_angle, gimbal_max, T_min, T_max, options)
arguments
    t
    x 
    u 
    glideslope_angle
    gimbal_max
    T_min
    T_max
    options.step = 4
    options.title = "State Trajectory"
end
%PLOT_6DOF_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

Nu = size(u, 2);

tiledlayout(1, 2, "TileSpacing","compact")

%% Trajectory Plot
nexttile
z_max = max(x(3, :));
x_lim = [min(x(1, :)) - z_max / 10, max(x(1, :)) + z_max / 10];
y_lim = [min(x(2, :)) - z_max / 10, max(x(2, :)) + z_max / 10];
z_lim = [min(x(3, :)) - z_max / 10, z_max * 1.1];

C_BN = angle2dcm(x(7, :), x(8, :), x(9, :), "XYX");

C_NB = pagetranspose(C_BN);
u_N = squeeze(pagemtimes(C_NB(:, :, 1:Nu), permute(reshape(u(1:3, :), [3, Nu, 1]), [1, 3, 2])));

plot3(x(1, :), x(2, :), x(3, :)); hold on
%line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle), 'Color', 'k', 'LineStyle', '--'); hold on
quiver3(x(1, 1:options.step:Nu), x(2, 1:options.step:Nu), x(3, 1:options.step:Nu), -u_N(1, 1:options.step:end), -u_N(2, 1:options.step:end), -u_N(3, 1:options.step:end), ShowArrowHead = "off", color = "r", AutoScaleFactor=0.4)
title(options.title)
%legend("", "Glideslope", "Thrust", Location="southoutside", Orientation="horizontal")
xlabel("X [km]")
ylabel("Y [km]")
xlabel("Z [km]")
axis equal
grid on
xlim(x_lim)
ylim(y_lim)
zlim(z_lim)

%% Thrust Plot
nexttile

cmap = colormap;
u_color = interp1(linspace(T_min / T_max*0, 1, size(cmap, 1)), cmap, u(4, :) .* exp(x(13, 1:Nu)) / T_max, "linear", "extrap");

gimbal_angle = acos(u(1, :) ./ u(4, :));
gimbal_clock = atan2(u(3, :), u(2, :));

for i = 1:(Nu - 1)
    polarplot(gimbal_clock(i:(i + 1)), rad2deg(gimbal_angle(i:(i + 1))), Color = u_color(i, :), LineWidth=1, Marker="o", MarkerFaceColor=u_color(i, :)); hold on
end
polarregion([-pi, pi], [0, rad2deg(gimbal_max)])
hold off

pax = gca;
pax.ThetaZeroLocation = "bottom";

cb = colorbar(pax, "southoutside",TickLabels=round(linspace(T_min / T_max, 1, 11) * 100, 2));
xlabel(cb, "Thrust [%]")

rmax = rad2deg(1.3 * gimbal_max);
rlim([0 rmax]);

title(sprintf("Control Trajectory\n"))

sgtitle(sprintf("6DoF Rocket Landing Trajectory"))

end

