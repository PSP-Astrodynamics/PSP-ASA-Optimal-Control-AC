function [] = plot_3DoFc_trajectory(t, x, u, glideslope_angle, gimbal_max, T_min, T_max, options)
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
%PLOT_3DOF_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

Nu = size(u, 2);

tiledlayout(1, 2, "TileSpacing","compact")

%% Trajectory Plot
nexttile
y_max = max(x(2, :));
x_lim = [min(x(1, :)) - y_max / 10, max(x(1, :)) + y_max / 10];
y_lim = [min(x(2, :)) - y_max / 10, y_max * 1.1];

C_BN_3D = make_R(x(5, :), 3);
C_BN = C_BN_3D(1:2, 1:2, :);
C_NB = pagetranspose(C_BN);
u_N = squeeze(pagemtimes(C_NB(:, :, 1:Nu), permute(reshape(u(1:2, :), [2, Nu, 1]), [1, 3, 2])));

plot(x(1, :), x(2, :)); hold on
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle), 'Color', 'k', 'LineStyle', '--'); hold on
quiver(x(1, 1:options.step:Nu), x(2, 1:options.step:Nu), -u_N(1, 1:options.step:end), -u_N(2, 1:options.step:end), ShowArrowHead = "off", color = "r", AutoScaleFactor=0.4)
title(options.title)
legend("", "Glideslope", "Thrust", Location="southoutside", Orientation="horizontal")
xlabel("X [km]")
ylabel("Y [km]")
axis equal
grid on
xlim(x_lim)
ylim(y_lim)

%% Thrust Plot
nexttile

cmap = colormap;
u_color = interp1(linspace(1, Nu, size(cmap, 1)), cmap, 1:Nu);

gimbal_angle = atan2(u(2, :), u(1, :));

for i = 1:(Nu - 1)
    polarplot(gimbal_angle(i:(i + 1)), u(3, i:(i + 1)) / T_max .* exp(x(7, i:(i + 1))), Color = u_color(i, :), LineWidth=0.5, Marker="o", MarkerFaceColor=u_color(i, :)); hold on
end
polarregion([-gimbal_max, gimbal_max], [T_min / T_max, 1])
hold off

pax = gca;
pax.ThetaZeroLocation = "bottom";

cb = colorbar(pax, "southoutside",TickLabels=[t(1), round(t(end) / 2), t(end)]);
xlabel(cb, "Time [s]")

rmax = 1.1;
rlim([0 rmax]);

thetalim(rad2deg(gimbal_max) * [-2, 2])
title(sprintf("Control Trajectory\n"))

sgtitle(sprintf("3DoF Rocket Landing Trajectory"))

end

