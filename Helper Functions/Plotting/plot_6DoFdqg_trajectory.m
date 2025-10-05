function [] = plot_6DoFdqg_trajectory(t, x, u, glideslope_angle, gimbal_max, T_min, T_max, options)
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
q = x(2:5, :);
r_N = 2 * eye(3, 4) * q_mul_array(x(6:9, :), q_conj(q));

nexttile
z_max = max(r_N(3, :));
x_lim = [min(r_N(1, :)) - z_max / 10, max(r_N(1, :)) + z_max / 10];
y_lim = [min(r_N(2, :)) - z_max / 10, max(r_N(2, :)) + z_max / 10];
z_lim = [min(r_N(3, :)) - z_max / 10, z_max * 1.1];

T = u(1, :);
gimbal_angle = u(2, :);
gimbal_clock = u(3, :);

u_B = [T .* sin(gimbal_angle) .* cos(gimbal_clock); T .* sin(gimbal_angle) .* sin(gimbal_clock); T .* cos(gimbal_angle)];

u_N = quat_rot_array(q_conj(q), u_B);

plot3(r_N(1, :), r_N(2, :), r_N(3, :)); hold on
%line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle), 'Color', 'k', 'LineStyle', '--'); hold on
quiver3(r_N(1, 1:options.step:Nu), r_N(2, 1:options.step:Nu), r_N(3, 1:options.step:Nu), -u_N(1, 1:options.step:end), -u_N(2, 1:options.step:end), -u_N(3, 1:options.step:end), ShowArrowHead = "off", color = "r", AutoScaleFactor=0.4)
title(options.title)
%legend("", "Glideslope", "Thrust", Location="southoutside", Orientation="horizontal")
xlabel("X [m]")
ylabel("Y [m]")
xlabel("Z [m]")
axis equal
grid on
xlim(x_lim)
ylim(y_lim)
zlim(z_lim)

%% Thrust Plot
nexttile

cmap = colormap;
u_color = interp1(linspace(T_min / T_max, 1, size(cmap, 1)), cmap, T / T_max, "linear", "extrap");

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

