function [] = plot_6DoF_MC_trajectories(t_mean, x_mean, t_fb, x_MC_fb, x_ref, t_k, P_k, t_no_fb, x_MC_no_fb, Pf, glideslope_angle, h_glideslope)
%PLOT_3DOF_MC_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here

m = size(x_MC_fb, 3);

figure
tiledlayout(2, 2, "TileSpacing","compact")

%% Feedback controlled

nexttile

[P_eigvecs, P_eigvals] = pageeig(P_k([1,3], [1,3], :)); % Looks more correct then projecting ellipsoid...

thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, numel(t_k)]);
for k = 1:numel(t_k)
    ellipse_3sigma(:, :, k) = x_ref([1,3], k) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
end

plot(squeeze(x_MC_fb(1, :, :)), squeeze(x_MC_fb(3, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(x_mean(1, :), x_mean(3, :), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
x_lim = xlim;
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle) - h_glideslope, 'Color', 'k', 'LineStyle', '--', "DisplayName", "Glideslope"); hold on
plot(squeeze(ellipse_3sigma(1, :, 1)), squeeze(ellipse_3sigma(2, :, 1)), Color = "k", DisplayName="Covariance"); hold on
plot(squeeze(ellipse_3sigma(1, :, 2:end)), squeeze(ellipse_3sigma(2, :, 2:end)), Color = "k", HandleVisibility='off'); hold off
title("XZ - With Optimized Feedback Policy")
xlabel("X [km]")
ylabel("Z [km]")
legend(location = "best")
axis equal
grid on


nexttile

[P_eigvecs, P_eigvals] = pageeig(P_k(2:3, 2:3, :)); % Looks more correct then projecting ellipsoid...

thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, numel(t_k)]);
for k = 1:numel(t_k)
    ellipse_3sigma(:, :, k) = x_ref(2:3, k) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
end

plot(squeeze(x_MC_fb(2, :, :)), squeeze(x_MC_fb(3, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(x_mean(2, :), x_mean(3, :), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
x_lim = xlim;
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle) - h_glideslope, 'Color', 'k', 'LineStyle', '--', "DisplayName", "Glideslope"); hold on
plot(squeeze(ellipse_3sigma(1, :, 1)), squeeze(ellipse_3sigma(2, :, 1)), Color = "k", DisplayName="Covariance"); hold on
plot(squeeze(ellipse_3sigma(1, :, 2:end)), squeeze(ellipse_3sigma(2, :, 2:end)), Color = "k", HandleVisibility='off'); hold off
title("YZ - With Optimized Feedback Policy")
xlabel("Y [km]")
ylabel("Z [km]")
axis equal
grid on



%% No FB
nexttile
plot(squeeze(x_MC_no_fb(1, :, :)), squeeze(x_MC_no_fb(3, :, :)), Color = [192, 192, 192] / 256); hold on
plot(x_mean(1, :), x_mean(3, :), Color = [30, 144, 255] / 256, LineWidth=1); hold on
x_lim = xlim;
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle) - h_glideslope, 'Color', 'k', 'LineStyle', '--', "DisplayName", "Glideslope"); hold off
title("XZ - Without Trajectory Corrections")
xlabel("X [km]")
ylabel("Z [km]")
grid on
axis equal

nexttile
plot(squeeze(x_MC_no_fb(2, :, :)), squeeze(x_MC_no_fb(3, :, :)), Color = [192, 192, 192] / 256); hold on
plot(x_mean(2, :), x_mean(3, :), Color = [30, 144, 255] / 256, LineWidth=1); hold on
x_lim = xlim;
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle) - h_glideslope, 'Color', 'k', 'LineStyle', '--', "DisplayName", "Glideslope"); hold off
title("YZ - Without Trajectory Corrections")
xlabel("Y [km]")
ylabel("Z [km]")
grid on
axis equal

sgtitle("Monte Carlo Simulation of 6DoF Rocket Landing with and without Optimized Feedback Policy")

% 
% % Plot zoomed in terminal section
% proj_Pf_r = project_ellipsoid(Pf, [1,2]);
% [Pf_eigvecs, Pf_eigvals] = eigs(proj_Pf_r);
% terminal_ellipse = x_ref(1:2, end) + Pf_eigvecs * [3 * sqrt(Pf_eigvals(1, 1)) * cos(thetas); 3 * sqrt(Pf_eigvals(2, 2)) * sin(thetas)];
% 
% axes('Position',[.3 .16 .09 .14])
% box on
% plot(squeeze(x_MC_fb(1, :, :)), squeeze(x_MC_fb(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
% plot(x_mean(1, :), x_mean(2, :), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
% plot(terminal_ellipse(1, :), terminal_ellipse(2, :), Color = "k", HandleVisibility='off', LineWidth=1, LineStyle="-."); hold on
% plot(squeeze(ellipse_3sigma(1, :, end)), squeeze(ellipse_3sigma(2, :, end)), Color = "k", DisplayName="Solution", LineWidth=0.7); hold on
% xlim(2 * [-max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2))), max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2)))])
% ylim(2 * [-max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2))), max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2)))])
% x_lim = xlim;
% line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glideslope_angle) - h_glideslope, 'Color', 'k', 'LineStyle', '--', "DisplayName", "Glideslope"); hold off
% grid on
% end

