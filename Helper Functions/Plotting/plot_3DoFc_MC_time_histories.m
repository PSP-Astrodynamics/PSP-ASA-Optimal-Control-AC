function [] = plot_3DoFc_MC_time_histories(t_mean, x_mean, u_mean, t_fb, x_MC_fb, u_MC_fb, t_k, X_k, S_k, t_no_fb, x_MC_no_fb, T_max, T_min, max_gimbal, use_legend)
%PLOT_3DOF_MC_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

m = size(x_MC_fb, 3);

P_k = pagemtimes(X_k, pagetranspose(X_k));

%% States
titles = ["X Position", "Y Position", "X Velocity", "Y Velocity", "Orientation", "Angular Velocity", "Mass"];
ylabels = ["r_x [km]", "r_y [km]", "v_x [m / s]", "v_y [m / s]", "\theta [deg]", "\omega [deg / s]", "m [kg]"];

ops = {@(x) x, @(x) x, @(x) 1000 * x, @(x) 1000 * x, @(x) rad2deg(x), @(x) rad2deg(x), @(x) exp(x)};

for x = 1:7
    %x_3sigbound = 3 * sqrt(squeeze(project_ellipsoid(P_k, x)))';
    x_3sigbound = 3 * sqrt(squeeze(P_k(x, x, :)))';
    x_3sigbound_cont = interp1(t_k, x_3sigbound, t_mean);
    
    figure
    tiledlayout(1, 2)
    
    nexttile
    for i = 1:m
        plot(t_fb, ops{x}(x_MC_fb(x, :, i)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
    end
    plot(t_mean, ops{x}(x_mean(x, :)), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
    plot(t_mean, ops{x}(x_mean(x, :) + x_3sigbound_cont), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
    plot(t_mean, ops{x}(x_mean(x, :) - x_3sigbound_cont), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, HandleVisibility='off'); hold off
    title(titles(x) + " vs Time with Optimized Feedback Policies")
    xlabel("Time [s]")
    ylabel(ylabels(x))
    if use_legend
        legend("Location","best")
    end
    grid on
    
    nexttile
    for i = 1:m
        plot(t_no_fb, ops{x}(x_MC_no_fb(x, :, i)), Color = [192, 192, 192] / 256); hold on
    end
    plot(t_mean, ops{x}(x_mean(x, :)), Color = "k",LineWidth=1); hold off
    title(titles(x) + " vs Time without Trajectory Corrections")
    xlabel("Time [s]")
    ylabel(ylabels(x))
    grid on
end

%% Thrust Magnitude
figure
%Pu_k = pagemtimes(S_k, pagetranspose(S_k));
%thrust_3sigbound = 3 * sqrt(squeeze(project_ellipsoid(Pu_k, 3)))';
nu = 2;
epsilon = 1e-3; % 99.9% 
S_k_norm = zeros([size(u_MC_fb, 2), 1]);
for j = 1:size(u_MC_fb, 2)
    S_k_norm(j) = norm(S_k(:, (tri(j - 1) + 1):tri(j)));
end
thrust_3sigbound = squeeze(sigma_mag_confidence(epsilon / 2, nu) * S_k_norm);

for i = 1:m
    stairs(t_fb(1:size(u_MC_fb, 2)), u_MC_fb(3, :, i) .* exp(x_MC_fb(7, 1:size(u_MC_fb, 2))), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

u_mean_full = interp1(t_k(1:size(u_mean, 2)), u_mean', t_k, "previous", "extrap")';
thrust_3sigbound_full = interp1(t_k(1:size(u_mean, 2)), thrust_3sigbound, t_k, "previous", "extrap");

stairs(t_k, vecnorm(u_mean_full(1:2, :),2,1) .* exp(x_mean(5, 1:size(u_mean_full, 2))), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(t_k, (vecnorm(u_mean_full(1:2, :),2,1) + thrust_3sigbound_full) .* exp(x_mean(7, 1:size(u_mean_full, 2))), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(t_k, (vecnorm(u_mean_full(1:2, :),2,1) - thrust_3sigbound_full) .* exp(x_mean(7, 1:size(u_mean_full, 2))), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1,HandleVisibility='off'); hold on
yline(T_max, LineWidth = 1, LineStyle="--", Color="k", DisplayName = "Constraint"); hold on
yline(T_min, LineWidth = 1, LineStyle="--", Color="k", HandleVisibility='off'); hold off
title("Thrust Magnitude vs Time with Optimized Feedback Policies")
legend(Location="southeast")
xlabel("Time [s]")
ylabel("Thrust [kN]")
grid on

%% Gimbal Angle
figure

epsilon = 1e-3; % 99.9% 
% Best guess at how to evaluate 99.9% bound... how to turn into actual
% angle bound??
u1_3sigbound = norminv(1 - epsilon / 2) * S_k_norm;
u2_3sigbound = sigma_mag_confidence(epsilon / 2, 1) * S_k_norm;

for i = 1:m
    stairs(t_fb(1:size(u_MC_fb, 2)), atan2d(u_MC_fb(2, :, i), u_MC_fb(1, :, i)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

u_mean_full = interp1(t_k(1:size(u_mean, 2)), u_mean', t_k, "previous", "extrap")';
u1_3sigbound_full = interp1(t_k(1:size(u_mean, 2)), u1_3sigbound, t_k, "previous", "extrap");
u2_3sigbound_full = interp1(t_k(1:size(u_mean, 2)), u2_3sigbound, t_k, "previous", "extrap");
gimbal_99p9bound_up = atan2d(u_mean_full(2,:) + u2_3sigbound_full .* sign(u_mean_full(2,:)), u_mean_full(1, :) - u1_3sigbound_full);
gimbal_99p9bound_dn = atan2d(u_mean_full(2,:) - u2_3sigbound_full .* sign(u_mean_full(2,:)), u_mean_full(1, :) + u1_3sigbound_full);
gimbal_99p9bound2 = acosd((u_mean_full(1, :) - u1_3sigbound_full) ./ (u_mean_full(3, :)));

gimbal_99p9bound = rad2deg((-(u_mean_full(2,:) .* u1_3sigbound_full) + (u_mean_full(1, :) .* u2_3sigbound_full)) ./ (u_mean_full(1, :) .^ 2 + u_mean_full(2, :) .^ 2));

ck = u_mean_full(3, :) * cos(max_gimbal) - (u_mean_full(1, :) - u1_3sigbound_full);


stairs(t_k, atan2d(u_mean_full(2, :), u_mean_full(1, :)), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(t_k, gimbal_99p9bound_up, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(t_k, gimbal_99p9bound_dn, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, HandleVisibility='off'); hold on
yline(rad2deg(max_gimbal), LineWidth = 1, LineStyle="--", Color="k", DisplayName = "Constraint"); hold on
yline(-rad2deg(max_gimbal), LineWidth = 1, LineStyle="--", Color="k", HandleVisibility='off'); hold off
title("Gimbal Angle vs Time with Optimized Feedback Policies")
legend(Location="southeast")
xlabel("Time [s]")
ylabel("Gimbal Angle [deg]")
grid on

 
end

