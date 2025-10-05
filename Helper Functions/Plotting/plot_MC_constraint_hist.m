function [] = plot_MC_constraint_hist(constraint_func, i, x_MC_fb, x_MC_no_fb, constraint_val, name, label, N)
%PLOT_MC_CONSTRAINT_HIST Summary of this function goes here
%   Detailed explanation goes here

tiledlayout(1,2,"TileSpacing","compact")

nexttile
histogram(constraint_func(x_MC_fb, i(1)), N);
xline(constraint_val, LineWidth=1, LineStyle="--");
xlabel(label)
title("With Optimized Feedback Policy")
grid on

nexttile
histogram(constraint_func(x_MC_no_fb, i(2)), N);
xline(constraint_val, LineWidth=1, LineStyle="--");
xlabel(label)
title("Without Trajectory Corrections")
grid on

sgtitle(name + " Constraint Satisfaction Histograms")

end

