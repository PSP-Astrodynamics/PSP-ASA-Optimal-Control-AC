%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Reference for MPC
ref = load("reference.mat").ref_sol;

N_horizon = 10;
t_horizon = linspace(4, 9, N_horizon); % MPC time horizon

x_ref_sol = interp1(ref.t_ref, ref.x_ref_sol', t_horizon)';
u_ref_sol = interp1(ref.t_ref, ref.u_ref_sol', t_horizon)';

% Set cost matrices for MPC
r_cost = 100 * ones([3, 1]);
v_cost = 100 * ones([3, 1]);
q_cost = 10 * ones([4, 1]);
omega_cost = 10 * ones([3, 1]);
m_cost = 1 * ones([1, 1]);
Q = diag([r_cost; v_cost; q_cost; omega_cost; m_cost].^2);

thrust_cost = 0.1 * ones([3, 1]);
torque_cost = 0.1 * ones([1, 1]);
R = diag([thrust_cost; torque_cost] .^ 2);

% Define perturbation to initial state
r_pert = [0.005; 0; 0];
v_pert = [0; 0; 0];
q_pert = [0; 0; 0; 0];
omega_pert = [0; 0; 0];
mass_pert = 0;

perturbation = [r_pert; v_pert; q_pert; omega_pert; mass_pert];

%% Initialize
% Vehicle Parameters
Isp = 225; % [s]
alpha = 1 / (9.81e-3 * Isp); % [s / km]
T_min = 6000e-3; % [kg km / s2]
T_max = 22500e-3; % [kg km / s2]
I = [19150; 13600; 13600] * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 0.25e-3; % [km] Distance from CoM to nozzle
m_dry = 2100; % [kg]
m_wet = 1150; % [kg]
m_0 = m_dry + m_wet;
gimbal_max = deg2rad(20); % [rad]
g = 1.62e-3;
time_min_max_thrust = 3; % [s] time to throttle from min to max thrust
max_gimbal_rate = 10; % [deg / s] max rate of change of gimbal angle

vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I);

% Problem Parameters
tf = t_horizon(end) - t_horizon(1); % [s]
N = N_horizon; % []
r_0 = [250; -100; 433] * 1e-3; % [km]
v_0 = [0; 0; -35] * 1e-3; % [km / s]
theta_0 = [deg2rad(0); deg2rad(90); deg2rad(0)]; % [rad]
initial_roll = deg2rad(0);
R_0 = make_R(initial_roll, 3) * angle2dcm(theta_0(1), theta_0(2), theta_0(3));
q_0 = qexp(RLog(R_0));
w_0 = deg2rad([0; 0; 0]); % [rad / s]
glideslope_angle_max = deg2rad(65); % [rad]

theta_f = [0; deg2rad(90); 0]; % [rad]
R_f = angle2dcm(theta_f(1), theta_f(2), theta_f(3));
q_f = qexp(RLog(R_f));

x_0 = x_ref_sol(:, 1) + perturbation;
x_f = [[0; 0; 30] * 1e-3; [0; 0; -1] * 1e-3; q_f; zeros(3, 1)];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

nx = 14; 
nu = 4;
np = 0;

% PTR algorithm parameters
ptr_ops.iter_max = 5;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 5e-4;
ptr_ops.w_vc = 1e3;
ptr_ops.w_tr = ones(1, Nu) * 4e-3;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 2e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;
ptr_ops.include_v_0 = false;
ptr_ops.include_u_0 = true;

scale = true;

scale_hint.x_max = [max(r_0) * ones([3, 1]); max(v_0) * ones([3, 1]); ones([4, 1]); max(w_0) * ones([3, 1]); m_0];
scale_hint.x_min = [-max(r_0) * ones([3, 1]); -max(v_0) * ones([3, 1]); -ones([4, 1]); -max(w_0) * ones([3, 1]); m_dry];
scale_hint.u_max = [T_max; sin(gimbal_max) * ones([2, 1]); pi / 4];
scale_hint.u_min = [T_min; -sin(gimbal_max) * ones([2, 1]); -pi / 4];
scale_hint.p_max = 60;
scale_hint.p_min = 20;

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsQuat6DoF_localrot_noumag(x, u, L, I, alpha, g);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = {1:N, @(t, x, u, p) norm(x(1:3)) - x(3) / cos(glideslope_angle_max)};
mass_constraint = {1:N, @(t, x, u, p) m_dry - x(14)};
angular_velocity_constraint = {1:N, @(t, x, u, p) norm(x(11:13), Inf) - norm([w_0; deg2rad(20)], Inf)};

state_convex_constraints = {glideslope_constraint, mass_constraint, angular_velocity_constraint};

% Convex control constraints
max_thrust_constraint = {1:N, @(t, x, u, p) norm(u(1:3)) - T_max};
%min_thrust_constraint = {1:N, @(t, x, u, p) T_min - u(4)};
max_gimbal_constraint = {1:N, @(t, x, u, p) norm(u(1:3)) - u(1) / cos(gimbal_max)};
%lcvx_thrust_constraint = {1:N, @(t, x, u, p) norm(u(1:3)) - u(4)}; 
max_vane_angle_constraint = {1:N, @(t, x, u, p) abs(u(4)) - deg2rad(10)};
control_convex_constraints = {max_gimbal_constraint,max_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
pitch_constraint = @(t, x, u, p) 1 - (2 * (x(7, :) .* x(9, :) - x(8, :) .* x(10, :))) .^ 2 - sin(deg2rad(45)) ^ 2;
pitch_constraint_linearized = {1:N, linearize_constraint(pitch_constraint, nx, nu, np, "x", 7:10)};
pitch_func = @(t, x, u, p) 2 * (x(7, :) .* x(9, :) - x(8, :) .* x(10, :));
pitch_func_linearized = linearize_constraint(pitch_func, nx, nu, np, "x", 7:10);
state_nonconvex_constraints = {pitch_constraint_linearized};

% Nonconvex control constraints
min_thrust_constraint = @(t, x, u, p) T_min ^ 2 - sum_square(u(1:3));
min_thrust_constraint_linearized = {1:N, linearize_constraint(min_thrust_constraint, nx, nu, np, "u", 1:3)};
max_thrust_rate_constraint = {1:(N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) abs((norm(u_ref(1:3, k + 1)) + u_ref(1:3, k + 1)' / norm(u_ref(1:3, k + 1)) * (u(1:3, k + 1) - u_ref(1:3, k + 1))) - (norm(u_ref(1:3, k)) + u_ref(1:3, k)' / norm(u_ref(1:3, k)) * (u(1:3, k) - u_ref(1:3, k)))) / delta_t - (T_max - T_min) / time_min_max_thrust};
max_gimbal_rate_constraint = {1:(N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) -(u_ref(1:3, k)' * u_ref(1:3, k + 1) + u_ref(1:3, k)' * (u(1:3, k + 1) - u_ref(1:3, k + 1)) + u_ref(1:3, k + 1)' * (u(1:3, k) - u_ref(1:3, k))) + (norm(u_ref(1:3, k + 1)) * norm(u_ref(1:3, k))  + norm(u_ref(1:3, k)) * u(1:3, k + 1)' / norm(u_ref(1:3, k + 1)) * (u(1:3, k + 1) - u_ref(1:3, k + 1)) + norm(u_ref(1:3, k + 1)) * u_ref(1:3, k)' / norm(u_ref(1:3, k)) * (u(1:3, k) - u_ref(1:3, k))) * cosd(max_gimbal_rate * delta_t)};
control_nonconvex_constraints = {min_thrust_constraint_linearized, max_thrust_rate_constraint, max_gimbal_rate_constraint};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

% Terminal boundary conditions
terminal_bc = @(x, p, x_ref, p_ref) zeros([nx, 1]);

%% Specify Objective
min_fuel_objective = @(x, u, p, k) einsum(@(k) quad_form(u(4, k) - u_ref_sol(1:3, k), R(1:3, 1:3)) + Q(4, 4) * abs(u(4, k) - u_ref_sol(4, k)), 1:N) + quad_form(x(:, N) - x_ref_sol(:, N), Q);

%% Create Guess
guess.x = x_ref_sol;
guess.u = u_ref_sol;
guess.p = [];

%% Construct Problem Object
prob_6DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, scale_hint = scale_hint, terminal_bc = terminal_bc, nonconvex_constraints = nonconvex_constraints, discretization_method = "error", N_sub = 1);

%%
Delta = calculate_defect(prob_6DoF, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess

[prob_6DoF, Delta_disc] = prob_6DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_6DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_6DoF.cont_prop(guess.u, guess.p);
% 
figure
comparison_plot_6DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Solve Problem with PTR
%ptr_ops.Delta_min = 5e-5;
ptr_sol_vc = ptr(prob_6DoF, ptr_ops);

ptr_sol = ptr_sol_vc;

if ~ptr_sol.converged
    ptr_sol.converged_i = ptr_ops.iter_max;
end

%%
figure
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_6DoF.objective(prob_6DoF.guess.x, prob_6DoF.guess.u, prob_6DoF.guess.p), [ptr_sol.info.J]]); hold on
%yline(CasADi_sol.objective); hold off
hold off
legend("PTR Iterations", "CasADi Solution")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
title("Stopping Criteria vs Iteration")
yscale("log")
grid on

nexttile
plot(0:ptr_sol.converged_i, squeeze(sum(vecnorm(ptr_sol.Delta(:, :, 1:(ptr_sol.converged_i + 1)), 2, 1))))
yscale("log")
title("Defect Norm vs Iteration")
grid on

%%
i = ptr_sol.converged_i;

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_6DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
ref = prob_6DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
save("refrence_path.mat", "ref")

figure
plot_6DoFq_trajectory(t_k, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
comparison_plot_6DoF_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, ["Guess", "Continuous Propagation", "Solution Output"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Continuous vs Discrete Propagation of Solution")

figure
comparison_plot_6DoFq_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

function [q] = qexp(tau)
    theta = norm(tau);
    u = tau / theta;
    q = [u * sin(theta / 2); cos(theta / 2)];
end

function [tau] = qLog(q)
    N = size(q, 2);
    for k = 1 : N
        w = q(4, k);
        v = q(1:3, k);
        w = w * sign(w);
        v = v * sign(w);
        tau(:, k) = 2 * v * atan2(norm(v), w) / norm(v);
    end
end

function [tau] = RLog(R)
    theta = acos((trace(R) - 1) / 2);
    u = vee(R - R') / (2 * sin(theta));

    tau = theta * u;

    function [tau] = vee(tau_hat)
        tau = [tau_hat(3, 2); tau_hat(1, 3); tau_hat(2, 1)];
    end
end

function [] = plot_basis(n_i, basis_name, line_style)
%PLOT_BASIS Summary of this function goes here
%   Detailed explanation goes here

% Plot the basis
plot_vec(n_i(:, 1), "b", line_style, "$\hat " + basis_name + "_1$"); hold on;
plot_vec(n_i(:, 2), "r", line_style, "$\hat " + basis_name + "_2$"); hold on;
plot_vec(n_i(:, 3), "g", line_style, "$\hat " + basis_name + "_3$");
legend(interpreter = "latex")
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])
xlabel("X")
ylabel("Y")
zlabel("Z")

end

function [] = plot_vec(vec, color, line_style, name)
%PLOT_VEC Summary of this function goes here
%   Detailed explanation goes here
quiver3([0],[0],[0],[vec(1)],[vec(2)],[vec(3)], Color = color, LineStyle = line_style, DisplayName = name);
end

function [R] = quat_rotmatrix(q)
    w = q(4);
    v = q(1:3);

    R = (w ^ 2 - v' * v) * eye(3) + 2 * v * v' + 2 * w * skew(v);
end
% %%
% theta_0 = [deg2rad(0); deg2rad(90); deg2rad(0)]; % [rad]
% R_0 = make_R(deg2rad(10), 3) * angle2dcm(theta_0(1), theta_0(2), theta_0(3));
% q_0 = qexp(RLog(R_0));
% w_0 = deg2rad([0; 0; 0]); % [rad / s]
% glideslope_angle_max = deg2rad(80); % [rad]
% 
% theta_f = [0; deg2rad(90); 0]; % [rad]
% R_f = angle2dcm(theta_f(1), theta_f(2), theta_f(3));
% q_f = qexp(RLog(R_f));
% 
% Rq_0 = quat_rotmatrix(q_0);
% Rq_f = quat_rotmatrix(q_f);
% % 
% figure
% plot_basis(R_0, "q_0", "-"); hold on
% % plot_basis(R_f, "q_f", "--"); hold off
% %%
% quat_rotmatrix(x_cont_sol(7:10, end))
% 
% %%
% theta_1 = zeros([numel(t_cont_sol), 1]);
% theta_2 = zeros([numel(t_cont_sol), 1]);
% theta_3 = zeros([numel(t_cont_sol), 1]);
% for n = 1 : numel(t_cont_sol)
%     R = quat_rotmatrix(x_cont_sol(7:10, n));
%     [theta_1(n), theta_2(n), theta_3(n)] = dcm2angle(R, "YXZ");
% end
% 
% figure
% plot(t_cont_sol, rad2deg(theta_1), DisplayName="Pitch"); hold on
% plot(t_cont_sol, rad2deg(theta_2), DisplayName="Yaw"); hold on
% plot(t_cont_sol, rad2deg(theta_3), DisplayName="Roll"); hold off
% grid on
% legend()
% %%
% figure
% plot(t_cont_sol, acosd(pagemtimes([0, 0, 1], quat_rot_array(x_cont_sol(7:10, :), repmat([1; 0; 0], 1, numel(t_cont_sol)))))); hold on
% plot(t_cont_sol, acosd(2 * (x_cont_sol(7, :) .* x_cont_sol(9, :) - x_cont_sol(8, :) .* x_cont_sol(10, :)))); hold off
% %%
% % figure
% % plot(t_cont_sol, 1 - (2 * (x_cont_sol(7, :) .* x_cont_sol(9, :) - x_cont_sol(8, :) .* x_cont_sol(10, :))) .^ 2); hold on
% % yline(sin(deg2rad(45)) ^ 2)
% % hold off
% %%
% q_sym = sym("q", [4, 1]);
% dot([0; 0; 1], quat_rot(q_sym, [1; 0; 0]));
% 
% %%
t_iters = {};
x_iters = {};
u_iters = {};
linestyle = string(1:ptr_sol.converged_i);
for i = 1:ptr_sol.converged_i
    t_iters{i} = t_k;
    x_iters{i} = ptr_sol.x(:, :, i);
    u_iters{i} = ptr_sol.u(:, :, i);
    linestyle(i) = "-";
end
figure
comparison_plot_6DoF_trajectory(x_iters, "iter " + string(1:ptr_sol.converged_i), glideslope_angle_max, linestyle = linestyle, title = "Solution vs Iteration")

figure 
comparison_plot_6DoFq_time_histories(t_iters, x_iters, u_iters, "iter " + string(1:ptr_sol.converged_i), linestyle = linestyle, title = "Solution vs Iteration")

% %%
% u_sol = ptr_sol.u(:, :, i);
% gimbal_rate_check = acosd(dot(u_cont_sol(1:3, 1: (end - 1)), u_cont_sol(1:3, 2:end)) ./ (vecnorm(u_cont_sol(1:3, 1 : (end - 1))) .* vecnorm(u_cont_sol(1:3, 2 : end)))) ./ diff(t_cont_sol)';
% gimbal_rate_check_disc = acosd(dot(u_sol(1:3, 1: (end - 1)), u_sol(1:3, 2:end)) ./ (vecnorm(u_sol(1:3, 1 : (end - 1))) .* vecnorm(u_sol(1:3, 2 : end)))) / delta_t;
% 
% figure
% plot(t_cont_sol(2:end), gimbal_rate_check, t_k(1:(end - 1)), gimbal_rate_check_disc); 
