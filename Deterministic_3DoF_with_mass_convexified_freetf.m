%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Vehicle Parameters
alpha = 0.5086; % [s / km]
m_dry = 1500; % [kg]
m_wet = 600; % [kg]
m_0 = m_dry + m_wet;
T_max = 3 * m_0 * 9.81e-3; % [kg km / s2]
T_min = 0.55 * T_max; % [kg km / s2]
I = 150000 * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 3e-3; % [km] Distance from CoM to nozzle
gimbal_max = deg2rad(6); % [rad]
 
vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I, alpha = alpha);
%vehicle_big_gimbs = Vehicle(m_dry, L, L * 3, deg2rad(20), T_min, T_max, I = I, alpha = alpha);

% Problem Parameters
tf = 25; % [s]
N = 25; % []
r_0 = [0; 4.6]; % [km]
theta_0 = deg2rad(120); % [rad]
v_0 = make_R2(-deg2rad(60)) * [0.306; 0]; % [km / s]
w_0 = deg2rad(0); % [rad / s]
glideslope_angle_max = deg2rad(50); % [rad]

x_0 = [r_0; v_0; theta_0; w_0; log(m_0)];
x_f = [zeros(2, 1); zeros(2, 1); pi / 2; 0];

tspan = [0, 1];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

nx = 7;
nu = 3;
np = 1;

initial_guess = "straight line"; % "CasADi" or "straight line"

% PTR algorithm parameters
ptr_ops.iter_max = 10;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e3;
ptr_ops.w_tr = ones(1, Nu) * 1e-3;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 3e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = true;

%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF_mass_convexified_freetf(t, x, u, p, vehicle.L, vehicle.I(2), vehicle.alpha);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = {1:N, @(t, x, u, p) norm(x(1:2)) - x(2) / cos(glideslope_angle_max)};
timef_min_constraint = {1, @(t, x, u, p) 20 - p(1)};
timef_max_constraint = {1, @(t, x, u, p) p(1) - 35};
state_convex_constraints = {glideslope_constraint, timef_min_constraint, timef_max_constraint};

z_lb = @(t, tf) log(m_0 - alpha * T_max * t * tf);
z_lb_k = z_lb(t_k, tf);

% Convex control constraints
max_gimbal_constraint = {1:N, @(t, x, u, p) u(3) - u(1) / cos(gimbal_max)};
lcvx_thrust_constraint = {1:N, @(t, x, u, p) norm(u(1:2)) - u(3)}; 
control_convex_constraints = {max_gimbal_constraint,lcvx_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
glideslope_STC_constraint = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) (norm(x(1:2)) - x(2) / cos(deg2rad(5))) * (x_ref(2, k) <= 0.5), ...
                                                       @(t, x, u, p, k) (norm(x(1:2)) - x(2) / cos(deg2rad(5))) * (x(2) <= 0.5)};
state_nonconvex_constraints = {glideslope_STC_constraint};

% Nonconvex control constraints
max_thrust_constraint = @(t, x, u, p) u(3) - T_max * exp(-z_lb(t, p(1))) * (1 - (x(7) - z_lb(t, p(1))));
max_thrust_constraint_linearized = {1:N, linearize_constraint(max_thrust_constraint, nx, nu, np, "p", 1), ...
                                         @(t, x, u, p, k) max_thrust_constraint(t, x, u, p)};
min_thrust_constraint = @(t, x, u, p) T_min * exp(-z_lb(t, p(1))) * (1 - (x(7) - z_lb(t, p(1))) + 0.5 * (x(7) - z_lb(t, p(1))) ^ 2) - u(3);
min_thrust_constraint_linearized = {1:N, linearize_constraint(min_thrust_constraint, nx, nu, np, "p", 1), ...
                                         @(t, x, u, p, k) min_thrust_constraint(t, x, u, p)};
control_nonconvex_constraints = {max_thrust_constraint_linearized, min_thrust_constraint_linearized};

nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

% Terminal boundary conditions
terminal_bc = @(x, u, p) [x(1:6, :) - x_f; 0];

%% Specify Objective
min_fuel_angular_velocity_objective = @(x, u, p) sum(u(3, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p, x_ref, u_ref, p_ref) -x(7);%sum(u(3, :)) * p_ref(1) / N + sum(u_ref(3, :)) / N * (p(1) - p_ref(1));
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p, x_ref, u_ref, p_ref) -x(7);%sum((u(3, 1:(end - 1)) + u(3, 2:end)) / 2) * p_ref(1) / N + sum((u_ref(3, 1:(end - 1)) + u_ref(3, 2:end)) / 2) * (p(1) - p_ref(1)) / N;
end

%% Create Guess
sl_guess = guess_3DoF(x_0(1:6), x_f + [0; 0; 0; 0; 0; 0], N, Nu, delta_t, vehicle);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(3, :) * delta_t), sum(sl_guess.u(3, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(3, :) * delta_t)];
end
sl_guess.x(7, :) = log(sl_guess.x(7, :));
sl_guess.u = sl_guess.u .* exp(-sl_guess.x(7, 1:Nu));
sl_guess.p = tf;

CasADi_sol = CasADi_solve_mass_convexified(x_0, sl_guess.x, sl_guess.u, vehicle, N, delta_t * tf, glideslope_angle_max);%

if initial_guess == "straight line"
    guess = sl_guess;
elseif initial_guess == "CasADi"
    guess = CasADi_sol;
end
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
guess.p = sl_guess.p;

t_k_scaled = t_k * tf;

figure
plot_3DoFc_trajectory(t_k_scaled, sl_guess.x, sl_guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)

figure
plot_3DoFc_time_histories(t_k_scaled, sl_guess.x, sl_guess.u)

figure
plot_3DoFc_trajectory(t_k_scaled, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
plot_3DoFc_time_histories(t_k_scaled, guess.x, guess.u)

%% Construct Problem Object
prob_3DoF = DeterministicProblem(x_0, x_f, N, u_hold, 1, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc, nonconvex_constraints = nonconvex_constraints, discretization_method = "error");

%% Test Scaling
% guess_scaled.x = prob_3DoF.scale_x(guess.x);
% guess_scaled.u = prob_3DoF.scale_u(guess.u);
% guess_scaled.p = prob_3DoF.scale_p(guess.p);
% 
% figure
% plot_3DoF_trajectory(t_k, guess_scaled.x, guess_scaled.u, glideslope_angle_max, gimbal_max, 0, 0)
% 
% figure
% plot_3DoF_time_histories(t_k, guess_scaled.x, guess_scaled.u)

%%
Delta = calculate_defect(prob_3DoF, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess

[prob_3DoF, Delta_disc] = prob_3DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_3DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.u, guess.p);

figure
comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

figure
comparison_plot_3DoFc_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Solve Problem with PTR
ptr_sol_vc = ptr(prob_3DoF, ptr_ops);

%%
% ptr_ops.w_vse = 1e4;
% ptr_ops.w_tr = 5e-2;
% ptr_ops.w_prime = 1e2;
% ptr_sol_vs = ptr_virtual_state(prob_3DoF, ptr_ops, "CVX");

%%
scvxstar_ops.D_x = 1e-1;
scvxstar_ops.D_u = 1e-1;
scvxstar_ops.D_p = 1e-1;
scvxstar_ops.opt_tol = ptr_ops.delta_tol;
scvxstar_ops.feas_tol = 5e-4;%ptr_ops.Delta_min;
scvxstar_ops.eta_0 = 1;
scvxstar_ops.eta_1 = 0.5;
scvxstar_ops.eta_2 = 0.1;
scvxstar_ops.alpha_1 = 1.1;
scvxstar_ops.alpha_2 = 4;
scvxstar_ops.beta = 2;
scvxstar_ops.gamma = 0.95;
scvxstar_ops.w_0 = 1e1;
scvxstar_ops.w_max = 1e6;
scvxstar_ops.r_0 = 1;
scvxstar_ops.r_min = 1e-8;
scvxstar_ops.r_max = 1;
scvxstar_ops.tau = 1.1;
scvxstar_ops.iter_max = ptr_ops.iter_max * 5;
scvxstar_ops.iter_min = 2;

scvxstar_sol = SCvx_star(prob_3DoF, scvxstar_ops, "CVX");

%%
ptr_sol = scvxstar_sol;

if ~ptr_sol.converged
    ptr_sol.converged_i = ptr_ops.iter_max;
end

%%
figure
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_3DoF.objective(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p), [ptr_sol.info.J]]); hold on
%yline(CasADi_sol.objective); hold off
legend("PTR Iterations", "CasADi Solution")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
yscale("log")
title("Stopping Criteria vs Iteration")
grid on

nexttile
plot(0:ptr_sol.converged_i, vecnorm(ptr_sol.Delta(:, 1:(ptr_sol.converged_i + 1)), 2, 1))
yscale("log")
title("Defect Norm vs Iteration")
grid on

%%
Js = [ptr_sol.info.J];
%FOH_Js = [FOH_ptr_sol.info.J];
(Js(end) - CasADi_sol.objective) / CasADi_sol.objective * 100
%(FOH_Js(end) - CasADi_sol.objective) / CasADi_sol.objective * 100

%%
figure
plot(0:ptr_sol.converged_i, squeeze(squeeze(sum(ptr_sol.u(3, :, 1:(ptr_sol.converged_i + 1)), 2)) .* ptr_sol.p(1, 1:(ptr_sol.converged_i + 1))' / N)); hold on
yline(CasADi_sol.objective); hold off
xlabel("PTR Iteration")
ylabel("\Delta V [km / s]")
legend("PTR Iterations", "CasADi Solution", location = "southeast")
grid on

%%
figure
plot(0:ptr_sol.converged_i, vecnorm(ptr_sol.Delta(:, 1:(ptr_sol.converged_i + 1)), 2, 1))
xlabel("PTR Iteration")
ylabel("Defect")
yscale("log")
xlim([0, ptr_sol.converged_i])
grid on

%%
i = ptr_sol.converged_i + 1;

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));

t_scaled = t_cont_sol * ptr_sol.p(:, i);
t_k_scaled = t_k * ptr_sol.p(:, i);

figure
plot_3DoFc_trajectory(t_k_scaled, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, gimbal_max, T_min, T_max, step = 1, title = "")
%%
figure
comparison_plot_3DoF_trajectory({guess.x, x_cont_sol, CasADi_sol.x}, ["Guess", "PTR", "CasADi"], glideslope_angle_max, linestyle = [":", "-", "-"], title = "")

%%
figure
comparison_plot_3DoFc_time_histories({t_k_scaled, t_scaled, t_k_scaled}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "")

%%
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
comparison_plot_3DoF_trajectory(x_iters, "iter " + string(1:ptr_sol.converged_i), glideslope_angle_max, linestyle = linestyle, title = "Solution vs Iteration")

figure 
comparison_plot_3DoFc_time_histories(t_iters, x_iters, u_iters, "iter " + string(1:ptr_sol.converged_i), linestyle = linestyle, title = "Solution vs Iteration")