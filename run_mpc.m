%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 25 January, 2026
% Description: 6DoF landing of rocket using MPC
% Most Recent Change: 27 January, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Isp_actual = 220;
I_actual = [19200; 13600; 13500] * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA 

ref = load("reference.mat").ref_sol;

N_horizon = 10;
MPC_time_horizon = 5; % [s] How far MPC looks into the future
iter_max = 5; % Max convex subproblems per MPC solve
MPC_timestep = 0.3; % [s] How often MPC is run
max_time = ref.t_ref(end) / 2;

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
r_pert = [0.005; 0.0; 0];
v_pert = [0; 0; 0];
q_pert = [0; 0; 0; 0];
omega_pert = [0; 0; 0];
mass_pert = 0;

initial_perturbation = [r_pert; v_pert; q_pert; omega_pert; mass_pert];

x_0 = ref.x_ref_sol(:, 1) + initial_perturbation;
u_0 = ref.u_ref_sol(:, 1);

% Disturbances
disturbances = @(t, x) [zeros([3, 1]); [1; 0; 0]; zeros([4, 1]); [0.1; 0; 0]; 0];

%% Create MPC Solver
t_horizon = linspace(0, MPC_time_horizon, N_horizon); % MPC time horizon
[prob_6DoF_true] = MPC_6DoF_reynolds_quat_func(ref, t_horizon, Q, R, x_0, u_0, Isp_actual, I_actual);
[prob_6DoF_MPC] = MPC_6DoF_reynolds_quat_func(ref, t_horizon, Q, R, x_0, u_0);

% PTR algorithm parameters
ptr_ops.iter_max = iter_max;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 5e-4;
ptr_ops.w_vc = 1e3;
ptr_ops.w_tr = ones(1, N_horizon) * 4e-3;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 2e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;
ptr_ops.include_v_0 = false;
ptr_ops.include_u_0 = true;

%% Simulate 
x_hist = x_0;
u_hist = u_0;

k = 1;
t_k = 0;
while t_k < max_time
    % Run MPC
    t_horizon = linspace(t_k, t_k + MPC_time_horizon, N_horizon); % MPC time horizon
    
    [x_sol, u_sol, ptr_sol] = MPC_solve(prob_6DoF_MPC, x_0, u_0, ref, t_horizon, Q, R, ptr_ops);

    % Simulate forward until next MPC call
    u_func = @(t) interp1(t_horizon - t_k, u_sol', t)';

    [t_cont, x_cont] = ode45(@(t, x) prob_6DoF_true.cont.f(t, x, u_func(t), []) + disturbances(t, x), [0, MPC_timestep], x_0, prob_6DoF_true.tolerances);
    x_cont = x_cont';

    u_cont = u_func(t_cont);

    % Iterate
    k = k + 1;
    t_k = t_k + MPC_timestep;

    x_0 = x_cont(:, end);
    u_0 = u_cont(:, end);

    x_ref_t_k = interp1(ref.t_ref, ref.x_ref_sol', t_k, "previous")';
    x_error = x_0 - x_ref_t_k;
    fprintf("Completed iteration: %d, r err: %.3f, v err: %.3f, q_err: %.3f, w_err: %.3f \n", k, norm(x_error(1:3)) * 1e3, norm(x_error(4:6)) * 1e3, norm([0; 0; 0; 1] - q_mul(x_0(7:10), q_conj(x_ref_t_k(7:10)))), norm(x_error(11:13)))

    x_hist = [x_hist, x_cont(:, 2:end)];
    u_hist = [u_hist, u_cont(:, 2:end)];
end

%% Plot

figure
comparison_plot_6DoF_trajectory({ref.x_ref_sol, x_hist}, ["Reference", "MPC"], deg2rad(65), linestyle = [":", "-"], title = "MPC Tracking Performance")


%%
function [x_sol, u_sol, ptr_sol] = MPC_solve(prob, x_0, u_0, ref, t_horizon, Q, R, ptr_ops)
    %% Update reference trajectory
    x_ref_sol = interp1(ref.t_ref, ref.x_ref_sol', t_horizon, "previous")';
    u_ref_sol = interp1(ref.t_ref, ref.u_ref_sol', t_horizon, "previous")';
    u_ref_sol(:, 1) = u_0; % Make control start at the current control (for rate constraints)

    guess.x = x_ref_sol;
    guess.u = u_ref_sol;
    guess.p = []; % No parameters

    %% Update Problem
    prob.x0 = x_0; % Update initial state
    prob.initial_bc = @(x, p) x - x_0; % Update initial boundary condition
    prob.guess = guess; % Update reference trajectory

    prob.objective = @(x, u, p, k) einsum(@(k) quad_form(u(4, k) - u_ref_sol(1:3, k), R(1:3, 1:3)) + Q(4, 4) * abs(u(4, k) - u_ref_sol(4, k)), 1:numel(t_horizon)) + quad_form(x(:, end) - x_ref_sol(:, end), Q);

    %% Solve Problem with PTR
    ptr_sol_vc = ptr(prob, ptr_ops);
    
    ptr_sol = ptr_sol_vc;
    
    if ~ptr_sol.converged
        ptr_sol.converged_i = ptr_ops.iter_max;
    end

    x_sol = ptr_sol.x(:, :, ptr_sol.converged_i);
    u_sol = ptr_sol.u(:, :, ptr_sol.converged_i);
end