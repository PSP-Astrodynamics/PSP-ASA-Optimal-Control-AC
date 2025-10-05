%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
% PSP AC Optimal Control
% Author: Travis Hastreiter 
% Created On: 11 July, 2025 
% Description: Test of error state MPC on Lie groups for tracking 
% Most Recent Change: 11 July, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ASSUMING CONSTANT MASS

%% Define Group
G = SE2_2();
n = G.dim;

% State: [psi; twist] where psi is group error in tangent space and twist 
% is the velocity in the body frame
psi_ind = 1:n;
twist_ind = n + (1:n);

% Outputs: Configuration error and its velocity
C = @(twist_d) [eye(n), zeros(n, n); -G.ad(twist_d), zeros(n, n)];
d = @(twist_d) [zeros(n, 1); twist_d];
y = @(x) C(x(twist_ind)) * x - d(x(twist_ind));

%% Initial Conditions
r_0 = [0; 0; 4.6];
v_0 = make_R(-deg2rad(60), 2) * [0.306; 0; 0];
q_0 = dcm2quat(make_R(deg2rad(120), 2));
w_0 = [0; 0; 0];

r_d_0 = [0; 0; 0];
v_d_0 = [0; 0; 0];
q_d_0 = dcm2quat(make_R(deg2rad(90), 2));
w_d_0 = [0; 0; 0];

G_0 = SE2_2(q_0, r_0);
G_d_0 = SE2_2(q_0_d, r_0_d);

Psi_0 = group_error(G_0, G_d_0);
psi_0 = G.Log(Psi_0);
Psi_0_lin = linearized_group_error(G, psi_0);
twist_0 = [w_0; v_0]; 
x_0 = [psi_0; twist_0];
twist_d_0 = [w_d_0; v_d_0];

%% Problem Parameters
N = 20;
tf = 35;
t_k = linspace(0, tf, N);

nx = 2 * n;
nr = 2;
nu = nr + 1; 

%% Define Control
I = 150000 * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 3e-3; % [km] Distance from CoM to nozzle
m_dry = 1500; % [kg]
m_wet = 600; % [kg]
m_0 = m_dry + m_wet;

% ZERO ORDER HOLD
T_max = 3 * 9.81e-3; % [kg km / s2]
T_min = 0.5 * 9.81e-3; % [kg km / s2]
gimbal_max = deg2rad(6); % [rad]

%% Get Dynamics
f = @(x, u) Euler_Poincare(G, x(twist_ind), u);

%% Define Objective
P = eye(nx);
Q = eye(nx);
R = eye(nu);
objective = @(x, u) quad_form(y(x(:, N)), P) + einsum(@(k) quad_form(x(:, k), Q) + quad_form(u(:, k), R), 1:(N - 1));

%% Linearize about Desired Pose
%A = ;
%B = ;
%c = ;

%% Set up QP
cvx_begin
    variable x
    variable u
    minimize( objective(x, u) )
    subject to
        % Dynamics
        for k = 1:(N - 1)
            x(:, k + 1) == A(:, :, k) * x(:, k) ...
                         + B(:, :, k) * u(:, k) ...
                         + c(:, k);
        end

        % Constraints
        norms(u(1:nr, :)) - u(nr + 1) <= 0; % Lcvx constraint
        u(nr + 1) - T_max / m_0 <= 0; % Max thrust constraint
        T_min / m_0 - u(nr + 1) <= 0; % Min thrust constraint
        u(nr + 1) - u(1) / cos(gimbal_max) <= 0; % Gimbal constraint

        % Boundary Conditions
        x(:, 1) == x_0;
cvx_end

%% Get Solution
x_sol = x;
u_sol = u;
y_sol = zeros(nx, N);
for k = 1:N
    y_sol(:, k) = y(x(:, k));
end

%% Plot Solution
