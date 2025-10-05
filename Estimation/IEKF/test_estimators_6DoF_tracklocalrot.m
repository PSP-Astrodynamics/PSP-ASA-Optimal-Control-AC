%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
% PSP AC Optimal Control
% Author: Travis Hastreiter 
% Created On: 21 July, 2025
% Description: Filtering test on simple 6DoF rocket while tracking local
% orientation (R from body to inertial is in state)
% Most Recent Change: 21 July, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
% - get working then show Adam on HA ASAP! + AC GNC
% - "estimate" mass
% - continuous versions
% - also propagate error
%   - plot to show that IEKF linear error propagation holds
% - sick plots
% - another version with stochastic perturbations and external forces for MPC
% - add right variations of ESKF and IEKF, should be worse since everything
%      is left invariant
% - add UKF

% Monte Carlo steps: (there is only 1 true vehicle trajectory)
% - Sample from initial distribution
% - Propagate to next measurement deterministically (because only 1 true vehicle trajectory)
%   - Propagation of true vehicle trajectory is stochastic
% - Update through measurement of true trajectory but each vehicle samples
% its own measurement
% - Analyze error evolution and final distribution

%% Initialize
% Rocket
% Vehicle Parameters
alpha = 0.5086; % [s / km]
m_dry = 1100; % [kg]
m_wet = 1000; % [kg]
m_0 = m_dry + m_wet;
g = 9.81e-3;
g_vec = -[0; 0; g];
T_max = 2 * m_0 * g; % [kg km / s2]
T_min = 0.55 * m_0 * g; % [kg km / s2]
I = [5000 * (1e-3) ^ 2; 15000 * (1e-3) ^ 2; 15000 * (1e-3) ^ 2]; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 1e-3; % [km] Distance from CoM to nozzle
gimbal_max = deg2rad(6); % [rad]
 
% Dynamics
f = @(t, x, u, p) SymDynamicsQuat6DoF_localrot(x, u, L, I, alpha, g);

% Measurements
seed = 42;

measure_magnetometer = false;
measure_barometer = true;
measure_gps = true;

% IMU just used for propagation (for right now)
imu_sample_rate = 0.1; % [s]
imu_acc_accel = 0.1 / 3; % [m / s2]
imu_acc_gyro = deg2rad(1) / 3; % [rad / s]
imu_acc_mag = 1e-0 / 3; % [micro Tesla]
imu_sensor = IMU(1 / imu_sample_rate, ...
                 accelparams("NoiseDensity",imu_acc_accel), ...
                 gyroparams("NoiseDensity",imu_acc_gyro), ...
                 magparams("NoiseDensity",imu_acc_mag), ...
                 seed = seed, scale = [1e-3 * ones([3, 1]); ones([3, 1]); ones([3, 1])]);
mag_H_EKF = mag_H_EKF_func();

barometer_sample_rate = 0.2; % [s]
baro_acc_press = 5 / 3; % [Pa]
baro_sensor = Barometer(1 / barometer_sample_rate, baro_acc_press, seed = seed, scale = 1e-3);

gps_sample_rate = 1; % [s]
gps_acc_horiz = 3 / 3; % [m]
gps_acc_vert = 3 / 3; % [m]
gps_sensor = GPS(1 / gps_sample_rate, gps_acc_horiz, gps_acc_vert, seed = seed, scale = 1e-3);

sensors = {baro_sensor, gps_sensor};

% Estimators
% Set up A, G, Q, etc for each filter - ideally based on IMU noise but probably
% inflated for Q_EKF and Q_ESKF?
Q_IMU = diag([imu_acc_gyro * ones([3, 1]); imu_acc_accel * ones([3, 1]) .* 1e-3; zeros([3, 1])] .^ 2);
Q_EKF = Q_IMU([7, 8, 9, 4, 5, 6, 1, 2, 3], [7, 8, 9, 4, 5, 6, 1, 2, 3]);
Q_ESKF = Q_IMU([7, 8, 9, 4, 5, 6, 1, 2, 3], [7, 8, 9, 4, 5, 6, 1, 2, 3]);
Q_IEKF = Q_IMU;

[A_EKF, G_EKF] = A_G_EKF_euler212(); % Both @(x, u)
[A_ESKF, G_ESKF] = A_G_LESKF(); % Left attitude error state EKF - should also try right

%% Load in reference trajectory and control
load("takeoffhoverdivert_3DoF.mat")

t_2D_opt = linspace(0, sol_struct.prob.tf, sol_struct.prob.N);
t_2D_cont = sol_struct.cont_prop_t;
% [r_x, r_y, v_x, v_y, theta, w, log(m)] 
x_2D_opt = sol_struct.opt_sol.x(:, :, sol_struct.opt_sol.converged_i);
x_2D_cont = sol_struct.cont_prop_x;
% 2D u is [a_x, a_y, |a|] in body frame (a is acceleration)
u_2D_opt = sol_struct.opt_sol.u(:, :, sol_struct.opt_sol.converged_i);
u_2D_cont = sol_struct.cont_prop_u;
% Convert 2D inputs to 3D - y is up, z is normal to plane of movement
% 3D u is [T_x, T_y, T_z. |T|, tau] in body frame where tau is roll moment (T is thrust)
u = [u_2D_cont(1, :); zeros([1, length(u_2D_cont)]); u_2D_cont(2, :); u_2D_cont(3, :); zeros([1, length(u_2D_cont)])] .* exp(x_2D_cont(7, :));
u_ref = @(t) interp1(t_2D_cont(1:size(u, 2)), u', t, "linear", "extrap")'; % There exists faster way to get u(t)...

%% Simulation
tf = sol_struct.prob.tf; % [s]
N = ceil(tf / imu_sample_rate) + 1;
t_k = (0 : (N - 1)) * imu_sample_rate;
tspan = [0, t_k(end)];
r_0 = [50; 0; 0] .* 1e-3; % [km]
theta_0 = [0; deg2rad(90); 0]; % [rad]
R_0 = angle2dcm(theta_0(1), theta_0(2), theta_0(3));
q_0 = qexp(RLog(R_0));
% 
% theta = -90;
% u_vec = [0; 1; 0];
% q_0 = [sind(theta / 2) * u_vec; cosd(theta / 2)];
v_0 = [0; 0; 0] .* 1e-3; % [km / s]
w_0 = [deg2rad(0); deg2rad(0); deg2rad(0)]; % [rad / s]

x_0 = [r_0; v_0; q_0; w_0; m_0];
x_f = [zeros(2, 1); zeros(2, 1); pi / 2; 0];

% DOESN'T INCLUDE MASS OR ANGUlAR VELOCITY
X_0 = SE2_3_RotationMatrix(quat_rotmatrix(x_0(7:10)), v_0, r_0);
% 
% r = eye(3) + sind(theta) * skew(u_vec) + (1 - cosd(theta)) * skew(u_vec) ^ 2
% rq = quat_rotmatrix(q_0)
% rqm = quat2dcm(q_0([4, 1, 2, 3])')
% 
% tau = RLog(rq);
% theta_c = rad2deg(norm(tau));
% u_vec_c = tau / deg2rad(theta_c);
% q_0_c = [sind(theta_c / 2) * u_vec_c; cosd(theta_c / 2)];
% 
% r = eye(3) + sind(theta_c) * skew(u_vec_c) + (1 - cosd(theta_c)) * skew(u_vec_c) ^ 2
% 
% plot_basis(eye(3), "L", "-");
% %plot_basis(r, "Ge", "--");
% %plot_basis(rq, "Gq", "-.");
% plot_vec(rq * u(1:3, 1) / u(4, 1), "k", ":", "u_1");
% axis equal
% 
% angle2dcm(theta_0(1), theta_0(2), theta_0(3)) * u(1:3, 1)
% quat_rot(x_0(7:10), u(1:3, 1))
% quat_rotmatrix(x_0(7:10)) * u(1:3, 1)
% X_0.R * u(1:3, 1)

nu = 4;
nr = 3;
nx = 14;

integration_tolerance = 1e-12;
tolerances = odeset(RelTol=integration_tolerance, AbsTol=integration_tolerance);

x_disc = zeros([nx, N]);
x_disc(:, 1) = x_0(:, 1);
xtilde_disc = x_disc;
xtilde_disc(:, 1) = x_0(:, 1);

Ptilde0_l = diag(([deg2rad(2), ...
                deg2rad(2), ...
                deg2rad(2), ...
                0.1 .* 1e-3, ...
                0.1 .* 1e-3, ...
                0.1 .* 1e-3, ...
                0.2 .* 1e-3, ...
                0.2 .* 1e-3, ...
                0.2 .* 1e-3] / 3) .^ 2); % In global frame

Ptilde0_r = X_0.Ad * Ptilde0_l * X_0.Ad';

Sigma_hat_EKF = zeros([9, 9, N]); % 2-1-2 Body Euler angle error covariance
Sigma_hat_EKF(:, :, 1) = Ptilde0_r([7, 8, 9, 4, 5, 6, 1, 2, 3], [7, 8, 9, 4, 5, 6, 1, 2, 3]); 
Sigma_hat_ESKF = zeros([9, 9, N]); % Left error covariance
Sigma_hat_ESKF(:, :, 1) = blkdiag(Ptilde0_r([7, 8, 9, 4, 5, 6], [7, 8, 9, 4, 5, 6]), Ptilde0_l([1, 2, 3], [1, 2, 3]));
Sigma_hat_IEKF = zeros([9, 9, N]); % Left invariant error covariance
Sigma_hat_IEKF(:, :, 1) = Ptilde0_l;

n_sample = 100;


%% Sample Initial Conditions
% using X_0.sample_right_gaussian because it does rotations fine :)
for k = 1 : n_sample
    X_0_sample(k) = sample_right_gaussian(X_0, Ptilde0_r);
end

x_0_sample = zeros([10, 1]);
x_0_sample(1:3) = X_0_sample(1).r;
x_0_sample(4:6) = X_0_sample(1).v;
x_0_sample(7:10) = qexp(RLog(X_0_sample(1).R)); %[0, 0, 0, -1; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0] * dcm2quat(X_0_sample(1).R)';

x_hat_EKF = zeros([9, N]);
x_hat_EKF(:, 1) = quatstate2eulerstate(x_0_sample(1:10));
x_hat_ESKF = zeros([10, N]);
x_hat_ESKF(:, 1) = x_0_sample(1:10);
X_hat_IEKF(1) = X_0_sample(1);

%% Simulation Loop
last_mag_reading_t = -1;
last_barometer_reading_t = -1;
last_gps_reading_t = -1;

t_prop = [];
x_prop = [];
x_true = [];
x_true = [x_true, x_0];

X_trues(1) = X_0;

imu_reading = zeros(9, N - 1); % [accel, gyro, mag]
baro_reading = zeros(1, N - 1);
gps_reading = zeros(3, N - 1);

twists(1) = posquat2twist(x_0, f(0, x_0, u_ref(0), []));

x_am1 = x_0;
for a = 2 : N
    t = t_k(a);
    %% Propagate Simulation
    % Integrate dynamics for true trajectory - deterministic for now
    [t_int, x_int] = ode45(@(t, x) f(t, x, u_ref(t), []), [t_k(a - 1), t_k(a)], x_am1, tolerances);
    x_a = x_int(end, :)';
    x_true = [x_true, x_a];

    t_prop = [t_prop; t_int];
    x_prop = [x_prop, x_int'];

    X_true = posquat2SE2_3_Rotm(x_a);
    X_trues(a) = X_true;
    twist = posquat2twist(x_a, f(t, x_a, u_ref(t), []));
    twists(a) = twist;

    % Sample measurements
    imu_reading(:, a - 1) = imu_sensor.measure(X_true, twist);
    if measure_magnetometer
        mag_reading = imu_reading(7:9, a - 1);
        last_mag_reading_t = t;
    end
    if t >= last_barometer_reading_t + barometer_sample_rate && measure_barometer
        baro_reading(:, a - 1) = baro_sensor.measure(X_true, twist);
        last_barometer_reading_t = t;
    end
    if t >= last_gps_reading_t + gps_sample_rate && measure_gps
        gps_reading(:, a - 1) = gps_sensor.measure(X_true, twist);
        last_gps_reading_t = t;
    end

    %% EKF
    % EKF Time Update
    [x_hat_plus_EKF, Sigma_hat_plus_EKF] = EKF_disc_time_update(x_hat_EKF(:, a - 1), Sigma_hat_EKF(:, :, a - 1), Q_EKF, imu_sample_rate, imu_reading(:, a - 1), g_vec, A_EKF, G_EKF, A_rot_euler212(), rotation_euler212());
   
    R_b_i = rotation_euler212();

    X_hat_plus_EKF = SE2_3_RotationMatrix(R_b_i(x_hat_plus_EKF(7:9), eye(3)), x_hat_plus_EKF(4:6), x_hat_plus_EKF(1:3));

    % EKF Measurement Update
    if last_mag_reading_t == t
        [x_hat_plus_EKF, Sigma_hat_plus_EKF] = EKF_measurement_update(x_hat_plus_EKF, Sigma_hat_plus_EKF, mag_reading, mag_h(x_hat_plus_EKF(7:9), imu_sensor.magnetic_field', @(x, p) R_b_i(x, eye(3))' * p), mag_H_EKF(x_hat_plus_EKF, imu_sensor.magnetic_field'), imu_sensor.mag_N());
        %display(mag_h(x_hat_plus_EKF(7:9), imu_sensor.magnetic_field', R_b_i) - mag_reading)
    end
    if last_barometer_reading_t == t
        [x_hat_plus_EKF, Sigma_hat_plus_EKF] = EKF_measurement_update(x_hat_plus_EKF, Sigma_hat_plus_EKF, baro_reading(:, a - 1), baro_h(x_hat_plus_EKF), baro_H_EKF(), baro_sensor.N(X_hat_plus_EKF));
    end
    if last_gps_reading_t == t
        [x_hat_plus_EKF, Sigma_hat_plus_EKF] = EKF_measurement_update(x_hat_plus_EKF, Sigma_hat_plus_EKF, gps_reading(:, a - 1), gps_h(x_hat_plus_EKF), gps_H_EKF(), gps_sensor.N);
    end

    x_hat_EKF(:, a) = x_hat_plus_EKF;
    Sigma_hat_EKF(:, :, a) = Sigma_hat_plus_EKF;

    %% Left QEKF (ESKF)
    % QEKF Time Update
    [x_hat_plus_ESKF, Sigma_hat_plus_ESKF] = ESKF_disc_time_update_localrot(x_hat_ESKF(:, a - 1), Sigma_hat_ESKF(:, :, a - 1), Q_ESKF, imu_sample_rate, imu_reading(:, a - 1), g_vec, A_ESKF, G_ESKF, x_a);
        
    X_hat_plus_ESKF = SE2_3_RotationMatrix(quat_rotmatrix(x_hat_plus_ESKF(7:10)), x_hat_plus_ESKF(4:6), x_hat_plus_ESKF(1:3));
 
    % % QEKF Measurement Update
    if last_mag_reading_t == t
        [x_hat_plus_ESKF, Sigma_hat_plus_ESKF] = ESKF_measurement_update(x_hat_plus_ESKF, Sigma_hat_plus_ESKF, mag_reading, mag_h(x_hat_plus_ESKF(7:10), imu_sensor.magnetic_field', @(q, p) quat_rot(q_conj(q), p)), mag_H_LESKF(x_hat_plus_ESKF, imu_sensor.magnetic_field'), imu_sensor.mag_N(), @(x1, x2) LESKF_add_posquat(x1, x2), [0; 0; 0; 0; 0; 0; 1; 1; 1], x_a);
        %[x_hat_plus_ESKF, Sigma_hat_plus_ESKF] = ESKF_measurement_update(x_hat_plus_ESKF, Sigma_hat_plus_ESKF, skew(mag_reading / norm(mag_reading)) * mag_h(x_hat_plus_ESKF(7:10), imu_sensor.magnetic_field' / norm(imu_sensor.magnetic_field), @(q, p) quat_rot(q_conj(q), p)), 0, mag2_H_LESKF(x_hat_plus_ESKF, mag_reading, imu_sensor.magnetic_field'), imu_sensor.mag_N(), @(x1, x2) LESKF_add_posquat(x1, x2), x_a);
        %[x_hat_plus_ESKF, Sigma_hat_plus_ESKF] = ESKF_measurement_update(x_hat_plus_ESKF, Sigma_hat_plus_ESKF, -2 * skew(mag_reading) * mag_h(x_hat_plus_ESKF(7:10), imu_sensor.magnetic_field', @(q, p) quat_rot(q_conj(q), p)), 0, mag3_H_LESKF(x_hat_plus_ESKF, mag_reading, imu_sensor.magnetic_field'), imu_sensor.mag_N(), @(x1, x2) LESKF_add_posquat(x1, x2), 0, x_a);
    end
    if last_barometer_reading_t == t
        [x_hat_plus_ESKF, Sigma_hat_plus_ESKF] = ESKF_measurement_update(x_hat_plus_ESKF, Sigma_hat_plus_ESKF, baro_reading(:, a - 1), baro_h(x_hat_plus_ESKF), baro_H_EKF(), baro_sensor.N(X_hat_plus_ESKF), @(x1, x2) LESKF_add_posquat(x1, x2), [0; 0; 1; 0; 0; 1; 0; 0; 1e-12], x_a);
    end
    if last_gps_reading_t == t
        [x_hat_plus_ESKF, Sigma_hat_plus_ESKF] = ESKF_measurement_update(x_hat_plus_ESKF, Sigma_hat_plus_ESKF, gps_reading(:, a - 1), gps_h(x_hat_plus_ESKF), gps_H_EKF(), gps_sensor.N, @(x1, x2) LESKF_add_posquat(x1, x2), [1; 1; 1; 1; 1; 1; 0; 0; 1e-12], x_a);
    end

    x_hat_ESKF(:, a) = x_hat_plus_ESKF;
    Sigma_hat_ESKF(:, :, a) = Sigma_hat_plus_ESKF;

    %% LIEKF

    % LIEKF Time Update
    [X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF] = LIEKF_disc_time_update(X_hat_IEKF(a - 1), Sigma_hat_IEKF(:, :, a - 1), Q_IEKF, imu_sample_rate, imu_reading(:, a - 1), g_vec, X_true);
    
    % % LIEKF Measurement Update
    if last_mag_reading_t == t
        [X_hat_plus_LIEKF_test, Sigma_hat_plus_LIEKF_test] = LIEKF_measurement_update(X_hat_plus_LIEKF, Sigma_hat_plus_ESKF, imu_sensor.mag_innovation(X_hat_plus_LIEKF, mag_reading), imu_sensor.mag_inv_covariance_innovation(X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF, imu_sensor.mag_H_l(X_hat_plus_LIEKF)), imu_sensor.mag_H_l(X_hat_plus_LIEKF), imu_sensor.mag_N(), X_true);
        H_l = imu_sensor.mag_H_l(X_hat_plus_LIEKF);
        v_l_linear = imu_sensor.mag_innovation(X_hat_plus_LIEKF, mag_reading);
        v_l = -2 * skew(mag_reading / norm(mag_reading)) * X_hat_plus_LIEKF.R' * imu_sensor.magnetic_field' / norm(imu_sensor.magnetic_field);
        S_inv = imu_sensor.mag_inv_covariance_innovation(X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF, H_l);
        L = Sigma_hat_plus_LIEKF * H_l' * S_inv;
        delta_X = L * v_l_linear; %[-v_l; zeros(6, 1)];
        X_hat_plus = X_hat_plus_LIEKF.rplus(-delta_X);
        
        
        X_hat_plus_LIEKF.R' * imu_sensor.magnetic_field'
        X_hat_plus.R' * imu_sensor.magnetic_field'
        X_true.R' * imu_sensor.magnetic_field'
        mag_reading

        err_l_before = norm(left_invariant_error(X_true, X_hat_plus_LIEKF));
        err_l_after = norm(left_invariant_error(X_true, X_hat_plus));
    end
    if last_barometer_reading_t == t
        [X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF] = LIEKF_measurement_update(X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF, baro_sensor.innovation(X_hat_plus_LIEKF, [x_a(1); x_a(2); baro_reading(:, a - 1)]), baro_sensor.inv_covariance_innovation(X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF, baro_sensor.H_l), baro_sensor.H_l, baro_sensor.N(X_hat_plus_LIEKF), X_true);
    end
    if last_gps_reading_t == t
        [X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF] = LIEKF_measurement_update(X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF, gps_sensor.innovation(X_hat_plus_LIEKF, gps_reading(:, a - 1)), gps_sensor.inv_covariance_innovation(X_hat_plus_LIEKF, Sigma_hat_plus_LIEKF, gps_sensor.H_l), gps_sensor.H_l, gps_sensor.N, X_true);
    end

    X_hat_IEKF(a) = X_hat_plus_LIEKF;
    Sigma_hat_IEKF(:, :, a) = Sigma_hat_plus_LIEKF;

    %%
    
    x_am1 = x_a;
end
%% Convert SE2_3 to quatstate
x_hat_IEKF = zeros([10, N]);
for k = 1 : N
    x_hat_IEKF(1:3, k) = X_hat_IEKF(k).r;
    x_hat_IEKF(4:6, k) = X_hat_IEKF(k).v;
    % Black magic
    %x_hat_IEKF(7:10, k) = [0, 0, 0, -1; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0] * dcm2quat(X_hat_IEKF(k).R)';
    % Correct
    x_hat_IEKF(7:10, k) = qexp(RLog(X_hat_IEKF(k).R));
    r = quat_rotmatrix(x_hat_IEKF(7:10, k));
end

%%
for k = 1 : N - 1
    mag_EKF(:, k) = R_b_i(x_hat_EKF(7:9, k + 1), eye(3)) * imu_reading(7:9, k);
    mag_ESKF(:, k) = quat_rot(x_hat_ESKF(7:10, k + 1), imu_reading(7:9, k));
    mag_IEKF(:, k) = quat_rot(x_hat_IEKF(7:10, k + 1), imu_reading(7:9, k));%X_hat_IEKF(k + 1).R * imu_reading(7:9, k);
    mag_IEKF_true(:, k) = X_trues(k + 1).R * imu_reading(7:9, k);
end

%%
figure
plot(t_k(1:N-1), imu_sensor.magnetic_field' - mag_EKF); hold on
plot(t_k(1:N-1), imu_sensor.magnetic_field' - mag_ESKF, LineStyle="--")
plot(t_k(1:N-1), imu_sensor.magnetic_field' - mag_IEKF, LineStyle="-.")
plot(t_k(1:N-1), imu_sensor.magnetic_field' - mag_IEKF_true, LineStyle=":")
grid on
hold off

%%
figure
tiledlayout(1, 2)
nexttile
plot([twists.a]');
grid on
nexttile
plot(imu_reading(1:3, :)')
grid on

% figure
% tiledlayout(1, 2)
% nexttile
% plot([twists.w]');
% grid on
% nexttile
% plot(imu_reading(4:6, :)')
% grid on

% %%
% figure
% plot(imu_sensor.magnetic_field - quat_rot_array(x_true(7:10, 2:end), imu_reading(7:9, :))')
% grid on

%%
% -2 * skew(mag_reading) * mag_h(x_hat_plus_ESKF(7:10), imu_sensor.magnetic_field', @(q, p) quat_rot(q_conj(q), p))
% mag3_H_LESKF(x_hat_plus_ESKF, mag_reading, imu_sensor.magnetic_field')


%%

figure
plot3(x_prop(1, :) .* 1e3, x_prop(2, :) .* 1e3, x_prop(3, :) .* 1e3, DisplayName="Prop"); hold on
plot3(x_2D_cont(1, :) .* 1e3, x_2D_cont(1, :) .* 0, x_2D_cont(2, :) .* 1e3, DisplayName="True");
%plot3(x_hat_EKF(1, :) .* 1e3, x_hat_EKF(2, :) .* 1e3, x_hat_EKF(3, :) .* 1e3, DisplayName="EKF");
%plot3(x_hat_ESKF(1, :) .* 1e3, x_hat_ESKF(2, :) .* 1e3, x_hat_ESKF(3, :) .* 1e3, DisplayName="ESKF"); 
%plot3(x_hat_IEKF(1, :) .* 1e3, x_hat_IEKF(2, :) .* 1e3, x_hat_IEKF(3, :) .* 1e3, DisplayName="IEKF"); 
xlim([-10, 60] * 1e-0)
ylim([-35, 35] * 1e-0)
zlim([0, 70] * 1e-0)
hold off
grid on
legend()


% 
% figure
% 
% tiledlayout(3, 2)
% 
% nexttile
% plot(t_prop, x_prop(1, :) * 1e3, DisplayName="x"); hold on
% plot(t_prop, x_prop(2, :) * 1e3, DisplayName="y")
% plot(t_prop, x_prop(3, :) * 1e3, DisplayName="z")
% hold off
% legend()
% grid on
% 
% nexttile
% plot(t_prop, x_prop(4, :) * 1e3, DisplayName="vx"); hold on
% plot(t_prop, x_prop(5, :) * 1e3, DisplayName="vy")
% plot(t_prop, x_prop(6, :) * 1e3, DisplayName="vz")
% hold off
% legend()
% grid on
% 
% nexttile
% plot(t_prop, x_prop(7, :), DisplayName="qx"); hold on
% plot(t_prop, x_prop(8, :), DisplayName="qy")
% plot(t_prop, x_prop(9, :), DisplayName="qz")
% plot(t_prop, x_prop(10, :), DisplayName="qw")
% hold off
% legend()
% grid on
% 
% nexttile
% plot(t_prop, x_prop(11, :), DisplayName="wx"); hold on
% plot(t_prop, x_prop(12, :), DisplayName="wy")
% plot(t_prop, x_prop(13, :), DisplayName="wz")
% hold off
% legend()
% grid on
% 
% nexttile
% plot(t_2D_cont, u(1, :), DisplayName="Tx"); hold on
% plot(t_2D_cont, u(2, :), DisplayName="Ty")
% plot(t_2D_cont, u(3, :), DisplayName="Tz")
% hold off
% legend()
% grid on
% 
% nexttile
% plot(t_prop, x_prop(14, :) * g)
% grid on
figure
x_prop_euler = quatstate2eulerstate(x_prop);
x_ESKF_euler = quatstate2eulerstate(x_hat_ESKF);
x_IEKF_euler = quatstate2eulerstate(x_hat_IEKF);
plot(t_prop, rad2deg(x_prop_euler(7, :)), DisplayName="e2_1"); hold on
plot(t_prop, rad2deg(x_prop_euler(8, :)), DisplayName="e1_2")
plot(t_prop, rad2deg(x_prop_euler(9, :)), DisplayName="e2_3")
plot(t_k, rad2deg(x_hat_EKF(7, :)), DisplayName="e2_1_EKF", LineStyle=":")
plot(t_k, rad2deg(x_hat_EKF(8, :)), DisplayName="e1_2_EKF", LineStyle=":")
plot(t_k, rad2deg(x_hat_EKF(9, :)), DisplayName="e2_3_EKF", LineStyle=":")
plot(t_k, rad2deg(x_ESKF_euler(7, :)), DisplayName="e2_1_ESKF", LineStyle="--")
plot(t_k, rad2deg(x_ESKF_euler(8, :)), DisplayName="e1_2_ESKF", LineStyle="--")
plot(t_k, rad2deg(x_ESKF_euler(9, :)), DisplayName="e2_3_ESKF", LineStyle="--")
plot(t_k, rad2deg(x_IEKF_euler(7, :)), DisplayName="e2_1_IEKF", LineStyle="-.")
plot(t_k, rad2deg(x_IEKF_euler(8, :)), DisplayName="e1_2_IEKF", LineStyle="-.")
plot(t_k, rad2deg(x_IEKF_euler(9, :)), DisplayName="e2_3_IEKF", LineStyle="-.")
hold off
legend()
grid on

%%
figure
tiledlayout(2, 2)
nexttile
plot(t_k, x_true(7:10, :) - x_hat_ESKF(7:10, :), DisplayName="ESKF"); hold on
plot(t_k, x_true(7:10, :) - x_hat_IEKF(7:10, :), LineStyle="-.", DisplayName = "IEKF");
hold off
grid on
legend()
xlabel("Time [s]")
ylabel("q component diff")
title("Component Difference")

nexttile
plot(t_prop, x_prop(7:10, :), DisplayName="true"); hold on
plot(t_k, x_hat_ESKF(7:10, :), DisplayName="ESKF", LineStyle=":");
plot(t_k, x_hat_IEKF(7:10, :), DisplayName="ISKF", LineStyle="-."); 
hold off
legend()
xlabel("Time [s]")
ylabel("q component")
title("Components")
grid on

nexttile
xi_theta_l_ESKF = rad2deg(quat_left_error(x_true(7:10, :), x_hat_ESKF(7:10, :)));
xi_theta_l_IEKF = rad2deg(quat_left_error(x_true(7:10, :), x_hat_IEKF(7:10, :)));
plot(t_k, xi_theta_l_ESKF, t_k, vecnorm(xi_theta_l_ESKF), DisplayName="ESKF"); hold on
plot(t_k, xi_theta_l_IEKF, t_k, vecnorm(xi_theta_l_IEKF), DisplayName="IEKF", LineStyle="-.");
hold off
grid on
xlabel("Time [s]")
ylabel("Error [deg]")
title("Left Error")
legend()

nexttile
xi_theta_r_ESKF = rad2deg(quat_right_error(x_hat_ESKF(7:10, :), x_true(7:10, :)));
xi_theta_r_IEKF = rad2deg(quat_right_error(x_hat_IEKF(7:10, :), x_true(7:10, :)));
plot(t_k, xi_theta_r_ESKF, t_k, vecnorm(xi_theta_r_ESKF), DisplayName="ESKF"); hold on
plot(t_k, xi_theta_r_IEKF, t_k, vecnorm(xi_theta_r_IEKF), DisplayName="IEKF", LineStyle="-."); 
hold off
grid on
xlabel("Time [s]")
ylabel("Error [deg]")
title("Right Error")
legend()

%%
figure
plot(t_prop, x_prop(1:3, :) * 1e3, DisplayName="true"); hold on
plot(t_k, x_hat_EKF(1:3, :) * 1e3, DisplayName="EKF", LineStyle=":");
plot(t_k, x_hat_ESKF(1:3, :) * 1e3, DisplayName="ESKF", LineStyle="-.");
plot(t_k, x_hat_IEKF(1:3, :) * 1e3, DisplayName="IEKF", LineStyle="--"); hold off
hold off
legend()
grid on
title("Position Comparison")
%%
figure
plot(t_prop, x_prop(4:6, :) * 1e3, DisplayName="true"); hold on
plot(t_k, x_hat_EKF(4:6, :) * 1e3, DisplayName="EKF", LineStyle=":");
plot(t_k, x_hat_ESKF(4:6, :) * 1e3, DisplayName="ESKF", LineStyle="-."); 
plot(t_k, x_hat_IEKF(4:6, :) * 1e3, DisplayName="IEKF", LineStyle="--"); 
hold off
legend()
grid on
title("Velocity Comparison")

%%
figure
tiledlayout(3, 3)
scale = [ones(6, 1) * 1e-3; pi / 180 * ones(3, 1)];
error_EKF = abs(quatstate2eulerstate(x_true) - x_hat_EKF);
for state = 1:9
    nexttile
    x_3sigbound = 3 * sqrt(squeeze(Sigma_hat_EKF(state, state, :)))' ./ scale(state);
    plot(t_k, x_3sigbound); hold on
    plot(t_k, error_EKF(state, :) ./ scale(state))
end
hold off
sgtitle("3-sigma Uncertainty Plots EKF")

figure
tiledlayout(3, 3)
scale = [ones(6, 1) * 1e-3; pi / 180 * ones(3, 1)];
error_ESKF = abs([x_true(1:6, :) - x_hat_ESKF(1:6, :); deg2rad(xi_theta_l_ESKF)]);
for state = 1:9
    nexttile
    x_3sigbound = 3 * sqrt(squeeze(Sigma_hat_ESKF(state, state, :)))' ./ scale(state);
    plot(t_k, x_3sigbound); hold on
    plot(t_k, error_ESKF(state, :) ./ scale(state))
end
hold off
sgtitle("3-sigma Uncertainty Plots ESKF")

figure
tiledlayout(3, 3)
scale = [ones(6, 1) * 1e-3; pi / 180 * ones(3, 1)];
error_IEKF = abs([quat_rot_array(x_true(7:10, :), x_hat_IEKF(1:3, :) - x_true(1:3, :)); quat_rot_array(x_true(7:10, :), x_hat_IEKF(4:6, :) - x_true(4:6, :)); deg2rad(xi_theta_l_IEKF)]);
Sigma_hat_IEKF_rearranged = Sigma_hat_IEKF([7, 8, 9, 4, 5, 6, 1, 2, 3], [7, 8, 9, 4, 5, 6, 1, 2, 3], :);
for state = 1:9
    nexttile
    x_3sigbound = 3 * sqrt(squeeze(Sigma_hat_IEKF_rearranged(state, state, :)))' ./ scale(state);
    plot(t_k, x_3sigbound); hold on
    plot(t_k, error_IEKF(state, :) ./ scale(state))
end
hold off
sgtitle("3-sigma Uncertainty Plots LIEKF (in body frame)")

figure
tiledlayout(3, 3)
scale = [ones(6, 1) * 1e-3; pi / 180 * ones(3, 1)];
error_IEKF = zeros([9, N]);
Sigma_hat_IEKF_r = Sigma_hat_IEKF;
for k = 1 : N
    error_IEKF(:, k) = abs([x_hat_IEKF(1:3, k) - X_hat_IEKF(k).R * X_trues(k).R' * x_true(1:3, k); x_hat_IEKF(4:6, k) - X_hat_IEKF(k).R * X_trues(k).R' * x_true(4:6, k); deg2rad(xi_theta_r_IEKF(:, k))]);

    Sigma_hat_IEKF_r(:, :, k) = X_hat_IEKF(k).Ad * Sigma_hat_IEKF(:, :, k) * X_hat_IEKF(k).Ad';
end
Sigma_hat_IEKF_rearranged = Sigma_hat_IEKF_r([7, 8, 9, 4, 5, 6, 1, 2, 3], [7, 8, 9, 4, 5, 6, 1, 2, 3], :);
for state = 1:9
    nexttile
    x_3sigbound = 3 * sqrt(squeeze(Sigma_hat_IEKF_rearranged(state, state, :)))' ./ scale(state);
    plot(t_k, x_3sigbound); hold on
    plot(t_k, error_IEKF(state, :) ./ scale(state))
end
hold off
sgtitle("3-sigma Uncertainty Plots LIEKF (in inertial frame)")

%% Plot
% trajectories

figure
plot3(x_prop(1, :) .* 1e3, x_prop(2, :) .* 1e3, x_prop(3, :) .* 1e3, DisplayName="True"); hold on
plot3(x_hat_EKF(1, :) .* 1e3, x_hat_EKF(2, :) .* 1e3, x_hat_EKF(3, :) .* 1e3, DisplayName="EKF");
plot3(x_hat_ESKF(1, :) .* 1e3, x_hat_ESKF(2, :) .* 1e3, x_hat_ESKF(3, :) .* 1e3, DisplayName="ESKF"); 
plot3(x_hat_IEKF(1, :) .* 1e3, x_hat_IEKF(2, :) .* 1e3, x_hat_IEKF(3, :) .* 1e3, DisplayName="IEKF"); 
xlim([-10, 60] * 1e-0)
ylim([-35, 35] * 1e-0)
zlim([0, 70] * 1e-0)
hold off
grid on
legend()
% all the uncertainties

% all the errors

%% Helper Functions
% Stuff to help plotting funky bananas and what not

function [X] = posquat2SE2_3_Rotm(x)
    q = x(7:10);
    X = SE2_3_RotationMatrix(quat_rotmatrix(q), x(4:6), x(1:3));
end

function [twist] = posquat2twist(x_a, xdot)
    twist.w = x_a(11:13);
    twist.a = xdot(4:6); 
    twist.v = xdot(1:3);
end

function [jac] = jacobian_func(f, p, wrt)
    t_sym = sym("t");
    x_sym = sym("x", [nx, 1]);
    u_sym = sym("u", [nu, 1]);
    
    if strcmp(wrt, "x")
        jac_func = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), x_sym(var_ind)),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
        jac = @(t, x, u) jac_func(t, x, u, p);
    elseif strcmp(wrt, "u")
        jac_func = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), u_sym(var_ind)),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
        jac = @(t, x, u) jac_func(t, x, u, p);
    end
end

function [R] = quat_rotmatrix(q)
    w = q(4);
    v = q(1:3);
    % x = v(1);
    % y = v(2);
    % z = v(3);
    % R = [w ^ 2 + x ^ 2 - y ^ 2 - z ^ 2, 2 * (x * y - w * z), 2 * (x * z + w * y);
    %      2 * (x * y + w * z), w ^ 2 - x ^ 2 + y ^ 2 - z ^ 2, 2 * (x * z - w * x);
    %      2 * (x * z - w * y), 2 * (y * z + w * x), w ^ 2 - x ^ 2 - y ^ 2 + z ^ 2];

    R = (w ^ 2 - v' * v) * eye(3) + 2 * v * v' + 2 * w * skew(v);
    %R = eye(3) + 2 * skew(v) * (w * eye(3) + skew(v));
    % R = [1 - 2 * y ^ 2 - 2 * z ^ 2, 2 * (x * y + z * w), 2 * (x * z - y * w);
    %      2 * (x * y - z * w), 1 - 2 * x ^ 2 - 2 * z ^ 2, 2 * (y * z + x * w);
    %      2 * (x * z + y * w), 2 * (y * z - x * w), 1 - 2 * x ^ 2 - 2 * y ^ 2];
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

function [xi_theta_l] = quat_left_error(q_true, q_est)
    % q is from inertial to body frame
    xi_theta_l = qLog(q_mul_array(q_true, q_conj(q_est)));
end

function [xi_theta_r] = quat_right_error(q_est, q_true)
    % q is from inertial to body frame
    xi_theta_r = qLog(q_mul_array(q_conj(q_est), q_true));
end

function [A_rot] = A_rot_euler212()
    b_inverse = @(theta) (1/sin(theta(2))) .* [0 sin(theta(3)) cos(theta(3));
                         0 sin(theta(2))*cos(theta(3)) -sin(theta(2))*sin(theta(3));
                         sin(theta(2)) -cos(theta(2))*sin(theta(3)) -cos(theta(2))*cos(theta(3))];
    A_rot = @(theta, w) -b_inverse(theta) * w;
end

function [rot_func] = rotation_euler212()
    Rxtheta1 = @(theta) [1 0 0; 0 cos(theta(1)) sin(theta(1)); 0 -sin(theta(1)) cos(theta(1))];
    Rxtheta2 = @(theta) [cos(theta(2)) 0 -sin(theta(2)); 0 1 0; sin(theta(2)) 0 cos(theta(2))];
    Rxtheta3 = @(theta) [1 0 0; 0 cos(theta(3)) sin(theta(3)); 0 -sin(theta(3)) cos(theta(3))];
 
    C_eb = @(theta) (Rxtheta3(theta) * Rxtheta2(theta) * Rxtheta1(theta));

    rot_func = @(theta, p) C_eb(theta) * p;
end

function [A, G] = A_G_EKF_euler212()
    % For position - velocity - euler angle state representation
    b_inverse = @(theta) (1/sin(theta(2))) .* [0 sin(theta(3)) cos(theta(3));
                         0 sin(theta(2))*cos(theta(3)) -sin(theta(2))*sin(theta(3));
                         sin(theta(2)) -cos(theta(2))*sin(theta(3)) -cos(theta(2))*cos(theta(3))];

    Rxtheta1 = @(theta) [1 0 0; 0 cos(theta(1)) sin(theta(1)); 0 -sin(theta(1)) cos(theta(1))];
    Rxtheta2 = @(theta) [cos(theta(2)) 0 -sin(theta(2)); 0 1 0; sin(theta(2)) 0 cos(theta(2))];
    Rxtheta3 = @(theta) [1 0 0; 0 cos(theta(3)) sin(theta(3)); 0 -sin(theta(3)) cos(theta(3))];
 
    C_eb = @(theta) (Rxtheta3(theta) * Rxtheta2(theta) * Rxtheta1(theta));
    R = @(theta) C_eb(theta);
    S = @(theta) -b_inverse(theta);

    theta_sym = sym("theta", [3, 1]);
    u_w = sym("u_w", [3, 1]);
    u_a = sym("u_a", [3, 1]);
    S_jac_func = matlabFunction(jacobian(S(theta_sym) * u_w, theta_sym),"Vars", [{theta_sym}, {u_w}]);
    S_jac = @(theta, u) S_jac_func(theta, u);

    R_jac_func = matlabFunction(jacobian(R(theta_sym) * u_a, theta_sym),"Vars", [{theta_sym}, {u_a}]);
    R_jac = @(theta, u) R_jac_func(theta, u);

    z33 = zeros(3, 3);
    A = @(x, u) [z33, eye(3), z33;
                 R_jac(x(7:9), u(1:3)), z33, z33;
                 S_jac(x(7:9), u(4:6)), z33, z33];

    % Don't know why G is negative in reference
    G = @(x, u) -blkdiag(R(x(7:9)), R(x(7:9)), S(x(7:9)));
end

function [A, G] = A_G_RESKF()
    % For right attitude error - attitude deviation in inertial frame
    % error state expressed as xi_i_theta = Log(q * q_inv)
    % For position - veloicty - quaternion state representation
 
    z33 = zeros(3, 3);
    A = @(x, u) [z33, eye(3), z33;
                 zeros(3, 6), -skew(quat_rot(x(7:10), u(1:3)));
                 zeros(3, 9)];

    G = @(x, u) kron(eye(3), quat_rotmatrix(x(7:10)));
end

function [A, G] = A_G_LESKF()
    % For left attitude error - attitude deviation in body frame
    % error state expressed as xi_b_theta = Log(q_inv * q)
    % For position - veloicty - quaternion state representation
    
    z33 = zeros(3, 3);
    A = @(x, u) [z33, eye(3), z33;
                 zeros(3, 6), -quat_rotmatrix(x(7:10)) * skew(u(1:3));
                 zeros(3, 6), -skew(u(4:6))];

    G = @(x, u) blkdiag(kron(eye(2), quat_rotmatrix(x(7:10))), eye(3));
end

function [eulerstate] = quatstate2eulerstate(quatstate)
    eulerstate = quatstate(1:9, :);

    quat = quatstate(7:10, :);
    for k = 1 : size(quatstate, 2)
        q = quat(:, k);
        [theta1, theta2, theta3] = dcm2angle(quat_rotmatrix(q), "XYX");
        eulerstate(7:9, k) = [theta1; theta2; theta3];
    end
end

function [h] = gps_h(x)
    h = x(1:3);
end

function [H] = gps_H_EKF()
    H = [eye(3), zeros(3, 6)];
end

function [h] = baro_h(x)
    h = x(3);
end

function [H] = baro_H_EKF()
    H = [zeros(1, 2), ones(1, 1), zeros(1, 6)];
end

function [h] = mag_h(x, B, rot_func)
    h = rot_func(x, B);
end

function [H_func] = mag_H_EKF_func()
    Rxtheta1 = @(theta) [1 0 0; 0 cos(theta(1)) sin(theta(1)); 0 -sin(theta(1)) cos(theta(1))];
    Rxtheta2 = @(theta) [cos(theta(2)) 0 -sin(theta(2)); 0 1 0; sin(theta(2)) 0 cos(theta(2))];
    Rxtheta3 = @(theta) [1 0 0; 0 cos(theta(3)) sin(theta(3)); 0 -sin(theta(3)) cos(theta(3))];
 
    C_eb = @(theta) (Rxtheta3(theta) * Rxtheta2(theta) * Rxtheta1(theta))';
    R = @(theta) C_eb(theta);

    theta_sym = sym("theta", [3, 1]);
    B_sym = sym("B", [3, 1]);
    R_jac_func = matlabFunction(jacobian(R(theta_sym)' * B_sym, theta_sym),"Vars", [{theta_sym}, {B_sym}]);
    R_jac = @(theta, B) R_jac_func(theta, B);

    H_func = @(x, B) [zeros(3, 6), R_jac(x(7:9), B)];
end

function [H] = mag_H_LESKF(x, B)
    % q must be from inertial to body
    H = [zeros(3, 6), -skew(quat_rot(q_conj(x(7:10)), B))];
end

function [H] = mag2_H_LESKF(x, b, B)
    % q must be from inertial to body
    H = [zeros(3, 6), -2 * b / norm(b) * quat_rot(q_conj(x(7:10)), B / norm(B))'];
end

function [H] = mag3_H_LESKF(x, b, B)
    % q must be from inertial to body
    H = [zeros(3, 6), -2 * skew(b) * skew(quat_rot(q_conj(x(7:10)), B))];
end

function [H] = mag_H_RESKF(x, B)
    % q must be from inertial to body
    H = [zeros(3, 6), -quat_rotmatrix(x(7:10))' * skew(B)];
end

function [x_composed] = LESKF_add_posquat(x1, x2)
    x_composed = [x1(1:6) + x2(1:6); q_mul(x1(7:10), qexp(x2(7:9)))];
end

function [q] = qexp(tau)
    theta = norm(tau);
    u = tau / theta;
    q = [u * sin(theta / 2); cos(theta / 2)];
end

function [] = plot_EKF_cov()
    %% States
    titles = ["X Position", "Y Position", "Z Position", "X Velocity", "Y Velocity", "Z Velocity", "Theta 1", "Theta 2", "Theta 3", "X Angular Velocity", "Y Angular Velocity", "Z Angular Velocity", "Mass"];
    ylabels = ["r_x [km]", "r_y [km]", "r_z [km]", "v_x [m / s]", "v_y [m / s]", "v_z [m / s]", "\theta_1 [deg]", "\theta_2 [deg]", "\theta_3 [deg]", "\omega_x [deg / s]", "\omega_y [deg / s]", "\omega_z [deg / s]", "m [kg]"];
    
    ops = {@(x) x, @(x) x, @(x) x, @(x) 1000 * x, @(x) 1000 * x, @(x) 1000 * x, @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) exp(x)};
    
    for x = 1:13
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


