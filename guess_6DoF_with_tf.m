function [x_guess, u_guess, tf_guess] = guess_6DoF_with_tf(x_0, x_f, steps, vehicle)
%STATE_GUESS Summary of this function goes here
%   Detailed explanation goes here

%x_initial (12, 1) double
%States: %States: x, y, z, x_dot, y_dot, z_dot, theta1, theta2, theta3, w1,
        %w2, w3. z is up
%x_final (12, 1) double = [0, 0, 0, 0, 0, 0, pi/2, 0, 0, 0, 0, 0]

%3dof states: % States: x, y, x_dot, y_dot, theta, theta_dot


g = 9.81; % [m/s2]

% Control guess
T_guess = 0.4; % [N] Initial thrust guess
u_guess = [0, 0, T_guess, 0] .* ones([steps, 4]);

% Parameter guess
v_0 = x_0(4:6);
v_f = x_f(4:6);
tf_guess = norm(v_f-v_0 + sqrt(2 * [0; 0; g] * x_0(3)), 2)/(g - T_guess * vehicle.max_thrust/vehicle.m);

% State guess
x_guess = zeros([steps, 12]);

% r_guess
x_guess(:, 1) = straight_line_interpolate(x_0(1), x_f(1), steps); %x
x_guess(:, 2) = straight_line_interpolate(x_0(2), x_f(2), steps); %y
x_guess(:, 3) = straight_line_interpolate(x_0(3), x_f(3), steps); %z

% v_guess
v_cst = (x_f(1:3) - x_0(1:3)) / tf_guess;
x_guess(:, 4:6) = repmat(v_cst, steps, 1); %check repmat argument

% θ1_guess 
x_guess(:, 7) = straight_line_interpolate(x_0(7), x_f(7), steps);

% θ2_guess 
x_guess(:, 8) = straight_line_interpolate(x_0(8), x_f(8), steps);

% θ3_guess 
x_guess(:, 9) = straight_line_interpolate(x_0(9), x_f(9), steps);

% ω1_guess
w1_cst = (x_0(10) - x_f(10)) / tf_guess;

% ω2_guess
w2_cst = (x_0(11) - x_f(11)) / tf_guess;

% ω3_guess
w3_cst = (x_0(12) - x_f(12)) / tf_guess;

x_guess(:, 10) = w1_cst;
x_guess(:, 11) = w2_cst;
x_guess(:, 12) = w3_cst;
end

