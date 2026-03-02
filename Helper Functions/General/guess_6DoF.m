function [guess] = guess_6DoF(x_0, x_f, N, Nu, t_step, vehicle)
%STATE_GUESS Summary of this function goes here
%   Detailed explanation goes here

g = 9.81; % [m/s2]

% Control guess
T_guess = 0.4; % [N] Initial thrust guess
u_guess = [T_guess; 0; 0; T_guess; 0] .* ones([5, Nu]);

% Parameter guess
v_0 = x_0(4:6);
v_f = x_f(4:6);
%tf_guess = norm(v_f-v_0 + sqrt(2 * [0; g] * x_0(2)), 2)/(g - T_guess * vehicle.max_thrust/vehicle.m);

% State guess
x_guess = zeros([numel(x_f), N]);

% r_guess
x_guess(1:3, :) = straight_line_interpolate(x_0(1:3), x_f(1:3), N);

% v_guess
v_cst = (x_f(1:3) - x_0(1:3)) / (N * t_step);
x_guess(4:6, :) = repmat(v_cst, 1, N);

% CHANGE TO SPHERICALLY INTERPOLATE ANGLES AND ANGULAR VELOCITIES!!
% θ_guess  
if numel(x_f) == 12
    q_start = quaternion(x_0(7:9)', "euler", "XYX", "frame");
    q_end = quaternion(x_f(7:9)', "euler", "XYX", "frame");
    quat2 = slerp(q_start, q_end, linspace(0, 1, N));
    x_guess(7:9, :) = euler(quat2, "XYX", "frame")';
    %x_guess(7:9, :) = straight_line_interpolate(x_0(7:9), x_f(7:9), N);
    
    % ω_guess
    w_cst = (x_0(7:9) - x_f(7:9)) / (N * t_step);
    
    x_guess(10:12, :) = repmat(w_cst(1:numel(10:12)), 1, N);
else
    q_0 = x_0(7:10);
    q_f = x_f(7:10);
    x_guess(7:10, :) = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0] * slerp(quaternion(q_0([4, 1, 2, 3])'), quaternion(q_f([4, 1, 2, 3])'), linspace(0, 1, N)).compact';

    % Leave ω zero
end

guess.x = x_guess;
guess.u = u_guess;
guess.p = [];

end