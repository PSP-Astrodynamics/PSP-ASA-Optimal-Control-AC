simout = Get_Sim_Output(out);

num_timestep = size(simout.utraj, 1);
num_optim = size(simout.utraj, 3);
op_time = 0.04;
dtime = linspace(0.2, 0.2, num_optim);
offset = op_time * (0:num_optim-1);
scaled = (dtime' * (0:num_timestep-1))';
solution_time = offset + scaled;
min_sol_time = min(solution_time,[],"all");
max_sol_time = max(solution_time,[],"all");
min_dtime = min(dtime,[],"all");
min_optime = min(op_time,[],"all");
k = linspace(min_sol_time, max_sol_time, num_optim);

q = linspace(0, max(simout.time), num_optim);
c = min_sol_time:op_time:max_sol_time;

[X, Y] = meshgrid(q,c);

scaled_values = zeros(numel(c), numel(q));
for i = 1:length(q) 
    scaled_values(:, i) = interp1(solution_time(:,i), simout.xtraj(:,2,i), c, 'linear', 0); 
end

scaled_first_axis = zeros(numel(c),numel(q));
for i = 1:length(q) 
    scaled_first_axis(:, i) = interp1(solution_time(:,i), simout.xtraj(:,1,i), c, 'linear', 0); 
end

scaled_third_axis = zeros(numel(c), numel(q));
for i = 1:length(q) 
    scaled_third_axis(:, i) = interp1(solution_time(:,i), simout.xtraj(:,3,i), c, 'linear', 0); 
end


% plots

% height
waterfall(X, Y, scaled_values);
xlabel("Simulation Time [s]")
ylabel("Solution Timestep [s]")
zlabel("Height [m]")
title("Height Trajectory vs Time [Fixed]")


% crossrange trajectory
waterfall(X, Y, scaled_first_axis);
% not sure if it should be transposed or not
ylabel("Simulation Time [s]")
xlabel("Solution Timestep [s]")
zlabel("Crossrange [m]")
title("Crossrange Trajectory vs Time [Fixed]")

% throttle
waterfall(X, Y, scaled_first_axis * 100);
ylabel("Simulation Time [s]")
xlabel("Solution Timestep")
zlabel("Throtte [%]")
title("Throttle Trajectory vs Time")

view([61.4 29.4])