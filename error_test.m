%%%%%
% Test for 

% initialize
glideslope_angle_max = deg2rad(65);


% 
test_1 = load("Results_ptr2.mat").results; 
test_2 = load("Results_Time_Ptr.mat").results; 

% extract quaternions
quat_1 = test_1.x(7:10, :);
quat_2 = test_2.x(7:10, :); 


error_mag = quat_error(quat_1, quat_2);
error_mag = rad2deg(error_mag); 

figure
plot(test_1.t, error_mag,'b')

figure
comparison_plot_6DoF_trajectory({test_1.x, test_2.x}, ["Test 1", "Test 2"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Test 1 and 2")

figure
comparison_plot_6DoFq_time_histories({test_1.t, test_2.t}, {test_1.x, test_2.x}, {test_1.u, test_2.u}, ["Test 1", "Test 2"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")
