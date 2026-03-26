%%%%%
% Test for errors 


% 
test_1 = load("Results_ptr2.mat").results; 
test_2 = load("Results_Time_Ptr.mat").results; 

% extract quaternions
quat_1 = test_1.x(7:10, :);
quat_2 = test_2.x(7:10, :); 


error_mag = quat_error(quat_1, quat_2);
error_mag = rad2deg(error_mag); 

plot(test_1.t, error_mag,'b')
