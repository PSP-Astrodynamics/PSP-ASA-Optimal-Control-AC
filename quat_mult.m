function[quat_prod] = quat_mult(test_quat, true_quat)

% MULTIPLY TRUE QUATERNION MATRICES- TEST FROM OUR ESTIMATOR AND REAL ONE

% EXTRACT SCALAR AND VECTOR 
w_test = test_quat(1,:);
w_true = true_quat(1,:); 
v_test = [test_quat(2,:); test_quat(3,:); test_quat(4,:)]; 
v_true = [true_quat(2,:); true_quat(3,:); true_quat(4,:)]; 


% quaternion multiplication 
quat_prod = [(w_test .* w_true) - (dot(v_test, v_true)); ...
             (w_test .* v_true) + (w_true.*v_test)+ cross(v_test,v_true)];



