function[error_pv] = pv_error(test_vec, true_vec)

% OBTAIN DIFFERENCE BETWEEN FINAL AND INITIAL POSITION 
error_vec = test_vec - true_vec; 

% obtain magnitude of position 
error_pv = norm(error_vec); 
