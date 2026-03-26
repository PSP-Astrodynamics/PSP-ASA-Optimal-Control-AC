function[log_q, prod_vec] = quat_orien_err(test_vec, true_vec)

% obtain inverse of true position 
inv_true = quat_inv(true_vec); 

% product 
prod_vec = quat_mult(inv_true,test_vec); % error quat

% logarithmic mapping
log_q = quat_log(prod_vec); 


