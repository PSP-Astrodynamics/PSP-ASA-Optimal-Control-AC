function[error_mag] = quat_error(test_vec, true_vec)

% THIS FUNCTION ESTIMATES THE ERROR 

% obtain inverse of true value
inv_q = quat_inv(true_vec); 

% obtain product
quat_prod = quat_mult(test_vec, inv_q);

% create an identity vector
ideal_q = eye(4,1); 

% obtain error vec
error_quat = quat_prod - ideal_q; 

% obtain error magnitude 
error_mag = vecnorm(error_quat); 

