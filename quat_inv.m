function[inv_quat] = quat_inv(q_vec)

% THIS FUNCTION OBTAINS THE INVERSE OF A QUATERNION VECTOR/ MATRIX

% OBTAIN CONJUGATE OF THE VECTOR
conj_q = [q_vec(1,:); (-1).*(q_vec(2,:)); (-1).*(q_vec(3,:)); (-1).*q_vec(4,:)]; 

% OBTAIN NORM 
norm_q = vecnorm(q_vec); 

% INVERSE
inv_quat = (conj_q) ./ (norm_q.^2); 
