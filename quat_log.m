function[orient_error] = quat_log(prod_vec)



% extract angle of rotation 
rot_0 = acos(prod_vec(1)) .* 2; 
s_0 = sin(rot_0 /2); % sine of rotation angle

% obtain rotation axis  
rot_ax = prod_vec(2:4) ./ s_0;

% error
orient_error = rot_0 .* rot_ax; 