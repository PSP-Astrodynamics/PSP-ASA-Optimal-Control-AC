function[exp_quat] = quat_exp(angle_vec)

% exponentiating the log quat 


loop_n = size(angle_vec, 2); 
exp_quat = zeros(4,loop_n); 


for k = 1:loop_n
    norm_angle = norm(angle_vec(:,k));

    if norm_angle > (1e-10)
        direction = angle_vec(:, k) ./ norm_angle;

    else
        direction = zeros(3,1);
    end

    exp_quat(:,k) = [direction .* sin(norm_angle./2) ; cos(norm_angle./2)];
end

 
