%% monte carlo 


%load
test_1 = load("Results_ptr2.mat").results; 

problem = load("prob_6DoF.mat").prob_6DoF;

% std dev
position_std = [10; 10; 10] ./(1000);   % km 
velocity_std = [10; 10; 10] ./1000;     % km/s 
orien_std = deg2rad([3; 3; 3]);      % std of orientation in radians
w_velo_std = deg2rad([1; 1; 1]);         % angular velo in rad

mass_std = 20; 

std_vec = [position_std; velocity_std; orien_std; w_velo_std; mass_std];

new_problem = problem;

num = 100;       % no. of iterations
noise = std_vec .* randn(13,num);   % no mean 

% 
new_x0 = add_noise(test_1.x(:,1), noise*1e-0);

time = linspace(0, problem.tf, 1000); 
x = zeros(14,numel(time), num); 

x_cell = {}; 

glideslope_angle_max = deg2rad(65);
% 
% for i = 1 : num 
% 
%     new_problem.x0 = new_x0(:,i) ;% new x0 with deviation
%     [~,x(:,:,i),~] = new_problem.cont_prop(test_1.u , [], tspan= time);
%     x_cell(end + 1) = {x(:,:,i)};
% 
% end 
% 

for i = 1 : num 
    new_problem.x0 = new_x0(:,i); % new x0 with deviation
    [~,x,~] = new_problem.cont_prop(test_1.u , []);
    x_cell(end + 1) = {x};
end 
%% 
figure
comparison_plot_6DoF_trajectory(x_cell, strings(1,num), glideslope_angle_max, linestyle = strings(1,num)+"-", title = "Monte Carlo test")

% function for adding noise 

function[sum_noise] = add_noise(mean, noise)

sum_noise = zeros(14, size(noise,2)); 

sum_noise([1:6, 11:14] , :) = mean([1:6, 11:14]) + noise([1:6, 10:13] , :); 

sum_noise([7:10],:) = quat_mult(repmat(mean(7:10),1,size(noise,2)), quat_exp(noise(7:9,:))); 
end