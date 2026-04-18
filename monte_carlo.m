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

std_vec = [];


new_problem = problem;

for i = 1 : 10
    new_problem.x0 = % new x0 with deviation
    [] = new_problem.cont_prop( , []);


end 