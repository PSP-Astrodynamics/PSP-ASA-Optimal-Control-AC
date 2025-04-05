% constants (same as 3DoF)
g = 9.81; % [m / s2]

syms mass L I1 I2 I3 max_thrust; 
c = [g; mass; L; I1; I2; I3];

% states

r = sym("r", [3;1]);
v = sym("v", [3,1]);
theta = sym("theta", [3;1]);
w = sym("w", [3;1]);

x = [r; v; theta; w];

% controls

norm_T = sym("T", [3;1]);
gamma = sym("gamma", 1);

u = [norm_T; gamma];

T = norm_T * max_thrust;

% calculate theta dot
b_inverse =  (1/sin(theta(2))) .* [0 sin(theta(3)) cos(theta(3));
     0 sin(theta(2))*cos(theta(3)) -sin(theta(2))*sin(theta(3));
     sin(theta(2)) -cos(theta(2))*sin(theta(3)) -cos(theta(2))*cos(theta(3))];
thetadot = b_inverse * w;

% r dot
r_dot = v;

Rxtheta1 = [1 0 0; 0 cos(theta(1)) sin(theta(1)); 0 -sin(theta(1)) cos(theta(1))];
Rxtheta2 = [cos(theta(2)) 0 -sin(theta(2)); 0 1 0; sin(theta(2)) 0 cos(theta(2))];
Rxtheta3 = [1 0 0; 0 cos(theta(3)) sin(theta(3)); 0 -sin(theta(3)) cos(theta(3))];
C_be = Rxtheta3 * Rxtheta2 * Rxtheta1;
C_eb = C_be.';
T_e = C_eb * T;
% v dot
v_dot = T_e / mass * cos(gamma) - [0; 0; g];
%
M = norm(T) .* sin(gamma) .* [1; 0; 0] + cross([-L; 0; 0], cos(gamma) * T);
% w dot
MOI = [I1; I2; I3];
w_dot = (M + (MOI([2; 3; 1]) - MOI([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ MOI([1; 2; 3;]);

x_dot = [r_dot; v_dot; thetadot; w_dot];
j_a = jacobian(x_dot, x);
j_b = jacobian(x_dot, u);

vars = [{x}; {u}; {mass; L; MOI; max_thrust}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","Dynamics Models/6DoF/SymDynamicsEuler6DoF","Vars",vars);

% Create equations of motion block for Simulink model
open("EoM_Euler6DoF.slx")
matlabFunctionBlock('EoM_Euler6DoF/SymDynamicsEuler6DoF',x_dot,'Vars',vars);
close_system("EoM_Euler6DoF.slx", 1)

% Create Jacobian functions for Kalman filter
matlabFunction(j_a,"File","Dynamics Models/6DoF/SymXJacobianEuler6DoF","Vars",vars);
matlabFunction(j_b,"File","Dynamics Models/6DoF/SymUJacobianEuler6DoF","Vars",vars);


%xdot_func6dof = matlabFunction(state);
%j_a_func6dof = matlabFunction(j_a);
%j_b_func6dof = matlabFunction(j_b);

%save("dynamics6dof.mat", "xdot_func6dof", "j_a_func6dof", "j_b_func6dof")
%save("SymDynamics6DoF.m", "xdot_func6dof")
%save("SymUJacobian6DoF.m", "j_b_func6dof")
%save("SymXJacobian6DoF.m", "j_a_func6dof")