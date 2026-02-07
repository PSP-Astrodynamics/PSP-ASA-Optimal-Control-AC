syms L m alpha g;
I = sym("I", [3; 1]);
c = [L; m; I; alpha; g];

% states

r = sym("r", [3;1]);
v = sym("v", [3,1]);
q = sym("theta", [4;1]);
w = sym("w", [3;1]);

x = [r; v; q; w];

% controls

T = sym("T", [3;1]);
tau = sym("tau", 1);

u = [T; tau];

thrust_mag = sqrt(T(1) ^ 2 + T(2) ^ 2 + T(3) ^ 2);

% calculate qdot
qdot = 1 / 2 * q_mul(q, [w; 0]);

% r dot
r_dot = v;

T_e = quat_rot(q, T);
% v dot
v_dot = T_e / m - [0; 0; g];
%
M = tau .* [1; 0; 0] + cross([-L; 0; 0], T);
% w dot
w_dot = (M + (I([2; 3; 1]) - I([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I([1; 2; 3;]);

x_dot = [r_dot; v_dot; qdot; w_dot];
%j_a = jacobian(x_dot, x);
%j_b = jacobian(x_dot, u);

vars = [{x}; {u}; {L; m; I; alpha; g}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","Dynamics Models/6DoF/SymDynamicsAstra.m","Vars",vars);

% % Create equations of motion block for Simulink model
% open("EoM_Euler6DoF.slx")
% matlabFunctionBlock('EoM_Euler6DoF/SymDynamicsEuler6DoF',x_dot,'Vars',vars);
% close_system("EoM_Euler6DoF.slx", 1)
% 
% % Create Jacobian functions for Kalman filter
% matlabFunction(j_a,"File","Dynamics Models/6DoF/SymXJacobianEuler6DoF","Vars",vars);
% matlabFunction(j_b,"File","Dynamics Models/6DoF/SymUJacobianEuler6DoF","Vars",vars);
% 

%xdot_func6dof = matlabFunction(state);
%j_a_func6dof = matlabFunction(j_a);
%j_b_func6dof = matlabFunction(j_b);

%save("dynamics6dof.mat", "xdot_func6dof", "j_a_func6dof", "j_b_func6dof")
%save("SymDynamics6DoF.m", "xdot_func6dof")
%save("SymUJacobian6DoF.m", "j_b_func6dof")
%save("SymXJacobian6DoF.m", "j_a_func6dof")