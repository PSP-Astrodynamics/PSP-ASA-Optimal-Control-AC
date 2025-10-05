%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF rocket landing dynamics with changing mass
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%g = 3.7114e-3; % [km / s2]
g = 9.81e-3; % [km / s2]

syms L I alpha;
t = sym("t");
r = sym("r", [2, 1]);
v = sym("v", [2, 1]);
theta = sym("theta", 1);
w = sym("w", 1);
z = sym("z", 1);
x = [r;v;theta;w;z];

thrust_accel = sym("thrust_accel", [2, 1]); % Thrust over mass
thrust_accel_mag = sym("thrust_accel_mag", 1); % Thrust over mass
u = [thrust_accel; thrust_accel_mag];

tf = sym("tf");

p = [tf];

m = exp(z);
T = m * thrust_accel;

rdot = v;
thetadot = w;
M = cross([-L; 0; 0], [T; 0]);
wdot = M(3) / I;

rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
a_T_e = rotationMatrix * thrust_accel;
vdot = a_T_e - [0; g];

zdot = -alpha * thrust_accel_mag;

xdot = [rdot; vdot; thetadot; wdot; zdot] * tf;

%j_a = jacobian(xdot, x);
%j_b = jacobian(xdot, u);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","Dynamics Models/3DoF/SymDynamics3DoF_mass_convexified_freetf","Vars", [{t}; {x}; {u}; {tf}; {L; I; alpha}]);

% Create equations of motion block for Simulink model
%matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions for Kalman filter
%matlabFunction(j_a,"File","3DoF/SymXJacobian3DoF","Vars",[x; u; mass; L; I]);
%matlabFunction(j_b,"File","3DoF/SymUJacobian3DoF","Vars",[x; u; mass; L; I]);