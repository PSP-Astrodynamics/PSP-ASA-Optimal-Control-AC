function j_a = SymXJacobian3DoF(r1,r2,v1,v2,theta1,w1,thrust1,gimbal1,mass,L,I,max_thrust)
%SymXJacobian3DoF
%    J_A = SymXJacobian3DoF(R1,R2,V1,V2,THETA1,W1,THRUST1,GIMBAL1,MASS,L,I,MAX_THRUST)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    16-Dec-2024 22:34:53

t2 = gimbal1+theta1;
t3 = 1.0./mass;
j_a = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,-max_thrust.*t3.*thrust1.*sin(t2),max_thrust.*t3.*thrust1.*cos(t2),0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],[6,6]);
end
