function [dtwist] = Euler_Poincare(G, twist, u, J_b)
%EULER_POINCARE Dynamics of twists (velocity in body frame)
%   The dynamical equation for the group in the body frame. Dual to the
%   Newton-Euler equations which are for dynamics in the space or inertial
%   frame.
%   G - group
%   twist - twist (velocity in body frame)
%   u - control (in body frame)
%   J_b - generalized inertia matrix


% coadjoint action is always negative transpose of adjoint?
ad_star_twist = -G.ad(twist)';

dtwist = J_b \ (ad_star_twist * J_b * twist + u);

end

