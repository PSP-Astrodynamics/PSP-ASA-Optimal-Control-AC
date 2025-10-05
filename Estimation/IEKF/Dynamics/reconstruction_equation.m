function [Xdot] = reconstruction_equation(G, X, twist)
%RECONSTRUCTION_EQUATION Group element differential equation
%   Differential equation on group used to reconstruct the group element
%   trajectory from the twist trajectory. Also works for reconstructing
%   desired pose trajectory X_d.
%   G - group
%   X - group element
%   twist - twist (velocity in body frame)

Xdot = X.compose(G.hat(twist)); 

end

