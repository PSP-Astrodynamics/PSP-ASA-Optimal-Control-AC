function [Psi] = group_error(X, X_d)
%GROUP_ERROR Calculate group error from desired pose
%   X - current pose group element
%   X_d - desired pose group element

X_d_inv = X_d.inv;
Psi = X_d_inv.compose(X);

end

