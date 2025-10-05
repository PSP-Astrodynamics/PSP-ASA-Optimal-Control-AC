function [A_proj] = project_ellipsoid(A_array, proj_plane)
%PROJECT_ELLIPSOID Summary of this function goes here
%   Algorithm from: 
% https://math.stackexchange.com/questions/2431159/how-to-obtain-the-equation-of-the-projection-shadow-of-an-ellipsoid-into-2d-plan

N = size(A_array, 3);
M = size(A_array, 4);

A_proj = zeros([numel(proj_plane), numel(proj_plane), N]);

for j = 1:M
    for i = 1:N
        A = A_array(:, :, i, j);
    
        J_mask = logical(blkdiag(ones(numel(proj_plane)), zeros(size(A, 1) - numel(proj_plane))));
        J_mask_shifted = circshift(J_mask,[proj_plane(1) + size(A_array, 1) - 1,proj_plane(1) + size(A_array, 1) - 1]);
        
        J = reshape(A(J_mask_shifted), [numel(proj_plane), numel(proj_plane)]);
        
        K_mask = logical(blkdiag(zeros(numel(proj_plane)), ones(size(A, 1) - numel(proj_plane))));
        
        K_mask_shifted = circshift(K_mask,[proj_plane(1) + size(A_array, 1) - 1,proj_plane(1) + size(A_array, 1) - 1]);
        
        K = reshape(A(K_mask_shifted), [(size(A, 1) - numel(proj_plane)), (size(A, 1) - numel(proj_plane))]);
        
        L_mask = zeros(size(A));
        L_mask((numel(proj_plane) + 1):end, 1:numel(proj_plane)) = 1;
        
        L_mask_shifted = logical(circshift(L_mask,[proj_plane(1) + size(A_array, 1) - 1, proj_plane(1) + size(A_array, 1) - 1]));
        
        L = reshape(A(L_mask_shifted), [(size(A, 1) - numel(proj_plane)), numel(proj_plane)]);
        
        A_proj(:, :, i, j) = J - L' * inv(K) * L;
    end
end
end

