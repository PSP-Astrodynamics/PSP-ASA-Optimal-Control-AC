function [A] = sqrtm_array(A)
%SQRTM_ARRAY Summary of this function goes here
%   Detailed explanation goes here
m = size(A, 3);

for i = 1:m
    A(:, :, i) = sqrtm(A(:, :, i));
end
end

