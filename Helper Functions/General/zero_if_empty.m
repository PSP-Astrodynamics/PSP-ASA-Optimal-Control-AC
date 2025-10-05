function [v] = zero_if_empty(v)
%ZERO_IF_EMPTY Summary of this function goes here
%   Detailed explanation goes here
new_size = size(v);
new_size(new_size == 0) = 1;
v = zeros(new_size);
end

