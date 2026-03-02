function [x_interp] = straight_line_interpolate(initial, final, N)
%STRAIGHT_LINE_INTERPOLATE Summary of this function goes here
%   Detailed explanation goes here

    x_interp = linspace(0, 1, N) .* (final - initial) + initial;
end

