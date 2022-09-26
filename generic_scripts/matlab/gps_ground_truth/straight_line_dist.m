function [distance] = straight_line_dist(x0,y0,z0,x1,y1,z1)
%STRAIGHT_LINE_DIST Summary of this function goes here
%   Detailed explanation goes here
distance = sqrt((x0-x1)^2 + (y0-y1)^2 + (z0-z1)^2);
end

