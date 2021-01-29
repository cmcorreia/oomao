function [cx, cy] = cog(a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dx = numel(a(:,1));
dy = numel(a(1,:));
[X, Y] = meshgrid(1:dx, 1:dy);

cx = sum(sum(a.*X))/sum(sum(a)) - 1;  % consistant with centroiding (index 0 Npixel - 1)
cy = sum(sum(a.*Y))/sum(sum(a)) - 1;  % consistant with centroiding (index 0 Npixel - 1)

end