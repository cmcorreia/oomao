function out = graduatedHeaviside(x,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%size(x)

dx = x(1,2)-x(1,1);
if dx==0
    dx = x(2,1)-x(1,1);
end

out = zeros(size(x));
out(-dx*n<=x) = 0.5;
out(x>dx*n) = 1;

end

