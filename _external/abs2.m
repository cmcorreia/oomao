function [output] = abs2(input)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if isreal(input)
    output = (abs(input)).^2;
else
    output = input .* conj(input);
end

end
