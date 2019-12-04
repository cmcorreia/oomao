function Y = heaviside(X)
%HEAVISIDE    Step function.
%    HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and .5 for X == 0.
%    HEAVISIDE(X) is not a function in the strict sense.
%    See also DIRAC.

%   Copyright 1993-2012 The MathWorks, Inc.

Y = zeros(size(X));
Y(X > 0) = 1;
Y(X == 0) = .5;