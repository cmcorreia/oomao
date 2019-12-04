function [Nj, Mj] = make_zernike_index(Nmodes)
% ------------HEADER-----------------
% Objective:
%   Output the radial order and azimuthal frequency indexes for up to a
%   given number of Nmodes
%INPUT:  
%   Nmodes       ::  The number of modes
% 
%OUTPUT:
%   Nj           ::  The radial order
%   Mj           ::  The azimuthal frequency
%
%Created by      ::  Carlos Correia 
%Creation date   ::  15 Mar 04
%Change Record:  ::  
% ------------HEADER END-----------------
Nj = 0;
n = 1;
while size(Nj,2) < Nmodes +1
    Nj = [Nj n*ones(1,n+1)];
    n = n+1;
end

even = [0 2 2 4 4 6 6 8 8 10 10 12 12 14 14 16 16 18 18 20 20 22 22 24 24 26 26 28 28 30 30 32 32 34 34 36 36 38 38 40 40 42 42 44 44 46 46 48 48 50 50 52 52 54 54 56 56 58 58 60 60 ];
odd  = [1 1 3 3 5 5 7 7 9 9 11 11 13 13 15 15 17 17 19 19 21 21 23 23 25 25 27 27 29 29 31 31 33 33 35 35 37 37 39 39 41 41 43 43 45 45 47 47 49 49 51 51 53 53 55 55 57 57 59 59 61 61];

Mj = 0;
isodd = 1;
n = 1;
while size(Mj,2) < Nmodes +1
    if isodd
        Mj = [Mj odd(1:n+1)];
        isodd = 0;
    else
        Mj = [Mj even(1:n+1)];
        isodd = 1;
    end
    n = n+1;
end