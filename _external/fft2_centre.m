function out = fft2_centre(x,dir,n1,n2)
%{
------------HEADER-----------------
Objective         ::  Compute centred 2-dimensional FFT

Comments:             

INPUT VARS
   x               :: The 2D object 
   dir             :: 1: direct;  2: inverse FFT
   n1, n2          :: dimensions (padding)

OUTPUT VARS
   out             :: centred FFT
Created by         :: C. Correia
Creation date      :: 02/11/2012

Change Record:     ::  
------------HEADER END----------------
%}

if nargin<3
    [n1,n2] = size(x);
else
    try
    [n1o,n2o] = size(x);
    x = padarray(x,[(n1-n1o)/2,(n2-n2o)/2],'both');
    %disp('Info.: centred padding of input function possible...')
    catch err
       disp('Info.: centred padding of input function not possible, check dimensions...')
        x(n1,n2) = 0;
    end
    [n1,n2] = size(x);
end

if dir == 1
    out = fftshift(fft2(fftshift(x)));
else
    out = fftshift(ifft2(fftshift(x)));
end