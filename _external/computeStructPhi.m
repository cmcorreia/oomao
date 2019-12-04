%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function

Comments:          :: Based on function compStructFunc.pro, form Jean-Pierre Veran

INPUT VARS
 phi               :: the phase matrix
 pup               :: pupil matrix (0-non illumunated, 1- illuminated), same size as phi

OUTPUT VARS
 dph               :: bi-dimensional phase structure function
Created by         :: C. Correia - Canada National Research Council
Creation date      :: 03/13/2011
                      
Change Record:     ::
------------HEADER END----------------
%}

function dph = computeStructPhi(phi,pup)

fftcorrel = @(x,y) ifft2(fft2(x).*conj(fft2(y)));

% Create phi and pup an a 2x larger grid
[nx ny] = size(phi);
phi2    = padarray(phi.*pup, [nx ny],0,'post');
pup2    = padarray(pup, [nx ny],0,'post');

dph     = zeros(nx,ny);
corp2w  = fftcorrel(phi2.^2,pup2);
corpp   = fftcorrel(phi2,phi2);
dph     = corp2w + fftsym(corp2w) - 2.*corpp;
corww   = round(abs(fftcorrel(pup2,pup2)));

% Compute mask of locations where there is overlap between shifted pupils and normilized by the number of points used to compute dphi
mask    = (corww > 0);
corww(corww <= 1) = 1;
dph     = dph .* mask ./ corww;

function b = fftsym(x)
[nx,ny]        = size(x);
b              = zeros(nx,ny);
b(1,1)         = x(1,1);
b(1,2:end)     = x(1,end:-1:2);
b(2:end,1)     = x(end:-1:2,1);
b(2:end,2:end) = x(end:-1:2,end:-1:2);



