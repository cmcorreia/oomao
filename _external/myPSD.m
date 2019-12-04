%{
------------HEADER-----------------
Objective         ::  Compute and plot the one-sided PSD of a signal, knowing the sampling rate

Comments:             

INPUT VARS
   y               :: The sampled data sequence
   Fs              :: Sampling frequency in Hz

OUTPUT VARS
   Y               :: The PSD in units^2/Hz
   f               :: The frequency vector
Created by         :: C. Correia
Creation date      :: 16/10/2010

Change Record:     ::  
------------HEADER END----------------
%}

function [Y,f] = myPSD(y,Fs,doplot)

if nargin < 3
    doplot = 1;
end
if size(y,2) > size(y,1);
    y = y';
end
y = y - mean(y);                       %the sequence must be zero-mean

L = length(y);
NFFT = L;%2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
df = gradient(f);

Y = Y.*conj(Y);
Y = 2*abs(Y(1:NFFT/2+1))./df'*L/NFFT;

% Plot single-sided amplitude spectrum.
if doplot
    loglog(f,Y) 
    title('Single-Sided Amplitude Spectrum of y(t)','fontsize',16)
    xlabel('Temporal frequency [Hz]','fontsize',14)
    ylabel('|Y(f)|^2/Hz','fontsize',14)
    box on
    axis tight
end
