
function [GM, PM] = mybode(nu,h_ol,margin)
% ------------HEADER-----------------
% Objective:
%   depict the bode plot
%INPUT:
%   nu           ::  frequency vector in Hz or m^-1 
%   h_ol         ::  open-loop TF
%   margin       ::  flag 1: compute gain and phase margins, 0 don't
%OUTPUT:
%
%Created by      ::  Carlos Correia 
%Creation date   ::  29 Jul 11
%Change Record:  :: 
%{
%*******************************************
%%        EXAMPLE
%*******************************************
T           = 0.001; Ts = T;
nu          = logspace(-3, log10(1/2/T),1000);
sNum        = sqrt(-1) * 2 * pi* nu;
TF_int = 1./(1-exp(-T*sNum));
TF_lag = exp(-1*T*sNum);
mybode(nu, 0.5*TF_int.*TF_lag.^2,1); % numerically evaluated integrator
% Compare to native matlab discrete-time TFs
z           = tf('z',T);
TFd_intI    = (1/(1-z^(-1)));
figure
margin(0.5*TFd_intI/z/z)
%COMMENT: GM AND PM ARE THE SAME HERE

% COMPARING INTEGRATORS AND MARGINS
figure
g1 = 1; 
sLap  = tf('s'); % Laplace 's' var
bode(g1/T/sLap, nu*2*pi) % cont-time integrator
hold on
bode(TFd_intI, nu*2*pi) % disc-time integrator
bode(tf([g1 0],[1 -1],T),nu*2*pi) % =c2d(1/s) 
mybode(nu, TF_int) % numerically evaluated integrator
fprintf('COMMENT:: The c2d and discrete-time controllers have same GAIN and PHASE as the numerically evaluated TFs.\n However the phase response differs from the cont-time integrator alone.\n')

% CLOSING THE LOOP
figure
bode(feedback(1,1/z/z*TFd_intI), nu*2*pi) % disc-time integrator
mybode(nu, 1./(1+TF_int.*TF_lag.^2)) % numerically evaluated integrator

%}
% ------------HEADER END-----------------

% parse input
if nargin <3
    margin = 0;
end
AMP_db  = 20*log10(abs(h_ol));
PHI_deg = unwrap(angle(h_ol))*180/pi;

res = abs(AMP_db);
idxP = find(res == min(res));
PM = PHI_deg(idxP) + 180;
if margin
    if nnz((PHI_deg - (-180) ) <0) > 0 %i.e. crossed the -180 level then need to eval GM
        
        res = abs(-180 - PHI_deg);
        %idxG = find(res == min(res));
        idxG = find(diff(res) == min(diff(res)));
        GM = -AMP_db(idxG);
        string = ['Bode plot. GM=' num2str(GM) 'dB(@' num2str(nu(idxG)) 'Hz), PM=' num2str(PM) 'deg (@' num2str(nu(idxP)) ,'Hz)'];
    else
        GM = inf;
        idxG = inf;
        string = ['Bode plot. GM=' num2str(GM) ', PM=' num2str(PM) 'deg (@' num2str(nu(idxP)) 'Hz)'];
    end
    disp(string)
else
    string= 'Bode diagram';
    idxG = inf;
end

figure
subplot(2,1,1)
semilogx(nu, AMP_db)
ylabel('Amplitude [dB]','fontsize',14)
title(string,'fontsize',14)
grid on
hold on
plot(nu, zeros(size(nu)),'k--')

plot(nu(idxP), AMP_db(idxP),'+')
if ~isinf(idxG) && margin
    plot(nu(idxG), AMP_db(idxG),'o')
end

subplot(2,1,2)
semilogx(nu, PHI_deg )
ylabel('phase [deg]','fontsize',14)
xlabel('temporal frequency [Hz]','fontsize',14)
grid on
hold on
plot(nu, -180*ones(size(nu)),'k--')
plot(nu(idxP), PHI_deg(idxP),'+')
if ~isinf(idxG)
plot(nu(idxG), PHI_deg(idxG),'o')
end


