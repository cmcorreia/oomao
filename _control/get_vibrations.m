function [outTilt, outTip] = get_vibrations(data, tel, atm, vib_central_freq, fMinRange, fMaxRange, fMaxRangeNoise)

%{
------------HEADER-----------------
Objective         ::  Compute the vibration identification model from Raven
                      OL WFS telemetry

INPUT VARS
                      
data              :: a data structure with 
                        - samples with nSlopes x nSamples telemetry data
                        - timestamp from which to compute the sampling freq
                        or, if algebraic, set the sampling freq directly
                      
tel               :: Telescope object
atm               :: atmosphere object
vib_central_freq  :: [2 x nVib] matrix with vibration central frequencies
fMinRange         :: lower frequency bound for the model fitting, default 0.1Hz
fMaxRange         :: upper frequency bound for the model fitting, default 50Hz
fMaxRangeNoise    :: noise frequency bound, frequency above which noise
                     \geq signal, default 50Hz


OUTPUT VARS       
outTilt, outTip    :: cells with 
                        1) PSD atmospheric model and telemetry
                        2) ... nVib+1) the 2nd-order vibration model

Created by         ::  Carlos Correia 
Creation date      ::  07/2014
Change Record:     ::  
------------HEADER END----------------
%}

%% Parse input
if ~exist('nSlopes','var')
    nSlopes = 144;      % number of valid slopes
end
nSamples = length(data.samples);

if numel(data.timestamp) > 1
    Ts =  double(data.timestamp(2) - data.timestamp(1))*1e-9;   % sample frequency in seconds
else
    Ts = data.timestamp;
end

f = linspace(1/(nSamples*Ts), 1/(2*Ts), nSamples/2);        % frequency vector in physical units (Hertz)

if ~exist('fMinRange','var')
    fMinRange = .1;         % lower frequency bound, user-defined
end
if ~exist('fMaxRange','var')
    fMaxRange = 50;         % upper frequency bound, user-defined
end
if ~exist('fMaxRangeNoise','var')
    fMaxRangeNoise = 50;    % noise frequency bound, frequency above which noise \geq signal
end


if ~exist('vib_central_freq','var')
    vib_central_freq = [5.52, 11.03, 0.68; 5.52, 6.143, 0.4725];
end
nVib = size(vib_central_freq,2);
whichModel = 'f';
if ~exist('whichModel','var')
    whichModel = 'AoA';     
end
                        % which is the underlying signal model, 
                        % 'AoA' - angle of arrival
                        % 'power law' - a power law \propto f^-2 (default)
                        % Zernike Tip/Tilt - not implemented yet

z = tf('z',Ts);         % Z-transform var

verbose = 1;            % display outputs

%% get averaged TT modes time sequence
if size(data.samples,1) > 1 %&& ~any(size(data.samples)) == 1
    Tilt = mean(data.samples(1:nSlopes/2,:));
    Tip = mean(data.samples(nSlopes/2+1:end,:));
end

% ccorreia overridden to accept a sequence of TT directly (Feb'16)
    Tilt = data.samples;
    Tip = Tilt;
    
%% Temporal PSD
PSDTilt = abs(fft(Tilt)).^2/length(Tilt)^2;
PSDTilt = PSDTilt/trapz(f, PSDTilt(1:end/2))*sum(PSDTilt);  % normalisation such that integral over time is variance of time sequence (Parseval theorem)
%trapz(f,PSDTilt(1:end/2)) % = var(Tilt)

subplot(2,1,1)
loglog(f, PSDTilt(1:nSamples/2))
hold on
legend('Tilt')
xlabel('temporal Frequency [Hz]')
ylabel('Arbitrary units')
title(' Tip and Tilt temporal PSDs from telemetry')
axis tight

PSDTip = abs(fft(Tip)).^2/length(Tilt)^2;                   
PSDTip = PSDTip/trapz(f, PSDTip(1:end/2))*sum(PSDTip);      % normalisation such that integral over time is variance of time sequence (Parseval theorem)
%trapz(f,PSDTip(1:end/2))% = var(Tip)

subplot(2,1,2)
loglog(f, PSDTip(1:nSamples/2),'r')
hold on
legend('Tip')
xlabel('temporal Frequency [Hz]')
ylabel('Arbitrary units')
title(' Tip and Tilt temporal PSDs from Raven telemetry')
axis tight

%% Model Noise
    fimaxNoise = find(f >= fMaxRangeNoise);
    Nfloor(1) = mean(PSDTilt(fimaxNoise(1):nSamples/2));
    Nfloor(2) = mean(PSDTip(fimaxNoise(1):nSamples/2));

%% Model PSD   

if strcmp(whichModel, 'AoA')
    
    % Angle-of-arrival temporal PSD
    %Parameters conversion for local variables
    
    
    lambda          = 500e-9;%
    arcsectorad     = 1./3600./180.*pi;
    %seeingref       = 1.00; %arcsec
    %r0              = 0.978*lambda/(seeingref*arcsectorad);
    cn2dh           = (atm.r0)^(-5./3.)*1/0.423*(lambda/2./pi)^2;
    
    params.D        = tel.D;        % telescope diameter
    params.r0       = atm.r0;       % coherence length
    %params.seeing   = seeingref;
    params.cn2dh    = cn2dh*[1];    % Turbulence Profile (Cn2*dh) [m^{1/3]
    params.lambda   = lambda;       % wavelength
    params.R        = params.D/2;   % telescope radius %deprecated in newer versions
    params.V        = [atm.layer.windSpeed];  % wind speed in m/s
    params.dky      = 0.05;         % frequency step for integration in perpendicular direction to the wind (if linear scale)
    params.dnu      = 0.05;         % frequency step for integration in paralel direction to the wind (if linear scale)
    params.Fmax     = f(end);       % Max integration frequency in 1/m
    params.Fmin     = f(1);         % Min integration frequency in 1/m
    params.L0       = [atm.L0];     % outer-scale of turbulence in meters (m)
    params.Nmodes   = 3;            % Number of modes; 1 for piston, 2 for piston and tip, 3 for piston, tip and tilt etc...
    params.Nlayers  = 1;            % Number of layers on which the temporal DSP are to be computed
    params.thetaW   = atm.layer.windDirection;  % wind direction
    params.type     = 'aoa';        % 'Zernike' for Zernikes, 'aoa' for angle-of-arrival
    
    try
       load AoAPSD.mat
       fprintf('ATTENTION: READING FILE FROM DISK. MAKE SURE PARAMS ARE CONSISTENT\n')
    catch
        [Faoa,nu] = make_multi_layered_tempo_spectra(params);
    end
    %load('temporalPSD_3layers.mat')
    %Faoa = Faoa*params.D^2/16; % to convert from angle in sky to displacement at the borders of the pupil
    PSDaoaX = zeros(1,length(Faoa));
    PSDaoaY = PSDaoaX;
    for l = 1:params.Nlayers
        PSDaoaX = PSDaoaX + cos(params.thetaW(l))^2*Faoa(1,:,l) + sin(params.thetaW(l))^2*Faoa(2,:,l);
        PSDaoaY = PSDaoaY + cos(pi/2-params.thetaW(l))^2*Faoa(1,:,l) + sin(pi/2-params.thetaW(l))^2*Faoa(2,:,l);
    end
    %loglog(nu, PSDaoa,'k--')
    
    Var_aoa_th = 0.97*6.88*params.r0^(-5/3)*params.D^(-1/3)*params.D^2;% Expression from Tallon, PhD Thesis, pg 188. Factor D^2 at the end is to convert back to displacement at the edges of the pupil
    AoA_var = sum(PSDaoaX.*nu*(log(nu(2)) - log(nu(1))));
    rad2toskyarcsec2_aoa =  (lambda/2./pi)^2./(arcsectorad)^2;
    AoA_var_arcsec2 = AoA_var * rad2toskyarcsec2_aoa;

    if verbose
        fprintf('AoA theoretical value: %1.4f arcseconds^2\n',AoA_var_arcsec2)
    end
    %to have AoA_var in Z2,3 equivalent, divide AoA_var_arcsec/rad2toskyarcsec2_zer
    
    modelPSDX = interp1(nu,PSDaoaX, f);
    modelPSDX(end) = 0;
    modelPSDY = interp1(nu,PSDaoaY, f);
    modelPSDY(end) = 0;
    fimin = find(f> fMinRange); fimin = fimin(1);
    fimax = find(f <= fMaxRange); fimax = fimax(end);


    
    %p = fminsearch(@(p) AoA_fit_parms(modelPSDX, modelPSDY, PSDTilt, PSDTip, p, f, fMinRange, fMaxRange), [15 10]);

else
    fimin = find(f> fMinRange); fimin = fimin(1);
    fimax = find(f <= fMaxRange); fimax = fimax(end);

    seq1 = PSDTilt(fimin:fimax);    
    p = fminsearch(@(p) fitPowerLaw(seq1-Nfloor(1),f(fimin:fimax),p), [-3 0]);
    modelPSDX = p(2)*f.^(p(1))+Nfloor(1);
    
    seq1 = PSDTip(fimin:fimax);
    p = fminsearch(@(p) fitPowerLaw(seq1-Nfloor(2),f(fimin:fimax),p), [-2 1]);
    modelPSDY = p(2)*f.^(p(1))+Nfloor(2);;
    
    modelPSDX = f.^(-2);
    modelPSDY = f.^(-2);
    
end


nMode = 2;
for kMode = 1:nMode
    if kMode == 1
        PSDmode = PSDTilt;
        seq1 = PSDmode(fimin:fimax);
        seq2 = modelPSDX(fimin:fimax);
        k(kMode) = fminbnd(@(k) sum(abs(seq1 - k*seq2).^2), 1e-10,1e4,optimset('TolX',1e-13,'Display','off'));
        modelPSDstruct{kMode} = k(kMode)*modelPSDX;
        subplot(2,1,kMode)
        loglog(f, modelPSDstruct{kMode},'k--')
        
    elseif kMode == 2
        PSDmode = PSDTip;
        seq1 = PSDmode(fimin:fimax);
        seq2 = modelPSDY(fimin:fimax);
        k(kMode) = fminbnd(@(k) sum(abs(seq1 - k*seq2).^2), 1e-10,1e4,optimset('TolX',1e-13,'Display','off'));
        modelPSDstruct{kMode} = k(kMode)*modelPSDY;
        subplot(2,1,kMode)
        loglog(f, modelPSDstruct{kMode},'k--')
    end
end

if verbose
    p = fminsearch(@(p) AoA_fit_parms(modelPSDX*k(1), modelPSDY*k(2), PSDTilt, PSDTip, p, f, fMinRange, fMaxRange), [15 10]);
    fprintf('Fitted windSpeed and windDirection: %f m/s, %f deg\n',p(1), p(2))
end

% save model and telemetry PSD
outTilt{1}.PSD_ATM_model = modelPSDstruct{1};
outTilt{1}.PSD_Telemetry = PSDTilt;
outTilt{1}. f = f;

outTip{1}.PSD_ATM_model = modelPSDstruct{2};
outTip{1}.PSD_Telemetry = PSDTip;
outTip{1}. f = f;

%% user info to help locate frequencies - not implemented yet
nVibIdent = size(vib_central_freq,2);
deltafexclude = 5*0.4;
df = mean(diff(f));
nfiexclude = floor(deltafexclude/df);
for kMode = 1:nMode
    if kMode == 1
        subplot(2,1,1)
        PSDmode = PSDTilt;
        
    else
        subplot(2,1,2)
        PSDmode = PSDTip;
    end
    PSDmode_LP = conv(PSDmode, [0.8 1 0.8],'same');

    TRPSDTilt = abs(PSDmode_LP(fimin:fimax) - modelPSDstruct{1}(fimin:fimax)).*f(fimin:fimax).^2;
    a = sort(TRPSDTilt,'descend');
    fc = fimin-1 + find(TRPSDTilt == a(1));
    vib(kMode,1) = f(fc);
    %TRPSDTilt(fc-nfiexclude:fc+nfiexclude) = 0;
    plot(f(fc),PSDmode(fc),'ro','markersize',10)
    for  kvib=2:nVibIdent
        TRPSDTilt(fc-nfiexclude-fimin:fc+nfiexclude-fimin) = 0;
        a = sort(TRPSDTilt,'descend');
        a(1)
        fc = fimin-1 + find(TRPSDTilt == a(1));
        vib(kMode,kvib) = f(fc)
        plot(f(fc),PSDmode(fc),'ro','markersize',10)
    end
end
fprintf('Automated central vibration frequency search: %2.2f\n',vib')
vib 

%vib_central_freq = vib;
%nVib = nVibIdent;
%% Vibration identification
for kMode = 1:nMode
    if kMode == 1
        PSDmode = PSDTilt;
        legend_strTilt = {'Tilt - telemetry','Tilt - model'};
    elseif kMode == 2
        PSDmode = PSDTip;
        legend_strTip = {'Tip - telemetry','Tip - model'};
    end
    
    % noise floor
    fimaxNoise = find(f >= fMaxRangeNoise);
    Nfloor(kMode) = abs(mean(PSDmode(fimaxNoise(1):nSamples/2) - modelPSDstruct{kMode}(fimaxNoise(1):nSamples/2)));
    if kMode == 1
        outTilt{1}.Nfloor = Nfloor(1);
    else
        outTip{1}.Nfloor = Nfloor(2);
    end
    model_atm_vib_acc = modelPSDstruct{kMode};
    for kVib = 1:nVib
        
        % fitting 1st vibration peak
        vibf = vib_central_freq(kMode,kVib);
        tic
        parms_vib{kVib} = fminsearch(@(x) QuadraticDifference(PSDmode(fimin:fimax), modelPSDstruct{kMode}(fimin:fimax)+Nfloor(kMode), vibf, x, f(fimin:fimax), Ts), [0.004*3, 0.5*6.5e-5]);
        toc
        
        % create model and plot
        csi = parms_vib{kVib}(1);
        w0v = vib_central_freq(kMode, kVib)*2*pi;
        
        a1v = 2*exp(-csi*w0v*Ts)*cos(w0v*sqrt(1-csi^2)*Ts);
        a2v = -exp(-2*csi*w0v*Ts);
        %e2 = (parms_vib{kVib}(2)/(1-a1v^2 - a2v^2 -2*a1v*a2v*a1v/(1-a2v)))^2;
        %betav = sqrt(e2);
        
        betav = parms_vib{kVib}(2);
        sig_model_var(kVib) = betav/(1-a1v^2 - a2v^2 -2*a1v*a2v*a1v/(1-a2v));
        
        model_vib = tf(betav, [1 -a1v -a2v], Ts);
        [mag_mv1, phase_mv] = bode(model_vib, 2*pi*f);
        mag_mv1 = reshape(mag_mv1, 1, length(mag_mv1));
        if verbose
            tt = 0:Ts:1000;
            uu = randn(size(tt));
            yd = lsim(model_vib, uu, tt);
            fprintf('Vibration model output variance: %1.4e,  numerical integration PSD: %1.4e\n', betav*sig_model_var(kVib), var(yd))
            trapz(f, mag_mv1.^2)*2*Ts; % = power (factor 2 is for two-sided PSD)
        end
        model_atm_vib_acc = model_atm_vib_acc + (mag_mv1).^2;
        
        subplot(2,1,kMode)
        loglog(f, model_atm_vib_acc + Nfloor(kMode),'k--', 'linewidth',2)
        if kMode == 1
            legend_strTilt =  [legend_strTilt, {['Vib1, f_0=', num2str(vib_central_freq(kMode,kVib)) 'Hz; FWHM =' num2str(2*csi) 'Hz;']}];
            legend(legend_strTilt)
        else
            legend_strTip =  [legend_strTip, {['Vib1, f_0=', num2str(vib_central_freq(kMode,kVib)) 'Hz; FWHM =' num2str(2*csi) 'Hz;']}];
            legend(legend_strTip)
        end
        legend('Location','SouthWest')
        % save outputs
        if kMode == 1
            outTilt{kVib+1}.parms_vib = parms_vib{kVib};
            outTilt{kVib+1}.model = model_vib;
            outTilt{kVib+1}.sig_model_var = sig_model_var(kVib);
            outTilt{kVib+1}.sig_var = betav*sig_model_var(kVib);
        else
            outTip{kVib+1}.parms_vib = parms_vib{kVib};
            outTip{kVib+1}.model = model_vib;
            outTip{kVib+1}.sig_model_var = sig_model_var(kVib);
            outTip{kVib+1}.sig_var = betav*sig_model_var(kVib);
        end
        
    end
end

function quadDiff = QuadraticDifference(telemPSD, modelPSD, model_vib_freq, model_parms, f, Ts)
 
 csi = model_parms(1);
 
% fit an AR2 model to the vibration
w0v = model_vib_freq*2*pi;

a1v = 2*exp(-csi*w0v*Ts)*cos(w0v*sqrt(1-csi^2)*Ts);
a2v = -exp(-2*csi*w0v*Ts);

betav = model_parms(2);

model_vib = tf(betav,[1 -a1v -a2v],Ts);
[mag_mv] = bode(model_vib, 2*pi*f);
mag_mv = reshape(mag_mv, 1, length(mag_mv));
quadDiff = sum(abs(telemPSD - modelPSD - abs(mag_mv).^2));

 function out = AoA_fit_parms(modelPSDX, modelPSDY, PSDTilt, PSDTip, p, f, fMinRange, fMaxRange)
 
 
 PSDaoaX = zeros(size(modelPSDX));
 PSDaoaY = zeros(size(modelPSDY));
 
 
 PSDaoaX = PSDaoaX + cos(p(2))^2*modelPSDX + sin(p(2)/2/pi)^2*modelPSDY;
 PSDaoaY = PSDaoaY + cos(pi/2 - p(2))^2*modelPSDX + sin(pi/2 - p(2)/2/pi)^2*modelPSDY;
 
 
 fimin = find(f> fMinRange); fimin = fimin(1);
 fimax = find(f <= fMaxRange); fimax = fimax(end);
 modelPSDX = interp1(f,PSDaoaX, f*p(1)/20);
 modelPSDY = interp1(f,PSDaoaY, f*p(1)/20);
 
 out = sum(abs(PSDTilt(fimin:fimax) - modelPSDX(fimin:fimax)).^2) + sum(abs(PSDTip(fimin:fimax) - modelPSDY(fimin:fimax)).^2);
