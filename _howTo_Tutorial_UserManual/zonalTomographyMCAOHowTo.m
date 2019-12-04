%% ADAPTIVE OPTICS ZONAL MCAO DEMO
% Demonstrate how to build a zonal tomographic reconstructor and the
% optimal projection of the reconstruction on multiples DM conjugated at
% different altitudes for a MCAO system
clear all;
close all;

%% Atmosphere definition

seeing          =  0.8;
zenithAngle     =  pi/6;
r0              =  0.976*0.5e-6/arcsec(seeing)*cos(zenithAngle)^(3/5);
L0              =  25;

% 10-layer atmosphere
fractionalR0    = [0.59 0.02 0.04 0.06 0.01 0.05 0.09 0.04 0.05 0.05];
altitude        = 1/cos(zenithAngle)*[30 140 281 562 1125 2250 4500 7750 11000 14000];
windSpeed       = [10 13 5 8 17 3 14 2 12 5];
windDirection   = 2*pi*rand(1,10);

% 3-layer atmosphere
% fractionalR0    = [0.7,0.25,0.05];
% altitude        = 1/cos(zenithAngle)*[0,4,10]*1000;
% windSpeed       = [5,10,20];
% windDirection   = [0,pi/4,pi];


atm = atmosphere(photometry.V0,r0,L0,...
    'fractionnalR0',fractionalR0,'altitude',altitude,...
    'windSpeed',windSpeed,'windDirection',windDirection);

%% Telescope
nL              = 32;                % number of lenslets
nPx             = 10;                % number of pixels per lenslet
nRes            = nL*nPx;            % resolution on the pupil plane (no of pixels)
D               = 8;                   % telescope primary mirror diameter
d               = D/nL;              % lenslet pitch
samplingFreq    = 1000;      % WFS sampling time
obstructionRatio= 0.15;  % central obscuration ratio
fieldOfViewInArcsec = 43; %fieldOfViewInArcsec

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);


%% Calibration sources
ngs = source;


%% Wavefront sensor
wfs             = shackHartmann(nL,nRes,0.75);
ngs.wavelength  = photometry.Na;

ngs = ngs.*tel*wfs;
wfs.INIT;
wfs.camera.photonNoise = false;
+wfs;
wfs.camera.frameListener.Enabled    = false;
wfs.slopesListener.Enabled          = false;

% adjust the optical gain such that CoG is calibrated
wfs.gainCalibration(tel,ngs);

%% STACK OF DMs
nDm             = 2;
dmPitch         = [0.25 0.4 0.4];
bif(1)          = gaussianInfluenceFunction(25/100,dmPitch(1));
bif(2)          = gaussianInfluenceFunction(25/100, dmPitch(2));
bif(3)          = gaussianInfluenceFunction(25/100, dmPitch(3));

% low-resolution influence functions to fit modes to the DM
bifLowRes(1)    = gaussianInfluenceFunction(25/100, dmPitch(1));
bifLowRes(2)    = gaussianInfluenceFunction(25/100,dmPitch(2));
bifLowRes(3)    = gaussianInfluenceFunction(25/100,dmPitch(3));

nActuators      = [33 25 25]; % actuators across the DM diameter

dmHeights       = [0 8e3 11e3];
validActuator   = {true(33) true(25) true(25)};

dmStack = multiDeformableMirror(dmHeights(1:nDm), nActuators(1:nDm) , tel, ...
    'modes',bif(1:nDm),...
    'validActuator',validActuator(1:nDm) );

nValidActuatorsStack = dmStack.nValidActuators;

dmStackLowRes = multiDeformableMirror(dmHeights(1:nDm), nActuators(1:nDm) , tel, ...
    'modes',bifLowRes(1:nDm),...
    'validActuator',validActuator(1:nDm),...
    'resolution',nActuators(1:nDm));



    
%% SOURCES

% Optimization sources for the MCAO fitting step
x           = linspace(-15,15,11);
[x y]       = meshgrid(x);
[theta rho] = cart2pol(x,y);
zenithOpt   = rho(:)*constants.arcsec2radian;
azimuthOpt  = theta;
optFitSrc   = source('zenith', zenithOpt, 'azimuth', azimuthOpt);

% LGS sources for tomography
radiusAst   = 15;
nAst        = 4;
thetaAst    = 0;
lgsHeight   = 90e3;
spotSize    = [];
flux        = 25e6; % 1000ph/frame/spup (40x40 spup, 1kHz, 0.2 side square)

% % LGS
lgsAst = laserGuideStar(d,tel.D, lgsHeight/cos(zenithAngle), ...
    spotSize, flux, [],'asterism',{[nAst,arcsec(radiusAst),thetaAst]}, 'wavelength', photometry.Na,'height',lgsHeight/cos(zenithAngle));


%% INTERACTION MATRIX
% Several interaction matrices are computed for each DM and guide-star
% location
lgsAst      = lgsAst.*tel*wfs;
calibDms    = calibrationMultiDm(dmStack,wfs,lgsAst,lgsAst(1).wavelength/100,25,'cond',1e2);

% create a multi-calibration vault which is a "concatenation" of calibrationVaults
multiDmCalibVault = multiCalibrationVault(calibDms);
%% MCAO TOMOGRAPHIC RECONSTRUCTOR
% sparse-based tomographic reconstructor that explitely estimates the phase
% at the layers 
lgsAst = lgsAst.*tel*wfs;
mmse = sparseLMMSE(wfs,tel,atm,lgsAst,...
    'iNoiseVar', 1/1e-14,...
    'layerEstimation',true);

%% LTAO RECONSTRUCTOR
%
% iterative Toeplitz reconstructor. Produces a phase estimate in teh pupil
% plane in the direction of the mmseStar
slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',[ngs ngs]);

% Alternative call
% noiseVar = 1e-14;
%         lgsAst_slmmse = slopesLinearMMSE(wfs,...
%             tel,...
%             atm,...
%             ast,...
%             'mmseStar',ngs,...
%             'noiseVar', noiseVar,...
%             'isTTRM', true,...
%             'isFocusRM', true,...
%             'covModel','FFT',...
%             'wavefrontMask', wfs.validActuator);

%% FITTING MATRIX
bifLowRes2(1)    = gaussianInfluenceFunction(25/100, dmPitch(1));
bifLowRes2(2)    = gaussianInfluenceFunction(25/100,dmPitch(2));
bifLowRes2(3)    = gaussianInfluenceFunction(25/100,dmPitch(3));

F = fittingMatrix(dmStack,mmse,optFitSrc,bifLowRes2(1:nDm));

%% SCIENCE CAMERAS

% Science Sources for performance estimation
xs = linspace(-15,15,3);
[Xs Ys] = meshgrid(xs);
[theta rho] = cart2pol(Xs,Ys);
zenithSci     = rho(:)*constants.arcsec2radian;
azimuthSci  = theta;

sciH = source('zenith', zenithSci, 'azimuth', azimuthSci,'wavelength',photometry.H);

nGs = length(lgsAst);
nSciH = length(sciH);
nNgs = length(ngs);


camH = imager('nyquistSampling',2,'fieldStopSize',75);
sciH = sciH.*tel*camH;
camH.referenceFrame = camH.frame;

tel = tel + atm;




%% Loop initialization
% reset(tel);
%% With noise on the WFS
% wfs.camera.readOutNoise = 0;
% wfs.camera.photonNoise = true;
% wfs.framePixelThreshold = 0;

lgsAst = lgsAst.*tel*dmStack*wfs;
flush(camH)
camH.frame          = camH.frame*0;
camH.clockRate      = 1;
exposureTime        = 1000;
camH.exposureTime   = exposureTime;
startDelay          = 10;
camH.startDelay     = startDelay;
flush(camH)
set(sciH,'logging',true)
nIt                 = camH.startDelay+camH.exposureTime;
dmStack.coefs       = 0;

camH.frameListener.Enabled = 1;

figure(31416)
imagesc(camH,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(sciH));
axis xy equal tight
colorbar


%% LTAO WITH POLC CONTROLLER
% As an initial test case, we set up a LTAO (or NTAO) to make sure
% everything is on track


FlowRes         = 2*bifLowRes(1).modes(wfs.validActuator,:);
iFlowRes        = pinv(full(FlowRes),1e-1);
dmStack.coefs   = 0;
gain_pol        = 0.5;

for k=1:nIt
    % Objects update
    +tel;
    +lgsAst;
    sciH=sciH.*tel*dmStack*camH;
    % Pseudo-open-loop controller
    dmStack.dms{1}.coefs = (1-gain_pol)*dmStack.dms{1}.coefs + ...
        gain_pol*iFlowRes*( slmmse*( bsxfun( @minus,wfs.slopes, calibDms{1}.D*dmStack.dms{1}.coefs ) ) );
    % Display
    set(h,'Cdata',catMeanRmPhase(sciH))
    drawnow
end
imagesc(camH)
set(h,'Cdata',catMeanRmPhase(sciH))

% PLOT RESULTS

var_wfe_ltaoH           = zeros(nIt,2,nSciH);
wfe_ltaoH               = zeros(nIt,2,nSciH) ;
marechalStrehl_ltaoH    = zeros(nSciH,1);
psfStrehl_ltaoH         = 1e2*camH.strehl;
for k = 1: nSciH
    var_wfe_ltaoH(:,:,k)        = reshape(sciH(k).phaseVar(1:nIt*2),2,[])';
    wfe_ltaoH(:,:,k)            = sqrt(var_wfe_ltaoH(:,:,k))/sciH(1).waveNumber*1e6;
    marechalStrehl_ltaoH(k,1)   = 1e2*exp(-mean(var_wfe_ltaoH(startDelay:end,2,k)));
end

% FOR 1000 iterations the camera SR and the Marechal approximation are
% relatively consistent (trends are equal), absolute values are not, with
% the Marechal being 10-20 % points below
figure(100)
bar([psfStrehl_ltaoH/10 marechalStrehl_ltaoH])
xlabel('#star position')
ylabel('Strehl-ratio, H-band')
text(2,90,'Blue: imager SR, yellow: Marechal approx')
set(gca,'fontSize',16)


%% PSEUDO OPEN-LOOP CONTROLLER IN MCAO MODE
%
% The loop is closed for one full exposure of the science camera.
%
gain_pol        = 0.5;
dmStack.coefs   = 0; %put commands to zero
slopesPsolMat  = zeros(wfs.nSlope,nGs);

for k=1:nIt
    tic
    fprintf('Frame %d/%d\n',k,nIt);
    
    +tel;
    +lgsAst;
    sciH=sciH.*tel*dmStack*camH;
    
    slopesPsolMat = multiDmCalibVault*dmStack; %
    
    dmStack.coefs = (1-gain_pol)*dmStack.coefs + ...
        gain_pol*F.commandMatrix*(mmse*bsxfun( @minus, wfs.slopes, slopesPsolMat ));
    
    % Display WFE
    set(h,'Cdata',catMeanRmPhase(sciH))
    drawnow
    toc
end

%% RESULTS

var_wfe_mcaoH = zeros(nIt,2,nSciH);
wfe_mcaoH =zeros(nIt,2,nSciH) ;
marechalStrehl_mcaoH = zeros(nSciH,1);
psfStrehl_mcaoH =1e2*camH.strehl;
for k = 1: nSciH
    var_wfe_mcaoH(:,:,k) = reshape(sciH(k).phaseVar(1:nIt*2),2,[])';
    wfe_mcaoH(:,:,k) = sqrt(var_wfe_mcaoH(:,:,k))/sciH(1).waveNumber*1e6;
    marechalStrehl_mcaoH(k,1) = 1e2*exp(-mean(var_wfe_mcaoH(startDelay:end,2,k)));
end
figure(110)
bar([psfStrehl_mcaoH marechalStrehl_mcaoH])
xlabel('#star position')
ylabel('Strehl-ratio, H-band')
text(2,90,'Blue: imager SR, yellow: Marechal approx')
set(gca,'fontSize',16)

%% MCAO WITH SPLIT-TOMOGRAPHY

% Create a LOWFS to measure TT and TTA modes
tel = tel - atm;
ngs.wavelength = photometry.H;
ttRebinFactor  = 1;
nResTT = 20;
LOWFS = shackHartmann(1,nResTT,1);
%LOWFS.camera.resolution = [10 10];
LOWFS.lenslets.nyquistSampling = 0.5;
%LOWFS.lenslets.fieldStopSize = 20;

LOWFS.camera.resolution = nRes/ttRebinFactor;
ngs = ngs.*tel*LOWFS;

LOWFS.INIT

+LOWFS;
figure
imagesc(LOWFS.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(LOWFS,'parent',subplot(3,2,[5,6]))

% LOWFS gain calibration
LOWFS.gainCalibration(tel,ngs)
%% Piston-tip-tilt removal
 
 
clear bif dmLowRes zern TT DMTTRem iF
for i = 1:dmStack.nDm
%     bifLowRes1{i} = gaussianInfluenceFunction(25/100,   dmPitch(i));
%     dmLowRes{i} = deformableMirror(dmStack.dms{i}.nActuator,'modes',bifLowRes1{i} ,...
%         'resolution',dmStack.dms{i}.nActuator,...
%         'validActuator',dmStack.dms{i}.validActuator);
%     Finfl = 2*bifLowRes1{i}.modes(dmStack.dms{i}.validActuator,:);
    Finfl = 2*bifLowRes(i).modes(dmStack.dms{i}.validActuator,:);
    iF{i} = pinv(full(Finfl),1e-1);
    zern  = zernike(2:3,'resolution',dmStack.dms{i}.nActuator, 'pupil',dmStack.dms{i}.validActuator);
    TT{i} = iF{i}*zern.modes(dmStack.dms{i}.validActuator,:);
    dmTTRem{i} = eye(dmStack.dms{i}.nValidActuator) - TT{i}*pinv(TT{i});
end
dmTTRemFull = blkdiag(dmTTRem{:});

t = ones(wfs.nValidLenslet,1);
TT = [t 0*t; 0*t,t];
slopeTTRem = eye(2*wfs.nValidLenslet) - TT*pinv(TT);

%% TTA :: set of modes that produce only tilt on the LGS measurements and are thus discarded

%See paper Correia et al, Increased sky coverage with optimal correction 
% of tilt and tilt-anisoplanatism modes in laser-guide-star multiconjugate
% adaptive optics, JOSA A (2013)

% TT NGS Sources
nTT         = 3;
ttRadius    = 20;
ttAst       = source('asterism',{[nTT,arcsec(ttRadius),0]},'wavelength',photometry.H,'magnitude',15);

%Zernike TT to Quad Zernike TTA
% Normalization factors
rn          = tel.D/tel.diameterAt(dmStack.dms{2}.zLocation); % ratio of meta-pupil at DM2 altitude to primary mirror diameter
rc          = 1-dmStack.dms{2}.zLocation/(lgsHeight/cos(zenithAngle)); % cone compression factor
rl          = rc*rn; 
% Zernike projection
maxRadialDegree = 2;
zern        = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nRes,'pupil',tel.pupil);
zern.D      = tel.D;
pistonRM    = 1;
[projAlpha] = analyticalSmallFootprintExpansion(zern,tel,ttAst,dmStack.zLocations,pistonRM);
projTheta   = projAlpha([1,2,6,7,11,12],:);
Pmu2phi = [[1 0 0 0 0];...
    [0 1 0 0 0];...
    [0 0 1 0 0];...
    [0 0 0 1 0];...
    [0 0 0 0 1];...
    [0 0 0 0 0];...
    [0 0 0 0 0];...
    [0 0 -1/rl^2 0 0];...
    [0 0 0 -1/rl^2 0];...
    [0 0 0 0 -1/rl^2];]
if nDm == 3
    Pmu2phi = [Pmu2phi;zeros(5)];
end

projNgs     = projTheta*Pmu2phi;
invProjNgs  = pinv(projNgs);
% Zernike to DMs Surface
clear zernTT TTQ
for kDm = 1:nDm
    zernTT{kDm} = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),...
        'resolution',dmStack.dms{kDm}.nActuator, 'pupil',dmStack.dms{kDm}.validActuator);
    TTQ{kDm} = iF{kDm}*zernTT{kDm}.modes(dmStack.dms{kDm}.validActuator,:);
end
TTQ = blkdiag(TTQ{:});

Z2C = TTQ*Pmu2phi*invProjNgs; % zernike coeffs to DM commands through TTA modes
% Zernike WFS calibration
zernT1 = zernike(2:3,'resolution',nRes, 'pupil',tel.pupil);
zernT1.lex = false;
zernT1.c = eye(zernT1.nMode)*ttAst(1).wavelength/2/pi;
ngs.wavelength = photometry.H;
ngs=ngs.*zernT1*LOWFS;
Z2S = LOWFS.slopes; % Z2S  zernike2slopes when zernike are in radians
S2Z = pinv(Z2S)*ttAst(1).wavelength/2/pi; % attention the S2Z uses convertion to meter since default Zernike units are meters
%

%% PSEUDO OPEN-LOOP CONTROLLER IN MCAO MODE
%
% The loop is closed for one full exposure of the science camera.
%
% set the atmosphere back in
tel = tel + atm;

%INTEGRATOR
gain_cl         = 0.5;
uNgs            = zeros(nValidActuatorsStack,1);
ttAst           = ttAst.*tel*dmStack*LOWFS;

gain_pol        = 0.5;
slopesPsolMat   = zeros(wfs.nSlope,nGs);
dmStack.coefs   = 0; %put commands to zero

slopesPsolCell  = cell(nDm,nGs);
slopesPsolMat   = zeros(wfs.nSlope,nGs);
uLgs            = zeros(nValidActuatorsStack,1);
uNgs            = zeros(nValidActuatorsStack,1);
zTT             = zernike(2:3);
for k=1:nIt
    tic
    fprintf('Frame %d/%d\n',k,nIt);
    
    +tel;
    +lgsAst;
    +ttAst;
    sciH=sciH.*tel*dmStack*camH;
    
    slopesPsolMat = multiDmCalibVault*dmStack; %
    
    uLgs = (1-gain_pol)*dmTTRemFull*uLgs + ...
        gain_pol*dmTTRemFull*F.commandMatrix*(mmse*(slopeTTRem*bsxfun( @minus, wfs.slopes, slopesPsolMat )));
        
    %Evaluation of TTA and corresponding DMs commands

    %zTT = zernike(2:3)\LOWFS;
    zCoeff = S2Z*LOWFS.slopes;
    
    % Add TT from TT stars, % Closed-loop controller
    uNgs = uNgs - gain_cl*Z2C*zCoeff(:);
    
    dmStack.coefs = uLgs - uNgs;
    % Display WFE
    set(h,'Cdata',catMeanRmPhase(sciH))
    drawnow
    toc
end

%% RESULTS

var_wfe_mcaoH_SplitTomo = zeros(nIt,2,nSciH);
wfe_mcaoH_SplitTomo =zeros(nIt,2,nSciH) ;
marechalStrehl_mcaoH_SplitTomo = zeros(nSciH,1);
psfStrehl_mcaoH_SplitTomo =1e2*camH.strehl;
for k = 1: nSciH
    var_wfe_mcaoH_SplitTomo(:,:,k) = reshape(sciH(k).phaseVar(1:nIt*2),2,[])';
    wfe_mcaoH_SplitTomo(:,:,k) = sqrt(var_wfe_mcaoH_SplitTomo(:,:,k))/sciH(1).waveNumber*1e6;
    marechalStrehl_mcaoH_SplitTomo(k,1) = 1e2*exp(-mean(var_wfe_mcaoH_SplitTomo(startDelay:end,2,k)));
end
figure
bar([marechalStrehl_mcaoH_SplitTomo, marechalStrehl_mcaoH])
legend('Full tomography','split tomography')
xlabel('star #')
ylabel('Strehl-ratio, H-band')


%% AGGREGATED RESULTS
figure
bar([psfStrehl_ltaoH, psfStrehl_mcaoH, psfStrehl_mcaoH_SplitTomo])
legend('LTAO','Full MCAO tomography','MCAO split tomography')
xlabel('star #')
ylabel('Strehl-ratio, H-band')


figure
subplot(1,3,1)
[C,h] = contourf(reshape(psfStrehl_ltaoH,3,3),5);
%clabel(C,h)
clabel(C,'FontSize',8,'Color','r','Rotation',0)
colorbar
title('LTAO')
subplot(1,3,2)
[C,h] = contourf(reshape(psfStrehl_mcaoH,3,3),5);
clabel(C,'FontSize',8,'Color','r','Rotation',0)
colorbar
title('Full LGS MCAO tomography')
subplot(1,3,3)
[C,h] = contourf(reshape(psfStrehl_mcaoH_SplitTomo,3,3),5)
clabel(C,'FontSize',8,'Color','r','Rotation',0)
colorbar
title('MCAO Split-tomography')