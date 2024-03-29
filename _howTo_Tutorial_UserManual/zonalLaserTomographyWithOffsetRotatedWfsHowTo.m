%% ADAPTIVE OPTICS ZONAL LTAO DEMO WITH WFS OFFSETS AND ROTATION
 
%{
This demo is based on the zonalLaserTomographyHowTo

We focus here on the use of the sparse tomographic reconstructor in the
presence of WFS offsets and rotations.
%}
%%
clear all;
close all;
%% SIMULATION OPTIONS (USER-DEFINED)
extendedFriedDmLayout   = 1; % use a stand-alone definition of valid actuators based on each IF influence within the pupil.
useWfsG = 1; % use a G-SHWFS in open-loop
useDoubleResolutionReconstruction = 1; % can only be used with the reconstructorType = 'sparse' or 'lmmse'
useExplicitReconFromSparse = 0;

if useDoubleResolutionReconstruction == 2
    os = 4;
elseif useDoubleResolutionReconstruction == 1
    os = 2;
elseif useDoubleResolutionReconstruction == 0
    os = 1;
end

%% SOURCES
% on--axis souce
ngs = source;
zenithAngle = 0;
%
%% ATMOSPHERE - overrides predifined ATMOSPHERE

r0              = 0.1587;   % coherence lenght in meters at 0.5microns
L0              = 100;      % Outer scale in meters


% Multi-layer atmosphere
fractionalR0    = [0.5, 0.05, 0.1, 0.13, 0.08, 0.04, 0.1];
altitude        = (0:2:12)*1e3;
windSpeed       = [10,5,20,5 10 10 10];
windDirection   = [0,pi/2,pi/4, -pi/2, pi/3, 2/5*pi, 7/3*pi];


% Multi-layer atmosphere
% fractionalR0   = [0.5,0.3,0.2];
% altitude        = [0e3,5e3,12e3];
% windSpeed        = [10,5,20];
% windDirection   = [0,pi/2,pi/4];


%Mono-layer atmosphere
% fractionalR0    = [1];
% altitude        = [1]*1e3;
% windSpeed       = [10];
% windDirection   = [0];


atm = atmosphere(photometry.V0,r0,L0,...
    'fractionnalR0',fractionalR0,'altitude',altitude,...
    'windSpeed',windSpeed,'windDirection',windDirection);


%% TELESCOPE - overrides predifined TELESCOPE
nOuterLenslet       = 4;
nInnerLenslet       = 16;
nL = nInnerLenslet + nOuterLenslet; % number of lenslets
nPx                 = 8;            % number of pixels per lenslet
nRes = nInnerLenslet*nPx;           % resolution on the pupil plane (no of pixels)
D                   = 8;            % telescope primary mirror diameter
d                   = D/nL;         % lenslet pitch
samplingFreq        = 500;          % WFS sampling time
obstructionRatio    = 0.3;          % central obscuration ratio
fieldOfViewInArcsec = 60;           %fieldOfViewInArcsec

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);


% enlarge the telescope computational grid to accommodate the WFS offsets
% hack to enlarge the telescope computational grid
p = utilities.piston((nL-nOuterLenslet)*nPx,nRes);
pin = utilities.piston(tel.obstructionRatio*nRes,nRes);
tel.pupil = p-pin;

tel.embed(nRes*nL/nInnerLenslet);
nRes = tel.resolution;

% %% TELESCOPE - overrides predifined TELESCOPE
% nSubapExtra = 4;
% nL   = 16+nSubapExtra;  % number of lenslets
% nPx  = 8;              % number of pixels per lenslet
% nRes = nL*nPx;          % resolution on the pupil plane (no of pixels)
% D    = 8*nL/(nL-nSubapExtra);   % telescope primary mirror diameter
% d    = D/nL;            % lenslet pitch
% samplingFreq        = 500;  % WFS sampling time
% obstructionRatio    = 0.3;  % central obscuration ratio
% fieldOfViewInArcsec = 60;   %fieldOfViewInArcsec
% 
% tel1 = telescope(D,'resolution',nRes,...
%     'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);
% 
% 
% % enlarge the telescope computational grid to accommodate the WFS offsets
% % hack to enlarge the telescope computational grid
% p = utilities.piston((nL-nSubapExtra)*nPx,nRes);
% pin = utilities.piston(tel1.obstructionRatio*nRes,nRes);
% tel1.pupil = p-pin;


%% Calibration sources
ngs = source;


%% Wavefront sensor
wfs                     = shackHartmann(nL,nRes,0.75);
ngs.wavelength          = photometry.Na;

ngs = ngs.*tel*wfs;

wfs.INIT;
wfs.camera.photonNoise = false;
+wfs;
wfs.camera.frameListener.Enabled    = false;
wfs.slopesListener.Enabled          = false;

% adjust the optical gain such that CoG is calibrated
wfs.gainCalibration(tel,ngs);


%% geomGrad
wfsG = geomGrad(wfs.validLenslet, tel.pupilLogical);

%% STACK OF DMs
nDm             = 1;
dmPitch         = [0.5];
dmCrossCouplingCoeff = 0.25;

bif(1)          = gaussianInfluenceFunction(dmCrossCouplingCoeff,dmPitch(1));


nActuators      = 21;%[nL+1]; % actuators across the DM diameter

dmHeights       = [0];

if extendedFriedDmLayout % create a scrap DM just for the definition if the extended valid actuator map
    % SCRAP DM
    bifa               = gaussianInfluenceFunction(dmCrossCouplingCoeff, dmPitch(1));
    dmScrap = deformableMirror(nActuators,'modes',bifa,...
        'resolution',tel.resolution,...
        'validActuator',true(nActuators), 'zLocation',0);
    validActMap = dmScrap.setValidActuator(tel.pupil, 0.1);
    clear dmScrap
else
    validActMap = wfs.validActuator;
end

dm = deformableMirror(nActuators,'modes',bif,...
    'resolution',tel.resolution,...
    'validActuator',validActMap, 'zLocation',0);

% low-resolution influence functions to fit modes to the DM
bifLowRes(1)    = gaussianInfluenceFunction(dmCrossCouplingCoeff, dmPitch(1));
bifLowRes2x(1)    = gaussianInfluenceFunction(dmCrossCouplingCoeff, dmPitch(1));


% WF reconstruction DMs
if useDoubleResolutionReconstruction
    if useDoubleResolutionReconstruction == 1
        dmLowRes2x = deformableMirror(nActuators,'modes',bifLowRes2x,...
            'resolution',2*nL+1,...
            'validActuator',validActMap, 'zLocation',0 );
        reconstructionGrid = wfs.reconstructionGrid(wfs,2);
    elseif useDoubleResolutionReconstruction == 2
        dmLowRes2x = deformableMirror(nActuators,'modes',bifLowRes2x,...
            'resolution',4*nL+1,...
            'validActuator',validActMap, 'zLocation',0 );
        [G, reconstructionGrid] = wfs.sparseGradientMatrixAmplitudeWeighted([],4);
    end
else
    dmLowRes = deformableMirror(nActuators,'modes',bifLowRes,...
        'resolution',nL+1,...
        'validActuator',validActMap, 'zLocation',0 );
    reconstructionGrid = wfs.reconstructionGrid(wfs,1);
end




%% SOURCES

% Optimization sources for the fitting step
x           = linspace(-15,15,11);
[x y]       = meshgrid(x);
[theta rho] = cart2pol(x,y);
zenithOpt   = rho(:)*constants.arcsec2radian;
azimuthOpt  = theta;
optFitSrc   = source('zenith', zenithOpt, 'azimuth', azimuthOpt);
optFitSrc   = source('zenith', 0, 'azimuth', 0);
% LGS sources for tomography
radiusAst   = 10;
nLgs        = 4;
thetaAst    = 0;
lgsHeight   = 90e3;
spotSize    = [];
flux        = 25e5; % 1000ph/frame/spup (40x40 spup, 1kHz, 0.2 side square)

% % LGS
lgsAst = laserGuideStar(d,tel.D, lgsHeight/cos(zenithAngle), ...
    spotSize, flux, [],'asterism',{[nLgs,arcsec(radiusAst),thetaAst]}, 'wavelength', photometry.Na,'height',lgsHeight/cos(zenithAngle));

theta = linspace(0,2*pi-2*pi/nLgs,nLgs);
for iGS = 1:nLgs
    lgsAst(iGS).viewPoint = D/2*[cos(theta(iGS)), sin(theta(iGS))];
end

%% OFFSETS AND ROTATION
% add a rotation in radians, one scalar per WFS
wfs.lenslets.rotation = -pi/4*ones(1,4);
wfs.lenslets.rotation = 1*[pi/10 -pi/2 pi pi*0];
% add an offset [x;y] per WFS in units of 1/D
wfs.lenslets.offset = 1*[1 -1 -1 1; 1 1 -1 -1]*0.25/nL;



%% INTERACTION MATRIX
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
ngs = source;
ngs.wavelength = photometry.Na;
ngs = ngs.*tel;


calibDm = calibration(dm,wfs,ngs,ngs.wavelength/8,nL+1,'cond',1e2);

%% The interaction matrix is replaced by the geometric IM
dm.coefs = eye(dm.nValidActuator)*ngs.wavelength/2/pi;
wfsG.offset = wfs.lenslets.offset;
wfsG.rotation = wfs.lenslets.rotation;

ngs = ngs.*tel*dm*wfsG;
Dgeom = wfsG.slopes*2*pi/ngs.wavelength;
%calibDm.D = Dgeom;

% %% SANITY-CHECK THAT THE WFS IMPLEMENT THE ROTATION + OFFSETS CORRECTLY
% wfsG.rotation = 0*-pi/4*ones(1,4);
% wfsG.offset = [1 -1 -1 1; 1 1 -1 -1]*0.25/nL;
% zern = zernike(5,tel.D,'resolution',tel.resolution,'pupil',ones(tel.resolution));
% %zern = zernike(5,tel.D,'resolution',tel.resolution,'pupil',abs(tel.pupil));
% zern.c = ngs.wavelength/2/pi;
% ngs = ngs.*zern*wfsG;
% figure(1)
% plot(wfsG.slopes)
% hold on
% 
% zern = zernike(6,tel.D,'resolution',tel.resolution,'pupil',ones(tel.resolution));
% %zern = zernike(5,tel.D,'resolution',tel.resolution,'pupil',abs(tel.pupil));
% zern.c = ngs.wavelength/2/pi;
% 
% wfsG.rotation = 1*-pi/4*ones(1,4);
% wfsG.offset = [1 -1 -1 1; 1 1 -1 -1]*0.25/nL;
% ngs = ngs.*zern*wfsG;
% figure(1)
% plot(wfsG.slopes)
% hold on
% 
% 
% % zern = zernike(6,tel.D,'resolution',tel.resolution,'pupil',abs(tel.pupil));
% % zern.c = ngs.wavelength/2/pi;
% % ngs = ngs.*tel*zern;
% % figure(2)
% % imagesc(ngs.meanRmPhase)

%% SPARSE MMSE
% sparse-based tomographic reconstructor that explitely estimates the phase
% at the layers

science = source('wavelength',photometry.H);


lgsAst = lgsAst.*tel*wfs;
%        if os == 1
sparseMmse = sparseLMMSE(wfs,tel,atm,lgsAst,...
    'iNoiseVar', 1/1e-14,...
    'layerEstimation',true,...
    'mmseStar',science,...
    'outputWavefrontMask', reconstructionGrid,'overSampling',os);%,'phaseCovariance','real');
%         else
%             % create an amplitude mask by down-sampling the telescope pupil to the WFS
%             % resolution
%             x = 1:nRes;
%             [Xi, Yi] = meshgrid(x,x);
%             x = 0.5:nPx/2:nRes+0.5;
%             [Xo, Yo] = meshgrid(x,x);
%             ampMask = interp2(Xi, Yi, tel.pupil, Xo, Yo);
%             ampMask(isnan(ampMask)) = 0;
%             %ampMask = logical(ampMask); % TODO: revisit this formulation. does not improve over using a binary mask (next line)
%             ampMask = logical(ones(2*nL+1));
%
%             science = source('wavelength',photometry.H);
%
%             lgsAst = lgsAst.*tel*wfs;
%             tstart = tic;
%             sparseMmse = sparseLMMSE(wfs,tel,atm,lgsAst,...
%                 'iNoiseVar', 1/1e-14,...
%                 'layerEstimation',true,...
%                 'mmseStar',science,...
%                 'outputWavefrontMask', reconstructionGrid,...
%                 'maskAmplitudeWeighted',ampMask);
%         end
if useExplicitReconFromSparse
    [Rec,Left,Right] = sparseMmse.getReconstructionMatrix;
end

%% SANITY-CHECK THAT THE SPARSE-LMMSE IMPLEMENT THE ROTATION + OFFSETS CORRECTLY
% %Mono-layer atmosphere
% atmModel = atmosphere(photometry.V0,r0,L0,...
%     'fractionnalR0',1,'altitude',0,...
%     'windSpeed',10,'windDirection',0);
% 
% 
% % add a rotation in radians, one scalar per WFS
% wfs.lenslets.rotation = 1*-pi/4*ones(1,4);
% % add an offset [x;y] per WFS in units of 1/D
% wfs.lenslets.offset = 1*[10 -1 -1 1; 
%                          -10 1 -1 -1]*0.25/nL;
% %wfs.lenslets.offset = 10*ones(2,4)*0.25/nL;
% 
% wfsG.offset = wfs.lenslets.offset;
% wfsG.rotation = wfs.lenslets.rotation;
% 
% 
% science = source('wavelength',photometry.H);
% 
% lgsAst = lgsAst.*tel*wfs;
% %        if os == 1
% sparseMmse = sparseLMMSE(wfs,tel,atmModel,lgsAst,...
%     'iNoiseVar', 1/1e-14,...
%     'layerEstimation',true,...
%     'mmseStar',science,...
%     'outputWavefrontMask', reconstructionGrid,'overSampling',os);%,'phaseCovariance','real');
% 
% zernMode = 5;
% zern = zernike(zernMode,tel.D,'resolution',41,'pupil',ones(41));
% zern.c = ngs.wavelength/2/pi;
% cross = zeros(41);
% cross(19:23,: ) = 1;
% cross(:, 19:23) = 1;
% zern.modes = cross(:);
% 
% 
% resTot = sparseMmse.H*zern.modes;
% ith = 1;
% 
% res = resTot((ith-1)*760 + (1:760));
% 
% res2D = zeros(41);
% res2D(reconstructionGrid) = res;
% figure(10)
% imagesc(res2D)
% 
% 
% 
% zernHR = zernike(zernMode,tel.D,'resolution',tel.resolution,'pupil',ones(tel.resolution));
% zernHR.c = ngs.wavelength/2/pi;
% 
% lgsAst = lgsAst.*tel*zernHR*wfsG;
% phiRes = (sparseMmse*wfsG.slopes);
% 
% phiRes = sparseMmse.Hss*(sparseMmse*wfsG.slopes);
% 
% rec2D = zeros(41);
% rec2D(reconstructionGrid) = phiRes;
% figure(20)
% imagesc(rec2D)
% 
% 
% figure(39)
% subplot(1,4,1)
% in = reshape(zern.modes,41,41).*reconstructionGrid;
% imagesc(in)
% title('input')
% subplot(1,4,2)
% imagesc(res2D)
% title('model prop through H')
% subplot(1,4,3)
% imagesc(rec2D)
% title('reconstructed')
% subplot(1,4,4)
% err = rec2D/ngs.wavelength*2*pi - in;
% sqrt(var(err(reconstructionGrid)))
% imagesc(err)
% title('err')

%% SCIENCE CAMERAS

% Science Sources for performance estimation
xs = linspace(-15,15,3);
[Xs Ys] = meshgrid(xs);
[theta rho] = cart2pol(Xs,Ys);
zenithSci     = rho(:)*constants.arcsec2radian;
azimuthSci  = theta;

sciH = source('zenith', zenithSci, 'azimuth', azimuthSci,'wavelength',photometry.H);
sciH = source('wavelength',photometry.H);

nGs = length(lgsAst);
nSciH = length(sciH);
nNgs = length(ngs);

tel = tel - atm;
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


if useWfsG
    lgsAst = lgsAst.*tel*wfsG*dm*wfs;
else
    lgsAst = lgsAst.*tel*dm*wfs;
end

flush(camH)
camH.frame          = camH.frame*0;
camH.clockRate      = 1;
exposureTime        = 200;
camH.exposureTime   = exposureTime;
startDelay          = 10;
camH.startDelay     = startDelay;
flush(camH)
set(sciH,'logging',true)
nIt                 = camH.startDelay+camH.exposureTime;
dm.coefs       = 0;

camH.frameListener.Enabled = 1;

figure(31416)
imagesc(camH,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(sciH));
axis xy equal tight
colorbar

%% LOOP OPTIONS
if useDoubleResolutionReconstruction
    FlowRes2x       = 2*bifLowRes2x(1).modes(reconstructionGrid,:);
    iFlowRes2x      = pinv(full(FlowRes2x),1e-3);
    fittingMatrix   = iFlowRes2x*sparseMmse.Hss;
else
    FlowRes         = 2*bifLowRes(1).modes(reconstructionGrid,:);
    iFlowRes        = pinv(full(FlowRes),1e-3);
    fittingMatrix   = iFlowRes*sparseMmse.Hss;
end

gain_pol        = 0.6;
gainLoop        = 0.35;


%% LTAO WITH POLC CONTROLLER
dm.coefs   = 0;
for k=1:nIt
    % Objects update
    +tel;
    +lgsAst;
    sciH=sciH.*tel*dm*camH;
    % Pseudo-open-loop controller
    if useWfsG
        dm.coefs = (1-gain_pol)*dm.coefs + ...
            gain_pol*fittingMatrix*( sparseMmse*( wfsG.slopes ) );
        if useExplicitReconFromSparse
            dm.coefs = (1-gain_pol)*dm.coefs + ...
                gain_pol*fittingMatrix*( Rec*( wfsG.slopes(:) ) );
        end
        
    else
        dm.coefs = (1-gain_pol)*dm.coefs + ...
            gain_pol*fittingMatrix*( sparseMmse*( bsxfun( @minus, wfs.slopes, calibDm.D*dm.coefs ) ) );
        if useExplicitReconFromSparse
            dm.coefs = (1-gain_pol)*dm.coefs + ...
                gain_pol*fittingMatrix*( Rec*( wfs.slopes(:) ) );
        end
    end
    
    
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

fprintf(['Simulation SR @' num2str(sciH.wavelength*1e6, '%2.2f') 'um is %2.2f \n'], psfStrehl_ltaoH)
fprintf(['Marechal SR @' num2str(sciH.wavelength*1e6, '%2.2f') 'um is %2.2f \n'], marechalStrehl_ltaoH)


