
%% ADAPTIVE OPTICS MODELING WITH OOMAO
% Demonstrate how to build a simple closed-loop single conjugated adaptive
% optics system with both a geometric diffractive Shack-Hartmann WFS or a
% Pyramid WFS

%% Choose which WFS mode to use
wfsModel = 'diff'; % Options: 'diff', 'geom'

%% Definition of the atmosphere


% The atmosphere class constructor has 2 required input:22
%
% * the wavelength [m]
% * the Fried parameter for the previously given wavelength [m]
%
% 1 optionnal input: [m]
%
% * the outer scale
%
% and  parameter/value pairs of optional inputs:
%
% * the altitudes of the turbulence layers  [m]
% * the fractionnal r0 which is the ratio of the r0 at altitude h on the
% integrated r0: $(r_0(h)/r_0)^{5/3}$
% * the wind speeds  [m/s]
% * the wind directions [rd]
%
% In the following the atmosphere is given for a r0=15cm in V band and an
% outer scale of 30m with 3 turbulence layers.
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.25,0.05],...
    'windSpeed',[5,10,20],...
    'windDirection',[0,pi/4,pi]);

%% Definition of the telescope
% The telescope class constructor has 1 required input:
%
% * the telescope diameter [m]
%
% 1 optionnal input:
%
% * the central obscuration ratio
%
% and  parameter/value pairs of optional inputs:
%
% * the field of view either in arcminute or arcsecond
% * the pupil sampling or resolution in pixels
% * the atmopheric layer motion sampling time [s]
nPx = 600;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/250);

%%
nOuterLenslet = 4;
nInnerLenslet = 10;
nLenslet = nInnerLenslet + nOuterLenslet;

tel.embed(nPx*nLenslet/nInnerLenslet);
nPx = nPx*nLenslet/nInnerLenslet;
%% Definition of a calibration source
% The source class constructor has parameters/value pairs of optional inputs:
%
% * the zenith angle [rd] (default: 0rd)
% * the azimuth angle [rd] (default: 0rd)
% * the wavelength [m] (default: V band)
% * the magnitude
%
% In the following, an on-axis natural guide star in V band is defined.
ngs = source('wavelength',photometry.V);

% and for later use a science object in H band is instantiated
science = source('wavelength',photometry.H);


%% Definition of the wavefront sensor
% Up to now, only the Shack--Hartmann WFS has been implemented in OOMAO.
% The shackHartmann class constructor has 2 required inputs:
%
% * the number of lenslet on one side of the square lenslet array
% * the resolution of the camera
%
% and  1 optionnal input:
%
% * the minimum light ratio that is the ratio between a partially
% illuminated subaperture and a fully illuminated aperture
d = tel.D/nLenslet;
wfs = shackHartmann(nLenslet,tel.resolution,0.5);
%wfs = pyramid(nLenslet,nPx,'modulation',5); % increase modulation to avoid loss of performance due to small linear range
%wfs = pyramid(nLenslet,nPx,'modulation',5, 'src', ngs,'tel',tel,'minLightRatio', 0.05);
%%
% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;
%%
% Selecting the subapertures illuminated at 75% or more with respect to a
% fully illuminated subaperture
% setValidLenslet(wfs)
% %%
% % A new frame read-out and slopes computing:
% +wfs;
% %%
% % Setting the reference slopes to the current slopes that corresponds to a
% % flat wavefront
% wfs.referenceSlopes = wfs.slopes;
wfs.INIT
%%
% A new frame read-out and slopes computing:
+wfs;
%%
% The WFS camera display:
figure
subplot(1,2,1)
imagesc(wfs.camera)
%%
% The WFS slopes display:
subplot(1,2,2)
slopesDisplay(wfs)

%% GAIN CALIBRATION
if ~strcmp(wfs.tag,'PYRAMID')
    wfs.gainCalibration(tel,ngs)
end
%% Combining the atmosphere and the telescope
tel = tel+atm;
figure(10)
imagesc(tel)
%% prefocal station with one natural source
offset = [-1.5;0];
pfs = preFocalStation(tel,ngs,offset,'rotation',180);

% propagate down to the WFS
ngs = ngs.*tel;
phaseBefore = ngs.phase;

ngs = ngs.*tel*pfs;
phaseAfter = ngs.phase;

figure(10)
imagesc([phaseBefore phaseAfter])



%% STACK OF DMs
nDm             = 1;
dmPitch         = d;
dmCrossCouplingCoeff = 0.25;

bif(1)          = gaussianInfluenceFunction(dmCrossCouplingCoeff,dmPitch(1));


nActuators      = [nLenslet+1]; % actuators across the DM diameter

dmHeights       = [0];

if 1 %extendedFried % create a scrap DM just for the definition if the extended valid actuator map
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

dmLowRes = deformableMirror(nActuators,'modes',bifLowRes,...
    'resolution',nLenslet+1,...
    'validActuator',validActMap, 'zLocation',0 );
reconstructionGrid = wfs.reconstructionGrid(wfs,1);


%%
big          = gaussianInfluenceFunction(dmCrossCouplingCoeff,dmPitch(1));

nActuator = nLenslet + 1;
dm = deformableMirror(nActuator,...
    'modes',big,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator);

%% CHECK OFFSETS APPLIED TO ATM PHASE AND DM PHASE ALTOGETHER
nParalelDmCommands = 2;
dm.coefs = randn(dm.nValidActuator,nParalelDmCommands)*ngs.wavelength/2/pi*10;
% propagate down to the WFS
ngs = ngs.*tel*dm;
phaseBefore = reshape(ngs.phase,nPx,nPx*nParalelDmCommands);

ngs = ngs.*tel*dm*pfs;
phaseAfter = reshape(ngs.phase,nPx,nPx*nParalelDmCommands);

figure(10)
imagesc([phaseBefore phaseAfter])


%% Interaction matrix: DM/WFS calibration
tel = tel - atm;

wfs.lenslets.offset = 0;
wfs.lenslets.rotation = 0;

ngs=ngs.*tel;
dmCalib = calibration(dm,wfs,ngs,ngs.wavelength/100,10);

%%
wfs.lenslets.offset = 0;
wfs.lenslets.rotation = 0;

wfs.lenslets.offset = -1/5*[d;d]/tel.D;
wfs.lenslets.rotation = pi/6;
ngs=ngs.*tel;
dmCalib = calibration(dm,wfs,ngs,ngs.wavelength/100,10);


dmCalib.nThresholded = 5;
if strcmp(wfsModel,'geom')
    commandMatrix = dmCalibG.M;
else
    commandMatrix = dmCalib.M;
end
%% The closed loop
% Combining the atmosphere and the telescope
tel = tel+atm;
figure(10)
imagesc(tel)
%%
% Resetting the DM command
dm.coefs = 0;
%%
% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;
%%
% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase;
%%
% Propagation to the WFS
ngs=ngs*dm*wfs;
%%
% Display of turbulence and residual phase
figure(11)
h = imagesc([turbPhase,ngs.meanRmPhase]);
axis equal tight
colorbar
snapnow
%%
zern = zernike(2:10,tel.D,'resolution',tel.resolution,'pupil',tel.pupil);
%%
% Closed loop integrator gain:
loopGain = 0.5;
gain_pol = 0.9;
%%
% closing the loop with an integrator controller 
nIteration = 200;
total  = zeros(1,nIteration);
residue = zeros(1,nIteration);
totalProj = zeros(zern.nMode,nIteration);
resProj = zeros(zern.nMode,nIteration);

tic
for kIteration=1:nIteration
    % Propagation throught the atmosphere to the telescope, +tel means that
    % all the layers move of one step based on the sampling time and the
    % wind vectors of the layers
    ngs=ngs.*+tel;
    % Saving the turbulence aberrated phase
    turbPhase = ngs.meanRmPhase;
    % Variance of the atmospheric wavefront
    total(kIteration) = var(ngs);
    % Projection onto Zernike Modes
    zern\turbPhase(:);
    totalProj(:,kIteration) = zern.c;
    % Propagation to the WFS
    if strcmp(wfsModel,'geom')
        ngs=ngs*dm*wfsG;
    else
        ngs=ngs*dm*wfs;
    end
    % Variance of the residual wavefront
    residue(kIteration) = var(ngs);
        % Projection onto Zernike Modes
    zern\ngs.meanRmPhase(:);
    resProj(:,kIteration) = zern.c;
    
    % Computing the DM residual coefficients
    if strcmp(wfsModel,'geom')
        residualDmCoefs = commandMatrix*wfsG.slopes;
    else
        residualDmCoefs = commandMatrix*wfs.slopes;
    end
    
    % Integrating the DM coefficients
    dm.coefs = dm.coefs - loopGain*residualDmCoefs;
    
    % POL 
    %dm.coefs = (1-gain_pol)*dm.coefs - gain_pol*commandMatrix*(wfs.slopes - dmCalib.D*dm.coefs);

    % Display of turbulence and residual phase
    set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
    drawnow
end
toc
snapnow
u = (0:nIteration-1).*tel.samplingTime;
atm.wavelength = ngs.wavelength;
%%
% Piston removed phase variance
totalTheory = phaseStats.zernikeResidualVariance(1,atm,tel);
atm.wavelength = photometry.V;
%%
% Phase variance to micron rms converter 
rmsMicron = @(x) 1e6*sqrt(x).*ngs.wavelength/2/pi;
figure(12)
plot(u,rmsMicron(total),u([1,end]),rmsMicron(totalTheory)*ones(1,2),u,rmsMicron(residue))
grid
legend('Full','Full (theory)','Residue','Location','Best')
xlabel('Time [s]')
ylabel('Wavefront rms [\mum]')

%% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN
varFit      = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
varAlias    = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
varTempo    = phaseStats.closedLoopVariance(atm, tel.samplingTime,0*tel.samplingTime,loopGain)*(atm.wavelength/science(1).wavelength)^2;
marechalStrehl_lsq_theoretical = 100*exp(-varFit-varAlias-varTempo)

figure(12)
text(0.5,1,['Theoretical SR (Marechal approx):' num2str(marechalStrehl_lsq_theoretical) '%'])
text(0.5,0.9,['Empirical SR (Marechal approx):' num2str(...
    100*exp(-(mean(rmsMicron(residue(50:end)))*1e-6/ngs.wavelength*2*pi)^2*(atm.wavelength/science(1).wavelength)^2)) '%'])
ylim([0 1])


%% pfs with a LGS

% LGS sources for tomography
zenithAngle = pi/6;
radiusAst   = 30;
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
    lgsAst(iGS).viewPoint = tel.D/2*[cos(theta(iGS)), sin(theta(iGS))];
end

offset = 0.5;
offset = [-1.9;0.6];

pfs = preFocalStation(tel,lgsAst,offset);


% propagate down to the WFS
tel = tel + atm;
lgsAst = lgsAst.*tel;
phaseBefore = catPhase(lgsAst);

lgsAst = lgsAst.*tel*pfs;
phaseAfter = catPhase(lgsAst);
figure(10)
imagesc([phaseBefore; phaseAfter])

    
    
%%
wfs.lenslets.offset = 0;
wfs.lenslets.rotation = 0;

wfs.lenslets.offset = -1/5*[d;d]/tel.D;
wfs.lenslets.rotation = pi/6;


%%
science = source('wavelength',photometry.H);


lgsAst = lgsAst.*tel*wfs;
%        if os == 1
sparseMmse = sparseLMMSE(wfs,tel,atm,lgsAst,...
    'iNoiseVar', 1/1e-14,...
    'layerEstimation',true,...
    'mmseStar',science,...
    'outputWavefrontMask', reconstructionGrid,'overSampling',1);%,'phaseCovariance','real');


%% SCIENCE CAMERAS

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

FlowRes         = 2*bifLowRes(1).modes(reconstructionGrid,:);
iFlowRes        = pinv(full(FlowRes),1e-3);

fittingMatrix   = iFlowRes*sparseMmse.Hss;

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
                dm.coefs = (1-gain_pol)*dm.coefs + ...
                    gain_pol*fittingMatrix*( sparseMmse*( bsxfun( @minus, wfs.slopes, dmCalib.D*dm.coefs ) ) );
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


