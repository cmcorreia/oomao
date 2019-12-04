%% ADAPTIVE OPTICS MODAL TOMOGRAPHY DEMO
% Demonstrate how to build a modal tomographic reconstructor and the
% optimal projection of the reconstruction on multiples DM conjugated at
% different altitudes for a MCAO system
clear all;
close all;

%% Atmosphere definition
% profile = fitsread('profil_turbulent_eso.fits');
% fractionalR0 = profile(:,4)/100;
seeing =  0.65; 
r0     =  0.976*0.5e-6/arcsec(seeing);
L0     =  25;
% fractionalR0 = fractionalR0/sum(fractionalR0);
% altitude = profile(1:35,2);
% windSpeed = profile(1:35,3);
% windDirection = profile(1:35,9);
fractionalR0   = [0.5 0.3 0.2];
altitude        = [0 4e3 12e3];
windSpeed        = [10 5 7];
windDirection   = [0, pi/2 , pi];
atm = atmosphere(photometry.V0,r0,L0,...
                         'fractionnalR0',fractionalR0,'altitude',altitude,...
                         'windSpeed',windSpeed,'windDirection',windDirection);


%% Telescope VLT
nL   = 40;                % number of lenslets
nPx  = 10;                % number of pixels per lenslet
nRes = nL*nPx;            % resolution on the pupil plane (no of pixels)
D    = 8;                 % telescope primary mirror diameter
d    = D/nL;              % lenslet pitch
samplingFreq = 1000;      % WFS sampling time
obstructionRatio =0.15;   % central obscuration ratio
fieldOfViewInArcsec = 90; %fieldOfViewInArcsec

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);

%% SOURCES
% Calibration sources
ngs = source;
zenithAngle = 0.;
% LGS sources
lgsHeight = 90e3;
nLgs = 4;
spotSize = [];
flux = 1e7;
ast =  ...%source('wavelength',photometry.Na);
 laserGuideStar(d,tel.D, lgsHeight, ...
    spotSize, flux, [],'asterism',{[nLgs,sqrt(2)*arcsec(15)/cos(zenithAngle),0]}, 'wavelength', photometry.Na,'height',lgsHeight);
% Science Sources
x = linspace(-20,20,2);
[x y] = meshgrid(x);
[theta rho] = cart2pol(x,y);
zenith = rho(:)*constants.arcsec2radian;
azimuth = theta;
sci = source('wavelength',photometry.R3,'zenith', zenith, 'azimuth', azimuth);
%sci = source('wavelength',photometry.R3);%,'zenith', zenith, 'azimuth', azimuth);

wpd = ones(size(sci));
gs = ast;
nGs = length(ast);
nSci = length(sci);
nNgs = length(ngs);


%% Wavefront sensor
wfs = shackHartmann(nL,nRes,0.75);
%setValidLenslet(wfs,utilities.piston(nPx))
ngs = ngs.*tel*wfs;
wfs.INIT;
wfs.camera.photonNoise = false;
+wfs;

%% Zernike definition and ZernikeWfs to Zernike Calibration (Dz)
maxRadialDegree = 9;
zern = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nRes,'pupil',tel.pupil);
zern.lex = false;
zern.c = eye(zern.nMode)*ngs.wavelength/2/pi;
ngs=ngs.*zern*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
Dz = z.c;
Dzs = wfs.slopes;
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
zernWfs = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);

%% With noise
ngs.wavelength = photometry.R;
ngs.magnitude = 0;
wfs.camera.readOutNoise = 0;
wfs.camera.photonNoise = true;
wfs.framePixelThreshold = 0;
ngs=ngs.*tel*wfs;

%% noise convariance matrix
nMeas = 250;
slopes = zeros(wfs.nSlope,nMeas);
parfor kMeas=1:nMeas
    +wfs
    slopes(:,kMeas) = wfs.slopes;
end
Cn = slopes*slopes'/nMeas;
figure(5)
subplot(1,2,1)
imagesc(Cn)
axis equal tight
colorbar
wfs.slopes = slopes;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
Czn = z.c*z.c'/nMeas;
subplot(1,2,2)
imagesc(Czn)
axis equal tight
colorbar

%% TOMOGRAPHY (Modal tomographic reconstructor)
tel = tel+atm;
figure
imagesc(tel,[ast])
nLayer = atm.nLayer;
altitudes = [atm.layer.altitude];
% Modal projector at the atmosphere layer altitudes in the guide star directions with piston removal
PistonRM = 1;
projalpha = cell(nLayer,nGs);
[projalpha] = analyticalSmallFootprintExpansion(zernWfs,tel,ast,atm,PistonRM);
% Modal projector at the atmosphere layer altitudes in the science directions with piston removal
projbeta = cell(nLayer,nSci);
[projbeta projbetaCell] = analyticalSmallFootprintExpansion(zernWfs,tel,sci,atm,PistonRM);

%Atmospheric zernike covariance matrix

CxxLayered = cell(nLayer,nLayer);
fr0 = [atm.layer.fractionnalR0];
altitudes = [atm.layer.altitude];
ws = [atm.layer.windSpeed];
wd = [atm.layer.windDirection];



for kLayer = 1:nLayer
    for kLayer1 = 1:nLayer
        if kLayer == kLayer1
            zernWfsi = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nRes,'pupil',tel.pupil);
            zernWfsi.D = tel.diameterAt(altitudes(kLayer));%/tel.D*2
            CxxLayered{kLayer,kLayer1} = phaseStats.zernikeCovariance(zernWfsi,atm)*fr0(kLayer);
        else
            CxxLayered{kLayer,kLayer1} = zeros(zernWfs.nMode);
        end
    end
end

CxxLayered = cell2mat(CxxLayered);
covPhiPup = projalpha*CxxLayered*projalpha';
CznAst = blkdiag( Czn  , Czn , Czn , Czn);
DzAst = blkdiag( Dz  , Dz , Dz  , Dz );

%Tomographic recontructor (GS --> Atmosphere layers)
%Least-Square
iMatLS = projbeta*pinv((projalpha)'*(projalpha),1e-5)*(projalpha)';

%MMSE
Sigma_thetaalpha = CxxLayered*projalpha';
iMatMMSE = Sigma_thetaalpha*pinv( covPhiPup + CznAst , 1e-4 );
%mmse = sparseLMMSE(wfs,tel,atm,ast,'mmseStar',sci,'isTTRM',true);

%% Phase recontruction on the meta pupil at atmosphere layer altitudes
% Propagation of the sources
sci = sci.*tel;
phase = [sci.meanRmPhase];
axis equal tight xy
colorbar

ast = ast.*tel*wfs;
figure
imagesc(tel)

% Phase reconstruction
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
zern.c = Dz\z.c;%(2:end,:);
zern.lex = true;
coeff =  iMatMMSE*zern.c(:);
ScienceCoeff = projbeta*coeff;
zernScience = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);
zernScience.c = ScienceCoeff;
PhaseSience = reshape(zernScience.p*reshape(zernScience.c,(zernModeMax-1),nSci),nRes,nRes*nSci);
figure
imagesc([PhaseSience])
axis equal tight xy
colorbar

coeffLayer = reshape(coeff,zernModeMax-1,nLayer);
coeffLayer1 = coeffLayer(:,1);
coeffLayer2 = coeffLayer(:,2);
coeffLayer3 = coeffLayer(:,3);
zernLayer1  = zernike(2:zernModeMax,'resolution',nRes);
zernLayer2  = zernike(2:zernModeMax,'resolution',floor(nRes*tel.diameterAt(altitude(2))/tel.D));
zernLayer3  = zernike(2:zernModeMax,'resolution',floor(nRes*tel.diameterAt(altitude(3))/tel.D));
zernLayer1.c = coeffLayer1;
zernLayer2.c = coeffLayer2;
zernLayer3.c = coeffLayer3;

phaseLayer1 = reshape(zernLayer1.p*zernLayer1.c,nRes,nRes);
phaseLayer2 = reshape(zernLayer2.p*zernLayer2.c,zernLayer2.resolution,zernLayer2.resolution);
phaseLayer3 = reshape(zernLayer3.p*zernLayer3.c,zernLayer3.resolution,zernLayer3.resolution);


figure
imagesc([phaseLayer1])
axis equal tight xy
colorbar
figure
imagesc([phaseLayer2])
axis equal tight xy
colorbar
figure
imagesc([phaseLayer3])
axis equal tight xy
colorbar

%% Optimal projection for MCAO on the metapupil at DMs altitudes
DMaltitudes = [0 12e3];
nDM = length(DMaltitudes);
projbetaDM = cell(nDM,nSci);
zernDM = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);

% DM projector in the science direction (at the DM altitudes)
[projbetaDM projbetaDMCell] = analyticalSmallFootprintExpansion(zernDM,tel,sci,DMaltitudes,PistonRM);

%  Optimal projection on DMs for MCAO
PbetaDMmean = zeros(zernDM.nMode*nDM,zernDM.nMode*nDM);
PbetaDMLmean = zeros(zernDM.nMode*nDM,zernDM.nMode*nLayer);
for kSci=1:nSci
   intDM = [projbetaDMCell{kSci,:}];
   intL  = [projbetaCell{kSci,:}];
   PbetaDMmean  = PbetaDMmean  +  intDM'*intDM;
   PbetaDMLmean = PbetaDMLmean +  intDM'*intL;
end
PbetaDMmean = PbetaDMmean/nSci;
PbetaDMLmean = PbetaDMLmean/nSci;

ProjOpt =   pinv(PbetaDMmean,1)*PbetaDMLmean;

%% Wavefront reconstruction on the metapupil at DM altitudes
DMCoeff = ProjOpt*coeff;
zernDM = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);
zernDM.c = DMCoeff;
PhaseDM = reshape(zernDM.p*reshape(zernDM.c,(zernModeMax-1),nDM),nRes,nRes*nDM);
figure
imagesc([PhaseDM])
axis equal tight xy
colorbar

%% Wavefront recontruction in the pupil for the Sciences directions through the DM
RecPhaseSciCoeff = projbetaDM*DMCoeff;
zernDM.c = RecPhaseSciCoeff;
PhaseSciDM = reshape(zernDM.p*reshape(zernDM.c,(zernModeMax-1),nSci),nRes,nRes*nSci);
figure
imagesc([PhaseSciDM])
axis equal tight xy
colorbar
figure
imagesc([phase])
axis equal tight xy
colorbar