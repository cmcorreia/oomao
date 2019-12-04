%% Add path
%
addpath('../../_libOomao/');


%% load test slope
load testSlope000.mat
S=data.slope*constants.radian2arcsec;


%% Load a parameter file
%
%param_makeTestSlope



%% Define a calibration source
%
calSource = source;



%% Create an atmosphere
%

atm = atmosphere(photometry.V0,...
    data.info.atm.r0,...
    data.info.atm.L0,...
    'fractionnalR0',data.info.atm.fractionalR0,...
    'altitude',data.info.atm.altitude,...
    'windSpeed',data.info.atm.windSpeed,...
    'windDirection',data.info.atm.windDirection);



%% Create telescope
%
% In OOMAO, a grid resolution of telescope (nRes) is related to
% a computation of Shack-Hartmann WFS image.
% The requirements are
%    nL : the number of lenslets on a telescope diameter
%    nF : FieldSpotSize for the shackHartmann object.
%    1. nRes < nF*nL
%    2. rem(nRes,nL)~=0
% Otherwise, you will get error in the sharckHaltmann image computation.
%
tel = telescope(data.info.tel.D,...
    'resolution',data.info.tel.nRes,...
    'obstructionRatio',data.info.tel.obstructionRatio,...
    'fieldOfViewInArcsec',data.info.tel.fieldOfViewInArcsec,...
    'samplingTime',data.info.simu.samplingTime);



%% Create NGS WFSs
%
% Size of diffraction-limited PSF is given by
%     diff = (wavelength of WFS)/(diameter of subaperture)
%
% FoV of WFS is defined as     
%     FoV = diff*(fieldSpotSize);
%
% Pixel scale of WFS is defiend as
%     nPx : number of pixel on diameter of a subaperture
%     pixelScale = FoV/nPx
%
% In the code, at first, the WFS image is created with a
% (nRes*NyquistSampling*2)x(nRes*NyquistSampling*2) grid.
% Then, central (nF*NyquistSampling*2)x(nRes*NyquistSampling*2) images
% are extracted and, finally the extracted image is binned to
% a nPx x nPx image.
% 
calSource.wavelength = data.info.ngs.wavelength;
for iGs = 1:data.info.ngs.nNGS
    ngsWfs(iGs) = shackHartmann(data.info.ngswfs.nL,...
        data.info.ngswfs.nPx*data.info.ngswfs.nL,...
        data.info.ngswfs.minLightRatio);
    ngsWfs(iGs).lenslets.nyquistSampling=data.info.ngswfs.nyquistSampling;
    ngsWfs(iGs).lenslets.fieldStopSize=data.info.ngswfs.fieldSpotSize;
    calSource = calSource.*tel*ngsWfs(iGs);
    ngsWfs(iGs).INIT;
    +ngsWfs(iGs);
    ngsWfs(iGs).camera.pixelScale = data.info.ngswfs.pixelScale;
    
    % if you need a threshold for CoG ...
    %{
    ngsWfs(iGs).framePixelThreshold = ngsThreshold;
    %}
end



%% Optical gain calibration for NGS WFS
% 

% Computing optical gain and linierlity
ngsWfs(1).pointingDirection = zeros(2,1);
calSource = calSource.*tel*ngsWfs(1);
tipStep = data.info.ngswfs.pixelScale/4;
nStep   = floor(data.info.ngswfs.FoV/tipStep/3);
sx      = zeros(1,nStep+1);
u       = 0:nStep;
ngsWfs(1).camera.frameListener.Enabled = false;
ngsWfs(1).slopesListener.Enabled = false;
ngsWfs(1).slopesUnits = 1; % slopes are given in pixel
warning('off','oomao:shackHartmann:relay')
for kStep=u
    calSource.zenith = -tipStep*kStep;
    +calSource;
    drawnow
    sx(kStep+1) = median(ngsWfs(1).slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec; % in arcsec
Ox_out = sx*data.info.ngswfs.pixelScale*constants.radian2arcsec; % in arcsec

% If you want to plot the slope linierlity...
%{
figure(200);
plot(Ox_in,Ox_out,'o');
axis equal;
xlim([0 tipStep*nStep]*constants.radian2arcsec);
xlabel('input [arcsec]');
ylim([0 tipStep*nStep]*constants.radian2arcsec);
ylabel('output [arcsec]');
title('Optical Gain');
hold on;
plot([0 tipStep*nStep]*constants.radian2arcsec,[0 tipStep*nStep]*constants.radian2arcsec,'-');
hold off;
drawnow;
%}

% Set optical gain & pixel scale
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
for iGs=1:data.info.ngs.nNGS
    ngsWfs(iGs).slopesUnits = 1/slopesLinCoef(1)*data.info.ngswfs.pixelScale; % slopes are given in radian (considering optical gain)
end
calSource.zenith = 0;
ngsWfs(1).pointingDirection = [];


%% Create NGSs
%
ngs = source('asterism',data.info.ngs.asterism,'magnitude',data.info.ngs.magnitude,...
    'wavelength', data.info.ngs.wavelength);

% If you want to test with laser guide stars... (but no elongation!!)
%{
ngs = source('asterism',ngsAsterism,'magnitude',ngsMagnitude,...
    'wavelength',ngsWavelength,'height',90000);
for iGs=1:nNGS
    ngs(iGs).nPhoton = ngs(iGs).nPhoton*90000.^2;
end
%}


%% Create turbulence plofiler

tP = newTurbulenceProfiler(tel.D,ngsWfs,ngs,'nLayer',length(data.info.atm.altitude));
%tP = newTurbulenceProfiler(D,ngsWfs,ngs,'atm',atm);


%% set flag

tP.isMap = true;
tP.isXY = false;
tP.flagTT = false;
tP.twoSteps = false;
tP.isFFT = false;
tP.rmDiagonal = true;
tP.isMask4Map = false;
tP.boxWidth4mask = 3;

%% set fit flag
tP.fith = true;
tP.fitr0 = true;
tP.fitL0 = true;

% if the number of layer is 3, you can do...
% tP.fith = true;
% tP.fith = [true false false];

%%

tP.maxIter = 50;
tP.tolFun = 1e-11;

%% set initial values
tP.h = [0 2000 12000];
tP.r0 = [.1 .1 .1];
tP.L0 = 100;

% if the number of layer is 3, you can do...
% tP.h = 1000;
% tP.h = [0 5000 10000];

init = tP.fitParam;


%% set lower and upper bounds
%tP.LBh = 3000;
%tP.LBh = 15000;
tP.LBfr0 = 0;
tP.UBfr0 = Inf;
tP.LBL0 = 5;
tP.UBL0 = 100;
%
% lb = tP.fitLB;
% ub = tP.fitUB;


%% fitting

%tP.fitting(S);
tP.fitting(S,'init',init);
%tP.fitting(S,'init',init,'LB',[],'UB',ub);


%% check fit results
fprintf('\n\n');
fprintf('---------input values---------\n');
for kLayer=1:length(data.info.atm.altitude)
   fprintf('h [%2d] : %5.0f m\n',kLayer,atm.layer(kLayer).altitude);
   fprintf('r0[%2d] : %.4f m\n',kLayer,atm.r0*atm.layer(kLayer).fractionnalR0^(-3/5));
   fprintf('L0[%2d] : %6.2f m\n',kLayer,atm.layer(kLayer).layeredL0);   
end
fprintf('\n\n');
fprintf('---------fit values---------\n');
for kLayer=1:length(data.info.atm.altitude)
   fprintf('h [%2d] : %5.0f m\n',kLayer,tP.h(kLayer));
   fprintf('r0[%2d] : %.4f m\n',kLayer,tP.r0(kLayer));
   fprintf('L0[%2d] : %6.2f m\n',kLayer,tP.L0(kLayer));
end
fprintf('\n\n');
