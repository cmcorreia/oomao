%% Add path
%
addpath('../_libOomao/');



%% Load a parameter file
%
param_makeTestSlope



%% Define a calibration source
%
calSource = source;



%% Create an atmosphere
%

atm = atmosphere(photometry.V0,...
    r0,...
    L0,...
    'fractionnalR0',fractionalR0,...
    'altitude',altitude,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);



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
tel = telescope(D,...
    'resolution',nRes,...
    'obstructionRatio',obstructionRatio,...
    'fieldOfViewInArcsec',fieldOfViewInArcsec,...
    'samplingTime',samplingTime);



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
calSource.wavelength = ngsWavelength;
for iGs = 1:nNGS
    ngsWfs(iGs) = shackHartmann(ngsnL,ngsnPx*ngsnL,ngsMinLightRatio);
    ngsWfs(iGs).lenslets.nyquistSampling=ngsNyquistSampling;
    ngsWfs(iGs).lenslets.fieldStopSize=ngsFieldSpotSize;
    calSource = calSource.*tel*ngsWfs(iGs);
    ngsWfs(iGs).INIT;
    +ngsWfs(iGs);
    ngsWfs(iGs).camera.pixelScale = ngsPixelScale;
    
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
tipStep = ngsPixelScale/4;
nStep   = floor(ngsFoV/tipStep/3);
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
Ox_out = sx*ngsPixelScale*constants.radian2arcsec; % in arcsec

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
for iGs=1:nNGS
    ngsWfs(iGs).slopesUnits = 1/slopesLinCoef(1)*ngsPixelScale; % slopes are given in radian (considering optical gain)
end
calSource.zenith = 0;
ngsWfs(1).pointingDirection = [];


%% Create NGSs
%
ngs = source('asterism',ngsAsterism,'magnitude',ngsMagnitude,...
    'wavelength', ngsWavelength);

% If you want to test with laser guide stars... (but no elongation!!)
%{
ngs = source('asterism',ngsAsterism,'magnitude',ngsMagnitude,...
    'wavelength',ngsWavelength,'height',90000);
for iGs=1:nNGS
    ngs(iGs).nPhoton = ngs(iGs).nPhoton*90000.^2;
end
%}



%% Tip-Tilt Removal for phase
zern = zernike(2:3,'resolution',tel.resolution, 'pupil',tel.pupil);
pTT = zern.modes(tel.pupil(:)==1,:);
pTTinv = inv(pTT'*pTT);
pTTremove = @(p) p-pTT*pTTinv*(pTT'*p);
pTTexstract = @(p) pTT*pTTinv*(pTT'*p);



%% Loop Setteing
tel = tel+atm;
for iGs = 1:nNGS
    ngs(iGs) = ngs(iGs).*tel*ngsWfs(iGs);
end
for iGs = 1:nNGS
    ngsWfs(iGs).camera.clockRate = frameTime/samplingTime;
end



%% AO LOOP
slopeStack = zeros(length(ngsWfs(1).slopes),nNGS,nIteration);
for iGs = 1:nNGS
    ngsWfs(iGs).camera.frameCount = 0;
end
for kIteration=1:nIteration
    fprintf('%d/%d\n',kIteration,nIteration);
    tic();
    
    % update atmosphere
    +tel;
    
    % propagation ngs~atm~tel~wfs
    for iGs = 1:nNGS
        +ngs(iGs);
    end
    
    for iGs = 1:nNGS;
        slopeStack(:,iGs,kIteration) = ngsWfs(iGs).slopes(:);
    end
    
    toc();
end

slopeStack=reshape(slopeStack,length(ngsWfs(1).slopes)*nNGS,nIteration);

data.slope = slopeStack;
data.info = rParam;

save testSlope001.mat data