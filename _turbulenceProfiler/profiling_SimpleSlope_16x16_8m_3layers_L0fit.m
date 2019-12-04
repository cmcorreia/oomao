%% INTRODUCTION
%
% This is a script of the turbulence profiling.
%


%% Add paths

% Add required paths

addpath ../../_libOomao/  % OOMAO codes
addpath ../../_libOomao/_turbulenceProfiler  % OOMAO codes
%addpath ../../_HARMONI/Templates/ % harmoni-Project-ESO-parameters
%addpath ../../_HARMONI/_inputData/ % 35-layers-atmosphere-models


%% PARAMETERS FILE

% Read HARMONI parameters from harmoniLtaoParameters.m if you want.

% harmoniLtaoParameters


%% CALIBRATION SOURCES

% on--axis souce
ngs = source;


%% ATMOSPHERE - overrides predifined TELESCOPE

% If you want to re-define atmospheric parameters here,

% Atmospheric parameters
r0 = 0.15; % coherence lenght in meters at 0.5microns
%L0 = 30;     % Outer scale in meters

% 3-layers atmosphere 
fractionnalR0   = [0.5,  0.3,  0.2]; 
altitude        = [0e3,  5e3, 12e3]; % altitude in km
windSpeed       = [10 ,    5,   20]; % wind speed in m/s
windDirection   = [0  , pi/2,   pi]; % wind direction in angle of radian

% 7-layers atmosphere 
%fractionnalR0   = [0.5960,  0.0963,  0.0325, 0.0372, 0.0869, 0.0684, 0.0826]; 
%altitude        = [0,  .5, 1, 2, 4, 8, 16]*1e3; % altitude in km
%windSpeed       = [7, 8, 9, 11.6, 19.7, 33, 7]; % wind speed in m/s
%windDirection   = [0, 45, 90, 135, 180, 225, 270]*pi/180; % wind direction in angle of radian

% Mono-layer atmosphere 
% fractionalR0    = [1];
% altitude        = [0e3]; % altitude in km
% windSpeed       = [8];   % wind speed in m/s
% windDirection   = [0];   % wind direction in angle of radian

% Create atmosphere object
atm = atmosphere(photometry.V0,... % wavelength = 500nm
    r0,...
    L0,...
    'fractionnalR0',fractionnalR0,...
    'altitude',altitude,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);


%% TELESCOPE - overrides predifined TELESCOPE

% If you want to re-define telescope parameters here,

% Telescope parameters
nL   = 16;                 % number of lenslets
nPx  = 20;                 % number of pixels per lenslet
nRes = nL*nPx;             % resolution on the pupil plane (no of pixels)
D    = 8;                  % telescope primary mirror diameter
d    = D/nL;               % lenslet pitch
samplingFreq = 100;        % WFS sampling time
startDelay = 100;
nFrame = 6000;
obstructionRatio=0.3;      % central obscuration ratio
fieldOfViewInArcsec=100; %fieldOfViewInArcsec
zenithAngle = 0;

% Create telescope object
tel = telescope(D,...
    'resolution',nRes,...
    'obstructionRatio',obstructionRatio,...
    'fieldOfViewInArcsec',fieldOfViewInArcsec,...
    'samplingTime',1/samplingFreq);


%% LGS WAVE-FRONT SENSOR

% Measurement is done in 589nm
ngs.wavelength = photometry.Na;

% Create shackHartmann object for a LGS WFS
wfs = shackHartmann(nL,nRes,0.85);
%wfs.camera.resolution = [160 160];

% WFS initialisation
ngs = ngs.*tel*wfs;
wfs.INIT
+wfs;

% Plots a WFS image
%figure(51);
%imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
%slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;


%% LGS WFS GAIN CALIBRATION 

% The WFS must be calibrated such as for 1rd of tip--tilt wavefront, it
% will measure a slopes of 1rd. To do so, the pointingDirection property in
% shackHartmann object should be set on-axis.

% Set on-axis
wfs.pointingDirection = zeros(2,1);

% Pixel scale for LGS WFS [radian]
FoV = wfs.lenslets.fieldStopSize*ngs.wavelength/d;
pixelScale = FoV/(wfs.camera.resolution(1)/nL);

% Optical gain calibration
tipStep = pixelScale/4;
nStep   = floor(FoV/tipStep/3);
sx      = zeros(1,nStep+1);
u       = 0:nStep;
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = pixelScale*sx*constants.radian2arcsec;

% Plot result
%plot(Ox_in,Ox_out);

% Set optical gain
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);

% back to default
ngs.zenith = 0;
wfs.pointingDirection = [];

% off automatic plot
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;


%% GENERATE ATMOSPHERE

% Connect atmosphere with telescope and generate atmospheric layers
%tel = tel + atm;


%% LGS SOURCES

% number of LGSs
nLGS = 2;

% multiple altitude vector (width 10km)
nLgsNaPoints = 1;
u = linspace(-1,1,nLgsNaPoints);
%lgsHeight = 1e3*(linspace(-5,5,nLgsNaPoints)+90);
lgsHeight = Inf;%1e3*90;

% Double Gaussian profile
%profileNa = (exp(-((u-0.35)*5).^2)+0.65*exp(-((u+0.45)*4).^2));


% Create enlongated source object
%lgsAst = laserGuideStar(tel.D/nL,...
%    tel.D,...
%    90e3/cos(zenithAngle),... % mean alitude
%    1 ,... % fwhmInArcsec
%    2e6,... % number of photon
%    profileNa,... % Sodium density profile
%    'asterism',{[nLGS,arcsec(20),0]},... % asterism for source object
%    'wavelength', photometry.Na,... % wavelength for source object
%    'height',lgsHeight); % height for source object

lgsAst = source('asterism',{[nLGS,arcsec(30),0]},...
    'wavelength', photometry.Na,...
    'height',lgsHeight);

% Equivalent guide star magnutude
ngs.wavelength = photometry.R;
ngs.magnitude = 12;
LGSnPhoton = ngs.nPhoton*2.4;

% Na profile
profileNa = exp(-((u)*2).^2);
profileNa = profileNa/sum(profileNa);
for kGS=1:nLgsNaPoints
    set(lgsAst(:,:,kGS),'nPhoton',LGSnPhoton*profileNa(kGS));
end

%{
% Kernel
ngs.wavelength = photometry.Na;
[x,y] = meshgrid(((0:nPx-1)-(nPx-1)/2)*ngs.wavelength/(2*d)); % Assume Nyquist
r = hypot(x,y);
fwhmLGS = arcsec(1);
sigLGS = fwhmLGS./(2.*sqrt(2*log(2)));
kernel = exp(-r.^2/(2*sigLGS.^2));
set(lgsAst,'extent',kernel/sum(kernel(:)));
%}

% Side-launch
theta = linspace(0,2*pi-2*pi/nLGS,nLGS);
for kGS=1:nLgsNaPoints
    for iGS = 1:nLGS
        lgsAst(1,iGS,kGS).viewPoint = D/2*[cos(theta(iGS)), sin(theta(iGS))];
    end
end


%%
flagTT = true;
S=cell(1,10);
gS=cell(1,10);
for k=1:10
    %load(sprintf('../../../data/slopes_16x16_8m_3layers_%02dm_%02d.mat',L0,k));
    load(sprintf('/result/yono/slopes_16x16_8m_3layers_%02dm_%02d.mat',L0,k));
    S{k} = data.slopes;
    S{k} = bsxfun(@minus, S{k}, mean(S{k},2));
    gS{k} = data.gslopes;
    gS{k} = bsxfun(@minus, gS{k}, mean(gS{k},2));
end

S = cell2mat(S);
gS = cell2mat(gS);


%% Create turbulence plofiler

N=2;
tPall = cell(N,1);

tPall{1} = newTurbulenceProfiler(tel.D,wfs,lgsAst,'nLayer',3,...
    'isMyLMA', true,...
    'isMap', true,...
    'isXY', false,...
    'flagTT',flagTT,...
    'twoSteps',false,...
    'isFFT',true,...
    'rmDiagonal',true,...
    'maxIter',40,...
    'tolFun',1e-15,...
    'isMask4Map',false,...
    'boxWidth4mask',0);

tPall{2} = newTurbulenceProfiler(tel.D,wfs,lgsAst,'nLayer',3,...
    'isMyLMA', true,...
    'isMap', true,...
    'isXY', false,...
    'flagTT',flagTT,...
    'twoSteps',false,...
    'isFFT',false,...
    'rmDiagonal',true,...
    'maxIter',40,...
    'tolFun',1e-15,...
    'isMask4Map',false,...
    'boxWidth4mask',0);

for k=1:N
    tPall{k}.fith = false;
    tPall{k}.fitr0 = true;
    tPall{k}.fitL0 = true;
    
    tPall{k}.h = data.info.altitude;
    tPall{k}.r0 =.2;
    tPall{k}.L0 = 30;
    
    init = tPall{k}.fitParam;
    
    tPall{k}.LBfr0 = 0;
    tPall{k}.UBfr0 = Inf;
    tPall{k}.LBL0 = 0;
    tPall{k}.UBL0 = Inf;
    
    tPall{k}.fitting(S*constants.radian2arcsec,'init',init);
end 

%%
Mobs=zeros(size(tPall{1}.mask));
Mmod1=zeros(size(tPall{1}.mask));
Mmod2=zeros(size(tPall{1}.mask));
Mobs(tPall{1}.mask)=tPall{1}.computeObsCov(S*constants.radian2arcsec);
Mmod1(tPall{1}.mask)=tPall{1}.computeModCov([]);
Mmod2(tPall{1}.mask)=tPall{2}.computeModCov([]);
a=[(-15:15)' Mobs(47,:)' Mmod1(47,:)' Mmod2(47,:)' Mobs(16,:)' Mmod1(16,:)' Mmod2(16,:)'];

if ~flagTT
    save('-ascii',sprintf('result/cross_profile_FFT_vs_Hudgin_%02d_L0fit.txt',L0),'a');
else
    save('-ascii',sprintf('result/cross_profile_FFT_vs_Hudgin_%02d_woTT_L0fit.txt',L0),'a');
end

%%
for k=1:10
    %load(sprintf('../../../data/slopes_16x16_8m_3layers_%02dm_%02d.mat',L0,k));
    load(sprintf('/result/yono/slopes_16x16_8m_3layers_%02dm_%02d.mat',L0,k));
    S = data.slopes;
    S = bsxfun(@minus, S, mean(S,2));
    gS = data.gslopes;
    gS = bsxfun(@minus, gS, mean(gS,2));

    tP{k,1} = newTurbulenceProfiler(tel.D,wfs,lgsAst,'nLayer',3,...
        'isMyLMA', true,...
        'isMap', true,...
        'isXY', false,...
        'flagTT',flagTT,...
        'twoSteps',false,...
        'isFFT',true,...
        'rmDiagonal',true,...
        'maxIter',40,...
        'tolFun',1e-15);
    
    tP{k,2} = newTurbulenceProfiler(tel.D,wfs,lgsAst,'nLayer',3,...
        'isMyLMA', true,...
        'isMap', true,...
        'isXY', false,...
        'flagTT',flagTT,...
        'twoSteps',false,...
        'isFFT',false,...
        'rmDiagonal',true,...
        'maxIter',40,...
        'tolFun',1e-15);
    
    for t=1:N
        tP{k,t}.fith = false;
        tP{k,t}.fitr0 = true;
        tP{k,t}.fitL0 = true;
        
        tP{k,t}.h = data.info.altitude;
        tP{k,t}.r0 =.2;
        tP{k,t}.L0 = 30;
        
        init = tP{k,t}.fitParam;
        
        tP{k,t}.LBfr0 = 0;
        tP{k,t}.UBfr0 = Inf;
        tP{k,t}.LBL0 = 0;
        tP{k,t}.UBL0 = Inf;
        
        tP{k,t}.fitting(S*constants.radian2arcsec,'init',init);
    end
end


%%
er0{1} = [tP{1,1}.r0 tP{2,1}.r0 tP{3,1}.r0 tP{4,1}.r0 tP{5,1}.r0,...
    tP{6,1}.r0 tP{7,1}.r0 tP{8,1}.r0 tP{9,1}.r0 tP{10,1}.r0];
err0{1} = std(er0{1},0,2);
mr0{1} = mean(er0{1},2);
er0{2} = [tP{1,2}.r0 tP{2,2}.r0 tP{3,2}.r0 tP{4,2}.r0 tP{5,2}.r0,...
    tP{6,2}.r0 tP{7,2}.r0 tP{8,2}.r0 tP{9,2}.r0 tP{10,2}.r0];
err0{2} = std(er0{2},0,2);
mr0{2} = mean(er0{2},2);

efr0{1} = [tP{1,1}.fr0 tP{2,1}.fr0 tP{3,1}.fr0 tP{4,1}.fr0 tP{5,1}.fr0,...
        tP{6,1}.fr0 tP{7,1}.fr0 tP{8,1}.fr0 tP{9,1}.fr0 tP{10,1}.fr0];
erfr0{1} = std(efr0{1},0,2);
mfr0{1}  = mean(efr0{1},2);
efr0{2} = [tP{1,2}.fr0 tP{2,2}.fr0 tP{3,2}.fr0 tP{4,2}.fr0 tP{5,2}.fr0,...
    tP{6,2}.fr0 tP{7,2}.fr0 tP{8,2}.fr0 tP{9,2}.fr0 tP{10,2}.fr0];
erfr0{2} = std(efr0{2},0,2);
mfr0{2}  = mean(efr0{2},2);

eL0{1} = [tP{1,1}.L0 tP{2,1}.L0 tP{3,1}.L0 tP{4,1}.L0 tP{5,1}.L0,...
    tP{6,1}.L0 tP{7,1}.L0 tP{8,1}.L0 tP{9,1}.L0 tP{10,1}.L0];
erL0{1} = std(eL0{1},0,2);
mL0{1}  = mean(eL0{1},2);
eL0{2} = [tP{1,2}.L0 tP{2,2}.L0 tP{3,2}.L0 tP{4,2}.L0 tP{5,2}.L0,...
    tP{6,2}.L0 tP{7,2}.L0 tP{8,2}.L0 tP{9,2}.L0 tP{10,2}.L0];
erL0{2} = std(eL0{2},0,2);
mL0{2}  = mean(eL0{2},2);


%% check fit results
inh = zeros(atm.nLayer,1);
inr0 = zeros(atm.nLayer,1);
inL0 = zeros(atm.nLayer,1);

for kLayer=1:length(data.info.altitude)
    inr0(kLayer) = atm.r0*atm.layer(kLayer).fractionnalR0^(-3/5);
    inL0(kLayer) = atm.layer(kLayer).layeredL0;
    inh(kLayer) = atm.layer(kLayer).altitude;
end

fp = fopen(sprintf('result/input_%02d.txt',L0),'w');
for kLayer=1:length(data.info.altitude)
    fprintf(fp,'%5.0f\t' , inh(kLayer));
    fprintf(fp,'%6.4f\t' , inr0(kLayer));
    fprintf(fp,'%7.4f\t' , inr0(kLayer)^(-5/3));
    fprintf(fp,'%6.4f\t' , inr0(kLayer)^(-5/3)/sum(inr0.^(-5/3)));
    fprintf(fp,'%6.2f\t' , inL0(kLayer));
    fprintf(fp,'\n');
end
fclose(fp);

fprintf('\n');

if ~flagTT
    fp = fopen(sprintf('result/output_FFT_vs_Hudgin_%02d_L0fit.txt',L0),'w');
else
    fp = fopen(sprintf('result/output_FFT_vs_Hudgin_%02d_woTT_L0fit.txt',L0),'w');
end

fprintf(fp,'# Map, with TT\n');
for kLayer=1:length(tPall{1}.h)
    for k=1:N
        fprintf(fp,'%5.0f\t' , tPall{k}.h(kLayer));
        fprintf(fp,'%6.4f\t' , tPall{k}.r0(kLayer));
        fprintf(fp,'%6.4f\t' , mr0{k}(kLayer));        
        fprintf(fp,'%6.4f\t' , err0{k}(kLayer));
        fprintf(fp,'%7.4f\t' , (inr0(kLayer)-tPall{k}.r0(kLayer)));
        fprintf(fp,'%7.4f\t' , tPall{k}.fr0(kLayer));
        fprintf(fp,'%6.4f\t' , mfr0{k}(kLayer));        
        fprintf(fp,'%7.4f\t' , erfr0{k}(kLayer));
        fprintf(fp,'%7.4f\t' , (inr0(kLayer)^(-5/3)-tPall{k}.fr0(kLayer)));
        fprintf(fp,'%6.4f\t' , tPall{k}.fr0(kLayer)/sum(tP{k}.fr0));
        fprintf(fp,'%6.2f\t' , tPall{k}.L0(kLayer));
        fprintf(fp,'%6.2f\t' , mL0{k}(kLayer));
        fprintf(fp,'%6.2f\t' , erL0{k}(kLayer));
        fprintf(fp,'%7.4f\t' , (inL0(kLayer)-tPall{k}.L0(kLayer)));
    end
    
    fprintf(fp,'\n');
end
if fp~=1
    fclose(fp);
end