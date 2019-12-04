clear all;
close all;

%Atmosphere parameters
r0            = 0.15;
L0            = 30;
fractionnalR0 = [0.75,0.25];
layerAlt      = [0,10].*1e3;
windSpeed     = [10,25];
windDirection = [0.,0.];
%System parameters
D             = 4.2;
nL            = 7;
d             = D/nL;
obs           = 0.2;
foV           = 300;
Fe            = 150; 
nPxGs         = 10;
nRes          = nPxGs*nL;
% --------------------- CREATING SUB-CLASSES ---------------------------%%
            
atm = atmosphere(photometry.V0,r0,L0,'fractionnalR0',fractionnalR0,'altitude',layerAlt,...
                'windSpeed',windSpeed,'windDirection',windDirection);         
            
tel = telescope(D,'resolution',nRes,'obstructionRatio',obs ,...
    'fieldOfViewInArcsec',foV ,'samplingTime',1/Fe);

%On-axis target
sci = source();

%Guide strs asterism
ast = source('asterism',{[3,1*constants.arcmin2radian,pi/3]});
%ast = source('zenith',10/constants.radian2arcsec);
%WFS type, defined only one is enough
wfs = shackHartmann(nL,nRes,0.75);
sci = sci.*tel*wfs;
wfs.INIT;
+wfs;
wfs.slopesUnits = calibrateWfsPixelScale(wfs,sci,tel.D,nPxGs,0);

% Creating a command matrix
bif                  = gaussianInfluenceFunction(0.35,d);
dm                   = deformableMirror(nL+1,'modes',bif,'resolution',nRes,'validActuator',wfs.validActuator);
dmCalib              = calibration(dm,wfs,sci,sci.wavelength/40);
dmCalib.cond         = 30;
ps                   = constants.radian2arcsec*wfs.slopesUnits*sci.wavelength/(2*d*wfs.lenslets.nyquistSampling);

%% Creating the profiler
p  = turbulenceProfiler(tel,atm,wfs,[sci,ast],'generatorFlag','simulation',...
    'fitr0',1,'fith',1,'fitL0',1,'nFrames',10000);
            
p.mod.Projector      = dmCalib.M;
p.mod.ProjectorUnits = 1e9/ps/sqrt(dm.nValidActuator);


%% LOOP ON FRAMES
nk = 2;
% Instantiation
wfeTomoRaw = zeros(nk,1);
wfeTomo0   = zeros(nk,1);
wfeTomo1   = zeros(nk,1);
wfeTomo2   = zeros(nk,1);
X1         = zeros(6,nk);
X2         = zeros(6,nk);
% Splitting into two parts
s          = p.sGen.slopes;
nf         = size(s,2);
sprof      = s(:,1:floor(end/2));
% Computing the covariance matrix of test
stest      = s(:,1+floor(end/2):end);
stest      = stest - repmat(mean(stest,2),1,size(stest,2));
Ctt        = stest*stest'/size(stest,2);
% Minimal raw tomographic error
R          = p.mod.getMMMSEreconstructor(Ctt);
Cee        = p.mod.getErrorCovarianceMatrix(Ctt,R);
wfeTomoMin = p.mod.getWaveFrontError(Cee);

% Minimal tomographic error with a fitting procedure
p.mod.atm  = atm;
p.mod      = p.mod.atmToInit();
p.mod      = p.mod.dataToObj(p.mod.parameters);
Css0       = p.mod.covarianceModel(0);
R0         = p.mod.getMMMSEreconstructor(Css0);
Cee        = p.mod.getErrorCovarianceMatrix(Ctt,R0);
wfeTomo0   = p.mod.getWaveFrontError(Cee);
    
%% 

for k=1:nk
    % ---------------------- Truncating the slopes ---------------------- %
    tk      = floor(k*nf/2/nk);
    sk      = sprof(:,1:tk);    
    % ------------- Getting the covariance matrix of slopes ------------- %
    sk      = sk - repmat(mean(sk,2),1,size(sk,2));
    Cssk    = sk*sk'./size(sk,2);
    p.Css   = Cssk;
    
    % Raw tomographic error with a fitting procedure
    Rraw            = p.mod.getMMMSEreconstructor(Cssk);
    Cee             = p.mod.getErrorCovarianceMatrix(Ctt,Rraw);
    wfeTomoRaw(k)   = p.mod.getWaveFrontError(Cee);    
    
    
    % -------- Performing the covariance fitting: 1 step procedure ------ %
    % Set the initial guess
    p.X0    = [20,1,0,3000,20,20];
    % Defining the fitting options        
    p       = p.doLsq1Step();   
    X1(:,k) = p.X(1:6);        
      
    % Error tomographic using the 1 step procedure
    p.mod.atm       = p.atmFit;
    p.mod           = p.mod.atmToInit();
    p.mod           = p.mod.dataToObj(p.mod.parameters);
    Css1            = p.mod.covarianceModel(0);
    R1              = p.mod.getMMMSEreconstructor(Css1);
    Cee             = p.mod.getErrorCovarianceMatrix(Ctt,R1);
    wfeTomo1(k)     = p.mod.getWaveFrontError(Cee);
    
    
    % ------ Performing the covariance fitting: 2 steps procedure ------- %
    % Set the initial guess
    p.X0    = [20,1,0,3000,20,20];
    p2      = p.doLsq2Steps();
    X2(:,k) = p2.X(1:6);   
          
    % Error tomographic using the 2 steps procedure
    p2.mod.atm      = p2.atmFit;   
    p2.mod          = p2.mod.atmToInit();
    p2.mod          = p2.mod.dataToObj(p2.mod.parameters);
    Css2            = p2.mod.covarianceModel(0);
    R2              = p2.mod.getMMMSEreconstructor(Css2);
    Cee             = p2.mod.getErrorCovarianceMatrix(Ctt,R2);
    wfeTomo2(k)     = p2.mod.getWaveFrontError(Cee);
end