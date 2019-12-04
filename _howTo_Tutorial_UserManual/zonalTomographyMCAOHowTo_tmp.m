%% Atmosphere definition

% profile = fitsread('profil_turbulent_eso.fits');
% fractionalR0 = profile(:,4)/100;
seeing =  0.8;
zenithAngle = pi/6;
r0     =  0.976*0.5e-6/arcsec(seeing)*cos(zenithAngle)^(3/5);
L0     =  25;
% fractionalR0 = fractionalR0/sum(fractionalR0);
% altitude = profile(1:35,2);
% windSpeed = profile(1:35,3);
% windDirection = profile(1:35,9);
fractionalR0      = [0.59 0.02 0.04 0.06 0.01 0.05 0.09 0.04 0.05 0.05];
altitude             = 1/cos(zenithAngle)*[30 140 281 562 1125 2250 4500 7750 11000 14000];
windSpeed        = [10 13 5 8 17 3 14 2 12 5];
windDirection   = 2*pi*rand(1,10);
atm = atmosphere(photometry.V0,r0,L0,...
                         'fractionnalR0',fractionalR0,'altitude',altitude,...
                         'windSpeed',windSpeed,'windDirection',windDirection);


%% Telescope
nL   = 32;                % number of lenslets
nPx  = 10;                % number of pixels per lenslet
nRes = nL*nPx;            % resolution on the pupil plane (no of pixels)
D    = 8;                   % telescope primary mirror diameter
d    = D/nL;              % lenslet pitch
samplingFreq = 1000;      % WFS sampling time
obstructionRatio = 0.15;  % central obscuration ratio
fieldOfViewInArcsec = 43; %fieldOfViewInArcsec

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);

%% SOURCES
% Calibration sources
ngs = source;

% Optimization source for the MCAO projector
x = linspace(-15,15,11);
[x y] = meshgrid(x);
[theta rho] = cart2pol(x,y);
zenithOpt = rho(:)*constants.arcsec2radian;
azimuthOpt = theta;
opt = source('zenith', zenithOpt, 'azimuth', azimuthOpt);

% LGS source for tomopgraphy
radiusAst = 15;
nAst = 4;
thetaAst = 0;
% LGS sources
lgsHeight = 90e3;
spotSize = [];
flux = 25e6; % 1000ph/frame/spup (40x40 spup, 1kHz, 0.2 side square)
flux2 = 25e5;% 100ph/fr/sap
% % LGS
ast = laserGuideStar(d,tel.D, lgsHeight/cos(zenithAngle), ...
    spotSize, flux, [],'asterism',{[nAst,arcsec(radiusAst),thetaAst]}, 'wavelength', photometry.Na,'height',lgsHeight/cos(zenithAngle));
ast1 = laserGuideStar(d,tel.D, lgsHeight/cos(zenithAngle), ...
    spotSize, flux, [], 'wavelength', photometry.Na,'height',lgsHeight/cos(zenithAngle));
%ast1 = source('wavelength',photometry.Na);
test = [ast];

% Science Sources for perf estimation
% rhoSci    = linspace(0,60,31);
xs = linspace(-15,15,11);
[xs ys] = meshgrid(xs);
[thetas rhos] = cart2pol(xs,ys);
zenithSci     = rhos(:)*constants.arcsec2radian;
azimuthSci  = thetas;
% azimuthSci1  = 0*ones(31,1);
% rho = cat(2,rhoSci,rhoSci);
% azim = cat(2,azimuthSci1,azimuthSci);
%zenithSci= rho*constants.arcsec2radian;

sciH = source('zenith', zenithSci, 'azimuth', azimuthSci,'wavelength',photometry.H);
sciR = source('zenith', zenithSci, 'azimuth', azimuthSci,'wavelength',photometry.R);

%wpd = ones(size(sci));
gs = ast;
nGs = length(test);
nSciR = length(sciR);
nSciH = length(sciH);
nNgs = length(ngs);
%% Wavefront sensor
wfs = shackHartmann(nL,nRes,0.75);
ngs.wavelength = photometry.Na;
%setValidLenslet(wfs,utilities.piston(nPx))
ngs = ngs.*tel*wfs;
wfs.INIT;
wfs.camera.photonNoise = false;
+wfs;
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
wfs.gainCalibration(tel,ngs);


%% DMs Stack
pitchs = [0.25 0.4];
bif  = gaussianInfluenceFunction(25/100,pitchs(1));
bif2 = gaussianInfluenceFunction(25/100, pitchs(2));
%bif3 = gaussianInfluenceFunction(25/100, pitchs(3));
bifLowRes = gaussianInfluenceFunction(25/100, pitchs(1));
bif2LowRes = gaussianInfluenceFunction(25/100,pitchs(2));
%bif3LowRes = gaussianInfluenceFunction(25/100,pitchs(3));
% dmStack = multiDeformableMirror([0], [41], tel, ...
%      'modes',[bif],...
%      'validActuator',{true(41)});
nActuators =  [33 25]; % 33 29 pour 60 et 33 37 pour 120

dmStack = multiDeformableMirror([0 8e3], nActuators , tel, ...
   'modes',[bif bif2],...
   'validActuator',{true(33) true(25)} );
nDm = dmStack.nDm;
nValidActuatorsStack = dmStack.nValidActuators;
% DM calibration
test = test.*tel*wfs;
calibDms = calibrationMultiDm(dmStack,wfs,test,test(1).wavelength/100,25,'cond',1e2);

% %% Tomographic reconstructor 
mmse = sparseLMMSE(wfs,tel,atm,test,...
'iNoiseVar', 1/1e-14,...
'layerEstimation',true);

% %% Fitting
F = fittingMatrix(dmStack,mmse,opt,[bifLowRes bif2LowRes]);
%% Science camera
camR = imager('nyquistSampling',2,'fieldStopSize',75);
sciR = sciR.*tel*camR;
camR.referenceFrame = camR.frame;

camH = imager('nyquistSampling',2,'fieldStopSize',75);
sciH = sciH.*tel*camH;
camH.referenceFrame = camH.frame;

tel = tel + atm;

%% Loop initialization
reset(tel);
%% With noise on the WFS
% wfs.camera.readOutNoise = 0;
% wfs.camera.photonNoise = true;
% wfs.framePixelThreshold = 0;

test = test.*tel*dmStack*wfs;
%sci = sci.*tel*dmStack*cam;
flush(camR)
flush(camH)
camR.frame = camR.frame*0;
camH.frame = camH.frame*0;
camR.clockRate    = 1;
camH.clockRate    = 1;
exposureTime = 10000;
camR.exposureTime = exposureTime;
camH.exposureTime = exposureTime;
startDelay = 99;
camR.startDelay   = startDelay;
camH.startDelay   = startDelay;
flush(camR)
flush(camH)
set(sciR,'logging',true)
set(sciR,'phaseVar',[])
set(sciH,'logging',true)
set(sciH,'phaseVar',[])
nIt = camR.startDelay+camR.exposureTime;
dmStack.dms{1}.coefs = 0;
dmStack.dms{2}.coefs = 0;
%dmStack.dms{3}.coefs = 0;

%% Display initialization
%f = figure;
% subplot(2,3,1);
% h1 = imagesc(zeros(nRes*nSci));
% axis xy equal tight
% title('Input WF');
% h1.Parent.XTickLabel = '';
% h1.Parent.YTickLabel = '';
% 
% subplot(2,3,2);
% h2 = imagesc(zeros(nRes*nSci));
% axis xy equal tight
% title('Residual WF');
% h2.Parent.XTickLabel = '';
% h2.Parent.YTickLabel = '';  
% % 
% subplot(2,3,3);
% h3 = imagesc(zeros(size(cam.referenceFrame)));
% title('PSF');
% axis equal tight xy
% h3.Parent.XTickLabel = '';
% h3.Parent.YTickLabel = '';
% 
% subplot(2,3,4);
% h4 = imagesc(zeros(size(dmStack.dms{1}.surface,1)));
% axis xy equal tight
% title('DM1');
% h4.Parent.XTickLabel = '';
% h4.Parent.YTickLabel = '';
% 
% 
% subplot(2,3,5);
% h5 = imagesc(zeros(size(dmStack.dms{2}.surface,1)));
% axis xy equal tight
% title('DM2');
% h5.Parent.XTickLabel = '';
% h5.Parent.YTickLabel = '';
% 
% 
% subplot(2,3,6);
% h6 = imagesc(zeros(size(dmStack.dms{3}.surface,1)));
% axis xy equal tight
% title('DM3');
% h6.Parent.XTickLabel = '';
% h6.Parent.YTickLabel = '';
% 
% % MCAO loop (PSOL)
% % PSEUDO OPEN-LOOP CONTROLLER
gain_pol = 0.5;
slopesPsolCell = cell(nDm,nGs);
slopesPsolMat  = zeros(wfs.nSlope,nGs);
CommandDMs     = zeros(nValidActuatorsStack,1);
%PhaseSci       = zeros(nIt,nRes,nRes*nSci);
%ResSci         = zeros(nIt,nRes,nRes*nSci);
surfDM         = cell(nDm,round(nIt/1000));
inWFR  = cell(round(nIt/1000),4);
inWFH  = cell(round(nIt/1000),4);
resWFR = cell(round(nIt/1000),4);
resWFH = cell(round(nIt/1000),4);
wfeH = zeros(camH.startDelay + camH.exposureTime,nSciH);
srwfeH =zeros(camH.startDelay + camH.exposureTime,nSciH) ;
wfeR = zeros(camR.startDelay + camR.exposureTime,nSciR);
srwfeR =zeros(camR.startDelay + camR.exposureTime,nSciR) ;
for k=1:nIt
    tic
    fprintf('Frame %d/%d\n',k,nIt);
    
    +tel;
    +test;
     sciH=sciH.*tel;
     sciR=sciR.*tel;
     %set(h1,'CData',sci.catMeanRmPhase/sci(1).waveNumber);
     if mod(k,1000)==0
         inWFR{round(k/1000),1} = catMeanRmPhase(sciR(61))/sciR(1).waveNumber;
         inWFR{round(k/1000),2} = catMeanRmPhase(sciR(83))/sciR(1).waveNumber;
         inWFR{round(k/1000),3} = catMeanRmPhase(sciR(94))/sciR(1).waveNumber;
         inWFR{round(k/1000),4} = catMeanRmPhase(sciR(116))/sciR(1).waveNumber;
         inWFH{round(k/1000),1} = catMeanRmPhase(sciH(61))/sciH(1).waveNumber;
         inWFH{round(k/1000),2} = catMeanRmPhase(sciH(83))/sciH(1).waveNumber;
         inWFH{round(k/1000),3} = catMeanRmPhase(sciH(94))/sciH(1).waveNumber;
         inWFH{round(k/1000),4} = catMeanRmPhase(sciH(116))/sciH(1).waveNumber;
         surfDM{1,round(k/1000)} = -2*dmStack.dms{1}.surface*2*pi/sciH(1).wavelength;
         surfDM{2,round(k/1000)} = -2*dmStack.dms{2}.surface*2*pi/sciH(1).wavelength;
         %surfDM{3,round(k/100)} = -2*dmStack.dms{3}.surface*2*pi/sciH(1).wavelength;
     end
    sciH=sciH*dmStack*camH;
    sciR=sciR*dmStack*camR;
    if mod(k,1000)==0
         resWFR{round(k/1000),1}  = catMeanRmPhase(sciR(61))/sciR(1).waveNumber;
         resWFR{round(k/1000),2}  = catMeanRmPhase(sciR(83))/sciR(1).waveNumber;
         resWFR{round(k/1000),3} = catMeanRmPhase(sciR(94))/sciR(1).waveNumber;
         resWFR{round(k/1000),4} = catMeanRmPhase(sciR(116))/sciR(1).waveNumber;
         resWFH{round(k/1000),1}  = catMeanRmPhase(sciH(61))/sciH(1).waveNumber;
         resWFH{round(k/1000),2}  = catMeanRmPhase(sciH(83))/sciH(1).waveNumber;
         resWFH{round(k/1000),3} = catMeanRmPhase(sciH(94))/sciH(1).waveNumber;
         resWFH{round(k/1000),4} = catMeanRmPhase(sciH(116))/sciH(1).waveNumber;
    end
         
    resWFHperf = reshape(catMeanRmPhase(sciH)/sciH(1).waveNumber,nRes*nRes,nSciH);
    wfeH(k,:) = std(resWFHperf(tel.pupilLogical(:),:),1)*1e+9;
    srwfeH(k,:) = 1e2*exp(-(wfeH(k,:)*sciH(1).waveNumber*1e-9).^2);
    resWFRperf = reshape(catMeanRmPhase(sciR)/sciR(1).waveNumber,nRes*nRes,nSciR);
    wfeR(k,:) = std(resWFRperf(tel.pupilLogical(:),:),1)*1e+9;
    srwfeR(k,:) = 1e2*exp(-(wfeR(k,:)*sciR(1).waveNumber*1e-9).^2);
    %set(h2,'CData',sci.catMeanRmPhase/sci(1).waveNumber);
    %ResSci(k,:,:) = [sci.catMeanRmPhase/sci(1).waveNumber];

    for kStar = 1: nGs
        for kDm = 1:nDm
            if kDm == 1
                slopesPsolCell{kDm,kStar} = calibDms{kDm,1}.D*dmStack.dms{kDm}.coefs;
            else
                slopesPsolCell{kDm,kStar} = calibDms{kDm,kStar}.D*dmStack.dms{kDm}.coefs;
            end
        end
        slopesPsolMat(:,kStar) = sum([slopesPsolCell{:,kStar}],2);
    end
    
    CommandDMs = (1-gain_pol)*CommandDMs + ...
        gain_pol*F.CommandMatrix*(mmse*bsxfun( @minus, wfs.slopes, slopesPsolMat ));
    
    dmStack.dms{1}.coefs = CommandDMs(1 : dmStack.dms{1}.nValidActuator,1);
    dmStack.dms{2}.coefs = CommandDMs(dmStack.dms{1}.nValidActuator +(1: dmStack.dms{2}.nValidActuator),1);
    %dmStack.dms{3}.coefs = CommandDMs(dmStack.dms{1}.nValidActuator+dmStack.dms{2}.nValidActuator +1 : end,1);
    
    %set(h4,'CData',-2*dmStack.dms{1}.surface*2*pi/sci(1).wavelength);
    %set(h5,'CData',-2*dmStack.dms{2}.surface*2*pi/sci(1).wavelength);
%      set(h6,'CData',-2*dmStack.dms{3}.surface*2*pi/sci(1).wavelength);
    

%     if cam.frameCount == 0 || isscalar(cam.frameBuffer)
%         set(h3,'CData',cam.frame);
%     else
%         set(h3,'CData',cam.frameBuffer);
%     end
%     
%     drawnow;
%     frame = getframe(f);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     
%     % Write to the GIF File
%     if k == 1
%         imwrite(imind,cm,'mcaoPsolConf1.gif','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'mcaoPsolConf1.gif','gif','WriteMode','append');
%     end
toc
end
var_wfe_mcaoH = zeros(nIt,2,nSciH);
wfe_mcaoH =zeros(nIt,2,nSciH) ;
marechalStrehl_mcaoH = zeros(nSciH,1);
psfStrehl_mcaoH =1e2*camH.strehl;
var_wfe_mcaoR = zeros(nIt,2,nSciR);
wfe_mcaoR =zeros(nIt,2,nSciR) ;
marechalStrehl_mcaoR = zeros(nSciR,1);
psfStrehl_mcaoR =1e2*camR.strehl;
for k = 1: nSciR
    var_wfe_mcaoR(:,:,k) = reshape(sciR(k).phaseVar(1:nIt*2),2,[])';
    wfe_mcaoR(:,:,k) = sqrt(var_wfe_mcaoR(:,:,k))/sciR(1).waveNumber*1e6;
    marechalStrehl_mcaoR(k,1) = 1e2*exp(-mean(var_wfe_mcaoR(startDelay:end,2,k)));
end
for k = 1: nSciH
    var_wfe_mcaoH(:,:,k) = reshape(sciH(k).phaseVar(1:nIt*2),2,[])';
    wfe_mcaoH(:,:,k) = sqrt(var_wfe_mcaoH(:,:,k))/sciH(1).waveNumber*1e6;
    marechalStrehl_mcaoH(k,1) = 1e2*exp(-mean(var_wfe_mcaoH(startDelay:end,2,k)));
end
result.atm.nLayer = atm.nLayer;
result.atm.fractionalR0 = fractionalR0;
result.atm.altitude = altitude;
result.atm.windSpeed = windSpeed;
result.atm.windDirection = windDirection;
result.atm.r0 = r0;
result.atm.L0 = L0;

result.tel.D = D;      % telescope primary mirror diameter
result.tel.nRes = nRes;  % resolution on the pupil plane (no of pixels)
result.tel.samplingFreq = samplingFreq; % WFS sampling time
result.tel.obstructionRatio = obstructionRatio;% central obscuration ratio
result.tel.fieldOfViewInArcsec = fieldOfViewInArcsec; %fieldOfViewInArcsec
result.tel.zenithAngle = zenithAngle;

result.wfs.nL = nL;
result.wfs.nPx =nPx;       % pixels/subap in phase space
result.wfs.d = d;

result.ast.nLGS = nGs;
result.ast.nAst = nAst;
result.ast.theta = thetaAst;
result.ast.astRadius = radiusAst;
result.ast.flux = flux;
result.ast.lgsHeight = lgsHeight;

result.opt.zenith    = zenithOpt;
result.opt.azimuth = azimuthOpt;

result.sci.zenith    = zenithSci;
result.sci.azimuth = azimuthSci;

result.dmStack.nDM = nDm;
result.dmStack.pitchs = pitchs;
result.dmStack.zLocations = dmStack.zLocations;
result.dmStack.nActuators = nActuators;
result.dmStack.Surface = surfDM;

result.camR.nyquistSampling = camR.imgLens.nyquistSampling;
result.camR.fieldStopSize       = camR.imgLens.fieldStopSize;
result.camR.resolution           = camR.resolution;
result.camR.wavelength         =[sciR(1).wavelength];

result.camH.nyquistSampling = camH.imgLens.nyquistSampling;
result.camH.fieldStopSize       = camH.imgLens.fieldStopSize;
result.camH.resolution           = camH.resolution;
result.camH.wavelength         =[sciH(1).wavelength];

result.loop.gain_pol = gain_pol;
result.loop.exposureTime = exposureTime ;
result.loop.startDelay = startDelay;

result.performanceR.wfe = wfeR;
result.performanceR.srwfe = srwfeR;
result.performanceR.var_wfe_mcao = var_wfe_mcaoR;
result.performanceR.wfe_mcao = wfe_mcaoR;
result.performanceR.marechalStrehl_mcao = marechalStrehl_mcaoR;
result.performanceR.psfStrehl_mcao = psfStrehl_mcaoR;
result.performanceR.psf = camR.frame;
result.performanceR.refpsf = camR.referenceFrame;
result.performanceR.inWF = inWFR;
result.performanceR.resWF = resWFR;

result.performanceH.wfe = wfeH;
result.performanceH.srwfe = srwfeH;
result.performanceH.var_wfe_mcao = var_wfe_mcaoH;
result.performanceH.wfe_mcao = wfe_mcaoH;
result.performanceH.marechalStrehl_mcao = marechalStrehl_mcaoH;
result.performanceH.psfStrehl_mcao = psfStrehl_mcaoH;
result.performanceH.psf = camH.frame;
result.performanceH.refpsf = camH.referenceFrame;
result.performanceH.inWF = inWFH;
result.performanceH.resWF = resWFH;

save('/result/ybrule/MAVIS/mcaoBaselineFov30ThierryCam300pix10s.mat','result','-v7.3')