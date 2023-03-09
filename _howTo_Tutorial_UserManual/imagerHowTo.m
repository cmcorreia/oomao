%% IMAGER HOWTO
% Demonstrate the use of the <matlab:doc('imager') detector> class
ccc
%% ATMOSPHERE
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.25,0.05],...
    'windSpeed',[5,10,20],...
    'windDirection',[0,pi/4,pi]);

%% TELESCOPE
D = 8;
nPx = 60;
tel = telescope(D,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/100);

%% SCIENCE CAMERA

% H-science camera
scienceH = source('wavelength',photometry.H);
camH = imager('diameter',tel.D, 'fieldStopSize',30,'nyquistSampling',2);
tel = tel - atm;
scienceH = scienceH.*tel*camH;
figure(31416)
imagesc(camH,'parent',subplot(2,1,1))

camH.referenceFrame = camH.frame;

%% EE in \lambda/D 
% Encircled (not ensquared) energy of an Airy pattern, from https://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/PSFtheory.pdf
% ? (?/D units) EE(?)
% 0.52    0.50
% 0.70    0.70
% 1.22    0.86
% 1.66    0.90

camH.eeWidth = camH.imgLens.nyquistSampling*2* [0.52 0.79 1.22 1.66]';
camH.flush
camH.ee(1,:)
%% ENSQUARED ENERGY
intBoxDiam = linspace(0,96,50); % in pixels
camH.eeWidth = intBoxDiam;
camH.flush

pixSizeInArcsec = scienceH.wavelength/tel.D/camH.imgLens.nyquistSampling/2*constants.radian2arcsec;
difLimInArcSec = scienceH.wavelength/tel.D*constants.radian2arcsec;

figure
plot(intBoxDiam*pixSizeInArcsec/difLimInArcSec, squeeze(camH.ee(1,:)),'displayName','numerical integration in OTF space')
title(['Ensquared energy, Airy function, sampling: \theta_{pix}=' num2str(camH.imgLens.nyquistSampling/2) '\lambda/D'])
xlabel('integration box diameter size, [\lambda/D units]')
set(gca,'FontSize',16)
grid on
hold on

pbaspect([1.618 1 1])
%% NUMERICAL INTEGRATION
nPts = size(camH.frame,1)/2;
totalSumRefImage = sum(camH.referenceFrame(:));
for d = 0:nPts-2
    roi = camH.referenceFrame(end/2-d:end/2+1+d,end/2-d:end/2+1+d);
    ee(d+2) = sum(roi(:))/totalSumRefImage;
end
plot((0:nPts-1)*pixSizeInArcsec/difLimInArcSec, ee,'--', 'displayName','numerical integration in PSF space') 
legend