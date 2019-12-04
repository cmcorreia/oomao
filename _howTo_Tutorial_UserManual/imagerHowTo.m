%% IMAGER HOWTO
% Demonstrate the use of the <matlab:doc('imager') detector> class

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
camH = imager('diameter',tel.D, 'fieldStopSize',30,'nyquistSampling',8*2);
tel = tel - atm;
scienceH = scienceH.*tel*camH;
figure(31416)
imagesc(camH,'parent',subplot(2,1,1))

camH.referenceFrame = camH.frame;

%% ENSQUARED ENERGY
nPts = 150;
pixSizeInArcsec = scienceH.wavelength/tel.D/camH.imgLens.nyquistSampling/2*constants.radian2arcsec;
camH.eeWidth = linspace(0,60*pixSizeInArcsec,nPts);% integration box size in arcsec
camH.flush
camH.ee(1) = 0;
xx = linspace(0,60*pixSizeInArcsec,nPts);

figure
difLimInArcSec = scienceH.wavelength/tel.D*constants.radian2arcsec;
plot(xx/difLimInArcSec, squeeze(camH.ee))
title('Ensquared energy, Airy function')
xlabel('integration box diameter size, [\lambda/D units]')
set(gca,'FontSize',13)
grid on
hold on
%% NUMERICAL INTEGRATION
totalSumRefImage = sum(camH.referenceFrame(:));
for d = 0:nPts-2
    roi = camH.referenceFrame(end/2-d:end/2+1+d,end/2-d:end/2+1+d);
    ee(d+2) = sum(roi(:))/totalSumRefImage;
end
plot(xx/difLimInArcSec*2, ee) % the factor of 2 here is due to the integration box being n x "2" x pixSizeInArcsec