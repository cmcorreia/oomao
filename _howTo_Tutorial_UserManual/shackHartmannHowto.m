%% SHACKHARTMANN HOWTO 
% Demonstrate the use of the <matlab:doc('shackHartmann') shackHartmann> class

ngs = source;
tel = telescope(8,'resolution',120);

%% A Nyquist-sampled SH :: case 1: pixelSize = lambda/D/2
% instantiation
wfs = shackHartmann(20,120,0.5);
% initialisation
ngs = ngs.*tel*wfs;
wfs.INIT


wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

ngs = ngs.*tel*wfs;
figure(212)
imagesc(wfs.camera)

% gain calibration :: default lambda/2/d units
wfs.gainCalibration(tel,ngs); 

binFactor = max(1,2*wfs.lenslets.fieldStopSize/wfs.lenslets.nyquistSampling/(wfs.camera.resolution(1)/wfs.lenslets.nLenslet));
lo2DInMas = ngs.wavelength/(2*tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
loDInMas = ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
pixelSizeInMas = lo2DInMas*binFactor;
fprintf('Pixel size is %4.0f mas\n' , pixelSizeInMas)
fprintf('Pixel size is %4.2f lambda/D \n' , pixelSizeInMas/loDInMas)

%alternative
wfs.lenslets.pixelScale(ngs,tel)
%% a random aberration
zern = zernike(5:6,'resolution',120,'pupil',tel.pupil);
zern.c = 20*(2*rand(zern.nMode,1)-1);
% reset(ngs)*tel*zern*wfs.lenslets;
% grabAndProcess(wfs)
ngs = ngs.*tel*zern*wfs;

%% a random aberration function
for k=1:100
    o = (k-1).*2*pi/99;
    zern.c = 10.*[cos(o);sin(o)];
    +ngs;
    drawnow
end


%% case of a subsampled SH :: case 2: pixelSize = lambda/D
wfs = shackHartmann(20,120,0.5);

% fieldStopSize is in units of the DL (=lambda/D). By choosing 6, 12 pixels 
% are used in the Fraunhofer propagation (because
% wfs.lenslets.nyquistSampling=1 by default). Since the detector only has 6
% pixels per subaperture, the pixelSize=lambda/D
 
wfs.lenslets.fieldStopSize = 6;
% initialisation
ngs = ngs.*tel*wfs;
wfs.INIT


wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

ngs = ngs.*tel*wfs;
figure(212)
imagesc(wfs.camera)

binFactor = max(1,2*wfs.lenslets.fieldStopSize/wfs.lenslets.nyquistSampling/(wfs.camera.resolution(1)/wfs.lenslets.nLenslet));
lo2DInMas = ngs.wavelength/(2*tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
loDInMas = ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
pixelSizeInMas = lo2DInMas*binFactor;
fprintf('Pixel size is %4.0f mas\n' , pixelSizeInMas)
fprintf('Pixel size is %4.2f lambda/D \n' , pixelSizeInMas/loDInMas)


%% case of a subsampled SH :: case 3: pixelSize = 2*lambda/D
wfs = shackHartmann(20,120,0.5);

% fieldStopSize is in units of the DL (=lambda/D). By choosing 6, 12 pixels 
% are used in the Fraunhofer propagation (because
% wfs.lenslets.nyquistSampling=1 by default). Since the detector only has 3
% pixels per subaperture, the pixelSize=2*lambda/D
wfs.lenslets.fieldStopSize = 6;

% camera resolution can be set directly in the wfs =
% shackHartmann(20,60,0.5); It controls the binning of the pixels
wfs.camera.resolution = 60;
% initialisation
ngs = ngs.*tel*wfs;
wfs.INIT


wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

ngs = ngs.*tel*wfs;
figure(212)
imagesc(wfs.camera)

binFactor = max(1,2*wfs.lenslets.fieldStopSize/wfs.lenslets.nyquistSampling/(wfs.camera.resolution(1)/wfs.lenslets.nLenslet));
lo2DInMas = ngs.wavelength/(2*tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
loDInMas = ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
pixelSizeInMas = lo2DInMas*binFactor;
fprintf('Pixel size is %4.0f mas\n' , pixelSizeInMas)
fprintf('Pixel size is %4.2f lambda/D \n' , pixelSizeInMas/loDInMas)


%% case of a subsampled SH :: case 4: quadCell 
wfs = shackHartmann(20,40,0.5);

% fieldStopSize is in units of the DL (=lambda/D). By choosing 6, 12 pixels 
% are used in the Fraunhofer propagation (because
% wfs.lenslets.nyquistSampling=1 by default). Since the detector only has 6
% pixels per subaperture, the pixelSize=lambda/D
wfs.lenslets.fieldStopSize = 6; % user can change this from 1 to 6 and check results

% camera resolution can be set directly in the wfs =
% shackHartmann(20,60,0.5); It controls the binning of the pixels
% wfs.camera.resolution = 40;

% initialisation
ngs = ngs.*tel*wfs;
wfs.INIT


wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

ngs = ngs.*tel*wfs;
figure(212)
imagesc(wfs.camera)

binFactor = max(1,2*wfs.lenslets.fieldStopSize/wfs.lenslets.nyquistSampling/(wfs.camera.resolution(1)/wfs.lenslets.nLenslet));
lo2DInMas = ngs.wavelength/(2*tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
loDInMas = ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
pixelSizeInMas = lo2DInMas*binFactor;
fprintf('Pixel size is %4.0f mas\n' , pixelSizeInMas)
fprintf('Pixel size is %4.2f lambda/D \n' , pixelSizeInMas/loDInMas)
