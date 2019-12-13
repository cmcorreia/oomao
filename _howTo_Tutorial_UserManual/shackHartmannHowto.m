%% SHACKHARTMANN HOWTO 
% Demonstrate the use of the <matlab:doc('shackHartmann') shackHartmann> class

%{
The detector pixel size can be adjusted considering the following
- The FoV of a sub-aperture is at most 2*nPixElecField*lo2D, where
    - nPixElecField = number of pixels per sub-aperture used to sample the
    electric field. This is tel/resolution/wfs.nLenslet
    - lo2D is half the diffraction limit of the sub-aperture
- the wfs.lenslets.fieldStopSize the lenslet field of view given in diffraction fwhm units
- wfs.lenslets.nyquistSampling gives the 2x number of pixels per lo2D
-wfs.camera.resolution bins the electric field pixels to a smaller number

Example 1: Detector with 8 50mas pixels in K band
    tel.D = 11.25;
    ngs.wavelength = photometry.Ks;
    lo2DInMas =
    ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas; %40mas
    wfs.lenslets.nyquistSampling = 2;% i.e. 4 pixels per loD -> 10mas/pixel
    wfs.lenslets.fieldStopSize = 10; % i.e there are 10*4=40 pixels in the output wave
    wfs.camera.resolution = [1 1]*8; % i.e. from 40 to 8 a factor 5 binning
    is obtained from 10mas pixels -> 50 mas pixels as desired

- wfs.lenslets.nyquistSampling: controls the angular size of the
electric-field pixel
- wfs.camera.resolution: controls the binning of the detector camera pixels
- wfs.lenslets.fieldStopSize: controls the field bounds passed to the
detector

            lo2DInMas = ngs.wavelength/(2*d)*constants.radian2mas;
            binFactor = 2*obj.lenslets.fieldStopSize*obj.lenslets.nyquistSampling/nPxDetector;
            detectorPixelSizeInMas = lo2DInMas/obj.lenslets.nyquistSampling*binFactor;



%}
%%
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
