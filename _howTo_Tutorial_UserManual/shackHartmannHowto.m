%% SHACKHARTMANN HOWTO 
% Demonstrate the use of the <matlab:doc('shackHartmann') shackHartmann> class

%{
The detector pixel size can be adjusted considering the following
- The FoV of a sub-aperture is at most nPixElecField*loD, where
    - nPixElecField = number of pixels per sub-aperture used to sample the
    electric field. This is tel.resolution/wfs.nLenslet
    - loD is the diffraction limit of the sub-aperture
- the wfs.lenslets.fieldStopSize the lenslet field of view given in diffraction fwhm units
- wfs.lenslets.nyquistSampling gives the 2x number of pixels per loD; ontrols the angular size of the electric-field pixel
- wfs.camera.resolution bins the electric field pixels to a smaller number; controls the binning of the detector camera pixels

Example 1: Detector with 8 50mas pixels in K band
    tel.D = 11.25;
    ngs.wavelength = photometry.Ks;
    wfs.lenslets.nLenslet = 1;
    lo2DInMas =
    ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas; %40mas 
    wfs.lenslets.nyquistSampling = 2;% i.e. 4 pixels per loD -> 10mas/pixel
    wfs.lenslets.fieldStopSize = 10; % i.e there are 10*4=40 pixels in the output wave
    wfs.camera.resolution = [1 1]*8; % i.e. from 40 to 8 a factor 5 binning
    is obtained from 10mas pixels -> 50 mas pixels as desired
%}
%%
ccc

ngs = source;
nL = 1;
res = 24;
tel = telescope(8,'resolution',24*nL);

%% A Nyquist-sampled SH :: case 1: pixelSize = lambda/D/2
% instantiation
wfs = shackHartmann(nL,24*nL,0.5);

% 
wfs.lenslets.nyquistSampling = 1;
wfs.lenslets.fieldStopSize = 24;
%wfs.camera.resolution = [48,48]*2;
fprintf('Max fieldStopSize is : %2.0f \n', 2*tel.resolution/(wfs.lenslets.nyquistSampling*2))

2*wfs.lenslets.nyquistSampling * wfs.lenslets.fieldStopSize
% initialisation
ngs = ngs.*tel*wfs;
wfs.INIT


wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

ngs = ngs.*tel*wfs;
figure(212)
imagesc(wfs.camera)

%%
% gain calibration :: default lambda/2/d units
wfs.gainCalibration(tel,ngs); 

binFactor = max(1,2*wfs.lenslets.fieldStopSize/wfs.lenslets.nyquistSampling/(wfs.camera.resolution(1)/wfs.lenslets.nLenslet));
lo2DInMas = ngs.wavelength/(2*tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
loDInMas = ngs.wavelength/(tel.D/wfs.lenslets.nLenslet)*constants.radian2mas;
pixelSizeInMas = lo2DInMas*binFactor;

% diplay lensletArray (electric field pixels)
fprintf('Pixel size is %4.0f mas\n' , pixelSizeInMas)
fprintf('Pixel size is %4.2f lambda/D \n' , pixelSizeInMas/loDInMas)
wfs.lenslets
%alternative
wfs.lenslets.pixelScale(ngs,tel)


%% a random aberration
zern = zernike(5:6,'resolution',tel.resolution,'pupil',tel.pupil);
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
% done in two ways
%   1) set wfs.lenslets.nyquistSampling = 0.5; if done alone, limits the
%   FoV to half the default value, generating a detector image with half
%   the number of pixels
%       1.1) enlarge the FoV to have the same output detector pixels: wfs.lenslets.fieldStopSize = 6
%   2) set wfs.lenslets.fieldStopSize = 6; Since the pixel size is 1/2 loD,
%   the total FoV=12loD. As the detector resolution is still 6, the
%   detector will bin the pixels 2x, creating the desired effect
%
%   NOTE HOWEVER, that these two options are not equivalent. In 1) the FFT
%   sampling is changed to a smaller number (i.e. the sky pixel size is actually smaller), in 2) the FFT is done at
%   Nyquist sampling and then binned down in the detector. This leads to
%   numerical differences as the two operations are not interchangeable
nr = 24;
wfs = shackHartmann(nL,nr*nL,0.5);
option = 2;
% fieldStopSize is in units of the DL (=lambda/D). By choosing 6, 12 pixels 
% are used in the Fraunhofer propagation (because
% wfs.lenslets.nyquistSampling=1 by default). Since the detector only has 6
% pixels per subaperture, the pixelSize=lambda/D
 
switch option
    case 1
        wfs.lenslets.fieldStopSize = nr;
        wfs.lenslets.nyquistSampling = 0.5;
        
    case 2
        wfs.lenslets.fieldStopSize = nr;
        wfs.camera.resolution = [nr nr];
end


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


% diplay lensletArray (electric field pixels)
fprintf('Subap FoV is %4.0f mas\n' , lo2DInMas*wfs.lenslets.fieldStopSize/wfs.lenslets.nyquistSampling)
fprintf('Pixel size is %4.0f mas\n' , pixelSizeInMas)
fprintf('Pixel size is %4.2f lambda/D \n' , pixelSizeInMas/loDInMas)
%alternative
wfs.lenslets.pixelScale(ngs,tel)
%% case of a subsampled SH :: case 3: pixelSize = 2*lambda/D
tel = telescope(8,'resolution',120);
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
