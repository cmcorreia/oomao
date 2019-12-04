
%% SIMULATION PARAMTERS
% number of iteration
rParam.simu.nIteration           = 30000;
% sampling time seconds
rParam.simu.samplingTime         = 0.004;
% frame Time
rParam.simu.frameTime            = 0.004;
% sets the global random state
rParam.simu.seed                 = 4987; 

% local variants
nIteration                       = rParam.simu.nIteration;
samplingTime                     = rParam.simu.samplingTime;
frameTime                        = rParam.simu.frameTime;
seed                             = rParam.simu.seed;
rng(seed);


%% ATMOSPHERIC PARAMTERS
% zenith angle
if ~exist('zenithAngle','var')
    zenithAngle = 0; % zenith angle in radians
end
rParam.atm.zenithAngle = zenithAngle;

% 7 Layer profile
rParam.atm.nLayersProfile       = 3;%7;
% fractional R0
rParam.atm.fractionalR0         = [.5 .25 .25];%[59.60 9.63 3.25 3.72 8.69 6.84 8.26]/100;
% height in meters
rParam.atm.altitude             = [0 5000 10000];%[0 0.5 1 2 4 8 16]*1000;
% windSpeeds in m/s
rParam.atm.windSpeed            = [5 10 20];%[7 8 9 11.6 19.7 33.0 7.0];
%windDirection in rad
rParam.atm.windDirection        = [0 90 180];%[0 45 90 135 180 225 270]*pi/180;
% coherence lenght in meters at 0.5microns
rParam.atm.r0                   = 0.156*cos(zenithAngle)^(3/5);       
% Outer scale in meters
rParam.atm.L0                   = 30;   
% Apparent altitue in meters
rParam.atm.altitude             = rParam.atm.altitude./cos(zenithAngle);
% Wavelength
rParam.atm.wavelength           = photometry.V0;

% local variants
nLayer                          = rParam.atm.nLayersProfile;
r0                              = rParam.atm.r0;
fractionalR0                    = rParam.atm.fractionalR0;
L0                              = rParam.atm.L0;
altitude                        = rParam.atm.altitude;
windSpeed                       = rParam.atm.windSpeed;
windDirection                   = rParam.atm.windDirection;
atmwavelength                   = rParam.atm.wavelength;


%% TELESCOPE
% telescope primary mirror diameter
rParam.tel.D                    = 8;
% central obscuration ratio
rParam.tel.obstructionRatio     = 0.2;
% fieldOfViewInArcsec
rParam.tel.fieldOfViewInArcsec  = 120;
% resolution on the pupil plane (no of pixels)
rParam.tel.nRes                 = 240;

% local variants
D                               = rParam.tel.D;
nRes                            = rParam.tel.nRes;
obstructionRatio                = rParam.tel.obstructionRatio;
fieldOfViewInArcsec             = rParam.tel.fieldOfViewInArcsec; 


%% Natural Guide Star
% guide star wavelength
rParam.ngs.wavelength           = photometry.R;
% guide star magnitude
rParam.ngs.magnitude           = 0;
% zenith angle of guide stars in radian
rParam.ngs.zenith               = 45*constants.arcsec2radian;
% azimuth angle of guide stars in radian
rParam.ngs.azimuth              = 0*pi/180;
% number of guide star
rParam.ngs.nNGS                 = 3;
% ngs asterism
rParam.ngs.asterism             = {[rParam.ngs.nNGS,rParam.ngs.zenith,rParam.ngs.azimuth]};


%local variants
ngsWavelength                   = rParam.ngs.wavelength;
ngsMagnitude                    = rParam.ngs.magnitude;
ngsZenith                       = rParam.ngs.zenith;
ngsAzimuth                      = rParam.ngs.azimuth;
nNGS                            = rParam.ngs.nNGS;
ngsAsterism                     = rParam.ngs.asterism;


%% NGS WAVE-FRONT SENSOR
% FoV = wavelength/dsub*fieldSpotSize
% pixel scale = wavelength/dsub*fieldSpotSize/nPx 
%
% number of lenslets
rParam.ngswfs.nL                = 16;         
% number of pixels per lenslet
rParam.ngswfs.nPx               = 30;
% field spot size
rParam.ngswfs.fieldSpotSize     = 15;
% nyquistSampling
rParam.ngswfs.nyquistSampling   = 4;
% minumum intensity ratio to a fully-illuminated sub-aperture to be considered in
rParam.ngswfs.minLightRatio     = 0.85;
% subaperture diameter
rParam.ngswfs.dsub              = D/rParam.ngswfs.nL;
% wfs read out noise
rParam.ngswfs.readOutNoise      = .2;
% ngs FoV
rParam.ngswfs.FoV               = ngsWavelength.wavelength/rParam.ngswfs.dsub*rParam.ngswfs.fieldSpotSize;
% ngs pixel scale
rParam.ngswfs.pixelScale        = ngsWavelength.wavelength/rParam.ngswfs.dsub*rParam.ngswfs.fieldSpotSize/rParam.ngswfs.nPx;
% ngs CoG Threshold
rParam.ngswfs.Threshold         = [0 0.1];

% local variants
ngsnL                           = rParam.ngswfs.nL;              
ngsnPx                          = rParam.ngswfs.nPx;   
ngsMinLightRatio                = rParam.ngswfs.minLightRatio;
ngsDsub                         = rParam.ngswfs.dsub;
ngsReadOutNoise                 = rParam.ngswfs.readOutNoise;
ngsFieldSpotSize                = rParam.ngswfs.fieldSpotSize;
ngsNyquistSampling              = rParam.ngswfs.nyquistSampling;
ngsPixelScale                   = rParam.ngswfs.pixelScale;
ngsFoV                          = rParam.ngswfs.FoV;
ngsThreshold                    = rParam.ngswfs.Threshold;


if nRes < ngsFieldSpotSize*ngsnL || rem(nRes,ngsnL)~=0
    error('Telescope Resolution Error');
end