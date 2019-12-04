function out = calibrateWfsPixelScale(wfs,gs,psInMas,nPx)


psInRad = psInMas/constants.radian2mas;
% Start with ngs on axis
zen0      = gs.zenith;
azi0      = gs.azimuth;
gs.zenith = 0;
gs.azimuth = 0;
% Set wfs on-axis
wfs.pointingDirection = zeros(2,1);

tipStep = psInRad/2;
nStep   = max(floor(nPx/3)*2,2);
sx      = zeros(1,nStep+1);
u       = 0:nStep;
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled       = false;

warning('off','oomao:shackHartmann:relay')
for kStep=u
    gs.zenith = -tipStep*kStep;
    +gs;
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end

warning('on','oomao:shackHartmann:relay')

% Convert input and output to arcsecs
Ox_in  = u*tipStep;
Ox_out = sx*psInRad;

% Calibrate units
slopesLinCoef = polyfit(Ox_in,Ox_out,1);

% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting pointing direction [].
gs.zenith  = zen0;
gs.azimuth = azi0;
+gs;
wfs.pointingDirection = [];
out = 1/slopesLinCoef(1) * wfs.slopesUnits;
end