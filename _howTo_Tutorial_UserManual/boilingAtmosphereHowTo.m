%This script provides an implementation of the boiling atmosphere model in
%Srinath et al, An Efficient Autoregressive Method for Generating Phase
%Screens with Frozen Flow and Turbulence in Optical Simulations, Opt Exp
%2015
% The spatial and temporal statistics are evaluated and compared to the
% default OOMAO implementation of frozen-flow phase-screens
%
% Spatial statistics are quite close to the theory with frozen atmosphere
% converging much faster than boiling atmosphere. The frozen atmosphere
% follows quite well the theoretical temporal PSDs from Conan95 whereas the
% boiling deviates from those. Although this is expected, the deviations
% are rather difficult to explain. TBD develop an analytic modelfor the
% temporal PSD of the boiling atmosphere...
%% INIT
nPixel  = 128;  % output phaseScreen size
D       = 4;    % physical dimensions of the phaseScrenn in m
T       = 1/100;% samplingTime in s
vx      = 10;  % windSpeed in x    
vy      = 10;  % windSpeed in y
alpha   = 0.999;% auto-regressive "forgetting" factor
u       = 1:nPixel;
%% ATMOSPHERE
% 3-layer
% atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',[0,4,10]*1e3,...
%     'fractionnalR0',[0.7,0.25,0.05],...
%     'windSpeed',[5,10,20],...
%     'windDirection',[0,pi/4,pi]);

%1-layer
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0]*1e3,...
    'fractionnalR0',[1],...
    'windSpeed',[10*sqrt(2)],...
    'windDirection',[pi/4]);

% telescope
tel = telescope(D,'resolution',nPixel,'samplingTime',T);

%% start atmosphere (frozen-flow)
tel = tel + atm;

%% start atmoshere (boiling)
nTimes = 4;
N = nTimes*nPixel;
L = (N-1)*D/(nPixel-1);

[fx,fy]  = freqspace(N,'meshgrid');
[~,fr]  = cart2pol(fx,fy);
fr  = fftshift(fr.*(N-1)/L./2);

[X,Y] = ndgrid(-N/2:N/2-1);
X = X/(nPixel*D);
Y = X';

spatialShiftWindPhasor = fftshift(exp(-1i*2*pi*T*(X*vx + Y*vy)));


% sqrt(psd) of full phase screen
spectrumVK  = phaseStats.spectrum(fr,atm);
psdRoot     = sqrt(spectrumVK); % Phase FT magnitude
[idx]       = find(fr==0);
psdRoot(idx)= 0;
fourierSampling = 1./L;
WNF         = fft2(randn(atm.rngStream,N))./N; % White noise filtering
map         = psdRoot.*WNF;
map         = real(ifft2(map).*fourierSampling).*N.^2;

imagesc(map(u,u))

%% ITERATE (JUST FOR VISUALISATION PURPOSES, OTHERWISE SKIP TO NEXT CELL)
% spectrum_old is the spectrum from the previous timestep, alpha is the boiling coefficient
% WNF         = fft2(randn(atm.rngStream,N))./N; % White noise filtering
% spectrum    = psdRoot.*WNF;
% alpha       = 1;% auto-regressive "forgetting" factor
% CLIM        = [-50 50];
% for kIter = 1:100
%     WNF         = fft2(randn(atm.rngStream,N))./N; % White noise filtering
%     spectrum    = alpha*spatialShiftWindPhasor.*spectrum + sqrt(1-alpha^2)*psdRoot.*WNF;
%     map         = real(ifft2(spectrum).*fourierSampling).*N.^2;
%     imagesc(map, CLIM), colorbar, drawnow
% end

%% CHECK STATS
zern = zernike(2:100,D,'resolution',nPixel);
src = source;

spectrum    = psdRoot.*WNF;
alpha       = 0.9;% auto-regressive "forgetting" factor
CLIM        = [-50 50];
for kIter = 1:10000
    WNF         = fft2(randn(atm.rngStream,N))./N; % White noise filtering
    spectrum    = alpha*spatialShiftWindPhasor.*spectrum + sqrt(1-alpha^2)*psdRoot.*WNF;
    map         = real(ifft2(spectrum).*fourierSampling).*N.^2;
    %imagesc(map, CLIM), colorbar, drawnow
    ps = map(u,u);
    zern\ps(:);
    boiling(:,kIter) = zern.c;
    +tel;
    src = src.*tel;
    zern\src;
    frozen(:,kIter) = zern.c*2*pi/src.wavelength;
    kIter
end


  
% SPATIAL SPECTRUM
semilogy(var(boiling'))
hold on
semilogy(var(frozen'))

semilogy(zernikeStats.variance(zern, atm, src),'k','linewidth',2)
title('Cheking Zernike polynomials spatial variance')
ylabel('[rad^2 @ 550nm]')
xlabel('zernike mode number j')
set(gca,'fontSize',16)
legend('Boiling atmosphere','Frozen-flow atmosphere','Theory')
% COMMENT: Agreement is quite good withas few as 100 timeSteps. With 1000
% the agreement is quasi-perfect. The stats converge towards the
% theoretical ones
%% video

h = figure(10);
axis tight manual square% this ensures that getframe() returns a consistent size
filename = 'testAnimatedBoiling.gif';
alpha = 0.995;
mask = utilities.piston(128);
idx = find(mask);
mask(mask == 0) = nan;
for kIter = 1:150
WNF         = fft2(randn(atm.rngStream,N))./N; % White noise filtering
    spectrum    = alpha*spatialShiftWindPhasor.*spectrum + sqrt(1-alpha^2)*psdRoot.*WNF;
    map         = real(ifft2(spectrum).*fourierSampling).*N.^2;
    %imagesc(map, CLIM), colorbar, drawnow
    ps = map(u,u);
    acc(:,:,kIter) = ps;
    imagesc(mask.*(ps-mean(ps(idx))),[-20 20]); colorbar, drawnow
          % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if kIter == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

%% TEMPORAL SPECTRUM
figure(2)
clf
jMode = 9;
for jMode = 6
nu = logspace(-1,log10(50),200);
zern1 = zernike(jMode,D,'resolution',nPixel);
[out] = lamTools.multiZernikeTemporalSpectra(nu,atm,zern1,tel);

[outBoiling] = lamTools.ZernikeTemporalSpectrum(nu,atm,tel, zern1.j,alpha,T);

fprintf('Zernike polynomial variance from temporal spectrum:%f, from spatial variance: %f \n', ...
    trapz(nu,out), zernikeStats.variance(zern1, atm, src))
loglog(nu, out)
hold on
%loglog(nu, outBoiling)
%[Y,f] = myPSD(boiling(jMode-1,:),1/T,0);
%fprintf('Zernike polynomial variance from EMPIRICAL temporal spectrum:%f\n', trapz(f,Y))

[psdj,nuj] = pwelch(boiling(jMode-1,:),[],[],[],T);
psdj = psdj*T^2;
nuj = nuj*2/T*1/2/T;
loglog(nuj, psdj)

[psdj,nuj] = pwelch(frozen(jMode-1,:),[],[],[],T);
psdj = psdj*T^2;
nuj = nuj*2/T*1/2/T;
loglog(nuj, psdj)
xlabel('temporal frequency,[Hz]')
ylabel('[rad^2/Hz]')
title(['Temporal PSD mode ' num2str(jMode)])
legend('Theory (frozen flow)','pwelch boiling','pwelch frozen')
set(gca,'fontSize',16)
SP = 0.3*(zern.n(jMode)+1)*atm.layer(1).windSpeed/tel.D;
line([SP SP],get(gca,'YLim'),'Color',[1 0 0])
drawnow
end