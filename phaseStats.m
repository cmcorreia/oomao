classdef phaseStats
    % Phase statistics static class
    
    methods (Static)
        
        function  out = variance(atm)
            %% VARIANCE Phase variance
            %
            % out = phaseStats.variance(atm) computes the phase variance
            % from an atmosphere object
            %
            % See also atmosphere
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            out   = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            layers = atm.layer;
            out = sum([layers.fractionnalR0]).*out;
        end
        
        function out = covariance(rho,atm)
            %% COVARIANCE Phase covariance
            %
            % out = phaseStats.covariance(rho,atm) computes the phase covariance from
            % the baseline rho and an atmosphere object
            %
            % See also atmosphere
             
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            out   = ones(size(rho)).*(24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            index         = rho~=0;
            u             = 2.*pi.*rho(index)./atm.L0;
            out(index) = cst.*u.^(5./6).*besselk(5./6,u);
            layers = atm.layer;
            out = sum([layers.fractionnalR0]).*out;
        end
        
        function out = covarianceXDerivative(rho,argRho,atm)
            %% COVARIANCE Phase covariance
            %
            % out = phaseStats.covariance(rho,atm) computes the phase covariance from
            % the baseline rho and an atmosphere object
            %
            % See also atmosphere
             
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            out   = zeros(size(rho));
            index         = rho~=0;
            u             = 2.*pi.*rho(index)./atm.L0;
            out(index) = -cst.*(1/6).*u.^(-1./6).*(2*pi/atm.L0).*cos(argRho(index)).*...
                (-5*besselk(5./6,u) + 3.*u.*(besselk(-1./6,u)+besselk(11./6,u)));
            layers = atm.layer;
            out = sum([layers.fractionnalR0]).*out;
        end
        
        function out = angularCovariance(theta,atm)
            %% ANGULARCOVARIANCE Phase angular covariance
            %
            % out = phaseStats.angularCovariance(theta,atm) computes the
            % phase angular covariance from the zenith angle theta and an
            % atmosphere object
            %
            % See also atmosphere
           
            out = zeros(size(theta));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                out = out + phaseStats.covariance(atmSlab.layer.altitude*tan(theta),atmSlab);
            end
        end        
        function out = angularStructureFunction(theta,atm)
            %% ANGULARSTRUCTUREFUNCTION Phase angular structure function
            %
            % out = phaseStats.angularStructureFunction(rho,atm) computes the
            % phase angular structure function from the zenith angle theta
            % and an atmosphere object
            %
            % See also atmosphere
           
            out = zeros(size(theta));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                out = out + 2.*( phaseStats.variance(atmSlab) - ...
                    phaseStats.covariance(atmSlab.layer.altitude*tan(theta),atmSlab) );
            end
       end        
        
        function out = temporalCovariance(tau,atm)
            %% TEMPORALCOVARIANCE Phase temporal covariance
            %
            % out = phaseStats.temporalCovariance(rho,atm) computes the
            % phase temporal covariance from the delay tau and an
            % atmosphere object
            %
            % See also atmosphere
           
            out = zeros(size(tau));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                out = out + phaseStats.covariance(atmSlab.layer.windSpeed*tau,atmSlab);
            end
        end        
        function out = temporalStructureFunction(tau,atm)
            %% TEMPORALSTRUCTUREFUNCTION Phase temporal structure function
            %
            % out = phaseStats.temporalStructureFunction(rho,atm) computes
            % the phase temporal structure function from the delay tau and
            % an atmosphere object
            %
            % See also atmosphere
           
            out = zeros(size(tau));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                out = out + 2.*( phaseStats.variance(atmSlab) - ...
                    phaseStats.covariance(atmSlab.layer.windSpeed*tau,atmSlab) );
            end
        end
        
        function L2 = sparseInverseCovarianceMatrix(layerGrid,pitch,atm)
            %% SPARSEINVERSECOVARIANCEMATRIX 
            %
            % out = sparseInverseCovarianceMatrix(gridMask,pitch,atm)
            % computes a sparse approximation of the inverse of the
            % covariance matrix from the pupil grid mask and the atmosphere
            %
            % See also atmosphere
            
%             nGrid = length(layerGrid);
%             rho = [0 1 nGrid]*pitch;
%             cov = phaseStats.covariance(atm,rho);
%             cov = cov.*[-4 2 2];
%             sigma = sum(cov);
%             L2 = cell(1,nGrid);
%             for kGrid=1:nGrid
%                 gridMask = layerGrid{kGrid};
%                 [n,m] = size(gridMask);
%                 if any([n,m])==1
%                     n = sqrt(numel(gridMask));
%                 end
%                 L = gallery('poisson',n);
%                 L(~gridMask,:) = [];
%                 L(:,~gridMask) = [];
%                 L = (L'*L);
%                 L2{kGrid} = L./(phaseStats.variance(slab(atm,kGrid))*sum(L(:)));
%             end

            nGrid = length(layerGrid);
            L2 = cell(1,nGrid);
            for kGrid=1:nGrid
                gridMask  = layerGrid{kGrid};
                n         = sqrt(numel(gridMask));
                e         = ones(n^2,1);
%                 e(~gridMask) = [];
                L         = spdiags([e e -4*e e e],[-n -1 0 1 n],n^2,n^2);
                L(~gridMask,:) = [];
                L(:,~gridMask) = [];
                L2{kGrid} = L'*L;
                [i,j]     = find(L2{kGrid});
                [x,y]     = meshgrid( (0:n-1)*pitch );
                z         = complex(x(:),y(:));
                psCov     = phaseStats.covariance(abs(z(i)-z(j)),slab(atm,kGrid));
%                 psCov     = psCov*atm.la
                alpha2    = sum(nonzeros(L2{kGrid}).*psCov);
                L2{kGrid} = L2{kGrid}*n^2/alpha2;
            end

        end
        
        function out = structureFunction(rho,atm)
            %% STRUCTUREFUNCTION Phase structure function
            %
            % out = phaseStats.structureFunction(rho,atm) computes the
            % phase structure function from the baseline rho and an
            % atmosphere object
            %
            % See also atmosphere
            
            if isinf(atm.L0)
                out   = zeros(size(rho));
                index = rho~=0;
                out(index) = 2.*(24.*gamma(6./5)./5).^(5./6).*(rho(index)./atm.r0).^(5./3);
            else
                out = 2.*(phaseStats.variance(atm)-phaseStats.covariance(rho,atm));
            end
        end
        
        function out = spectrum(f,atm)
            %% SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            out = out.*(f.^2 + 1./atm.L0.^2).^(-11./6);
            layers = atm.layer;
            out = sum([layers.fractionnalR0]).*out;
        end
        function S = ampSpectrum(f,atm)
            %% SPECTRUM Amplitude power spectrum density
            %
            % out = phaseStats.ampSpectrum(f,atm) computes the amplitude power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            [~,Y] = phaseStats.fresnel(f, atm);
            S = phaseStats.spectrum(f,atm).*Y;
        end
        
        function [O, S] = fresnelChromaticity(f, atm, wvl1, wvl2)
            % OPD and AMP power-spectral densities caused by Fresnel
            % chromaticity
            [o,a] = phaseStats.chromaticity(f, atm, wvl1, wvl2);
            [X, Y] = phaseStats.fresnel(f,atm);
            O = o./X.*phaseStats.spectrum(f,atm);
            S = a./Y.*phaseStats.ampSpectrum(f,atm);
        end
        
        function [O,S] = fresnel(f, atm)
            %% Fresnel phase and amplitude factors when propagating through a stratified turbulence
            O = zeros(size(f));
            for kLayer = 1:atm.nLayer
                O = O + atm.layer(kLayer).fractionnalR0.*cos(pi*atm.altitudeAtZenith(kLayer)*f.^2*atm.wavelength).^2;
            end
            S = 1-O;
        end
        function [O,S] = chromaticity(f, atm, wvl1, wvl2)
            %% Fresnel phase and amplitude differential chromaticity when propagating through a stratified turbulence
            O = zeros(size(f));
            S = O;
            for kLayer = 1:atm.nLayer
                O = O + atm.layer(kLayer).fractionnalR0.*(cos(pi*atm.altitudeAtZenith(kLayer)*f.^2*wvl1) - ...
                    cos(pi*atm.altitudeAtZenith(kLayer)*f.^2*wvl2)).^2;
                S = S + atm.layer(kLayer).fractionnalR0.*(sin(pi*atm.altitudeAtZenith(kLayer)*f.^2*wvl1) - ...
                    sin(pi*atm.altitudeAtZenith(kLayer)*f.^2*wvl2)).^2;
            end
        end
        
        function [O, S] = dispersionDisplacement(f, atm, theta)
            %% Differential refraction of the atmosphere that induces a
            % beam-shift thus becomes akin to an anisoplanatic error
            
            src = source('zenith',theta);
            zLayer = [atm.layer.altitude];
            fr0    = [atm.layer.fractionnalR0];
            A = zeros(size(f));
            for kLayer=1:atm.nLayer
                red = 2*pi*zLayer(kLayer)*...
                    ( f*src.directionVector(1) );
                A  = A + fr0(kLayer)*exp(1i*red);
            end
            O = (2-A-conj(A)).*phaseStats.spectrum(f, atm);
            S = (2-A-conj(A)).*phaseStats.ampSpectrum(f, atm);
        end
        
        function [O, S] = differentialRefraction(f, atm, wvl1, wvl2)
            %% OPD PSD due to differential atmospheric refraction
            % f - frequency vector
            % atm - atm object
            % wvl1 and wvl2 in meters
            wvl1InMicrons = wvl1*1e6;
            wvl2InMicrons = wvl2*1e6;
            n1 = phaseStats.nstp(wvl1InMicrons);
            n2 = phaseStats.nstp(wvl2InMicrons);
            O = phaseStats.spectrum(f, atm).*...
                ((n1 - n2)/n2)^2;
            S = phaseStats.ampSpectrum(f, atm).*...
                ((n1 - n2)/n2)^2;
            % The TFusco06 implementation in Eq 12 & 13 give the same result
            % as above
            %            nAverage = (n(wvl1InMicrons) + n(wvl2InMicrons))/2;
            %             DeltaWvl1 = abs(n(wvl1InMicrons) - nAverage);
            %             DeltaWvl2 = abs(n(wvl2InMicrons) - nAverage);
            %             alpha = DeltaWvl2/DeltaWvl1;
            %            alpha = n(wvl1InMicrons)/n(wvl2InMicrons);
            %            out = ((alpha-1)/alpha)^2*phaseStats.spectrum(f, atm);
            
        end
        %% correction chromatism = differential refraction
        function [O, S] = correctionChromatism(f, atm, wvl1, wvl2)
            [O, S] = phaseStats.differentialRefraction(f, atm, wvl1, wvl2);
        end
        
        %% refractive index fluctuations -  n at standard temperature and pressure
        function n = nstp(wvl, T, P)
            % Edlen 1966, Hardy Eq. 3.16, Fusco06 Eq. 8, Guyon 05 Eq. 28
            %n     = @(x) ...%1.0 ...
            %    + 8.34213e-5 + 0.0240603/(130-x^(-2))+ 0.00015997/(38.9-x^(-2)); % Hardy Eq. (3.16). Inex of refraction (and not refractivity) variations at standard temperature and pressure
            %
            if nargin == 1
                T = 12; %�C
                P = 1000; % bar
            end
            nstp = @(wvl, t, p) 1e-8*(8342.13+2406030./(130.0-wvl.^(-2))+15997./(38.9-wvl.^(-2))).*0.0010413*p/(1+0.003671*t); % Hardy Eq. (3.16) for standard Temperature (�C) and pressure (mbar)
            n = nstp(wvl,T,P);
        end
        %%
        function out = AoASpectrum(f,o,atm,tel)
              %% AOASPECTRUM Angle-of-Arrival temporal power spectrum density for a round aperture of diameter D (tel.D)
             fx = f.*cos(o);
             out = phaseStats.spectrum( f , atm ) .* abs(2*1i*pi*fx.*besselj(1,pi*tel.D*f)./(pi*tel.D*f)).^2;
         end
           function out = temporalAoASpectrum(nu,atm,tel,unit)
            %% TEMPORALSPECTRUM Phase temporal power spectrum density
            %
            % out = phaseStats.temporalAoASpectrum(nu,atm) computes the
            % angle-of-arrival single-axis temporal power spectrum density in rad^2 from the temporal
            % frequency nu and an atmosphere object 
            %
            % See also atmosphere and zernike
            %
            % EXAMPLE: Comparison of AoA to Z2. From numerical calculations, 
            % the AoA variance (for a squared aperture) is 0.162 * (d/r0)^(5/3) 
            % whereas the Z2 tilt variance is 0.172 * (d/r0)^(5/3) 
            
            %             tel = telescope(10,...
            %                 'fieldOfViewInArcMin',0,...
            %                 'resolution',60,...
            %                 'samplingTime',1/250);
            %
            %             atm = atmosphere(photometry.V0,0.15,30,...
            %                 'altitude',4e3,...
            %                 'fractionnalR0',1,...
            %                 'windSpeed',15,...
            %                 'windDirection',0);
            %
            %             nu = logspace(-3,log10(1000),300);
            %
            %             AoA = phaseStats.temporalAoASpectrum(nu,atm,tel); % AoA in rad
            %             fprintf('mas rms on-sky for a %d m telescope: %f mas rms\n',tel.D,sqrt(trapz(nu,AoA))*180/pi*3600*1000*atm.wavelength/2/pi)
            %
            %             zern = zernike(2:2,'resolution',tel.resolution,'pupil',tel.pupil,'D',tel.D);
            %             [out, outLayered] = zernikeStats.multiZernikeTemporalSpectra(nu,atm,zern,tel);
            %             rad2mas = 4/(tel.D)*180/pi*3600*1000*atm.wavelength/2/pi;
            %             fprintf('mas rms on-sky for a %d m telescope: %f mas rms\n',tel.D,sqrt(trapz(nu,out))*rad2mas)
            
            out = zeros(size(nu));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                [vx,vy] = pol2cart(atmSlab.layer.windDirection,atmSlab.layer.windSpeed);
                for k=1:numel(nu)
                    if abs(vx)>10*eps(atmSlab.layer.windSpeed)
                        out(k) = out(k) + quadgk( @integrandFy , -Inf, Inf);
                    else
                        out(k) = out(k) + quadgk( @integrandFx , -Inf, Inf);
                    end
                end
            end
            out = out*4; % factor here accounts for double-sided spectrum for a single-axis
            if nargin > 3
                switch unit
                    case 'rad'
                    case 'am' %arcmin^2
                        out = out*(180/pi*60*atm.wavelength/2/pi)^2;
                    case 'as' %arcsec^2
                        out = out*(180/pi*3600*atm.wavelength/2/pi)^2;
                    case 'mas' %mas^2
                        out = out*(180/pi*3600*1000*atm.wavelength/2/pi)^2;
                    otherwise
                end
            end
            function int = integrandFy(fy)
                fx = (nu(k) -fy*vy)/vx;
                int = 2*abs(phaseStats.AoASpectrum( hypot(fx,fy) ,atan2(fy,fx), atmSlab,tel )/vx);
            end
            
            function int = integrandFx(fx)
                fy = (nu(k) -fx*vx)/vy;
                int = 2*abs(phaseStats.AoASpectrum( hypot(fx,fy) ,atan2(fy,fx), atmSlab,tel )/vy);
            end
           end   
        
        function out = temporalSpectrum(nu,atm)
            %% TEMPORALSPECTRUM Phase temporal power spectrum density
            %
            % out = phaseStats.temporalSpectrum(nu,atm) computes the
            % phase temporal power spectrum density from the temporal
            % frequency nu, an atmosphere object and a zernike object
            %
            % See also atmosphere and zernike
            
            out = zeros(size(nu));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                [vx,vy] = pol2cart(atmSlab.layer.windDirection,atmSlab.layer.windSpeed);
                for k=1:numel(nu)
                    if vx>eps(atmSlab.layer.windSpeed)
                        out(k) = out(k) + quadgk( @integrandFy , -Inf, Inf);
                    else
                        out(k) = out(k) + quadgk( @integrandFx , -Inf, Inf);
                    end
                end
            end
            
            function int = integrandFy(fy)
                fx = (nu(k) -fy*vy)/vx;
                int = phaseStats.spectrum( hypot(fx,fy) , atmSlab )/vx;
            end
            
            function int = integrandFx(fx)
                fy = (nu(k) -fx*vx)/vy;
                int = phaseStats.spectrum( hypot(fx,fy) , atmSlab )/vy;
            end
        end    
        
        function out = closedLoopVariance(atm,T,tau,gain)
            %% SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            s = @(x) 2*1i*pi*x;
            z = @(x) exp(s(x)*T);
            
            G = @(x) ((1-exp(-s(x)*T))./(s(x)*T)).^2.*...
                exp(-tau*s(x)).*...
                gain./(1-exp(-s(x)*T));
            E = @(x) abs(1./(1+G(x)));
            
            figure
            nu = logspace(-2,log10(2/T),101);
            subplot(1,2,1)
            loglog(nu,abs(E(nu)).^2)
            xlabel('Hz')
            subplot(1,2,2)
            loglog(nu,phaseStats.temporalSpectrum(nu,atm),...
                nu,phaseStats.temporalSpectrum(nu,atm).*abs(E(nu)).^2)
            xlabel('Hz')
            drawnow
            
            out = 2*quadgk( @(nu) phaseStats.temporalSpectrum(nu,atm).*abs(E(nu)).^2 , 0 , Inf);

        end        
        
        function out = symSpectrum(symf)
            syms r0 L0
            out = (24*gamma(sym(6/5))/5)^(5./6)*...
                (gamma(sym(11/6))^2/(2*pi^(11/3)))*r0^(5/3);
            out = out*L0^(11/3)*( (symf*L0)^2 + 1 )^(-11/6);
        end
        
        function out = otf(rho,atm)
            %% OTF Phase optical transfert function
            %
            % out = phaseStats.otf(rho) compute the phase optical transfert function
            % from the baseline rho and an atmosphere object
            
            out = exp(-0.5*phaseStats.structureFunction(rho,atm));
        end
        
        function out = psf(f,atm)
            %% PSF Phase point spread function
            %
            % out = phaseStats.psf(f,atm) compute the phase point spread
            % from the frequency f and an atmosphere object

            fun = @(u) 2.*pi.*quadgk(@(v) psfHankelIntegrandNested(v,u),0,Inf);
            out = arrayfun( fun, f);
            function y = psfHankelIntegrandNested(x,freq)
                y = x.*besselj(0,2.*pi.*x.*freq).*phaseStats.otf(x,atm);
            end
        end
        
        function out = covarianceMatrix(varargin)
            %% COVARIANCEMATRIX Phase covariance matrix
            %
            % out = phaseStats.covarianceMatrix(rho1,atm) Computes the phase
            % auto-covariance matrix from the vector rho1 and an atmosphere
            % object
            %
            % out = phaseStats.covarianceMatrix(rho1,rho2,atm) Computes the phase
            % cross-covariance matrix from the vectors rho1 and rho2 and an
            % atmosphere object
            % Examples:
            % covariance matrix on a 1 metre square grid sampled on 16
            % pixels :
            % [x,y] = meshgrid((0:15)*1/15);
            % g = phaseStats.covarianceMatrix(complex(x,y),atm);
            % imagesc(g), axis square, colorbar
            % covariance matrix on a 1 meter square grid sampled on 16
            % pixels with the same grid but displaced of 1 meter:
            % [x,y] = meshgrid((0:15)*1/15);
            % z = complex(x,y);
            % g = phaseStats.covarianceMatrix(z,z+1,atm);
            % imagesc(g), axis square, colorbar
            %
            % See also atmosphere
            
            error(nargchk(2,3,nargin))
            rho1 = varargin{1}(:);
            if nargin==2
                atm  = varargin{2};
                rho  = abs(bsxfun(@minus,rho1,rho1.'));
            else
                rho2 = varargin{2}(:);
                atm  = varargin{3};
                rho  = abs(bsxfun(@minus,rho1,rho2.'));
            end
            [nRho,mRho] = size(rho);
            blockSize = 5000;
            if max(nRho,mRho)>blockSize % Memory gentle
                l  = floor(nRho/blockSize);
                le = nRho - l*blockSize;
                p  = floor(mRho/blockSize);
                pe = mRho - p*blockSize;
                fprintf(' @(phaseStats.covarianceMatrix)> Memory gentle! (%dX%d blocks)\n',l+1,p+1)
                rho = mat2cell( rho, ...
                    [ones(1,l)*blockSize, le],...
                    [ones(1,p)*blockSize, pe]);
%                 out = cell2mat( ...
%                     cellfun(@(x) phaseStats.covariance(x,atm), rho, 'UniformOutput', false) );
                out = cell(size(rho));
                for ii=1:numel(rho)
                    out{ii} = phaseStats.covariance(rho{ii},atm);
                end
                out = cell2mat(out);
            else % Memory intensive
                L0r0ratio= (atm.L0./atm.r0).^(5./3);
                cst      = (24.*gamma(6./5)./5).^(5./6).*...
                    (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                    L0r0ratio;
                out   = ones(size(rho),class(rho)).*(24.*gamma(6./5)./5).^(5./6).*...
                    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
                index         = rho~=0;
                u             = 2.*pi.*rho(index)./atm.L0;
                out(index) = cst.*u.^(5./6).*besselk(5./6,u);
                layers = atm.layer;
                out = sum([layers.fractionnalR0]).*out;
            end
        end
        
        function out = covarianceToeplitzMatrix(atm,z1,varargin)
            %% COVARIANCETOEPLITZMATRIX Phase covariance matrix
            %
            % out = phaseStats.covarianceMatrix(atm,rho1) Computes the phase
            % auto-covariance matrix from the vector rho1 and an atmosphere
            % object
            %
            % out = phaseStats.covarianceMatrix(atm,rho1,rho2) Computes the phase
            % cross-covariance matrix from the vectors rho1 and rho2 and an
            % atmosphere object
            %
            % out = phaseStats.covarianceMatrix(...,'mask',pupil) 
            %
            % Examples:
            % covariance matrix on a 1 metre square grid sampled on 16
            % pixels :
            % [x,y] = meshgrid((0:15)*1/15);
            % g = phaseStats.covarianceMatrix(complex(x,y),atm);
            % imagesc(g), axis square, colorbar
            % covariance matrix on a 1 meter square grid sampled on 16
            % pixels with the same grid but displaced of 1 meter:
            % [x,y] = meshgrid((0:15)*1/15);
            % z = complex(x,y);
            % g = phaseStats.covarianceMatrix(z,z+1,atm);
            % imagesc(g), axis square, colorbar
            %
            % See also atmosphere
                        
            inputs = inputParser; 
            inputs.addRequired('atm',@isobject);
            inputs.addRequired('z1',@isnumeric);
            inputs.addOptional('z2',z1,@isnumeric);
            inputs.addParamValue('mask',[],@islogical);
            inputs.parse(atm,z1,varargin{:});
            
            atm  = inputs.Results.atm;
            z1   = inputs.Results.z1;
            z2   = inputs.Results.z2;
            mask = inputs.Results.mask;
            
            [nz,mz] = size( z1 );
            
            % First Row
            r  =  phaseStats.covariance( abs(z2-z1(1)) , atm );
            r   = mat2cell( r  , nz , ones(mz,1));
            % First column in first blocks fow
            c  =  phaseStats.covariance( abs( bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' ) ) , atm );
            c   = mat2cell( c , nz , ones(mz,1) );
            % First block rows
            rr   = cellfun( @(x,y) localToeplitz(x,y) , c , r , 'UniformOutput' , false);
            
            % First Column
            c  =  phaseStats.covariance( abs(z1-z2(1)) , atm );
            c   = mat2cell( c  , nz , ones(mz,1));
            % First row in first blocks column
            r  =  phaseStats.covariance( abs( bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' ) ) , atm );
            r   = mat2cell( r , nz , ones(mz,1) );
            % First blocks column
            cc   = cellfun( @(x,y) localToeplitz(x,y) , c , r , 'UniformOutput' , false);
           
            out = cell2mat( localToeplitz(cc,rr) );
            out(~mask,:) = [];
            out(:,~mask) = [];
            
            function t = localToeplitz(c,r)
                % this version works on numeric vector as well as on cells vector
                r = r(:);                               % force column structure
                p = length(r);
                m = length(c);
                x = [r(p:-1:2) ; c(:)];                 % build vector of user data
                cidx = uint16(0:m-1)';
                ridx = uint16(p:-1:1);
                subscripts = cidx(:,ones(p,1)) + ridx(ones(m,1),:);  % Toeplitz subscripts
                t = x(subscripts);                                   % actual data
            end
            
        end

                function varargout = spatioAngularIrregularCovarianceMatrix(sampling, range,wfs,atm,srcAC,varargin)
            %% SPATIOANGULARCOVARIANCEMATRIX Phase spatio-angular covariance meta matrix
            %
            % S = spatioAngularCovarianceMatrix(sampling,range,atm,src1)
            % computes the spatio-angular auto-correlation meta-matrix of the
            % wavefront between all the sources srcAC. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % C = spatioAngularCovarianceMatrix(sampling,range,atm,src1,src2)
            % computes the spatio-angular cross-correlation meta-matrix of the
            % wavefront between all src2 and src1. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % [S,C] = spatioAngularCovarianceMatrix(...) computes both
            % auto- and cross-correlation meta-matrix
            %
            % ... = spatioAngularCovarianceMatrix(...,'mask',mask) restrict
            % the sampling within the mask
            
            inputs = inputParser;
                        inputs.addRequired('sampling',@isnumeric);
            inputs.addRequired('range',@isnumeric);
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann'));
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('srcAC',@(x) isa(x,'source'));
            inputs.addOptional('srcCC',[],@(x) isa(x,'source'));
            inputs.addOptional('tipTilt',false,@islogical);
            inputs.addParameter('mask',true(sampling),@islogical);
            inputs.addParameter('lag',0,@isnumeric);
            inputs.addParameter('xyOutput',[],@isnumeric);
            inputs.parse(sampling,range,wfs,atm,srcAC,varargin{:});
            
            m_srcCC = inputs.Results.srcCC;
            tipTilt = inputs.Results.tipTilt;
            m_mask  = inputs.Results.mask;
            m_tau   = inputs.Results.lag;
            xyOutput= inputs.Results.xyOutput;
                
            %[m_x,m_y] = meshgrid( linspace(-1,1,sampling)*range/2 );
            m_nGs   = length(srcAC);
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            m_cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            m_cstL0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            m_cstr0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            m_L0      = atm.L0;
            
            m_nLayer   = atm.nLayer;
            layers   = atm.layer;
            m_altitude = [layers.altitude];
            m_fr0      = [layers.fractionnalR0];
            [m_windVx,m_windVy] = pol2cart([layers.windDirection],[layers.windSpeed]);
            m_srcACdirectionVector = cat(2,srcAC.directionVector);
            m_srcACheight          = [srcAC.height];
            
            if nargout==2
                
                tic
                varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cstL0,m_cst);
                toc
                
                m_xAC = m_x(m_mask);
                m_yAC = m_y(m_mask);
                if isempty(xyOutput)
                    m_xCC = m_xAC;
                    m_yCC = m_yAC;
                else
                    m_xCC = xyOutput(:,1);
                    m_yCC = xyOutput(:,2);
                end
                
                tic                
                varargout{2} = crossCorrelation(sampling,range,wfs,m_mask,m_srcCC,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cstL0,m_cst,m_windVx,m_windVy,m_tau);
                toc
                
            else
                
                if isempty(m_srcCC)
                    
                    varargout{1} = autoCorrelation(sampling,range,wfs,m_mask,...
                        m_nGs,m_srcACdirectionVector,m_srcACheight,...
                        m_nLayer,m_altitude,m_fr0,...
                        m_L0,m_cstL0,m_cst);
                    
                else
                    
                    if tipTilt
                        varargout{1} = crossCorrelationTT(m_srcCC,m_x,m_y,m_mask,...
                            m_nGs,m_srcACdirectionVector,m_srcACheight,...
                            m_nLayer,m_altitude,m_fr0,...
                            m_L0,m_cstr0,range);
                    else
                        varargout{1} = crossCorrelation(sampling,range,wfs,m_mask,m_srcCC,...
                            m_nGs,m_srcACdirectionVector,m_srcACheight,...
                            m_nLayer,m_altitude,m_fr0,...
                            m_L0,m_cstL0,m_cst,m_windVx,m_windVy,m_tau);
                    end
                    
                end
                
            end
            
            
            function C = crossCorrelation(sampling,range,wfs,mask,srcCC,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstL0,cst,windVx,windVy,tau)
                    
                fprintf(' -->> Cross-correlation meta-matrix!\n')
                
                nSs = length(srcCC);
                srcCCdirectionVector = cat(2,srcCC.directionVector);
                srcCCheight          = [srcCC.height];
                C = cellfun( @(x) zeros(sampling,sampling) , cell(nSs,nGs) , 'UniformOutput' , false );
                %cstL0CC = cstL0*ones(length(xCC),length(xAC));
                
                for k=1:nSs*nGs
                    
                    [kGs,iGs] = ind2sub([nSs,nGs],k);
                    buf = 0;
                    
                        [x1,y1] = meshgrid( linspace(-1,1,sampling)*range/2 );
                        [x1,y1] = tools.rotateDM(x1(:),y1(:),wfs.lenslets.rotation(iGs)*180/pi);
                        x1 = x1 - wfs.lenslets.offset(1,iGs)*range;
                        y1 = y1 - wfs.lenslets.offset(2,iGs)*range;
                        x1 = reshape(x1, sampling, sampling);
                        y1 = reshape(y1, sampling, sampling);
                        
                        [x2,y2] = meshgrid( linspace(-1,1,sampling)*range/2 );
                        %[x2,y2] = tools.rotateDM(x2(:),y2(:),wfs.lenslets.rotation(kGs)*180/pi);
                        %x2 = x2 + wfs.lenslets.offset(1,kGs)*range;
                        %y2 = y2 + wfs.lenslets.offset(1,kGs)*range;
                        %x2 = reshape(x2, sampling, sampling);
                        %y2 = reshape(y2, sampling, sampling);
                        
                        for kLayer=1:nLayer
                            
                            %                         atmSlab = slab(atm,kLayer);
                            
                            beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(iGs);
                            iZ = complex( x1*scale + beta(1) , y1*scale + beta(2) );
                            
                            beta = srcCCdirectionVector(:,kGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcCCheight(kGs);
                            jZ  = complex( x2*scale + beta(1) , y2*scale + beta(2) );
                            
                            out = phaseStats.covarianceMatrix(jZ, iZ, slab(atm,kLayer));
                            out(~mask,:) = [];
                            out(:,~mask) = [];
                            
                            buf = buf + out;
                            
                        end
                    
                    C{k} = buf;
                    
                end
                
                    buf = C;
                    C = cell(nSs,1);
                    for k=1:nSs
                        C{k} = cell2mat(buf(k,:));
                    end
                
            end
            
            function S = autoCorrelation(sampling,range,wfs,mask,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstL0,cst)
                    
                    fprintf(' -->> Auto-correlation meta-matrix!\n')
                    
                    kGs = reshape( triu( reshape(1:nGs^2,nGs,nGs) , 1) , 1, []);
                    kGs(1) = 1;
                    kGs(kGs==0) = [];
                    S = cellfun( @(x) zeros(sum(mask(:))) , cell(1,length(kGs)) , 'UniformOutput' , false );
                    
                    for k=1:length(kGs)
                        
                        [iGs,jGs] = ind2sub( [nGs,nGs] , kGs(k) );
                        buf = 0;
                        
                        [x1,y1] = meshgrid( linspace(-1,1,sampling)*range/2 );
                        [x1,y1] = tools.rotateDM(x1(:),y1(:),wfs.lenslets.rotation(iGs)*180/pi);
                        x1 = x1 - wfs.lenslets.offset(1,iGs)*range;
                        y1 = y1 - wfs.lenslets.offset(2,iGs)*range;
                        x1 = reshape(x1, sampling, sampling);
                        y1 = reshape(y1, sampling, sampling);
                        
                        [x2,y2] = meshgrid( linspace(-1,1,sampling)*range/2 );
                        [x2,y2] = tools.rotateDM(x2(:),y2(:),wfs.lenslets.rotation(jGs)*180/pi);
                        x2 = x2 - wfs.lenslets.offset(1,jGs)*range;
                        y2 = y2 - wfs.lenslets.offset(2,jGs)*range;
                        x2 = reshape(x2, sampling, sampling);
                        y2 = reshape(y2, sampling, sampling);
                        
                        for kLayer=1:nLayer
                            
                            %                         atmSlab = slab(atm,kLayer);
                            
                            beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(iGs);
                            iZ = complex( x1*scale + beta(1) , y1*scale + beta(2) );
                            
                            beta = srcACdirectionVector(:,jGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(jGs);
                            jZ  = complex( x2*scale + beta(1) , y2*scale + beta(2) );
                            
                            out = phaseStats.covarianceMatrix(iZ, jZ, slab(atm,kLayer));
                            out(~mask,:) = [];
                            out(:,~mask) = [];
                            
                            buf = buf + out;
                            
                        end
                        
                        S{k} = buf;
                        
                    end
                    
                    buf = S;
                    S = cellfun( @(x) zeros(sum(mask(:))) , cell(nGs) , 'UniformOutput' , false );
                    S(kGs) = buf;
                    S(1:nGs+1:nGs^2) = S(1,1);
                    S = cell2mat(S);
                    S = triu(S,1)+triu(S)';

            end
            
            
        end

        function varargout = spatioAngularCovarianceMatrix(sampling,range,atm,srcAC,varargin)
            %% SPATIOANGULARCOVARIANCEMATRIX Phase spatio-angular covariance meta matrix
            %
            % S = spatioAngularCovarianceMatrix(sampling,range,atm,src1)
            % computes the spatio-angular auto-correlation meta-matrix of the
            % wavefront between all the sources srcAC. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % C = spatioAngularCovarianceMatrix(sampling,range,atm,src1,src2)
            % computes the spatio-angular cross-correlation meta-matrix of the
            % wavefront between all src2 and src1. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % [S,C] = spatioAngularCovarianceMatrix(...) computes both
            % auto- and cross-correlation meta-matrix
            %
            % ... = spatioAngularCovarianceMatrix(...,'mask',mask) restrict
            % the sampling within the mask
            
            inputs = inputParser;
            inputs.addRequired('sampling',@isnumeric);
            inputs.addRequired('range',@isnumeric);
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('srcAC',@(x) isa(x,'source'));
            inputs.addOptional('srcCC',[],@(x) isa(x,'source'));
            inputs.addOptional('tipTilt',false,@islogical);
            inputs.addParameter('mask',true(sampling),@islogical);
            inputs.addParameter('lag',0,@isnumeric);
            inputs.addParameter('xyOutput',[],@isnumeric);
            inputs.parse(sampling,range,atm,srcAC,varargin{:});
            
            m_srcCC = inputs.Results.srcCC;
            tipTilt = inputs.Results.tipTilt;
            m_mask  = inputs.Results.mask;
            m_tau   = inputs.Results.lag;
            xyOutput= inputs.Results.xyOutput;
                
            [m_x,m_y] = meshgrid( linspace(-1,1,sampling)*range/2 );
            m_nGs   = length(srcAC);
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            m_cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            m_cstL0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            m_cstr0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            m_L0      = atm.L0;
            
            m_nLayer   = atm.nLayer;
            layers   = atm.layer;
            m_altitude = [layers.altitude];
            m_fr0      = [layers.fractionnalR0];
            [m_windVx,m_windVy] = pol2cart([layers.windDirection],[layers.windSpeed]);
            m_srcACdirectionVector = cat(2,srcAC.directionVector);
            m_srcACheight          = [srcAC.height];
            
            if nargout==2
                
                tic
                varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cstL0,m_cst);
                toc
                
                m_xAC = m_x(m_mask);
                m_yAC = m_y(m_mask);
                if isempty(xyOutput)
                    m_xCC = m_xAC;
                    m_yCC = m_yAC;
                else
                    m_xCC = xyOutput(:,1);
                    m_yCC = xyOutput(:,2);
                end
                
                tic                
                varargout{2} = crossCorrelation(m_srcCC,m_xAC,m_yAC,m_xCC,m_yCC,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cstL0,m_cst,m_windVx,m_windVy,m_tau);
                toc
                
            else
                
                if isempty(m_srcCC)
                    
                    varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                        m_nGs,m_srcACdirectionVector,m_srcACheight,...
                        m_nLayer,m_altitude,m_fr0,...
                        m_L0,m_cstL0,m_cst);
                    
                else
                    
                    if tipTilt
                        varargout{1} = crossCorrelationTT(m_srcCC,m_x,m_y,m_mask,...
                            m_nGs,m_srcACdirectionVector,m_srcACheight,...
                            m_nLayer,m_altitude,m_fr0,...
                            m_L0,m_cstr0,range);
                    else
                        m_xAC = m_x(m_mask);
                        m_yAC = m_y(m_mask);
                        if isempty(xyOutput)
                            m_xCC = m_xAC;
                            m_yCC = m_yAC;
                        else
                            m_xCC = xyOutput(:,1);
                            m_yCC = xyOutput(:,2);
                        end
                        varargout{1} = crossCorrelation(m_srcCC,m_xAC,m_yAC,m_xCC,m_yCC,...
                            m_nGs,m_srcACdirectionVector,m_srcACheight,...
                            m_nLayer,m_altitude,m_fr0,...
                            m_L0,m_cstL0,m_cst,m_windVx,m_windVy,m_tau);
                    end
                    
                end
                
            end
            
            function C = crossCorrelationTT(srcTT,x,y,mask,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstr0,D)
                    
                fprintf(' -->> Cross-correlation meta-matrix (phase/tip-tilt)!\n')
                
                nSs = length(srcTT);
                srcTTdirectionVector = cat(2,srcTT.directionVector);
                C = cellfun( @(x) zeros(sum(mask(:))) , cell(nSs,nGs) , 'UniformOutput' , false );
                x = x(mask);
                y = y(mask);
                f02 = 1./L0.^2;
                
                for k=1:nSs*nGs
                    
                    [kSs,iGs] = ind2sub([nSs,nGs],k);
                    buf = 0;
                    
                    for kLayer=1:nLayer
                        
                        beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(iGs);
                        iZ = complex( x*scale + beta(1) , y*scale + beta(2) );
                        
                        betaTt = srcTTdirectionVector(:,kSs)*altitude(kLayer);
                        zTt = iZ.' - complex(betaTt(1),betaTt(2));
                        rTt = abs(zTt);
                        oTt = angle(zTt);
                        
                        Inm = @(r) quadgk( @(f) (f.^2 + f02).^(-11./6).*...
                            besselj(2,pi.*f.*D).*...
                            besselj(1,2.*pi.*f.*r),0,Inf);
                        out = fr0(kLayer)*cstr0.*arrayfun(Inm,rTt).*8/D/scale;

                        buf = buf + ...
                            [ out.*cos(oTt) ; out.*sin(oTt) ];
                        
                    end
                    
                    C{k} = buf;
                    
                end
                
                    buf = C;
                    C = cell(nSs,1);
                    for k=1:nSs
                        C{k} = cell2mat(buf(k,:));
                    end
                
            end
            
            
            function C = crossCorrelation(srcCC,xAC,yAC,xCC,yCC,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstL0,cst,windVx,windVy,tau)
                    
                fprintf(' -->> Cross-correlation meta-matrix!\n')
                
                nSs = length(srcCC);
                srcCCdirectionVector = cat(2,srcCC.directionVector);
                srcCCheight          = [srcCC.height];
                C = cellfun( @(x) zeros(length(xCC),length(xAC)) , cell(nSs,nGs) , 'UniformOutput' , false );
                cstL0CC = cstL0*ones(length(xCC),length(xAC));
                
                for k=1:nSs*nGs
                    
                    [kSs,iGs] = ind2sub([nSs,nGs],k);
                    buf = 0;
                    
                    for kLayer=1:nLayer
                        
                        
                        beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(iGs);
                        iZ = complex( xAC*scale + beta(1) , yAC*scale + beta(2) );
                        
                        betaSs = srcCCdirectionVector(:,kSs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcCCheight(kSs);
                        zSs = complex( ...
                            xCC*scale + betaSs(1) + windVx(kLayer)*tau, ...
                            yCC*scale + betaSs(2) + windVy(kLayer)*tau );
                        
                        rho   = abs(bsxfun(@minus,zSs,iZ.'));
                        out   = cstL0CC;
                        index = rho~=0;
                        u          = 2.*pi.*rho(index)./L0;
                        out(index) = cst.*u.^(5./6).*besselk(5./6,u);
                        
                        buf = buf + fr0(kLayer)*out;
                        
                    end
                    
                    C{k} = buf;
                    
                end
                
                    buf = C;
                    C = cell(nSs,1);
                    for k=1:nSs
                        C{k} = cell2mat(buf(k,:));
                    end
                
            end
            
            function S = autoCorrelation(x,y,mask,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstL0,cst)
                    
                    fprintf(' -->> Auto-correlation meta-matrix!\n')
                    
                    kGs = reshape( triu( reshape(1:nGs^2,nGs,nGs) , 1) , 1, []);
                    kGs(1) = 1;
                    kGs(kGs==0) = [];
                    S = cellfun( @(x) zeros(sum(mask(:))) , cell(1,length(kGs)) , 'UniformOutput' , false );
                    cstL0AC = cstL0*ones(sampling);
                    
                    for k=1:length(kGs)
                        
                        [iGs,jGs] = ind2sub( [nGs,nGs] , kGs(k) );
                        buf = 0;
                        
                        for kLayer=1:nLayer
                            
                            %                         atmSlab = slab(atm,kLayer);
                            
                            beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(iGs);
                            iZ = complex( x*scale + beta(1) , y*scale + beta(2) );
                            
                            beta = srcACdirectionVector(:,jGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(jGs);
                            jZ  = complex( x*scale + beta(1) , y*scale + beta(2) );
                            
                            z1 = iZ;
                            z2 = jZ;
                            [nz,mz] = size( z1 );
                            % First Row
                            %                         r  =  phaseStats.covariance( abs(z2-z1(1)) , atm );
                            rho = abs(z2-z1(1));
                            r   = cstL0AC;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            r(index) = cst.*u.^(5./6).*besselk(5./6,u);
                            
                            r   = mat2cell( fr0(kLayer)*r  , nz , ones(mz,1));
                            % First column in first blocks fow
                            %                         c  =  phaseStats.covariance( abs( bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' ) ) , atm );
                            rho = abs( bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' ) );
                            c   = cstL0AC;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            c(index) = cst.*u.^(5./6).*besselk(5./6,u);
                            
                            c   = mat2cell( fr0(kLayer)*c , nz , ones(mz,1) );
                            % First block rows
                            rr   = cellfun( @(x,y) myToeplitz(x,y) , c , r , 'UniformOutput' , false);
                            
                            % First Column
                            %                         c  =  phaseStats.covariance( abs(z1-z2(1)) , atm );
                            rho = abs(z1-z2(1));
                            c   = cstL0AC;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            c(index) = cst.*u.^(5./6).*besselk(5./6,u);
                            
                            c   = mat2cell( fr0(kLayer)*c  , nz , ones(mz,1));
                            % First row in first blocks column
                            %                         r  =  phaseStats.covariance( abs( bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' ) ) , atm );
                            rho = abs( bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' ) );
                            r   = cstL0AC;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            r(index) = cst.*u.^(5./6).*besselk(5./6,u);
                            
                            r   = mat2cell( fr0(kLayer)*r , nz , ones(mz,1) );
                            % First blocks column
                            cc   = cellfun( @(x,y) myToeplitz(x,y) , c , r , 'UniformOutput' , false);
                            
                            out = cell2mat( myToeplitz(cc,rr) );
                            out(~mask,:) = [];
                            out(:,~mask) = [];
                            
                            buf = buf + out;
                            
                        end
                        
                        S{k} = buf;
                        
                    end
                    
                    buf = S;
                    S = cellfun( @(x) zeros(sum(mask(:))) , cell(nGs) , 'UniformOutput' , false );
                    S(kGs) = buf;
                    S(1:nGs+1:nGs^2) = S(1,1);
                    S = cell2mat(S);
                    S = triu(S,1)+triu(S)';

            end
            
            
        end

        function out = upwardJitter(L,D,atm,zSrc)
            %% UPWARDJITTER Upward propagated beam motion
            %
            % varJit = upwardJitter(L,D,atm) computes the variance of the
            % upward propagated beam motion at a target at distance L from
            % source, the beam is lauched from a telescope of diameter D
            % and is propagated through the atmosphere object atm
            %
            % varJit = upwardJitter(L,D,atm,zSrc) computes the variance of
            % the upward propagated beam motion at a target at distance L
            % from source at zSrc, the beam is lauched from a telescope of
            % diameter D and is propagated through the atmosphere object
            % atm
            
            if nargin<4
                zSrc = L;
            end
            k0  = 2*pi/atm.wavelength;
            fun = @(f,z,atmSlab) f.*phaseStats.spectrum(f,atmSlab).*...
                (L-z).^2.*...
                (16/(k0.*(1-z./zSrc).*D)).^2.*...
                (besselj(2,2*pi*f*(1-z./zSrc)*D/2)./(2*pi*f*(1-z./zSrc)*D/2)).^2;
            out = 0;
            for kLayer=1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                z = atmSlab.layer.altitude;
                out = out + 2.*pi.*quadgk( @(x) fun(x,z,atmSlab) , 0 , Inf);
            end
        end
        function out = zernikeVariance(zern,atm)
            %% ZERNIKEVARIANCE Zernike coefficients variance
            %
            % out = variance(modes,atmosphere) computes the
            % variance of Zernike coefficients from the modes and the
            % atmosphere object
            %
            % out = variance(zernike,atmosphere) computes the
            % variance of Zernike coefficients from the Zernike polynomials
            % object and the atmosphere object
            %
            % Example:
            % atm = atmosphere(photometry.V,0.15,30);
            % tel = telescope(10);
            % modes = 1:15;
            % figure
            % semilogy(modes,phaseStats.zernikeVariance(modes,atm),'--.')
            % xlabel('Zernike modes')
            % ylabel('Variance [rd^2]')
            %
            % See also zernike, atmosphere
            
            
            if ~isa(zern,'zernike')
                zern = zernike(zern);
            end
            r0 = atm.r0;
            L0 = atm.L0;
            D  = zern.D;
            jv = zern.j;
            nv = zern.n;
            mv = zern.m;
            nv0   = nv;
            index = diff(nv)~=0;
            jv = [jv(index) jv(end)];
            mv = [mv(index) mv(end)];
            nv = [nv(index) nv(end)];
            nf  = length(nv);
            out = zeros(length(zern.j),1);
            
            for cpt = 1:nf
                
                j = jv(cpt);
                n = nv(cpt);
                m = mv(cpt);
                
                out(nv0==n,1) = zernCovCoef(r0,L0,D,j,j,n,m,n,m);
                
            end
            function out = zernCovCoef(r0,L0,D,i,j,ni,mi,nj,mj)
                if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                    if L0==Inf
                        if i==1 && j==1
                            out = Inf;
                        else
                            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                                (D./r0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                                newGamma(-5./6+(ni+nj)./2,...
                                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                        end
                    else
                        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                            (L0./r0).^(5./3).*(L0./D).^2.*...
                            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*D./L0);
                    end
                else
                    out = 0;
                end
                function out = newGamma(a,b)
                    % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                    % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                    % out = newGamma(a,b)
                    
                    out = prod(gamma(a))./prod(gamma(b));
                end
                function out = UnParamEx4q2(mu,alpha,beta,p,a)
                    % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                    % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                    % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                    % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                    % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                    
                    a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                    b1 = [1+alpha+beta 1+alpha 1+beta];
                    a2 = [(1-mu)./2+p 1+p p];
                    b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                    
                    out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                        newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                        pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                        newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                        pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                    function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                        % POCHAMMERSERIES Computes power series in Pochammer notation
                        % pochammerSeries(p,q,a,b,z)
                        % pochammerSeries(p,q,a,b,z,tol)
                        % pochammerSeries(p,q,a,b,z,[],nmax)
                        % pochammerSeries(p,q,a,b,z,tol,nmax)
                        
                        if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                            
                            if p==length(a) && q==length(b)
                                
                                switch nargin
                                    case 6
                                        nmax = 1e3;
                                    case 7
                                        if isempty(tol)
                                            tol = 1e-6;
                                        end
                                    otherwise
                                        tol = 1e-6;
                                        nmax = 1e3;
                                end
                                
                                out = zeros(size(z));
                                
                                indz = find(z==0);
                                if ~isempty(indz)
                                    out(indz) = 1;
                                end
                                
                                indnz = find(z~=0);
                                if ~isempty(indnz)
                                    z = z(indnz);
                                    ck = 1;
                                    step = Inf;
                                    k = 0;
                                    som = ck;
                                    while (k<=nmax) && (step>tol)
                                        ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                        step = abs(abs(ck)-abs(ckp1));
                                        som = som + ckp1;
                                        k = k+1;
                                        ck = ckp1;
                                    end
                                    if step>tol
                                        warning('pochammerSeries','Maximum iteration reached before convergence')
                                    end
                                    out(indnz) = som;
                                end
                                
                            else
                                error('p and q must be the same length than vectors a and b, respectively')
                                
                            end
                            
                        else
                            error('This generalized hypergeometric function doesn''t converge')
                        end
                    end
                end
                
            end
        end
        
        function out = zernikeCovariance(zern,atm)
            %% ZERNIKECOVARIANCE Zernike coefficients covariance
            %
            % out = phaseStats.zernikeCovariance(modes,atmosphere)
            % computes the covariance matrix of Zernike coefficients from
            % the modes and the atmosphere object
            %
            % out = phaseStats.zernikeCovariance(zernike,atmosphere) computes the
            % covariance matrix of Zernike coefficients from the Zernike
            % polynomials object and the atmosphere object
            %
            % Example:
            % atm = atmosphere(photometry.V,0.15,30);
            % tel = telescope(10);
            % modes = 1:15;
            % figure
            % spy(phaseStats.zernikeCovariance(modes,atm,tel))
            %
            % See also zernike, atmosphere
            
            
            if ~isa(zern,'zernike')
                zern = zernike(zern);
            end
            r0 = atm.r0;
            L0 = atm.L0;
            D  = zern.D;
            % -- covariance --
            n = zern.n;
            m = zern.m;
            [i,j]   = meshgrid(zern.j);
            [ni,nj] = meshgrid(n);
            [mi,mj] = meshgrid(m);
            index   = (mi==mj) & (rem(abs(i-j),2)==0 | ((mi==0) & (mj==0)));
            index   = logical(index - diag(diag(index)));
            i       = i(index);
            j       = j(index);
            ni      = ni(index);
            nj      = nj(index);
            mi      = mi(index);
            mj      = mj(index);
            
            nf       = sum(index(:));
            covValue = zeros(1,nf);
            for cpt = 1:nf
                covValue(cpt) = zernCovCoef(r0,L0,D,i(cpt),j(cpt),ni(cpt),mi(cpt),nj(cpt),mj(cpt));
            end
            out          = zeros(length(zern.j));
            out(index)   = covValue;
            % -- variance --
            jv = zern.j;
            nv = zern.n;
            mv = zern.m;
            nv0   = nv;
            index = diff(nv)~=0;
            jv = [jv(index) jv(end)];
            mv = [mv(index) mv(end)];
            nv = [nv(index) nv(end)];
            nf  = length(nv);
            variance = zeros(length(zern.j),1);
            
            for cpt = 1:nf
                
                j = jv(cpt);
                n = nv(cpt);
                m = mv(cpt);
                
                variance(nv0==n,1) = zernCovCoef(r0,L0,D,j,j,n,m,n,m);
                
            end
            % -- covariance + variance --
            out = out + diag(variance);
            function out = zernCovCoef(r0,L0,D,i,j,ni,mi,nj,mj)
                if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                    if L0==Inf
                        if i==1 && j==1
                            out = Inf;
                        else
                            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                                (D./r0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                                newGamma(-5./6+(ni+nj)./2,...
                                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                        end
                    else
                        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                            (L0./r0).^(5./3).*(L0./D).^2.*...
                            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*D./L0);
                    end
                else
                    out = 0;
                end
                function out = newGamma(a,b)
                    % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                    % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                    % out = newGamma(a,b)
                    
                    out = prod(gamma(a))./prod(gamma(b));
                end
                function out = UnParamEx4q2(mu,alpha,beta,p,a)
                    % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                    % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                    % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                    % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                    % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                    
                    a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                    b1 = [1+alpha+beta 1+alpha 1+beta];
                    a2 = [(1-mu)./2+p 1+p p];
                    b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                    
                    out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                        newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                        pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                        newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                        pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                    function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                        % POCHAMMERSERIES Computes power series in Pochammer notation
                        % pochammerSeries(p,q,a,b,z)
                        % pochammerSeries(p,q,a,b,z,tol)
                        % pochammerSeries(p,q,a,b,z,[],nmax)
                        % pochammerSeries(p,q,a,b,z,tol,nmax)
                        
                        if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                            
                            if p==length(a) && q==length(b)
                                
                                switch nargin
                                    case 6
                                        nmax = 1e3;
                                    case 7
                                        if isempty(tol)
                                            tol = 1e-6;
                                        end
                                    otherwise
                                        tol = 1e-6;
                                        nmax = 1e3;
                                end
                                
                                out = zeros(size(z));
                                
                                indz = find(z==0);
                                if ~isempty(indz)
                                    out(indz) = 1;
                                end
                                
                                indnz = find(z~=0);
                                if ~isempty(indnz)
                                    z = z(indnz);
                                    ck = 1;
                                    step = Inf;
                                    k = 0;
                                    som = ck;
                                    while (k<=nmax) && (step>tol)
                                        ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                        step = abs(abs(ck)-abs(ckp1));
                                        som = som + ckp1;
                                        k = k+1;
                                        ck = ckp1;
                                    end
                                    if step>tol
                                        warning('pochammerSeries','Maximum iteration reached before convergence')
                                    end
                                    out(indnz) = som;
                                end
                                
                            else
                                error('p and q must be the same length than vectors a and b, respectively')
                                
                            end
                            
                        else
                            error('This generalized hypergeometric function doesn''t converge')
                        end
                    end
                end
                
            end
        end
        
        function out = zernikeResidualVariance(N,atm,tel)
            %% ZERNIKERESIDUALVARIANCE
            %
            % out = zernikeResidualVariance(N,atm,tel)
            %
            % See also: atmosphere, telescope
            
            zern = zernike(1:N,tel.D,'logging',false);
            r0 = atm.r0;
            L0 = atm.L0;
            D  = tel.D;
            aiVar = phaseStats.zernikeVariance(zern,atm);
            atmVar = phaseStats.variance(atm);
            if isinf(L0)
                Delta1 = -(2.*gamma(11./6).^2./pi.^1.5).*(24.*gamma(6./5)./5).^(5./6).*...
                    (D./r0).^(5./3).*newGamma([-5./6,7./3],[23./6,17./6]);
                out = Delta1 - sum(aiVar(2:end));
            else
                out = atmVar - sum(aiVar);
            end
        end
        
        function aiaj = zernikeAngularCovariance(zern,atm,src,optSrc)
            %% ZERNIKEANGULARCOVARIANCE Zernike coefficients angular covariance
            %
            % aiaj = zernikeAngularCovariance(zern,atm,src) computes
            % the covariance matrix between Zernike coefficients of Zernike
            % polynomials zern corresponding to wavefront propagating from
            % two sources src(1) and src(2) through the atmosphere atm
            %atm
            % See also zernike, atmosphere, source
            
            nGs = numel(src);
            if nGs>2 % then its a meta-matrix
                disp(' @(phaseStats.zernikeAngularCovariance)> META-MATRIX:')
                if nargin<4 % a correlation meta-matrix
                    disp(' @(phaseStats.zernikeAngularCovariance)> AUTO CORRELATION META-MATRIX:')
                    iSrc = src;
                    jSrc = src;
                    mGs = nGs;
                    aiaj = cell(nGs,mGs);
%                     for iGs = 1:nGs
%                         fprintf(' @(Data covariance)> ');
%                         gsCurrent = iSrc(iGs);
%                         for jGs = iGs:mGs
%                             fprintf('gs#%d/gs#%d - ',iGs,jGs)
%                             aiaj{iGs,jGs} = phaseStats.zernikeAngularCovariance(zern,atm,[gsCurrent,jSrc(jGs)]);
%                         end
%                         fprintf('\b\b\b\n')
%                     end
                    nmGs  = [nGs mGs];
                    index = 1:nGs*mGs;
                    mask  = triu(true(nGs,mGs));
                    index = index(mask);
                    buffer = cell(1,length(index));
                    for kGs = 1:length(index)
                        ijGs = index(kGs);
                        [iGs,jGs] = ind2sub(nmGs,ijGs);
%                         fprintf(' @(phaseStats.zernikeAngularCovariance)> gs#%d/gs#%d \n',iGs,jGs);
                        buffer{kGs} = phaseStats.zernikeAngularCovariance(zern,atm,[iSrc(iGs),jSrc(jGs)]);
                    end
                    aiaj(mask) = buffer;
                    index = cellfun(@isempty,aiaj);
                    % address cells colums-wise, otherwise aiaj fails to be
                    % symmetric
                    [r,c] = find(triu(~index,1));
                    [r,s] = sort(r);
                    blocks = sub2ind(size(index),r,c(s));
                    aiaj(index) = cellfun(@transpose,aiaj(blocks),'UniformOutput',false);
                    
                    %index = cellfun(@isempty,aiaj);
                    %aiaj(index) = cellfun(@transpose,aiaj(triu(~index,1)),'UniformOutput',false);
                    
                else % a cross-correlation meta-matrix
%                     disp('@(phaseStats.zernikeAngularCovariance)> CROSS CORRELATION META-MATRIX:')
                    iSrc = src;
                    jSrc = optSrc;
                    mGs = numel(jSrc);
                    aiaj = cell(nGs,mGs);
%                     for iGs = 1:nGs
%                         fprintf(' @(phaseStats.zernikeAngularCovariance)> ');
%                         gsCurrent = iSrc(iGs);
%                         for jGs = 1:mGs
%                             fprintf('gs#%d/gs#%d - ',iGs,jGs);
%                             aiaj{iGs,jGs} = phaseStats.zernikeAngularCovariance(zern,atm,[gsCurrent,jSrc(jGs)]);
%                         end
%                         fprintf('\b\b\b\n')
%                     end
                    nmGs  = [nGs mGs];
                    for kGs = 1:nGs*mGs
                        [iGs,jGs] = ind2sub(nmGs,kGs);
%                         fprintf(' @(phaseStats.zernikeAngularCovariance)> gs#%d/gs#%d \n',iGs,jGs);
                        aiaj{kGs} = phaseStats.zernikeAngularCovariance(zern,atm,[iSrc(iGs),jSrc(jGs)]);
                    end
                end
%                 aiaj = cell2mat(aiaj);
            else
                if src(1)==src(2) && all( isinf( [src(1).height src(2).height] ) )
                    aiaj = phaseStats.zernikeCovariance(zern,atm);
                else
                    R   = zern.R;
                    zs1 = src(1).height;
                    zs2 = src(2).height;
                    xSrc = tan(src(1).zenith).*cos(src(1).azimuth) - ...
                        tan(src(2).zenith).*cos(src(2).azimuth);
                    ySrc = tan(src(1).zenith).*sin(src(1).azimuth) - ...
                        tan(src(2).zenith).*sin(src(2).azimuth);
                    rhoSrcLayer = hypot(xSrc,ySrc);
                    thetaSrcLayer = atan2(ySrc,xSrc);
                    nMode = length(zern.j);
                    znmj = mat2cell([zern.j;zern.n;zern.m],3,ones(1,nMode));
                    znmj = repmat(znmj,nMode,1);
                    znmi = znmj';
                    psdCst = (24.*gamma(6./5)./5).^(5./6).*...
                        (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                        atm.r0.^(-5./3);
                    if all( isinf( [zs1 zs2] ) ) % NGS CASE
                        a1l     = R*ones(1,atm.nLayer);
                        a2l     = R*ones(1,atm.nLayer);
                        denom   = pi.*a1l.*a2l;
                        sl      = [atm.layer.altitude]'.*rhoSrcLayer;
                        fr0     = [atm.layer.fractionnalR0]';
                        aiajFun = @ (znmi,znmj) ...
                            quadgk(@(x) integrand(x,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)), ...
                            0, Inf, 'AbsTol',1e-3, 'RelTol',1e-2);
%                         n = 201;
%                         r = linspace(0,20,n);
%                         r(1) = 1e-6;
%                         aiajFun = @ (znmi,znmj) ...
%                             trapz(r,integrandNgs(r,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)));                       
                    else % LGS CASE (TO DO: optimize for LGS as for NGS)
                        a1l   = zeros(1,atm.nLayer);
                        a2l   = zeros(1,atm.nLayer);
                        denom = zeros(1,atm.nLayer);
                        sl    = zeros(1,atm.nLayer);
                        for lLayer=1:atm.nLayer
                            a1l(lLayer) = R.*(1 - atm.layer(lLayer).altitude./zs1);
                            a2l(lLayer) = R.*(1 - atm.layer(lLayer).altitude./zs2);
                            denom(lLayer) = pi.*a1l(lLayer).*a2l(lLayer);
                            sl(lLayer) = atm.layer(lLayer).altitude.*rhoSrcLayer;
                        end
                        aiajFun = @ (znmi,znmj) ...
                            quadgk(@(x) integrand(x,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)), ...
                            0, Inf, 'AbsTol',1e-3, 'RelTol',1e-2);
%                         n = 201;
%                         r = linspace(0,20,n);
%                         r(1) = 1e-6;
%                         aiajFun = @ (znmi,znmj) ...
%                             trapz(r,integrandNgs(r,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)));                       
                    end
                    aiaj = zeros(nMode);
                    index = triu(true(nMode));
                    %                 tic
                    aiaj(index) = cellfun(aiajFun,znmj(index),znmi(index));
                    %                 toc
                    aiaj = aiaj + triu(aiaj,1)';
                    aiaj = bsxfun(@times,aiaj,(-1).^zern.m');
                    %                     aiaj = cellfun(aiajFun,znmj,znmi);
                end
            end
            function out = integrand(x,zi,ni,mi,zj,nj,mj)
                krkr_mi = mi==0;
                krkr_mj = mj==0;
                out = 0;
                factor1 = sqrt((ni+1)*(nj+1)).*...
                    (-1).^(0.5.*(ni+nj)).*...
                    2.^(1-0.5.*(krkr_mi+krkr_mj));%.*...
                %(-1).^mj;
                factor2 = (-1).^(1.5*(mi+mj)).*...
                    cos( ...
                    (mi+mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) + ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                factor3 = (-1).^(1.5*abs(mi-mj)).*...
                    cos( ...
                    (mi-mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) - ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                for kLayer=1:atm.nLayer
%                     a1l = R.*(1 - atm.layer(kLayer).altitude./zs1);
%                     a2l = R.*(1 - atm.layer(kLayer).altitude./zs2);
%                     denom = pi.*a1l.*a2l;
%                     sl = atm.layer(kLayer).altitude.*rhoSrcLayer;
                    red1 = a1l(kLayer).*x;
                    red2 = a2l(kLayer).*x;
                    red = sl(kLayer).*x;
%                     phasePSD = phaseStats.spectrum(0.5*x/pi,atmLayers{kLayer});
                    f = 0.5*x/pi;
                    %phasePSD = atm.layer(kLayer).fractionnalR0.*...
                    %    psdCst.*(f.^2 + 1./atm.L0.^2).^(-11./6);
                    phasePSD = atm.layer(kLayer).fractionnalR0.*...
                        psdCst.*(f.^2 + 1./atm.layer(kLayer).layeredL0.^2).^(-11./6);
%                     phasePSD = phaseStats.spectrum(0.5*x/pi,atm.slab(lLayer));
                    besselsRadialOrder = besselj(ni+1,red1).*besselj(nj+1,red2);
                    tripleBessel1 = besselj(mi+mj,red);
                    tripleBessel2 = besselj(abs(mi-mj),red);
%                     besselsRadialOrder = besselmx('J',ni+1,red1).*besselmx('J',nj+1,red2);
%                     tripleBessel1 = besselmx('J',mi+mj,red);
%                     tripleBessel2 = besselmx('J',abs(mi-mj),red);
                    out = out +  (factor1./denom(kLayer)).*(phasePSD./x).*...
                        ( factor2.*tripleBessel1 + factor3.*tripleBessel2 ).*...
                        besselsRadialOrder;
                end
                out = real(out);
            end
            function out = integrandNgs(x,zi,ni,mi,zj,nj,mj)
                mipmj = mi+mj;
                mimmj = abs(mi-mj);
                krkr_mi = mi==0;
                krkr_mj = mj==0;
                factor1 = sqrt((ni+1)*(nj+1)).*...
                    (-1).^(0.5.*(ni+nj)).*...
                    2.^(1-0.5.*(krkr_mi+krkr_mj));%.*...
                %(-1).^mj;
                factor1 = factor1./denom;
                factor2 = (-1).^(1.5*(mi+mj)).*...
                    cos( ...
                    (mi+mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) + ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                factor2 = factor2*factor1;
                factor3 = (-1).^(1.5*abs(mi-mj)).*...
                    cos( ...
                    (mi-mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) - ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                factor3 = factor3*factor1;
                red1 = a1l.*x;
                red2 = a2l.*x;
                besselsRadialOrder = besselj(ni+1,red1).*besselj(nj+1,red2);
                f = 0.5*x/pi;
                phasePSD = ...
                    psdCst.*(f.^2 + 1./atm.L0.^2).^(-11./6)./x;
                red = sl*x;
                tripleBessel1 = besselj(mipmj,red);
                tripleBessel2 = besselj(mimmj,red);
                out = (fr0*phasePSD).*...
                    ( factor2.*tripleBessel1 + factor3.*tripleBessel2 );
                out = sum( bsxfun( @times , out , besselsRadialOrder ) );
                out = real(out);
            end
        end
                
        function out = zern_aiaj(zi,ni,mi,zj,nj,mj,atm,tel)
            out = 0;
            if (ni==0) && (nj==0)
                out = Inf;
            elseif rem(abs(zi-zj),2)==0 || (mi==0 && mj==0)
                red = pi*tel.D;
                out = 2/(pi*tel.R^2)*...
                    sqrt((ni+1)*(nj+1))*(-1)^((ni+nj-mi-mj)/2)*(mi==mj).*...
                    quadgk( @(x) phaseStats.spectrum(x,atm).*...
                    besselj(ni+1,red.*x).*besselj(nj+1,red.*x)./x,0,Inf);
            end
        end
        
        function out = zern_phiai(r,o,j,n,m,atm,tel)
            Inm = @(r) quadgk( @(f) phaseStats.spectrum(f,atm).*...
                besselj(n+1,pi.*f.*tel.D).*...
                besselj(m,2.*pi.*f.*r),0,Inf);
            out = arrayfun(Inm,r);
            krkr = m~=0;
            out = out.*2.*sqrt(n+1).*2.^(0.5*krkr).*(-1).^((n-m*krkr)/2).*...
                cos(m.*o+pi.*krkr.*((-1).^j-1)./4)./tel.R;
        end
        
        
        function map = zern_resid_var(N,atm,tel,r,o)
            zern = zernike(1:N);
            
            mapFun = @(rr,oo) quadgk(@(f)integrand(f,rr,oo),0,Inf);
            tic
            map = arrayfun( mapFun , r , o );
            toc
            function out = integrand(f,rr,oo)
                cum_i = 0;
                for i=1:N
                    cum_ip = 0;
                    for ip=1:N
                        flag = (zern.m(i)==zern.m(ip)) && ...
                            ( rem( abs( zern.j(i)-zern.j(ip) ),2 )==0 || ...
                            ( (zern.m(i)==0) && (zern.m(ip)==0) ) );
                        if flag
                            cum_ip = cum_ip + ...
                                sqrt(zern.n(ip)+1).*...
                                (-1).^((zern.n(ip)-zern.m(ip))/2).*...
                                besselj(zern.n(ip)+1,pi.*tel.D.*f);
                            %{
                                zernike.fun(zern.j(ip),zern.n(ip),zern.m(ip),rr,oo).*...
                                %}

                        end
                    end
                    cum_ip = cum_ip./(pi.*f.*tel.R);
                    deltaNot = 1 - zern.m(i)==0;
                    alpha = pi.*((-1).^zern.j(i)-1).*deltaNot./4;
%                     cum_ip = cum_ip - ...
%                         2.*(-1).^(zern.m(i).*(zern.m(i)==0)/2).*...
%                         2.^(deltaNot/2).*...
%                         cos(zern.m(i).*oo + alpha).*...
%                         besselj(zern.m(i),2.*pi.*f.*rr);
                    cum_i = cum_i + ...
                        sqrt(zern.n(i)+1).*...
                        (-1).^((zern.n(i)-zern.m(i))/2).*...
                        besselj(zern.n(i)+1,pi.*tel.D.*f).*...
                        2.*cum_ip./tel.R;
                        %{
                        zernike.fun(zern.j(i),zern.n(i),zern.m(i),rr,oo).*...
                        %}
                end
%                 out = phaseStats.spectrum(f,atm).*(2.*pi.*f + cum_i);
                out = phaseStats.spectrum(f,atm).*(cum_i);
            end
            
        end
        
        function varargout = spatioAngularAngleOfArrivalCovarianceMatrix(sampling,range,D,atm,srcAC,varargin)
            %% SPATIOANGULARCOVARIANCEMATRIX Phase spatio-angular covariance meta matrix
            %
            % S = spatioAngularCovarianceMatrix(sampling,range,atm,src1)
            % computes the spatio-angular auto-correlation meta-matrix of the
            % wavefront between all the sources srcAC. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % C = spatioAngularCovarianceMatrix(sampling,range,atm,src1,src2)
            % computes the spatio-angular cross-correlation meta-matrix of the
            % wavefront between all src2 and src1. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % [S,C] = spatioAngularCovarianceMatrix(...) computes both
            % auto- and cross-correlation meta-matrix
            %
            % ... = spatioAngularCovarianceMatrix(...,'mask',mask) restrict
            % the sampling within the mask
            
            inputs = inputParser;
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('srcAC',@(x) isa(x,'source'));
            inputs.addOptional('srcCC',[],@(x) isa(x,'source'));
            inputs.addParamValue('slopesAxis',0,@isnumeric);
            inputs.addParamValue('phaseSlopesCC',false,@islogical);
            inputs.addParamValue('mask',true(sampling),@islogical);
            inputs.addParamValue('lag',0,@isnumeric);
            inputs.addParamValue('xyOutput',[],@isnumeric);
            inputs.parse(atm,srcAC,varargin{:});
            
            m_srcCC = inputs.Results.srcCC;
            m_mask  = inputs.Results.mask;
            m_tau   = inputs.Results.lag;
            xyOutput= inputs.Results.xyOutput;
            
            nCC = 0;
            if inputs.Results.phaseSlopesCC
                nCC = 1;
            end
                
            [m_x,m_y] = meshgrid( ((0:sampling-1)-(sampling-1)*0.5)*D );
            m_nGs   = length(srcAC);
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            m_cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            m_cstL0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            m_cstr0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            m_L0      = atm.L0;
            m_lambda  = atm.wavelength;
            m_cst_xy = 0.0719*m_lambda.^2.*atm.r0.^(-5/3).*m_L0^(11/3);
            m_cst_xyCC = 0.0719*m_lambda.^(2-nCC).*atm.r0.^(-5/3).*m_L0^(11/3);
            
            slopesAxis = inputs.Results.slopesAxis;

            m_nLayer   = atm.nLayer;
            layers   = atm.layer;
            m_altitude = [layers.altitude];
            m_fr0      = [layers.fractionnalR0];
            [m_windVx,m_windVy] = pol2cart([layers.windDirection],[layers.windSpeed]);
            m_srcACdirectionVector = cat(2,srcAC.directionVector);
            m_srcACheight          = [srcAC.height];
            
            if nargout==2
                
                tic
                varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cst_xy);
                toc
                
                m_xAC = m_x(m_mask);
                m_yAC = m_y(m_mask);
                if isempty(xyOutput)
                    m_xCC = m_xAC;
                    m_yCC = m_yAC;
                else
                    m_xCC = xyOutput(:,1);
                    m_yCC = xyOutput(:,2);
                end
                
                tic                
                varargout{2} = crossCorrelation(m_srcCC,m_xAC,m_yAC,m_xCC,m_yCC,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cst_xyCC,m_windVx,m_windVy,m_tau);
                toc
                
            else
                
                if isempty(m_srcCC)
                    
                    varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                        m_nGs,m_srcACdirectionVector,m_srcACheight,...
                        m_nLayer,m_altitude,m_fr0,...
                        m_L0,m_cst_xy);
                    
                else
                    
                    m_xAC = m_x(m_mask);
                    m_yAC = m_y(m_mask);
                    if isempty(xyOutput)
                        m_xCC = m_xAC;
                        m_yCC = m_yAC;
                    else
                        m_xCC = xyOutput(:,1);
                        m_yCC = xyOutput(:,2);
                    end
                    varargout{1} = crossCorrelation(m_srcCC,m_xAC,m_yAC,m_xCC,m_yCC,...
                        m_nGs,m_srcACdirectionVector,m_srcACheight,...
                        m_nLayer,m_altitude,m_fr0,...
                        m_L0,m_cst_xyCC,m_windVx,m_windVy,m_tau);
                    
                end
                
            end
            
            function C = crossCorrelation(srcCC,xAC,yAC,xCC,yCC,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cst_xy,windVx,windVy,tau)
                    
                fprintf(' -->> Cross-correlation meta-matrix!\n')
                
                nSs = length(srcCC);
                srcCCdirectionVector = cat(2,srcCC.directionVector);
                srcCCheight          = [srcCC.height];
                C = cellfun( @(x) zeros(length(xCC),length(xAC)) , cell(nSs,nGs) , 'UniformOutput' , false );
                f0 = 1/L0;
                nmax = 10;

                var_xy = cst_xy.*4*(pi*D).^-(4-nCC).*...
                    tools.oneParameterExample3(2-nCC,1,2,11/6,pi.*D*f0,nmax)*ones(length(xCC),length(xAC));

                for k=1:nSs*nGs
                    
                    [kSs,iGs] = ind2sub([nSs,nGs],k);
                    buf = 0;
                    
                    for kLayer=1:nLayer
                        
                        
                        beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(iGs);
                        iZ = complex( xAC*scale + beta(1) , yAC*scale + beta(2) );
                        
                        betaSs = srcCCdirectionVector(:,kSs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcCCheight(kSs);
                        zSs = complex( ...
                            xCC*scale + betaSs(1) + windVx(kLayer)*tau, ...
                            yCC*scale + betaSs(2) + windVy(kLayer)*tau );
                        
                        rho   = bsxfun(@minus,zSs,iZ.');
                        arg_rho = angle(rho);
                        rho   = abs(rho);
                        out   = var_xy;
                        index = rho~=0;
                        u          = 2.*pi.*rho(index)./L0;
                        out(index) = cst_xy.*(u*L0).^-(4-nCC).*...
                            ( tools.oneParameterExample2(4-nCC,0,2,11/6,u,nmax) - ...
                            cos(2*(arg_rho(index)-slopesAxis)).*...
                            tools.oneParameterExample2(4-nCC,2,2,11/6,u,nmax) );
                        
                        buf = buf + fr0(kLayer)*out;
                        
                    end
                    
                    C{k} = buf;
                    
                end
                                
                    buf = C;
                    C = cell(nSs,1);
                    for k=1:nSs
                        C{k} = cell2mat(buf(k,:));
                    end
                
            end
            
            function S = autoCorrelation(x,y,mask,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cst_xy)
                    
                    fprintf(' -->> Auto-correlation meta-matrix!!!!\n')
                    
                    kGs = reshape( triu( reshape(1:nGs^2,nGs,nGs) , 1) , 1, []);
                    kGs(1) = 1;
                    kGs(kGs==0) = [];
                    S = cellfun( @(x) zeros(sum(mask(:))) , cell(1,length(kGs)) , 'UniformOutput' , false );
%                     cstL0AC = cstL0*ones(sampling);
                    f0 = 1./m_L0;
                    nmax = 20;
                    
                    var_xy = cst_xy.*4*(pi*D).^-4.*...
                        tools.oneParameterExample3(2,1,2,11/6,pi.*D*f0,nmax)*ones(sampling);
                    
                    for k=1:length(kGs)
                        
                        [iGs,jGs] = ind2sub( [nGs,nGs] , kGs(k) );
                        buf = 0;
                        
                        for kLayer=1:nLayer
                            
                            %                         atmSlab = slab(atm,kLayer);
                            
                            beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(iGs);
                            iZ = complex( x*scale + beta(1) , y*scale + beta(2) );
                            
                            beta = srcACdirectionVector(:,jGs)*altitude(kLayer);
                            scale = 1 - altitude(kLayer)/srcACheight(jGs);
                            jZ  = complex( x*scale + beta(1) , y*scale + beta(2) );
                            
                            z1 = iZ;
                            z2 = jZ;
                            [nz,mz] = size( z1 );
                                                        
                            % First Row
                            %                         r  =  phaseStats.covariance( abs(z2-z1(1)) , atm );
                            rho = z2-z1(1);
                            arg_rho = angle(rho);
                            rho   = abs(rho);
                            r   = var_xy;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            r(index) = cst_xy.*(u.*L0).^-4.*...
                                ( tools.oneParameterExample2(4,0,2,11/6,u,nmax) ...
                                - cos(2*(arg_rho(index)-slopesAxis)).*tools.oneParameterExample2(4,2,2,11/6,u,nmax) );
                            
                            r   = mat2cell( fr0(kLayer)*r  , nz , ones(mz,1));
                            % First column in first blocks fow
                            %                         c  =  phaseStats.covariance( abs( bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' ) ) , atm );
                            rho = bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' );
                            arg_rho = angle(rho);
                            rho   = abs(rho);
                            c   = var_xy;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            c(index) = cst_xy.*(u.*L0).^-4.*...
                                ( tools.oneParameterExample2(4,0,2,11/6,u,nmax) ...
                                - cos(2*(arg_rho(index)-slopesAxis)).*tools.oneParameterExample2(4,2,2,11/6,u,nmax) );
                            
                            c   = mat2cell( fr0(kLayer)*c , nz , ones(mz,1) );
                            % First block rows
                            rr   = cellfun( @(x,y) myToeplitz(x,y) , c , r , 'UniformOutput' , false);
                            
                            % First Column
                            %                         c  =  phaseStats.covariance( abs(z1-z2(1)) , atm );
                            rho = z1-z2(1);
                            arg_rho = angle(rho);
                            rho   = abs(rho);
                            c   = var_xy;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            c(index) = cst_xy.*(u.*L0).^-4.*...
                                ( tools.oneParameterExample2(4,0,2,11/6,u,nmax) ...
                                - cos(2*(arg_rho(index)-slopesAxis)).*tools.oneParameterExample2(4,2,2,11/6,u,nmax) );
                            
                            c   = mat2cell( fr0(kLayer)*c  , nz , ones(mz,1));
                            % First row in first blocks column
                            %                         r  =  phaseStats.covariance( abs( bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' ) ) , atm );
                            rho = bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' );
                            arg_rho = angle(rho);
                            rho   = abs(rho);
                            r   = var_xy;
                            index         = rho~=0;
                            u             = 2.*pi.*rho(index)./L0;
                            r(index) = cst_xy.*(u.*L0).^-4.*...
                                ( tools.oneParameterExample2(4,0,2,11/6,u,nmax) ...
                                - cos(2*(arg_rho(index)-slopesAxis)).*tools.oneParameterExample2(4,2,2,11/6,u,nmax) );
                            
                            r   = mat2cell( fr0(kLayer)*r , nz , ones(mz,1) );
                            % First blocks column
                            cc   = cellfun( @(x,y) myToeplitz(x,y) , c , r , 'UniformOutput' , false);
                            
                            out = cell2mat( myToeplitz(cc,rr) );
                            out(~mask,:) = [];
                            out(:,~mask) = [];
                            
                            buf = buf + out;
                            
                        end
                        
                        S{k} = buf;
                        
                    end
                    
                    buf = S;
                    S = cellfun( @(x) zeros(sum(mask(:))) , cell(nGs) , 'UniformOutput' , false );
                    S(kGs) = buf;
                    S(1:nGs+1:nGs^2) = S(1,1);
                    S = cell2mat(S);
                    S = triu(S,1)+triu(S)';

            end
            
            
        end

        
    end
    
end