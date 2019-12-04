classdef lamTools %< eltTools
    % lamTools with shared functionality :: static class
    
    methods (Static)
        
        % REPLACED BY THE ONE IN psfrTools :: ccorreia 18 April 2019
        %% Estimate compressed profile using mean-weighted compression
%         function [Cn2eq altEq] = eqLayers(Cn2, altitudes, nEqLayers, power)
%             %{
%             Cn2         ::  The input Cn2 profile (vector)
%             altitudes   ::  The input altitudes (vector)
%             nEqLayers   ::  The number of output equivalent layers (scalar)
%             power       ::  the exponent of the turbulence (default 5/3)
%             
%             See: Saxenhuber17: Comparison of methods for the reduction of
%             reconstructed layers in atmospheric tomography, App Op, Vol. 56, No. 10 / April 1 2017
%             %}
%             nCn2 = numel(Cn2);
%             nAltitudes = numel(altitudes);
%             if nargin ~= 4
%                 power = 5/3;
%             end
%             
%             % if nargin ~= 5
%             nSlab = floor(round(nCn2)/fix(nEqLayers));
%             ppp = 1;
%             posSlab =  round((linspace(0, nEqLayers-1, nEqLayers))*nSlab)+1;
%             for iii = 1:nEqLayers-1
%                 if posSlab(iii) >= posSlab(iii+1)
%                     posSlab(iii+1) = posSlab(iii)+1;
%                 end
%             end
%             posSlab = [posSlab, nAltitudes+1];
%             Cn2eq = zeros(1,nEqLayers);
%             altEq = zeros(1,nEqLayers);
%             for ii = 1:nEqLayers
%                 Cn2eq(ii) = sum(Cn2(posSlab(ii):posSlab(ii+1)-1));
%                 altEq(ii) =  (sum(altitudes(posSlab(ii):posSlab(ii+1)-1) .^ (power) .* Cn2(posSlab(ii):posSlab(ii+1)-1))/Cn2eq(ii)) .^ (1./power);
%             end
%         end
%         
        
        %%
        % PSD: Power Spectral Density
        % Calculates PSD for circular aperture (i.e. not spiders, no central obscuration).
        % Uses a window for measuring PSD of variable [input] defined over 2D grid
        % (i.e. to smooth edges of the data) and data goes to the edges of
        % the 2D grid.
        function [psd, psdCirc, kx] = powerSpecrum(input, m, n)
            % Check that phase is square
            nRes = size(input, 1);
            w = fourierReconstructor.window('hann', nRes);
            psd = (abs(fftshift(fft2(input.*w))).^2).*sum(w(:).^2);
            
            % Calculer circular PSD
            if size(psd, 1) == size(psd, 2); % Input matrix must be square.
                if nargin < 2
                    [m,~] = size(psd);
                    centerm = ceil(m/2+1); %matrix is square , so use m or n
                    centern = centerm;
                else
                    centerm = n;
                    centern = m;
                    %[m,n] = size(psd);
                end
                [psdCirc, kx] = radial(psd, centerm, centern);
            else
                psdCirc = 0;
                kx = 0;
                disp('Not a square array, cannot calculate circular PSD');
            end
        end
        
        %% distances
        function rho = dist(l,u, nPts)
            [x,y] = meshgrid(linspace(l,u,nPts));
            z = complex(x,y);
            rho  = abs(bsxfun(@minus,z(:),z(:).'));
        end
        
        %% crops the central portion of a square matrix, zero pads a square matrix to extend it
        function out = crop(input, ncrop)
            nFrames = size(input,3);
            
            if isscalar(ncrop)
                ncrop = [ncrop ncrop];
            end
            
            dim = size(input);
            out = zeros(ncrop(1), ncrop(2), nFrames);
            for iFrame = 1:nFrames
                if all(ncrop<dim)
                    deb = round((dim - ncrop) / 2 + 1);
                    fin = round((dim + ncrop) / 2);
                    out(:,:,iFrame) = input(deb(1):fin(1), deb(2):fin(2), iFrame);
                else
                    deb = round((ncrop-dim) / 2 + 1);
                    fin = round((ncrop+dim) / 2);
                    out(deb(1):fin(1), deb(1):fin(1), iFrame) = input(:,:,iFrame);
                end
            end
        end
        
        %% Centre of Gravity
        function [c1, c2] = cog(image, threshold)
            %   usage: [cx, cy] = cog(image, threshold)
            % Noah Schwartz, February 2012
            
            % Make sure image is in Double
            image = double(image);
            
            %% Check for Threshold
            if nargin == 1
            else
                ind = image < threshold;
                image(ind) = 0.0;
            end
            
            % Calculate Centre of Gravity. Matlab index is [row column]
            % So c1 is row (i.e. y) and c2 is column (i.e. x)
            [rc, cc] = ndgrid(1:size(image, 1), 1:size(image, 2));
            Mt = sum(image(:));
            c1 = sum(image(:) .* rc(:))/Mt;
            c2 = sum(image(:) .* cc(:))/Mt;
        end
        
        %% subpixel shift by known amount
        function out = shift(im,x,y)
            %             subx = x-round(x);
            %             suby = y-round(y);
            %             im = circshift(im, round([x, y]));
            %             im = conv2(im, [subx, 1-subx], 'same');
            %             out = conv2(im, [suby, 1-suby], 'same');
            [nPx, nPy] = size(im);
            [X, Y] = meshgrid(1:nPx, 1:nPy);
            out = interp2(X, Y, im, X-x, Y-y, 'cubic', 0);
            
            %             [nPx, nPy] = size(im);
            %             amp = fftshift(abs(fft2(im,nPx*2,nPy*2)));
            %             [X, Y] = meshgrid(linspace(-1, 1, nPx*2), linspace(-1, 1, nPy*2));
            %             eField = amp .* exp(1i * (-X*pi*x-Y*pi*y));
            %             out = fftshift(abs(fft2(eField)));
            %             out = out(nPx/2+1:nPx/2+nPx, nPy/2+1:nPy/2+nPy);
            %             out = out*sum(im(:))/sum(out(:));
        end
        
        %% Invert a sparse measurement noise covariance matrix to use in a tomographic reconstructor
        function iCn = invertMeasurementNoiseCovarianceMatrix(input)
            if iscell(input)
                nGs = size(input,1);
                out = zeros([size(input{1}), nGs]);
                for iGs = 1:nGs
                    out(:,:,iGs) = input{iGs};
                end
            else
                nGs = size(input,3);
                out = input;
            end
            iCn = cell(nGs,1);
            nSlopes = size(out,1);
            for iGs = 1:nGs
                Cn = out(:,:,iGs);
                
                % Extract diagonal and cross-terms -> create sparse matrix
                noiseCovarDiag = diag(Cn);
                noiseCovarDiagP1 = diag(Cn,size(Cn,1)/2);
                B = zeros(size(Cn,1),3);
                B(:,1) = noiseCovarDiag;
                B(1:1:end/2,2) = noiseCovarDiagP1;
                B(end/2+1:1:end,3) = noiseCovarDiagP1;
                CnE = spdiags(B,[0,-nSlopes/2,nSlopes/2],nSlopes,nSlopes);
                iCn{iGs} = pinv(full(CnE));
                %iCn{iGs}(abs(iCn{iGs})< 1e-10) = 0;
                
                % extract only meaningful values on the main and two
                % cross-covar diagonals
                B(:,1) = diag(iCn{iGs});
                B(1:1:end/2,2) = diag(iCn{iGs},size(Cn,1)/2);
                B(end/2+1:1:end,3) = diag(iCn{iGs},size(Cn,1)/2);
                iCn{iGs} = spdiags(B,[0,-nSlopes/2,nSlopes/2],nSlopes,nSlopes);
            end
            iCn = blkdiag(iCn{:});
        end
        
        %% Rotate the DM actuator position by rotAngle in radian
        function [pxx, pyy] = rotateDM(px, py, rotAngle)
            % function [pxx, pyy] = rotateDM(px,py, rotAngle)
            % This function rotate the DM actuator positions.
            %
            % px (pxx)  :: The old (new) X actuator positions.
            % py (pyy)  :: The odl (new) Y actuator positions.
            % rotAngle  :: The rotation angle in radian.
            %
            % Created   :: N. Schwartz, Dec 2016
            
            % EXAMPLE:
            % px = real(bifM4.actuatorCoord);
            % py = imag(bifM4.actuatorCoord);
            % [pxx, pyy] = lamTools.rotateDM(px, py, rotAngle);
            % bifM4.actuatorCoord = pxx + 1j*pyy;
            pxx = px*cos(rotAngle) - py*sin(rotAngle);
            pyy = py*cos(rotAngle) + px*sin(rotAngle);
        end
        %         %% rotate DM coordinates
        %         function [px, py] = rotateDM(px,py, theta)
        %                 res = [cos(theta), sin(theta); -sin(theta), cos(theta)] * [px;py]';
        %                 px = res(1,:)';
        %                 py = res(2,:)';
        %         end
        
        
        %% Select valid detector pixel
        function [val, cog] = selectValidDetectorPixel(frame, minLightRatio)
            % Noah Schwartz, March 2018
            % ALGORITHM:
            %    - Estimates the shift of the 4 PYR pupils.
            %    - Shifts all the pupils back to the centre.
            %    - Add the 4 pupils together.
            %    - Calulate validDetectorPixels based on minLightRatio.
            % KNOW ISSUES:
            %    - It assumes that the centre is the best position, but it
            %    may be the case that a common centre, that is not the
            %    theoritical center, can be found.
            %    - Doesn't deal with binning at the moment...
            
            % Split the detector in 4
            %             frame = utilities.binning(frame, size(wfs.frameCalibration)/wfs.binning);
            no2 = size(frame,1)/2;
            val = mat2cell(frame, [no2 no2], [no2,no2]);
            frame = mat2cell(frame, [no2 no2], [no2,no2]);
            totalIntensity4QBefore = frame{1} + frame{2} + frame{3} + frame{4};
            
            for ii=1:4
                val{ii} = frame{ii} > minLightRatio*max(max(frame{ii}));
            end
            
            % Sum the pupils and select new valid ones
            val = val{1} + val{2} + val{3} + val{4};
            
            %             % Intersection
            %             val(val <4) = 0; val(val >=1) = 1;
            %             val = logical(val);
            
            % Union
            val(val <1) = 0; val(val >=1) = 1;
            val = logical(val);
            
            %             val1 = val; val1(val1 <1) = 0; val1(val1 >=1) = 1;
            %             val2=  val; val2(val2 <2) = 0; val2(val2 >=1) = 1;
            %             val3 = val; val3(val3 <3) = 0; val3(val3 >=1) = 1;
            %             val4 = val; val4(val4 <4) = 0; val4(val4 >=1) = 1;
            
            % Calculate pupil CoG and shift back to center
            cent = (no2+1)/2;
            cog = zeros(4,2);
            for ii=1:4
                [cog(ii,1), cog(ii,2)] = NScTools.cog(frame{ii});
                cog(ii,1) = cog(ii,1) - cent;
                cog(ii,2) = cog(ii,2) - cent;
                frame{ii} = lamTools.shift(frame{ii}, cog(ii,1), cog(ii,2));
            end
            
            %             % Sum the 4 shifted pupils and select valid pixels
            %             totalIntensity4Q = frame{1} + frame{2} + frame{3} + frame{4};
            %             val = totalIntensity4Q > minLightRatio*max(totalIntensity4Q(:));
            %             valBefore = totalIntensity4QBefore > minLightRatio*max(totalIntensity4QBefore(:));
        end
        
        
        %% saveIFCube
        function saveIFCube(dmModes, pupil, directory, fileSuffix)
            % function saveIFCube(dmModes, pupil, directory, fileSuffix)
            % This function saves the influence functions (IFs) in a data
            % cube, and saves the pupil. This function is useful to then
            % create the Karhunen-Loeve (see IDL function).
            %
            % dmModes     :: Data cube of NxMxk where k is the number of images/IFs/Modes.
            % pupil       :: The telescope pupil.
            % directory   :: The directory of where the data is saved.
            % fileSuffix  :: The 2 saved files will have this suffix.
            %
            %
            % Created     :: N. Schwartz, Dec 2016
            
            % EXAMPLE:
            %    dmModes = dm.modes.modes;
            %    pupil = tel.pupil;
            %    directory = '/result/SCAO-H/INPUT/KarhunenLoeve';
            %    fileSuffix = strcat(num2str(tel.resolution), 'pix_', num2str(size(dmModes, 2)),'modes')
            %    lamTools.saveIFCube(dmModes, pupil, directory, fileSuffix);
            
            % Reshape data into 3D data cube
            m = size(dmModes, 2);
            nRes = size(pupil,1);
            tmp = dmModes;
            if size(dmModes,3) ==1; tmp = (reshape(full(dmModes),nRes, nRes, m)); end;
            
            % Save IF cube
            filename = strcat(directory, '/IF_', fileSuffix,'.fits');
            fitswrite(tmp, filename);
            
            % Save pupil
            filename = strcat(directory, '/pupilMask_', fileSuffix,'.fits');
            fitswrite(double(pupil), filename);
        end %End of saveIFCube
        
        
        
        %% applyPetals
        function applyPetals(dm,actuatorCoord,petals,rotAngle)
            idx = lamTools.whichPetal(actuatorCoord, rotAngle,4);    % Find where the actuators belong (i.e. which petal/segment)
            for kActuator=1:dm.nActuator        % Apply petals to actuators influence functions
                dm.modes.modes(:,kActuator) = dm.modes.modes(:,kActuator).*petals(:,idx(kActuator));
            end
        end
        
        %% wraptopi
        function a = wraptopi(a, a_center )
            % function a = wrap(a,a_center)
            %
            % Wraps angles to a range of 2*pi.
            % Inverse of Matlab's "unwrap", and better than wrapToPi ( which has
            % redundant [-pi,pi])
            % Optional input "a_center" defines the center angle.  Default is 0, giving
            % angles from (-pi,pi], chosen to match angle(complex(-1,0)).  Maximum
            % possible value is pi.
            
            % T.Hilmer, UH
            % 2010.10.18 version 2
            %   removed code from version 1. Have not bug-checked second input
            %   "a_center"
            
            if nargin < 2, a_center = 0; end
            
            % new way
            a = mod(a,2*pi); % [0 2pi)
            
            % shift
            j = a > pi - a_center;
            a(j) = a(j) - 2*pi;
            j = a < a_center - pi;
            a(j) = a(j) + 2*pi;
        end
        
        %% multipleSpectra
        function [out, outLayered] = multiZernikeTemporalSpectra(nu,atm,zern,tel,beta,T)
            nMode = zern.nMode;
            out = zeros(nMode, length(nu));
            outLayered = zeros(nMode, length(nu),atm.nLayer);
            if nargin < 5
                beta = [];
                T = [];
            end
            if atm.nLayer > 1
                parfor kmode=zern.j
                    [out(kmode,:),outLayered(kmode,:,:)]  = lamTools.zernikeTemporalSpectrum(nu,atm,tel,kmode,beta,T);
                end
                out = out(zern.j,:); % keep only the requested modes
                outLayered = outLayered(zern.j,:,:);
            else
                parfor kmode=zern.j
                    [out(kmode,:)]  = lamTools.zernikeTemporalSpectrum(nu,atm,tel,kmode,beta, T);
                end
                out = out(zern.j,:); % keep only the requested modes
            end
        end
                %%
        function [out, outLayered] = zernikeTemporalSpectrum(nu,atm,tel,kmode,beta,T)
            %% zernikeTemporalSpectrummcomputes the Phase power spectrum density
            % two-sided instead of one-sided
            %
            % IT NOW PROVIDES THE SAME RESULTS OF
            % zernikeStats.temporalSpectrum after a fix was done on that
            % function to avoid adding negative spectra values and therefore bias the computation
            % C. Correia May 2019
            %
            
            outLayered = zeros(size(nu,2),atm.nLayer);
            for kLayer = 1:atm.nLayer
                nPx = size(atm.layer(kLayer).phase,1);
                pupil = logical(utilities.piston(nPx));
                D = tel.diameterAt(atm.layer(kLayer).altitude);
                
                zern = zernike(kmode,'resolution',nPx,'pupil',pupil,'D',D);
                
                atmSlab = slab(atm,kLayer);
                [vx,vy] = pol2cart(atmSlab.layer.windDirection,atmSlab.layer.windSpeed);
                if nargin > 4 && ~isempty(beta)%boiling atmosphere
                    for k=1:numel(nu)
                        outLayered(k,kLayer) = outLayered(k,kLayer) + integral2( @integrandFxFy , -Inf, Inf, -Inf, Inf);
                    end
                else
                    for k=1:numel(nu)
                        if abs(vx)>eps(atmSlab.layer.windSpeed)
                            outLayered(k,kLayer) = outLayered(k,kLayer) + quadgk( @integrandFy , -Inf, Inf);
                        else
                            outLayered(k,kLayer) = outLayered(k,kLayer) + quadgk( @integrandFx , -Inf, Inf);
                        end
                    end
                end
                %                 if vx == 0;
                %                     signvx = 1;
                %                 else
                %                     signvx = sign(vx);
                %                 end
                %                 if vy == 0;
                %                     signvy = 1;
                %                 else
                %                     signvy = sign(vy);
                %                 end
                %                 out(:,kLayer) = out(:,kLayer)*signvx*signvy;
                %%%%
                % --- apply wind direction using Conan95, Eq. 30
            end
            outLayered = 2*abs(outLayered); % x2 to consider two-sided instead of one-sided
            
            if atm.nLayer >1
                out = sum(outLayered,2);
            else
                out = outLayered;
            end
            
            
            function int = integrandFy(fy)
                fx = (nu(k) -fy*vy)/vx;
                %fx = (nu(k))/vx;
                int = zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vx;
            end
            
            function int = integrandFx(fx)
                fy = (nu(k) -fx*vx)/vy;
                int = zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vy;
            end
            
            function int = integrandFxFy(fx,fy)
                %factor = 2*log(beta)./(T*(4*pi^2*(nu(k) - fx*vx-fy*vy)+log(beta)^2/T^2));
                %int = factor.*zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern);
                %fy = (nu(k) -fx*vx)/vy;
                factor = 2*log(beta)./(T*(4*pi^2*(nu(k) -fx*vy)+log(beta)^2/T^2));
                int = factor.*zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vy;
            end
        end
        %% type I OMGI optimiser
        function g = omgi(nu, PSD, samplingFreq, pureDelay, noiseVar,verbose)
            if ~exist('verbose','var')
                verbose = 0;
            end
            T = 1/samplingFreq;
            s = @(x) 2*1i*pi*x;
            
            if size(PSD,1) > size(PSD,2)
                PSD = PSD';
            end
            
            % TFs
            hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T);
            hDac = hWfs; % Different meaning, same value
            hLag = @(x) exp(-pureDelay*s(x));
            hInt = @(x) 1./(1-exp(-s(x)*T));
            hDm  = @(x) 1;
            G = @(x,g) hWfs(x) .* hDac(x) .* hLag(x) .* g .* hInt(x) .* hDm(x);
            Gn = @(x,g) hDac(x) .* hLag(x) .* g .* hInt(x) .* hDm(x);
            % rejection transfer function
            E = @(x,g) abs(1./(1+G(x,g)));
            %Noise transfer
            hN = @(x,g) Gn(x,g) .* E(x,g);%./hWfs(x);
            %myfun = @(g) sum(hN(nu)*noiseVar + RTF*PSD);
            g = fminsearch(@(g)lamTools.integrate(g,nu,E,hN,PSD,noiseVar,T), 0.1);
            
            if verbose
                gains = linspace(0.001,0.5,100);
                for kGain = 1:length(gains)
                    outN(kGain) = trapz(nu, abs(hN(nu,gains(kGain))).^2*noiseVar*2*T);
                    outS(kGain) = trapz(nu, abs(E(nu,gains(kGain))).^2.*PSD);
                end
                figure
                semilogy(gains, outS,'r')
                hold on
                semilogy(gains, outN,'b')
                semilogy(gains, outS + outN,'k')
                xlabel('gain')
                ylabel('residual')
                title('single integrator optimal modal gian')
                legend('signal residual','noise residual','total residual')
                outNopt = trapz(nu, abs(hN(nu,g)).^2*noiseVar*2*T);
                outSopt = trapz(nu, abs(E(nu,g)).^2.*PSD);
                plot(g,outNopt + outSopt,'ro')
            end
        end
        
        function out = integrate(g,nu,E,hN, PSD, noiseVar,T)
            outN = trapz(nu, abs(hN(nu,g)).^2*noiseVar*2*T);
            outS = trapz(nu, abs(E(nu,g)).^2.*PSD);
            out = outS + outN;
        end
        
        %%
        function out = doubleIntParamOptim(nu, PSD, samplingFreq, pureDelay, varNoise)
            %% DOUBLE INTEGRATOR PARAMETER OPTIMISATION
            
            % out = doubleIntParamOptim(nu, PSD, samplingFreq, pureDelay,
            % varNoise) computes the double integrator gain and lead-filter
            % parameters for a 45 phase margin stability.
            % nu            :: the temporal frequency vector
            % PSD           :: the temporal PSD vector in units^2
            % samplingFreq  :: the loop sampling frequency in Hz
            % pureDelay     :: the loop pure delay
            % varNoise      :: noise variance in units^2
            % Created       :: C. Correia, Dec'15
            % Comment: based on do_optim_typeII from HIA and TMT,
            % developed originally by JPVeran
            
            T = 1/samplingFreq;
            s = @(x) 2*1i*pi*x;
            
            g0 = 1e-3;
            varTotalOld = inf;
            
            % OPEN-LOOP TRANSFER FUNCTION
            
            %             G = @(x) ((1-exp(-s(x)*T))./(s(x)*T)).^2.*...   % hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T); hDac = hWfs; % Different meaning, same value
            %                 exp(-tau*s(x)).*...                         % hLag = @(x) exp(-tau*s(x));
            %                 (1./(1-exp(-s(x)*T)).^2).*...               % hInt = @(x) 1./(1-exp(-s(x)*T)); squared for the double integrator
            %                 1;                                          % hDm  = @(x) 1;
            
            hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T);
            hDac = hWfs; % Different meaning, same value
            hLag = @(x) exp(-pureDelay*s(x));
            hInt = @(x) 1./(1-exp(-s(x)*T));
            hDm  = @(x) 1;
            G = @(x) hWfs(x) .* hDac(x) .* hLag(x) .* hInt(x) .* hInt(x) .* hDm(x);
            gcc = [];
            
            while 1
                % Open loop gain without phase lead
                %hOl = @(x) g0*G(x);
                
                [phaseMargin, fcross] = calcPhiMargin(g0*G(nu), nu);
                
                % Phase lead needed
                additionalPhaseNeeded=pi/4-phaseMargin;
                
                % calculating the lead filter parameters to achieve the required phase lead
                a=(1-sin(additionalPhaseNeeded))/(1+sin(additionalPhaseNeeded));
                f0=fcross*sqrt(a); %fprintf('fs=%g, fcross=%f\n',fs,fcross);
                Tlead=1/(2*pi*f0);
                
                % gain 1/g created by phase lead filter
                g=sqrt(a);
                
                % complete Hol. g is here to adjust g0 to remove the scaling caused by lead filter.
                %hOl1 = @(x) g * hOl(x) .* (1+Tlead*s(x))./(1+a*Tlead*s(x));
                hOl1 = @(x) g * g0*G(x) .* (1+Tlead*s(x))./(1+a*Tlead*s(x));
                % Rejection transfer function
                E = @(x) abs(1./(1+hOl1(x)));
                
                % Closed-loop transfer function
                %Ecl = @(x) hOl1(x) .* E(x);
                %Noise transfer
                hN = @(x) hOl1(x) .* E(x);%./hWfs(x);
                
                % RESIDUAL SIGNAL VARIANCE
                %rms=trapz(nus, PSD.*abs(Hrej.*Hrej));
                %varSignal = 2*quadgk( @(nu) PSD.*abs(E(nu)).^2 , 0 , Inf);
                varSignal = trapz(nu, PSD.*abs(E(nu)).^2);
                % NOISE GAIN
                %noiseGain=(trapz(nus, abs(Hn.*Hn))./trapz(nus, ones(length(nus),1)));
                %noiseGain = 2*quadgk( @(nu) PSD.*abs(hN(nu)).^2 , 0 , Inf)
                noiseGain = trapz(nu, abs(hN(nu)).^2)/(1/2/T);
                %Tot Error.
                varPropNoise=noiseGain*varNoise;
                varTotal=varSignal+varPropNoise;
                
                if varTotal < varTotalOld
                    %increase the gain.
                    g0=g0*1.01;
                    gcc = [gcc g0];
                    varTotalOld = varTotal;
                    out{1} = g0*g;
                    out{2} = a;
                    out{3} = Tlead;
                    out{4} = varSignal;
                    out{5} = varPropNoise;
                    out{6} = noiseGain;
                    out{7} = G;
                    out{8} = hN;
                    
                else
                    break
                end
            end
            function [margin, fcross]=calcPhiMargin(Hol, nus)
                %Computes Phase margin of the open loop transfer function.
                %Defined as the phase of Hol when abs(Hol) crosses 1.
                ind=abs(Hol)>1;
                indl=find(ind);
                indl=indl(end);
                if indl==length(Hol)
                    ph2=angle(Hol(end));
                    fcross=nus(end);
                else
                    abs1=log10(abs(Hol(indl:indl+1)));%interpolate the logrithm which is close to linear.
                    frac1=abs1(1)/(abs1(1)-abs1(2));%linear interpolation.
                    ph1=(angle(Hol(indl:indl+1)));
                    %Reduce distance
                    diff=(ph1(2)-ph1(1))/(2*pi);
                    diff=diff-fix(diff);
                    ph1(2)=ph1(1)+2*pi*diff;
                    
                    ph2=ph1(1)-(ph1(1)-ph1(2))*frac1;
                    
                    nu1=nus(indl:indl+1);
                    fcross=nu1(1)-(nu1(1)-nu1(2))*frac1;
                end
                %bring [-pi;pi] to [-2pi;0];
                normph=ph2/(2*pi);
                normph=normph-floor(normph)-1;
                ph2=normph*2*pi;
                if ph2>0
                    ph2=ph2-2*pi;
                end
                margin=pi+ph2;
            end
            
        end
        
        
        %%
        function plotDMWFSgrid(tel,dm,wfs)
            
            figure;
            
            % plot dm actuator
            plot(real(dm.modes.actuatorCoord),imag(dm.modes.actuatorCoord),'.k');
            
            %
            n = size(wfs.validLenslet,1);
            [sx,sy] = meshgrid(linspace(-1,1,n)*(tel.D/2-tel.D/n/2));
            hold on;
            plot(sx(wfs.validLenslet),sy(wfs.validLenslet),'or','MarkerSize',1);
            plot(sx(~wfs.validLenslet),sy(~wfs.validLenslet),'ob','MarkerSize',1);
            hold off;
            
        end
    end
end