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
        

        
        %% Rotate the DM actuator position by rotAngle in radian
        function [pxx, pyy] = rotateDM(px, py, rotAnglInRadians)
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
            pxx = px*cos(rotAnglInRadians) - py*sin(rotAnglInRadians);
            pyy = py*cos(rotAnglInRadians) + px*sin(rotAnglInRadians);
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