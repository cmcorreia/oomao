classdef pyramid < handle
    % PYRAMID Creates a pyramid WFS object
    %
    % pwfs = pyramid(nLenslet,nPx) creates a pyramid WFS object equivalent
    % to a (nLenslet X nLenslet) Shack-Hartmann WFS with an incoming
    % wavefront with a (nPx X nPx) resolution
    %
    % pwfs = pyramid(...,'modulation',m) modulates the pyramid with an
    % amplitude m lambda/D
    %
    % pwfs = pyramid(...,'binning',n) bins the detector by factor n before
    % computing the slopes
    
    % pwfs = pyramid(...,'obstructionRatio',r) creates a
    % validDetectorPixels with a central obstruction with ratio r
    
    % pwfs = pyramid(...,'src',source,'tel','telescope,'minLightRatio',
    % minLightRatio) creates an object with validDetectorPixels defined
    % over pixels with intensity higher than minLightRatio of the max, the
    % latter computed for a very large mudulation (default 25 \lambda/D)
    properties
        
        nLenslet                    % # of lenslet
        resolution;                 % resolution of the single quadrant detector images
        camera;                     % the detector used at the end of the chain
        referenceSlopesMap;         % the slopes map of reference
        rooftop;                    % pixels for the rooftop ([nx ny]) in pixels
        framePixelThreshold = -inf; % frame pixel threshold
        obstructionRatio;           % central obstrction ratio [0...1]
        slopesDisplayHandle;        % slopes display handle
        slopesListener;             % slope listener
        intensityDisplayHandle      % intensity display handle
        pupil;                      % User defined pupil for calibration purposes. (empty if we use the default)
        isInternalCalibrationPWFS;  % flag to indicate whether this is a temporary object for internal calibration only
        src;                        % source object for the calibration
        tel;                        % telescope object for the calibration
        object;                     % concatenation of sources sampling a user-defined object
        modulationCalibPwfs;        % modulation in lambda/D units for the internal pyramid WFS created only for choosing the valid pixels by applying a threshold (default 0.9) over the max intensity detected
        minLightRatio;              % the minimum amount of light per lenslet defined as a ratio to the maximum when a "large" modulation value is set
        tag = 'PYRAMID';            % tag
        slopesUnits=1;              % normalisation factor to have the 1-to-1 gain btw input and output
        wvlRange=0;                 % vector of wavelenghts for broadband WFSensing in absolute or relative units
        viewFocalPlane;             % Flag to view the psf at the focal plane of the Pyramid (will increase comutation time)
        flatField = 0;              % camera flat field
        pixelGains = 1;             % camera pixel gains
        extendedObject = 0;         % flag for extended object
        fftPhasor = 1;              % pre-computation of the fftPhasor for the extended object case
        useLinearIntensityLightMap=0;% use the physical-optics pyramidTransform model or instead a linear approximation (cf O. Fauvarque)
        graphicalDisplay = 0;       % sets the display ON or OFF (useful for cluster execution)
    end
    properties(GetAccess = 'public', SetAccess = 'private')
        validActuator;              % map of validActuators [default at corners of sub-apertures]
        lightMap;                   % map of the light passed through the pyramid
        slopesMap;                  % the slopes as read by processing the sensor
        wave;                       % incoming wave, a square complex matrix
        validSlopes;                % logical index indicating the validSlopes
        validDetectorPixels;        % logical index indicating the validDetectorPixels
        validLenslet;               % logical index corresponding to valid pixel points but
        % over dimensions of nLensletxnLenslet
        % (for compatability with shackHartmann
        % functions)
        pyrMask;                    % pyramid face transmitance and phase
        fpIm = 0;                   % image of the focal plane at the tip of the pyramid
        frameCalibration            % detector pixel intensity for the calibration WFS with larger modulation
    end
    
    properties(SetAccess=private,SetObservable=true)
        slopes;                     % the slopes as read by processing the sensor
    end
    
    properties (Dependent)
        modulation;                 % amplitude of the modulation in \lambda/D units, a decimal number
        alpha;                      % angle of incoming rays
        c;                          % Multiplicative term of the Nyquist sampling (which is 2 by default). With  c = 2 (default) the quadrants are spread by a factor 2. End result :: c=2 -> 4pixels across \lambda/D
        binning;                    % binning factor, default 1
    end
    
    properties (Dependent,SetAccess=private)
        nSlope;                     % # of slopes
    end
    
    properties (Access=private)
        
        p_modulation;
        p_alpha;
        p_c;
        p_binning = 1;              % binning factor, default 1
        isInitialized;
        u;
        pxSide;
        o;
        r;
        nTheta;
        q;
        I4Q4;                       % intensity map on the 4 pyramid quadrants before detector binning
        log;
    end
    
    %%
    methods
        %% Constructor
        function pwfs = pyramid(nLenslet,resolution,varargin)
            p = inputParser;
            addRequired(p,'nLenslet',@isnumeric);
            addRequired(p,'resolution',@isnumeric);
            addParameter(p,'modulation',0,@isnumeric);
            addParameter(p,'rooftop',[0 0],@isnumeric);
            addParameter(p,'binning',1,@isnumeric);
            addParameter(p,'graphicalDisplay',1,@isnumeric);
            addParameter(p,'c',2,@isnumeric);
            addParameter(p,'alpha',pi/2,@isnumeric);
            addParameter(p,'obstructionRatio',0,@isnumeric);
            addParameter(p,'pupil',[],@isnumeric);
            addParameter(p,'isInternalCalibrationPWFS',0,@isnumeric);
            addParameter(p,'minLightRatio',0.9,@isnumeric);
            %addParameter(p,'src',0,@isnumeric);
            addOptional(p,'tel', telescope(-1), @(x) isa(x,'telescopeAbstract') );
            addParameter(p,'src', source, @(x) isa(x,'source') );
            addParameter(p,'modulationCalibPwfs',50,@isnumeric);
            addParameter(p,'wvlRange',0,@isnumeric);
            addParameter(p,'viewFocalPlane',0)
            addParameter(p,'object',source, @(x) isa(x,'source') );
            addParameter(p,'extendedObject',0,@isnumeric);
            addParameter(p,'useLinearIntensityLightMap',0,@isnumeric);
            parse(p,nLenslet,resolution,varargin{:})
            
            pwfs.isInternalCalibrationPWFS= p.Results.isInternalCalibrationPWFS;
            pwfs.tel                    = p.Results.tel;
            pwfs.src                    = p.Results.src;
            pwfs.object                 = p.Results.object;
            pwfs.modulationCalibPwfs    = p.Results.modulationCalibPwfs;
            pwfs.minLightRatio          = p.Results.minLightRatio;
            pwfs.nLenslet               = p.Results.nLenslet;
            pwfs.rooftop                = p.Results.rooftop;
            pwfs.resolution             = p.Results.resolution;
            pwfs.p_alpha                = p.Results.alpha;
            pwfs.obstructionRatio       = p.Results.obstructionRatio;
            pwfs.p_binning              = p.Results.binning;
            pwfs.graphicalDisplay       = p.Results.graphicalDisplay;
            pwfs.c                      = p.Results.c; %triggers the setting of the validDetectorPixelss and validSlopes
            pwfs.modulation             = p.Results.modulation;
            pwfs.camera                 = detector(2*pwfs.c*nLenslet);
            pwfs.binning                = p.Results.binning;
            pwfs.validSlopes            = [flip(flip(pwfs.validDetectorPixels,1),2) flip(flip(pwfs.validDetectorPixels,1),2)];
            pwfs.pupil                  = p.Results.pupil;
            pwfs.wvlRange               = p.Results.wvlRange;
            pwfs.viewFocalPlane         = p.Results.viewFocalPlane;
            pwfs.extendedObject         = p.Results.extendedObject;
            pwfs.useLinearIntensityLightMap = p.Results.useLinearIntensityLightMap;
            
            
            if pwfs.extendedObject % CHECK for object size versus FOV on the pyramid
                sizeObjectMax = max([pwfs.object.zenith]) * 2; % maximum size of object in radians on sky
                fieldOfView   = size(pwfs.pyrMask, 1) * pwfs.src.wavelength / pwfs.tel.D / 2.; % FOV in the pyramid mask focal plane, assumed shannon sampling
                if fieldOfView < sizeObjectMax
                    disp('###########################################################################')
                    disp('WARNING - EXTENDED OBJECT SIZE LARGER THAN FOV - CASE NOT HANDLED - STOPPING')
                    disp('###########################################################################')
                    keyboard
                end
            end
            
            
            if isempty(pwfs.pupil) && ~isempty(pwfs.tel)
                pwfs.pupil = pwfs.tel.pupil;
            end
            
            pwfs.camera.frameGrabber    = pwfs; %
            
            pwfs.referenceSlopesMap = zeros(pwfs.nLenslet*pwfs.c/pwfs.binning,pwfs.nLenslet*2*pwfs.c/pwfs.binning);
            
            pwfs.slopesListener     = addlistener(pwfs,'slopes','PostSet',...
                @(src,evnt) pwfs.slopesDisplay );
            pwfs.slopesListener.Enabled = false;
            pwfs.log                = logBook.checkIn(pwfs);
        end
        
        %% Destructor
        function delete(pwfs)
            if isvalid(pwfs.camera)
                delete(pwfs.camera)
            end
            checkout(pwfs.log,pwfs);
        end
        
        %% SETS AND GETS
        %% Get nSlope
        function out = get.nSlope(pwfs)
            out = sum(pwfs.validSlopes(:));
        end
        
        %% Get and Set modulation
        function out = get.modulation(pwfs)
            out = pwfs.p_modulation;
        end
        function set.modulation(pwfs,val)
            pwfs.p_modulation = val;
            [uu,vv] = ndgrid((0:(pwfs.pxSide-1))./pwfs.pxSide);
            [pwfs.o,pwfs.r] = cart2pol(uu,vv);
            pwfs.nTheta = round(2*pi*pwfs.c*pwfs.modulation);
            setModulationData(pwfs,1)
        end
        
        function setModulationData(pwfs,nWave,nTheta)
            if ~exist('nTheta','var')
                nTheta = pwfs.nTheta;
            end
            pwfs.q      = zeros(pwfs.pxSide,pwfs.pxSide,nWave);
            pwfs.I4Q4   = zeros(pwfs.pxSide,pwfs.pxSide,...
                nWave,nTheta); %
            if size(pwfs.object,2) > 1 && size(pwfs.fftPhasor, 3) ~= size(pwfs.object,2)  % user-defined object
                makeFftPhasor(pwfs);
            end
        end
        %% Get and Set validDetectorPixels
        function setvalidDetectorPixels(pwfs,val)
            if nargin == 2 % user custom value
                pwfs.validDetectorPixels = logical(val);
                pwfs.validSlopes = [pwfs.validDetectorPixels pwfs.validDetectorPixels];
            elseif nargin == 1 && pwfs.tel.D < 0
                pwfs.validDetectorPixels = utilities.piston((pwfs.nLenslet)/pwfs.p_binning,...
                    pwfs.nLenslet*pwfs.c/pwfs.p_binning,'type','logical','xOffset',0,'yOffset',0);
                if pwfs.obstructionRatio
                    centralHole = utilities.piston(floor(pwfs.obstructionRatio*(pwfs.nLenslet)/pwfs.p_binning),...
                        pwfs.nLenslet*pwfs.c/pwfs.p_binning,'type','logical','xOffset',0,'yOffset',0);
                else
                    centralHole = 0;
                end
                pwfs.validDetectorPixels = logical(pwfs.validDetectorPixels - centralHole);
            elseif ~pwfs.isInternalCalibrationPWFS && pwfs.tel.D > 0    % if a telescope and a source are defined a criterion with the
                                                                        % minLightRation is applied
                pwfs_ = pyramid(pwfs.nLenslet,pwfs.resolution,'modulation',pwfs.modulationCalibPwfs,'alpha',pwfs.alpha,'c',pwfs.c,...
                    'isInternalCalibrationPWFS',1,...
                    'obstructionRatio',pwfs.obstructionRatio,...
                    'binning',pwfs.p_binning,...
                    'graphicalDisplay',pwfs.graphicalDisplay);
                pwfs.src = pwfs.src .* pwfs.tel * pwfs_; % propagate through
                
                if pwfs_.graphicalDisplay == 1 
                    figure
                    subplot(1,3,1)
                    imagesc(pwfs_.camera)
                end

                frame   = pwfs_.camera.frame;
                frame   = utilities.binning(frame,size(pwfs_.camera.frame)/pwfs_.binning);
                no2     = size(frame,1)/2;
                frame   = mat2cell(frame, [no2 no2], [no2,no2]);
                totalIntensity4Q = frame{1} + frame{2} + frame{3} + frame{4};
                
                if pwfs_.graphicalDisplay == 1 
                    subplot(1,3,2)
                    imagesc(pwfs.validDetectorPixels)
                    axis square
                end
                
                pwfs.validDetectorPixels = totalIntensity4Q>pwfs.minLightRatio*max(totalIntensity4Q(:));

                if pwfs_.graphicalDisplay == 1
                    subplot(1,3,3)
                    imagesc(pwfs.validDetectorPixels)
                    axis square
                end
                
                pwfs.frameCalibration = pwfs_.camera.frame;
                
                %idx = 1:length(pwfs.validDetectorPixels);
                %diff = length(idx)-pwfs.nLenslet;
                %idx = idx(diff/2+1:end-diff/2);
                %pwfs.setValidLenslet(pwfs.validDetectorPixels(idx,idx));
            end
        end
        
        %% Get and Set validLenslet
        function setValidLenslet(pwfs,val,nL)
            
            % Sets indices validLenslet for compatibility with other
            % functions for shackHartmann etc.  This gives the valid
            % lenslet/pixel pupil over a grid of nLenslet X nLenslet
            % dimensions.
            if nargin == 3 %select only the nLenslet central pixels
                idx = 1:length(pwfs.validDetectorPixels);
                diff = length(idx)-nL/pwfs.binning;
                idx = idx(diff/2+1:end-diff/2);
                pwfs.validLenslet = pwfs.validDetectorPixels(idx,idx);
            elseif nargin==2
                
                pwfs.validLenslet = logical(val);
                
            else % Modification : -1 and +1 to remove edge lenslets
                pwfs.validLenslet = utilities.piston(pwfs.nLenslet/pwfs.p_binning,...
                    'type','logical','xOffset',0,'yOffset',0);
                if pwfs.obstructionRatio
                    centralHole = utilities.piston(floor(pwfs.obstructionRatio*pwfs.nLenslet/pwfs.p_binning),...
                        pwfs.nLenslet/pwfs.p_binning,'type','logical','xOffset',0,'yOffset',0);
                else
                    centralHole = 0;
                end
                pwfs.validLenslet = logical(pwfs.validLenslet - centralHole);
            end
            
            %             if ~pwfs.isInternalCalibrationPWFS && pwfs.tel.D > 0
            %
            %                 idx = round((pwfs.c-1)/2*pwfs.nLenslet)+1:round((pwfs.c-1)/2*pwfs.nLenslet)+pwfs.nLenslet;
            %                 pwfs.validLenslet = pwfs.validDetectorPixels(idx,idx);
            %
            %             end
            pwfs.validSlopes = [pwfs.validLenslet pwfs.validLenslet];
        end
        %% Get and Set alpha
        function out = get.alpha(pwfs)
            out = pwfs.p_alpha;
        end
        function set.alpha(pwfs,val)
            pwfs.p_alpha = val;
            makePyrMask(pwfs);
        end
        
        %% Get and Set c
        function out = get.c(pwfs)
            out = pwfs.p_c;
        end
        function set.c(pwfs,val)
            pwfs.p_c = val;
            if rem(pwfs.nLenslet*pwfs.c,1) > 0 % size in pixels is not an integer
                fprintf('WARNING: Pupil separation changed from: %f, to: %f\n',pwfs.p_c,ceil(pwfs.nLenslet*pwfs.c)/pwfs.nLenslet)
                pwfs.p_c = ceil(pwfs.nLenslet*pwfs.c)/pwfs.nLenslet;
            end
            %pwfs.u = fix(1+pwfs.resolution*(2*pwfs.c-1)/2:pwfs.resolution*(2*pwfs.c+1)/2);
            pwfs.u = fix(2+pwfs.resolution*(2*pwfs.c-1)/2:pwfs.resolution*(2*pwfs.c+1)/2+1); % added one pixel shift to this index such that the intensities from the 4 quadrants when added together are symmetic!!! (14-Nov-2017)
            
            pwfs.pxSide = pwfs.resolution*2*pwfs.p_c;
            makePyrMask(pwfs);
            setvalidDetectorPixels(pwfs);
            setValidLenslet(pwfs, pwfs.validDetectorPixels,pwfs.nLenslet);
            
        end
        
        %% Get and Set binning
        function out = get.binning(pwfs)
            out = pwfs.p_binning;
        end
        function set.binning(pwfs,val)
            pwfs.p_binning = val;
            pwfs.isInitialized = false;
            setvalidDetectorPixels(pwfs);
            pwfs.validSlopes = [pwfs.validDetectorPixels pwfs.validDetectorPixels];
            pwfs.referenceSlopesMap = zeros(pwfs.nLenslet*pwfs.c/pwfs.p_binning,pwfs.nLenslet*2*pwfs.c/pwfs.p_binning);
        end
        
        %% Get valid actuators
        function val = get.validActuator(obj)
            nElements            = 2*obj.nLenslet+1; % Linear number of lenslet+actuator
            validLensletActuator = zeros(nElements);
            index                = 2:2:nElements; % Lenslet index
            validLensletActuator(index,index) = flipud(fliplr(obj.validLenslet));
            for xLenslet = index
                for yLenslet = index
                    if validLensletActuator(xLenslet,yLenslet)==1
                        xActuatorIndice = [xLenslet-1,xLenslet-1,...
                            xLenslet+1,xLenslet+1];
                        yActuatorIndice = [yLenslet-1,yLenslet+1,...
                            yLenslet+1,yLenslet-1];
                        validLensletActuator(xActuatorIndice,yActuatorIndice) = 1;
                    end
                end
            end
            index = 1:2:nElements; % Actuator index
            val   = logical(validLensletActuator(index,index));
        end
        %         %% Get valid actuators
        %         function val = get.validActuatorOld(obj)
        %             nElements            = 2*obj.nLenslet+1; % Linear number of lenslet+actuator
        %             validLensletActuator = zeros(nElements);
        %             index                = 2:2:nElements; % Lenslet index
        %             idx0 = ceil(obj.nLenslet*(obj.c-1)/2);
        %             w = idx0+1:idx0+obj.nLenslet;
        %             validLensletActuator(index,index) = obj.validDetectorPixels(w,w);
        %             for xLenslet = index
        %                 for yLenslet = index
        %                     if validLensletActuator(xLenslet,yLenslet)==1
        %                         xActuatorIndice = [xLenslet-1,xLenslet-1,...
        %                             xLenslet+1,xLenslet+1];
        %                         yActuatorIndice = [yLenslet-1,yLenslet+1,...
        %                             yLenslet+1,yLenslet-1];
        %                         validLensletActuator(xActuatorIndice,yActuatorIndice) = 1;
        %                     end
        %                 end
        %             end
        %             index = 1:2:nElements; % Actuator index
        %             val   = logical(validLensletActuator(index,index));
        %         end
        %% Set wvlRange
        function set.wvlRange(pwfs,val)
            pwfs.wvlRange = val;
        end
        
         %% Get flux
        function out = flux(pwfs)
            n = size(pwfs.camera.frame,1)/pwfs.binning;
            frame = pwfs.camera.frame;
            frame = (frame - pwfs.flatField)./pwfs.pixelGains;
            
            if isfinite(pwfs.framePixelThreshold)
                frame = frame - pwfs.framePixelThreshold;
            end
            frame(frame<0) = 0;
            I4Q = utilities.binning(frame,size(pwfs.camera.frame)/pwfs.binning);             % binning
            I4Q = reshape(I4Q,n,n,[]);
            
            n = size(I4Q,1);
            nFrames = size(I4Q,3);
            if nFrames > 1
                im = mat2cell(I4Q,[n/2 n/2],[n/2 n/2],nFrames);
            else
                im = mat2cell(I4Q,[n/2 n/2],[n/2 n/2]);
            end
            
            I = (im{4} + im{4} + im{4} + im{4});      %
            out = sum(I(pwfs.validDetectorPixels));

        end
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(pwfs, src)% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            %% RELAY pyramid to source relay
            wvl = [src.wavelength];
            if  numel(unique(wvl)) > 1 && any(pwfs.wvlRange) == 0
                if numel(unique(wvl)) > 1 %multi wavelength case
                    pwfs.wvlRange = wvl;
                end
            end
            if pwfs.useLinearIntensityLightMap
                linearSlopes(pwfs,src.catWave);
            else
                pyramidTransform(pwfs,src.catWave);
            end
            grab(pwfs.camera,pwfs.lightMap);% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            dataProcessing(pwfs);
        end
        % wfs=pyramid(nLenslet,nPx,'modulation',modul);
        function varargout = slopesDisplay(pwfs,varargin)
            %% SLOPESDISPLAY pyramide slopes display
            %% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            % slopesDisplay(obj) displays image plot of the slopes
            %
            % slopesDisplay(obj,'PropertyName',PropertyValue) displays% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            % image plot of the slopes and set the properties of the
            % graphics object quiver
            %% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            % h = slopesDisplay(obj,...) returns the graphics handle
            n = pwfs.nLenslet*pwfs.c/pwfs.p_binning;
            data = reshape( pwfs.slopesMap , n , 2*n, []);
            data = bsxfun( @times, data, pwfs.validSlopes);
            data = reshape( data, n , []);
            if ishandle(pwfs.slopesDisplayHandle)
                set(pwfs.slopesDisplayHandle,...
                    'Cdata',data,varargin{:})
            else
                pwfs.slopesDisplayHandle = imagesc(data,varargin{:});
                ax = gca;
                pos = get(ax,'position');
                axis xy equal tight
                ylabel(colorbar('location','EastOutside'),'Pixel')
                set(ax,'position',pos)
            end
            if nargout>1
                varargout{1} = pwfs.slopesDisplayHandle;
            end
        end
        
        function varargout = intensityDisplay(pwfs,varargin)
            %% INTENSITYDISPLAY pyramid lenslet intensity display
            %
            % intensityDisplay(obj) displays the intensity of the lenslets
            %
            % intensityDisplay(obj,'PropertyName',PropertyValue) displays
            % the intensity of the lenslets and set the properties of the
            % graphics object imagesc
            %
            % h = intensityDisplay(obj,...) returns the graphics handle
            %
            % See also: imagesc
            n = size(pwfs.camera.frame,1);
            frame = pwfs.camera.frame(1:n/2,:) + pwfs.camera.frame(1+n/2:n,:);
            frame = reshape(frame',size(frame,1)*2,[]);
            n = size(frame,1);
            intensity = frame(1:n/2,:) + frame(1+n/2:n,:);
            if ishandle(pwfs.intensityDisplayHandle)
                set(pwfs.intensityDisplayHandle,...
                    'CData',intensity,...
                    varargin{:})
            else
                pwfs.intensityDisplayHandle = ...
                    imagesc(intensity,varargin{:});
                ax = gca;
                pos = get(ax,'position');
                axis equal tight xy
                colorbar('location','EastOutside')
                set(ax,'position',pos)
            end
            if nargout>0
                varargout{1} = obj.intensityDisplayHandle;
            end
        end
        
        %% INITIALISATION :: referecence slopes and gain calibration
        function INIT(pwfs)
            %% INIT pyramid inialization
            %
            % pwfs.INIT sets the reference slopes, maps and gains
            
            pwfs.referenceSlopesMap = pwfs.slopesMap + pwfs.referenceSlopesMap;
            pwfs.slopesUnits = 1;
            pwfs.slopesMap = pwfs.slopesMap*0;
            pwfs.slopes = pwfs.slopes*0;
            gainCalibration(pwfs)
            pwfs.isInitialized = true;
        end
        
        %This function updates the detector with whatever light there
        %is in the pyramid, and then procedes to process the output of
        %the detector in order to obtain the slopes
        function uplus(pwfs)
            %% UPLUS + Update operator
            %
            % +obj grabs a frame and computes the slopes
            %
            % obj = +obj returns the pyramid object
            
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
            
        end
        
        %This function computes the slopes corresponding to a flat
        %incoming light wave.tel
        function dataProcessing(pwfs)
            %% DATAPROCESSING Processing a SH-WFS detector frame
            %tel
            % dataProcessing(obj) computes the WFS slopes
            
            n = size(pwfs.camera.frame,1)/pwfs.binning;
            frame = pwfs.camera.frame;
            frame = (frame - pwfs.flatField)./pwfs.pixelGains;
            
            if isfinite(pwfs.framePixelThreshold)
                frame = frame - pwfs.framePixelThreshold;
            end
            frame(frame<0) = 0;
            I4Q = utilities.binning(frame,size(pwfs.camera.frame)/pwfs.binning);             % binning
            I4Q = reshape(I4Q,n,n,[]);
            
            n = size(I4Q,1);
            nFrames = size(I4Q,3);
            if nFrames > 1
                im = mat2cell(I4Q,[n/2 n/2],[n/2 n/2],nFrames);
            else
                im = mat2cell(I4Q,[n/2 n/2],[n/2 n/2]);
            end
            computeSlopes(pwfs, im{1}, im{2}, im{4}, im{3});% modif
        end
        
        %% gain calibration
        function gainCalibration(pwfs,dSubap,onskyUnits)
            if  nargin == 1
                dSubap = pwfs.tel.D/pwfs.nLenslet;
            end
            pwfs.slopesUnits = 1; % reset gains
            
            nPx = pwfs.resolution;
            %ngs = source('wavelength',photometry.R);
            ngs = pwfs.src;
            if ~isempty(pwfs.pupil)
                zer = zernike(3,'resolution', nPx,'pupil',pwfs.pupil);
            else
                zer = zernike(3,'resolution', nPx);
            end
            if ~isempty(pwfs.tel)
                tel_ = pwfs.tel;
            else
                tel_ = telescope(1, 'resolution',pwfs.resolution, 'obstructionRatio',pwfs.obstructionRatio);
            end
            if ~isempty(pwfs.pupil)
                tel_.pupil = pwfs.pupil;
            end
            sx = zeros(1,5);
            sy = zeros(1,5);
            for i = 1:5
                zer.c = (i-3)*0.001/(ngs.waveNumber);
                ngs = ngs.*zer*tel_*pwfs;
                sx(i) = mean(pwfs.slopes(1:end/2));
                sy(i) = mean(pwfs.slopes(end/2+1:end));
            end
            if  dSubap < 0
                Ox_in = 4*(-2:2)*0.001;%/ngs.waveNumber;
            else
                pixelScale = ngs.wavelength/(2*dSubap);
                D = dSubap*pwfs.nLenslet;
                %if D ~=tel_.D
                %    'WARMING: TELESCOPE SIZE AND PIXELSIZE ARE INCOMPATIBLE'
                %end
                % something weird as unit! Left here for backwards
                % compatibility reasons
                Ox_in = 4*(-2:2)*0.001/ngs.waveNumber/D/pixelScale;
                
                if nargin == 3
                    if onskyUnits == 1
                        %actual phase-difference at the edge of the pupil in
                        %radians
                        Ox_in = 4/D*ngs.waveNumber*[(-2:2)*0.001/ngs.waveNumber]; %
                    end
                end

            end
            Ox_out = sy;
            slopesLinCoef = polyfit(Ox_in,Ox_out,1);
            pwfs.slopesUnits = 1/slopesLinCoef(1);
            
            zer.c = 0;
            ngs = ngs.*zer*tel_*pwfs;
        end
        %% Hilbert transform approximation for the slopes
        function hilbertSlopes(pwfs,telescope,source)
            P=telescope.pupil;
            phi=angle(source.wave);
            %            Sx = ???PHx (P??) + P??Hx (P) ??? Hxy (P)Hy (P??) + Hxy (P??)Hy (P),
            %            Sy = ???PHy (P??) + P??Hy (P) ??? Hxy (P)Hx (P??) + Hxy (P??)Hx (P).
            %            Sx = -P*imag(hilbert(P*phi)) + P*phi*imag(hilbert(P)) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(transpose(P*phi))) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(transpose(P))) ;
            %            Sy = -P*imag(hilbert(transpose(P*phi))) + P*phi*imag(hilbert(transpose(P))) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(P*phi)) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(P)) ;
            %Sx = -P.*imag(hilbert(P.*phi)) + P.*phi.*imag(hilbert(P)) - imag(hilbert(transpose(imag(hilbert(P))))).*imag(hilbert(transpose(P.*phi))) + imag(hilbert(transpose(imag(hilbert(P.*phi))))).*imag(hilbert(transpose(P))) ;
            %Sy = -P*imag(hilbert(transpose(P*phi))) + P*phi*imag(hilbert(transpose(P))) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(P*phi)) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(P)) ;
            myHilbertX = @(x) imag(hilbert(x));
            myHilbertY = @(x) imag(hilbert(x')');
            myHilbertXY = @(x) imag(hilbert(imag(hilbert(x')')));
            Sx = -P.*myHilbertX(P.*phi) + P.*phi.*myHilbertX(P) - myHilbertXY(P).*myHilbertY(P.*phi) + myHilbertXY(P.*phi).*myHilbertY(P);
            Sy = -P.*myHilbertY(P.*phi) + P.*phi.*myHilbertY(P) - myHilbertXY(P).*myHilbertX(P.*phi) + myHilbertXY(P.*phi).*myHilbertX(P) ;
            pwfs.slopes=[Sx,Sy];
            
        end
        
        %% O. Fauvarque's meta-intensity linear model
        function linearIntensity = linearSlopes(pwfs,wave,object)
            
            [n1,n2,n3] = size(wave);
            nWave = n2*n3/n1;
            
            fft2c       = @(x) fft2_centre(x,1);
            ifft2c      = @(x) (fft2_centre(x,-1));
                        
            nPxOut      = size(pwfs.pyrMask,1);
            nPx         = size(pwfs.pupil,1);    
            
            if isempty(pwfs.pupil)
                disp('CANNOT CONTINUE, WFS NEEDS A PUPIL\n')
            end
            % modulation signal
            pupilExtended = padarray(pwfs.pupil,[ nPxOut-nPx, nPxOut - nPx]/2,'both')/sum(pwfs.pupil(:));
            if pwfs.modulation == 0
                pupil2 = utilities.piston(2*4*pwfs.modulation+2, nPxOut);
            else
                pupil2 = utilities.piston(2*4*pwfs.modulation+1, nPxOut);
            end
            pupil1 = utilities.piston(2*4*pwfs.modulation-1, nPxOut);
            modulationPupil = pupil2 - pupil1;
            modulationPupil = modulationPupil/sum(modulationPupil(:));
            if nargin == 2
                weightSignal = ifft2c(pupilExtended.*fft2c(modulationPupil));% /nPxOut^2;
            elseif nargin == 3
                weightSignal = ifft2c(pupilExtended.*fft2c(conv2(modulationPupil,object,'same')));% /nPxOut^2;
            end
            weightSignal = weightSignal*nnz(pupilExtended); %sum(weightSignal(:)) = 1
            
            
            % FFT of the transparency mask
            mHat = fftshift(fft2(fftshift(pwfs.pyrMask)));
            
            % impulse response functions
            %IR = 2*imag(conj(mHat).*conv2(mHat,pwfs.pupil,'same'));
            IRm = -2*imag(conj(mHat).*fft2c((pwfs.pyrMask).*conj(weightSignal))); % this conj here is needed to improve match btw slopes from linear and full physical optics model
            
            linearIntensity = zeros(nPxOut, nPxOut, nWave);
            for z = 1:nWave
                linearIntensity(:,:,z) = conv2(IRm,-angle(wave(:,:,z)),'same');
            end
            
            pwfs.lightMap = linearIntensity;
            grab(pwfs.camera,pwfs.lightMap);% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            dataProcessing(pwfs);

        end
    end
    
    methods (Access=private)
        
        
        %% Propagate the wavefront transformed by the pyram WFS to the detector
        function pwfs = pyramidTransform(pwfs,wave)
            p = gcp('nocreate');
            [n1,n2,n3] = size(wave);
            nWave = n2*n3/n1;
            [n4,n5,n6] = size(pwfs.fftPhasor);
            nPointsObject = n6;
            if size(pwfs.q,3)~=nWave
                if isempty(p)
                    setModulationData(pwfs,nWave,1)
                else
                    setModulationData(pwfs,nWave,pwfs.nTheta)
                end
            else
                if isempty(p)
                    setModulationData(pwfs,nWave,1)
                else
                    setModulationData(pwfs,nWave,pwfs.nTheta)
                end
            end
            pwfs.q(pwfs.u,pwfs.u,:) = reshape(wave, n1,[],nWave);
            pwfs.I4Q4 = zeros(size(pwfs.I4Q4));
            ratioWvl = [pwfs.wvlRange]/pwfs.wvlRange(ceil(end/2));
            
            if pwfs.modulation>0
                if isempty(p)
                    % 1/ NO GPU ACCELERATION
                    if gpuDeviceCount == 0 %|| numel(pwfs.q) > 150000000 % if the required memory is too large to use a gpu...
                        tic
                        pwfs.fpIm = 0;
                        for kTheta = 1:pwfs.nTheta
                            theta = (kTheta-1)*2*pi/pwfs.nTheta;
                            
                            fftPhasor_ = exp(-1i.*pi.*4*pwfs.modulation*...
                                pwfs.c*pwfs.r.*cos(pwfs.o+theta));
                            
                            buf = bsxfun(@times,pwfs.q,fftPhasor_);
                            buf = bsxfun(@times,fft2(buf),pwfs.pyrMask);
                            if pwfs.viewFocalPlane
                                pwfs.fpIm = pwfs.fpIm + fftshift(sum(abs(buf).^2,3));
                                %imagesc(pwfs.fpIm);
                                %drawnow
                            end
                            pwfs.I4Q4(:,:,:) = bsxfun(@plus,pwfs.I4Q4(:,:,:), abs(fft2(buf)).^2);
                            if pwfs.viewFocalPlane
                                pwfs.fpIm = pwfs.fpIm + fftshift(sum(abs(buf).^2,3));
                                imagesc(pwfs.fpIm);
                                drawnow
                                hold on
                            end
                        end
                        I4Q = pwfs.I4Q4/pwfs.nTheta;
                        
                        toc
                    else
                        % 2/ WITH GPU ACCELERATION
                        tic
                        qGpu = gpuArray(pwfs.q);
                        pyrMaskGpu = gpuArray(pwfs.pyrMask);
                        I4Q4Gpu = gpuArray(zeros(size(pwfs.I4Q4)));
                        
                        %                         For test : construction of Fresnel Phasor
                        %                         lambda= pwfs.src.wavelength;  % wavelenght [??m]
                        %                         nPts = size(pwfs.q,1);
                        %                         freqx = linspace(-10,10,nPts);% setup frequency axis [1/??m]
                        %                         w=10;    % width of slit is 2*w [??m]
                        %                         u0=zeros(nPts);    % field at z=0
                        %                         a0=zeros(nPts);    % angular spectrum at z=0
                        %                         H=zeros(nPts);    % transfer function
                        %                         az=zeros(nPts);    % angular spectrum at z=z
                        %                         uz=zeros(nPts);    % field at z=z
                        %                         tab_freqx = ones(nPts,1) * freqx;
                        %                         tab_freqy = tab_freqx';
                        %                         phasor = sqrt(1-(lambda.*tab_freqx).^2-(lambda.*tab_freqy).^2) ;
                        %                         z=0;    % distance [??m]
                        %                         fresnelFourier = gpuArray(fftshift(exp(1i * 2 * pi * (z/lambda) * phasor)));
                        %pwfs.fpIm = 0;
                        if pwfs.extendedObject == 0 % Case point-source object
                            
                            for kTheta = 1:pwfs.nTheta
                                theta = (kTheta-1)*2*pi/pwfs.nTheta;
                                if length(unique(ratioWvl)) > 1
                                    fftPhasor_ = gpuArray(zeros([size(pwfs.o) nWave]));
                                    for i = 1:nWave
                                        fftPhasor_(:,:,i) = exp(-1i.*pi.*4*pwfs.modulation/ratioWvl(i)*...
                                            pwfs.c*pwfs.r.*cos(pwfs.o+theta));
                                    end
                                else
                                    fftPhasor_ = gpuArray( exp(-1i.*pi.*4*pwfs.modulation*...
                                        pwfs.c*pwfs.r.*cos(pwfs.o+theta))); % ATTENTION :: BUG :: should not use pwfs.c here ccorreia 12/04/2017
                                end
                                buf = bsxfun(@times,qGpu,fftPhasor_);
                                if length(unique(ratioWvl)) > 1
                                    buf = bsxfun(@times,myFft2mtx(buf,size(pwfs.I4Q4,1),ratioWvl),pyrMaskGpu)/length(ratioWvl);
                                else
                                    buf = bsxfun(@times,fft2(buf),pyrMaskGpu);
                                end
                                if pwfs.viewFocalPlane
                                    pwfs.fpIm = pwfs.fpIm + fftshift(sum(abs(buf).^2,3));
                                    imagesc(pwfs.fpIm);
                                    drawnow
                                    hold on
                                end
                                %if z == 0
                                I4Q4Gpu = bsxfun(@plus,I4Q4Gpu, abs(fft2(buf)).^2);
                                %else
                                %    psi_d = (fft2(buf)); % Champ dans le plan pupille d??tecteur
                                %    % JFS tests for Fresnel Propag
                                %    psi_d_fresnel = ifft2(bsxfun(@times,fft2(psi_d), fresnelFourier)); % Champ propag??
                                %    I4Q4Gpu = bsxfun(@plus,I4Q4Gpu, abs(psi_d_fresnel).^2);
                                %end
                            end
                            
                            
                            I4Q = gather(I4Q4Gpu)/pwfs.nTheta;
                            %clear I4Q4Gpu pyrMaskGpu qGpu fftPhaser buf; % reset gpuMemory
                            toc
                            
                        else % user defined modulation signal
                            % In this case, the modulation signal is
                            % defined by wfs.object
                            % (2) for zenith position, (3) for azimuth
                            % position and (6) for flux value
                            % wfs.fftPhasor contains the TT phasor for each
                            % of the modulation points
                            for kMode = 1:nWave
                                buf = bsxfun(@times,qGpu(:,:,kMode),pwfs.fftPhasor); % multiplication by phasor => adds a tip-tilt in a pupil plane
                                
                                buf = bsxfun(@times,fft2(buf),pyrMaskGpu);
                                amplitudeSource = sqrt(reshape([pwfs.object.nPhoton] * pwfs.tel.area * pwfs.tel.samplingTime, [1 1 pwfs.object(1).nSrc]));
                                buf = bsxfun(@times, buf, amplitudeSource);
                                
                                if pwfs.viewFocalPlane % visualisation of the pyramid focal plane
                                    pwfs.fpIm = (sum(abs(buf).^2,3));
                                    figure(6666)
                                    imagesc(fftshift(pwfs.fpIm)); title('Image in the pyramid plane') ; colorbar
                                    drawnow
                                end
                                I4Q4Gpu = sum(abs(fft2(buf)).^2,3);
                                I4Q(:,:,kMode) = gather(I4Q4Gpu);
                            end
                        end
                    end
                else %parallel implementation
                    % 1/ NO GPU ACCELERATION
                    
                    if pwfs.extendedObject == 0
                        tic
                        p_nTheta = pwfs.nTheta;
                        modulation = pwfs.modulation;
                        p_r = pwfs.r;
                        p_o = pwfs.o;
                        p_c = pwfs.c;
                        p_q = pwfs.q;
                        p_pyrMask = pwfs.pyrMask;
                        parfor kTheta = 1:pwfs.nTheta
                            theta = (kTheta-1)*2*pi/p_nTheta;
                            fftPhasor_ = exp(-1i.*pi.*4*modulation*...
                                p_c*p_r.*cos(p_o+theta));
                            buf = bsxfun(@times,p_q,fftPhasor_);
                            buf = bsxfun(@times,fft2(buf),p_pyrMask);
                            p_I4Q4(:,:,:,kTheta) = abs(fft2(buf)).^2;
                        end
                        I4Q = sum(p_I4Q4,4)/pwfs.nTheta;
                        toc
                    else % user defined modulation signal, simulating an extended object
                        % In this case, the modulation signal is
                        % defined by wfs.object
                        % (2) for zenith position, (3) for azimuth
                        % position and (6) for flux value
                        % wfs.fftPhasor contains the TT phasor for each
                        % of the modulation points
                        for kMode = 1:nWave
                            buf = bsxfun(@times,pwfs.q,pwfs.fftPhasor); % multiplication by phasor => adds a tip-tilt in a pupil plane
                            pyrMask_ = pwfs.pyrMask;
                            
                            %                             parfor kObject = 1:size(pwfs.object,3)
                            %                                 buf(:,:,kObject) = bsxfun(@times,fft2(buf(:,:,kObject)),pyrMask_);
                            %                             end
                            buf = bsxfun(@times,fft2(buf),pyrMask_);
                            
                                if pwfs.viewFocalPlane % visualisation of the pyramid focal plane
                                    pwfs.fpIm = (sum(abs(buf).^2,3));
                                    figure(6666)
                                    imagesc(fftshift(pwfs.fpIm)); title('Image in the pyramid plane') ; colorbar
                                    drawnow
                            end
                            I4Q4Gpu = sum(abs(fft2(buf)).^2,3);
                            I4Q(:,:,kMode) = gather(I4Q4Gpu);
                        end
                    end
                    
                    % 2/ WITH GPU ACCELERATION
                    %                     tic
                    %                     p_nTheta = pwfs.nTheta;
                    %                     modulation = pwfs.modulation;
                    %                     p_r = pwfs.r;
                    %                     p_o = pwfs.o;
                    %                     p_c = pwfs.c;
                    %                     p_q = gpuArray(pwfs.q);
                    %                     p_pyrMask = gpuArray(pwfs.pyrMask);
                    %                     p_I4Q4Gpu = gpuArray(zeros(size(pwfs.I4Q4)));
                    %                     parfor kTheta = 1:pwfs.nTheta
                    %                         theta = (kTheta-1)*2*pi/p_nTheta;
                    %                         fftPhasor = gpuArray(exp(-1i.*pi.*4*modulation*...
                    %                             p_c*p_r.*cos(p_o+theta)));
                    %                         buf = bsxfun(@times,p_q,fftPhasor);
                    %                         buf = bsxfun(@times,fft2(buf),p_pyrMask);
                    %                         p_I4Q4Gpu(:,:,:,kTheta) = abs(fft2(buf)).^2;
                    %                     end
                    %                     I4Q = sum(gather(p_I4Q4Gpu),4);
                    %                     toc
                    
                end
            else
                
                
                I4Q = bsxfun(@times,fft2(pwfs.q),pwfs.pyrMask);
                if pwfs.viewFocalPlane
                    pwfs.fpIm = fftshift(sum(abs(I4Q).^2,3));
                    %imagesc(pwfs.fpIm);
                    %drawnow
                end
                I4Q = abs(fft2(I4Q)).^2;
            end
            pwfs.lightMap = reshape(I4Q,pwfs.pxSide,...
                pwfs.pxSide*nWave);
        end
        %% make the pyramid phase mask
        function makePyrMask(pwfs)
            % NEW CODE TO BE ADOPTED SOON IN REPLACEMENT OF THE REST OF THE
            % FUNCTION makePyrMask
            % m = mask(4,pwfs.pxSide,'alpha',pwfs.alpha);
            % pwfs.pyrMask   = m.theMask;
            if isscalar(pwfs.alpha)
                realAlpha = ones(1,4)*pwfs.alpha;
            else
                realAlpha = pwfs.alpha;
            end
            % Complex alpha for PYR surface slope error in 2 directions...
            imagAlpha = imag(realAlpha);
            realAlpha = real(realAlpha);
            
            nx = pwfs.rooftop(1);
            ny = pwfs.rooftop(2);
            
            [fx,fy] = freqspace(pwfs.pxSide,'meshgrid');
            fx = fx.*floor(pwfs.pxSide/2);
            fy = fy.*floor(pwfs.pxSide/2);
            %pym = zeros(pxSide);
            
            % pyramid face transmitance and phase for fx>=0 & fy>=0
            mask  = graduatedHeaviside(fx,nx).*graduatedHeaviside(fy,nx);
            phase = -realAlpha(1).*(fx+fy) + -imagAlpha(1).*(-fx+fy);
            pym   = mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx>=0 & fy<=0
            mask  = graduatedHeaviside(fx,ny).*graduatedHeaviside(-fy,-ny);
            phase = -realAlpha(2).*(fx-fy) + -imagAlpha(2).*(fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx<=0 & fy<=0
            mask  = graduatedHeaviside(-fx,-nx).*graduatedHeaviside(-fy,-nx);
            phase = realAlpha(3).*(fx+fy) + -imagAlpha(3).*(fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx<=0 & fy>=0
            mask  = graduatedHeaviside(-fx,-ny).*graduatedHeaviside(fy,ny);
            phase = -realAlpha(4).*(-fx+fy) + -imagAlpha(4).*(-fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            pwfs.pyrMask   = fftshift(pym./sum(abs(pym(:))));
        end
        
        %% make the pyramid fftPhasor for a user-defined extended object
        function makeFftPhasor(pwfs)
                nSrc = size(pwfs.object,2);
                %             pwfs.fftPhasor = gpuArray(zeros([size(pwfs.pyrMask) nSrc]));
                cpuFftPhasor = zeros([size(pwfs.pyrMask) nSrc]);
                wvlOverD = pwfs.object(1).wavelength/pwfs.tel.D;
                'Creation of FftPhasor in GPU'
                'Object size : '
                nSrc
                for kSrc = 1:nSrc
                    %                 pwfs.fftPhasor(:,:,kSrc) = exp(-1i.*pi.*4*pwfs.object(kSrc).zenith/wvlOverD*...
                    %                     pwfs.c*pwfs.r.*cos(pwfs.o+pwfs.object(kSrc).azimuth));
                    cpuFftPhasor(:,:,kSrc) = exp(-1i.*pi.*4*pwfs.object(kSrc).zenith/wvlOverD*...
                                        pwfs.c*pwfs.r.*cos(pwfs.o+pwfs.object(kSrc).azimuth));
                end
                pwfs.fftPhasor = cpuFftPhasor;
                'Done !'
                clear cpuFftPhasor;
        end
        %% Compute slopes from 4 intensity maps
        function computeSlopes(pwfs, I1, I2, I3, I4)
            
            % normalisation options
            %    1) normalisation pixel-wise by the intensity
            I = (I1+I2+I3+I4);      %
            
            %    2) normalisation by the integrated flux
            I2D = utilities.toggleFrame(I,2);
            I = mean(I2D(pwfs.validDetectorPixels))*ones(size(I));
            
            SyMap = (I1-I2+I4-I3)./I;
            SxMap = (I1-I4+I2-I3)./I;
            SxMap = flip(flip(SxMap,1),2);
            SyMap = flip(flip(SyMap,1),2);
            pwfs.slopesMap = bsxfun(@minus,[SxMap,SyMap],pwfs.referenceSlopesMap);
            
            [n1,n2,n3] = size(pwfs.slopesMap);
            slopesMap_ = reshape(pwfs.slopesMap,n1*n2,n3);
            pwfs.slopes = slopesMap_(pwfs.validSlopes(:),:)*pwfs.slopesUnits;
            pwfs.slopesMap = reshape( slopesMap_ , n1, n2*n3);
        end
        
    end
end

