classdef shackHartmann < hgsetget
    % SHACKHARTMANN Create a Shack-Hartmann object
    %
    % obj = shackHartmann(nLenslet,detectorResolution) creates a
    % Shack-Hartmann object with a (nLenslet X nLenslet) lenslet array and
    % a detector with a detectorResolution resolution
    %
    % obj = shackHartmann(nLenslet,detectorResolution,validLenslet) creates
    % a Shack-Hartmann object with a (nLenslet X nLenslet) lenslet array, a
    % detector with a detectorResolution resolution and a logical mask of
    % size nLenslet setting the location of the valid lenslets inside the
    % lenslet array
    %
    % See also: lensletArray, detector, source, lensletArrayHowto,
    % detectorHowtot
    
    properties
        % lenslet array object
        lenslets;
        % detector object
        camera;
        % camera flat field
        flatField = 0;
        % camera pixel gains
        pixelGains = 1;
        % use quad-cell
        quadCell = false;
        % use center of gravity
        centroiding = true;
        % use matched filter
        matchedFilter = false;
        % use correlation
        correlation = false;
        % use brightes pixel
        brightestPixel = false;
        % centroiding mask
        centroidingMask = 1;
        % timer
        paceMaker;
        % slopes display handle
        slopesDisplayHandle;
        % slope listener
        slopesListener;
        % intensity display handle
        intensityDisplayHandle;
        % intensity listener
        intensityListener;
        % noise display handle
        noiseDisplayHandle;
        % frame pixel threshold
        framePixelThreshold = -inf;
        % slopes units (default:1 is pixel)
        slopesUnits = 1;
        % wavefront units (default:1 is slopes pixel!)
        wavefrontUnits = 1;
        % zernike to slopes conversion matrix
        zern2slopes;
        % wavefront sensor tag
        tag = 'SHACK-HARTMANN';
        % if true the mean of the slopes are removed
        rmMeanSlopes = false;
        % mean slopes
        meanSlopes;
        % handle to the function processing the lenslet intensity
        intensityFunction = @sum;
        % inverse of the sparse gradient matrix
        iG;
        % finiteDifferenceWavefront listener
        finiteDifferenceWavefrontListener
        % finiteDifferenceWavefront handle
        finiteDifferenceWavefrontHandle
        % zernCoefs listener
        zernCoefsListener
        % zernCoefs handle
        zernCoefsHandle
        % sub-aperture dependent elongation convolution kernel
        spotsLgsSrcKernel = [];
        % sub-aperture dependent elongation convolution kernel in Fourier
        % domain
        fftSpotsLgsSrcKernel = [];
        % correlantion reference frame
        correlationRefFrame = [];
        % correlantion reference frame Fourier Transform
        correlationRefFT = [];
        % matched filters for each subaperture
        matchedFilterR = [];
        % number of brightest pixels
        nBrightestPixels = 20;
        % sets the display ON or OFF (useful for cluster execution)
        graphicalDisplay = 0;
        % rotated subaperture frame
        rotatedFrame = false;
        % rotation angles (subaperture dependant)
        theta = 0;
    end
    
    properties (SetAccess=private)
        directionVector;
    end
    
    properties (SetObservable=true)
        % measurements
        slopes=0;
    end
    
    properties (Dependent)
        % valid lenslet mask
        validLenslet;
        % measurements reference
        referenceSlopes;
        % pointing direction [zenith,azimuth]
        pointingDirection;
    end
    
    properties (Dependent, SetAccess=private)
        % number of valid lenslet
        nValidLenslet;
        % number of slopes
        nSlope;
        % intensity in each lenslet
        lensletIntensity;
        % valid actuatord
        validActuator;
        % zernike coefficients
        zernCoefs;
        % X slopes map
        xSlopesMap
        % Y slopes map
        ySlopesMap
        % wavefront estimate from finite difference in wavelength unit
        finiteDifferenceWavefront
    end
    
    properties (Access=protected)
        %         p_slopes;
        p_referenceSlopes=0;
        p_validLenslet;
        p_pointingDirection;
        % index array to reshape a detector frame into a matrix with one
        % raster imagelet per column
        indexRasterLenslet = NaN;
        % lenslet centers
        lensletCenterX;
        lensletCenterY;
        log;
        quadCellX = [1 ;  1 ; -1 ; -1];
        quadCellY = [1 ; -1 ;  1 ; -1];
        meanProjection;
        spotTrail = zeros(2,10);
        p_finiteDifferenceWavefront
    end
    
    methods
        
        %% Constructor
        function obj = shackHartmann(nLenslet,detectorResolution,minLightRatio)
            if nargin>1
                narginchk(1, 4)
                obj.lenslets = lensletArray(nLenslet);
                obj.camera   = detector(detectorResolution);
                if detectorResolution==2
                    obj.quadCell = true;
                    obj.centroiding = false;
                end
                obj.lenslets.nLensletWavePx = ...
                    detectorResolution/nLenslet;
                if nargin>2
                    obj.lenslets.minLightRatio = minLightRatio;
                else
                    obj.lenslets.minLightRatio = 0;
                end
                obj.validLenslet = true(nLenslet);
                obj.camera.frameGrabber ...
                    = obj.lenslets;
                obj.referenceSlopes = zeros(obj.nValidLenslet*2,1);
                obj.p_referenceSlopes = ...
                    repmat(obj.p_referenceSlopes,obj.lenslets.nArray,1);
                
                %             % intensity listener (BROKEN: shackhartmann is not deleted after a clear)
                obj.intensityListener = addlistener(obj.camera,'frame','PostSet',...
                    @(src,evnt) intensityDisplay(obj) );
                obj.intensityListener.Enabled = false;
                
                obj.finiteDifferenceWavefrontListener = addlistener(obj,...
                    'slopes','PostSet',...
                    @(src,evnt) wavefrontDisplay(obj) );
                obj.finiteDifferenceWavefrontListener.Enabled = false;
                
                obj.zernCoefsListener = addlistener(obj,...
                    'slopes','PostSet',...
                    @(src,evnt) bar(obj) );
                obj.zernCoefsListener.Enabled = false;
                
                % Timer settings
                obj.paceMaker = timer;
                obj.paceMaker.name = 'Shack-Hartmann Wavefront Sensor';
                obj.paceMaker.TimerFcn = {@timerCallBack, obj};%(BROKEN: shackhartmann is not deleted after a clear)
                obj.paceMaker.ExecutionMode = 'FixedSpacing';
                obj.paceMaker.BusyMode = 'drop';
                obj.paceMaker.Period = 3;
                obj.paceMaker.ErrorFcn = 'disp('' @detector: frame rate too high!'')';
                %             function timerCallBack( timerObj, event, a)
                %                 %                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
                %                 a.grabAndProcess
                %             end
                display(obj)
                obj.log = logBook.checkIn(obj);
            end
            setSlopesListener(obj)
            function timerCallBack( timerObj, event, a)
                %                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
                %             a.grabAndProcess
                uplus(a)
            end
        end
        
        %% Destructor
        function delete(obj)
            if isvalid(obj.slopesListener)
                delete(obj.slopesListener)
            end
            if isvalid(obj.intensityListener)
                delete(obj.intensityListener)
            end
            if isvalid(obj.paceMaker)
                if strcmp(obj.paceMaker.Running,'on')
                    stop(obj.paceMaker)
                end
                delete(obj.paceMaker)
            end
            if ishandle(obj.slopesDisplayHandle)
                delete(obj.slopesDisplayHandle)
            end
            if ishandle(obj.intensityDisplayHandle)
                delete(obj.intensityDisplayHandle)
            end
            if isvalid(obj.lenslets)
                delete(obj.lenslets)
            end
            if isvalid(obj.camera)
                delete(obj.camera)
            end
            if ~isempty(obj.log)
                checkOut(obj.log,obj);
            end
        end
        
        function display(obj,ngs,tel)
            %% DISPLAY Display object informations
            %
            % display(obj) prints information about the Shack-Hartmann
            % wavefront sensor object
            %
            % display(obj,ngs,tel) provides information about the lenslet (numerical) and detector
            % pixel scale
            
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' Shack-Hartmann wavefront sensor: \n  . %d lenslets total on the pupil\n  . %d pixels per lenslet \n',...
                obj.nValidLenslet,obj.camera.resolution(1)/obj.lenslets.nLenslet)
            if nargin > 2
                d   = tel.D/size(obj.validLenslet,1);
                ps = obj.lenslets.pixelScale(ngs,tel);
                fprintf('   . The current lenslet pixel-scale is %4.2f mas\n',ps.convert('mas')) % sampling of the electic field
                binFactor = 2*obj.lenslets.fieldStopSize/obj.lenslets.nyquistSampling/(obj.camera.resolution(1)/obj.lenslets.nLenslet);
                lo2DInMas = ngs.wavelength/(2*d)*constants.radian2mas;
                fprintf('   . The current detector pixel-scale is %4.2f mas\n',lo2DInMas*binFactor); % sky-angle of the detector pixels after binning
            end
            
            if isinf(obj.framePixelThreshold)
                algoProp = ', no thresholding!';
            else
                algoProp = sprintf(', pixel threshold: %d\n',obj.framePixelThreshold);
            end
            
            algo = {'quadCell','centroiding','matchedFilter','correlation'};
            algoTF = [obj.quadCell,obj.centroiding,obj.matchedFilter,obj.correlation];
            fprintf('  . spot algorithm: %s%s\n',algo{algoTF},algoProp);
            fprintf('----------------------------------------------------\n')
            display(obj.lenslets)
            display(obj.camera)
            
        end
        
        %% Get and Set pointing direction
        function pointingDirection = get.pointingDirection(obj)
            pointingDirection = obj.p_pointingDirection;
        end
        function set.pointingDirection(obj,val)
            obj.p_pointingDirection = val;
            if ~isempty(obj.pointingDirection)
                obj.directionVector = [...
                    tan(obj.p_pointingDirection(1,:)).*cos(obj.p_pointingDirection(2,:));...
                    tan(obj.p_pointingDirection(1,:)).*sin(obj.p_pointingDirection(2,:));...
                    ones(1,size(obj.p_pointingDirection,2))];
            end
        end
        
        function obj = saveobj(obj)
            %% SAVEOBJ
            delete(obj.slopesListener)
            add(obj.log,obj,'Save!')
        end
        
        function INIT(obj)
            %% INIT WFS initialization
            %
            % obj.INIT computes the valid lenslet and set the reference
            % slopes based on the last measurements
            
            add(obj.log,obj,'Setting the valid lenslet and the reference slopes!')
            setValidLenslet(obj);
            obj.referenceSlopes = obj.slopes;
        end
        function lgsKernelInit
            
        end
        %         %% Get and Set slopes
        %         function slopes = get.slopes(obj)
        %             slopes = obj.p_slopes;
        %         end
        %         function set.slopes(obj,val)
        %             obj.p_slopes = val;
        %         end
        
        %% Get and Set valid lenslets
        function validLenslet = get.validLenslet(obj)
            validLenslet = obj.p_validLenslet;
        end
        function set.validLenslet(obj,val)
            obj.p_validLenslet = logical(val);
            index = ~[obj.p_validLenslet(:);obj.p_validLenslet(:)];
            obj.p_referenceSlopes(index) = [];
            obj.slopes(index) = [];
            obj.meanProjection = [ ...
                ones(obj.nValidLenslet,1)  zeros(obj.nValidLenslet,1)
                zeros(obj.nValidLenslet,1) ones(obj.nValidLenslet,1) ];
        end
        
        %% Get number of valid lenslet
        function nValidLenslet = get.nValidLenslet(obj)
            nValidLenslet = sum(obj.validLenslet(:));
        end
        
        %% Get number of slopes
        function nSlope = get.nSlope(obj)
            nSlope = obj.nValidLenslet*2;
        end
        
        %% Get X slopes map
        function out = get.xSlopesMap(obj)
            out = zeros(obj.lenslets.nLenslet);
            out(obj.validLenslet) = obj.slopes(1:end/2);
        end
        
        %% Get Y slopes map
        function out = get.ySlopesMap(obj)
            out = zeros(obj.lenslets.nLenslet);
            out(obj.validLenslet) = obj.slopes(1+end/2:end);
        end
        
        %% Get valid actuators
        function val = get.validActuator(obj)
            nElements            = 2*obj.lenslets.nLenslet+1; % Linear number of lenslet+actuator
            validLensletActuator = zeros(nElements);
            index                = 2:2:nElements; % Lenslet index
            validLensletActuator(index,index) = obj.validLenslet;
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
        
        %% Get/Set the reference spots and update spots location display if
        %% there is one
        function val = get.referenceSlopes(obj)
            val = obj.p_referenceSlopes;
        end
        function set.referenceSlopes(obj,val)
            obj.slopes = obj.slopes + obj.p_referenceSlopes;
            obj.p_referenceSlopes = val;
            obj.slopes = obj.slopes - obj.p_referenceSlopes;
            %             if ishandle(obj.slopesDisplayHandle)
            %                 hc = get(obj.slopesDisplayHandle,'children');
            %                 u = obj.p_referenceSlopes(1:end/2)+obj.lensletCenterX;
            %                 v = obj.p_referenceSlopes(1+end/2:end)+obj.lensletCenterY;
            %                 set(hc(2),'xData',u,'yData',v)
            %             end
        end
        
        %% Get the zernike coeficients
        function val = get.zernCoefs(obj)
            val = obj.zern2slopes'*obj.slopes;
            val = val*obj.wavefrontUnits;
            val(1,:) = []; % piston=0 removed
        end
        
        %% Computes the intensity in each lenslet
        function lensletIntensity = get.lensletIntensity(obj)
            if isempty(obj.camera.frame)
                lensletIntensity = [];
            else
                [nPx,mPx]  = size(obj.camera.frame);
                nLensletArray = obj.lenslets.nArray;
                nPxLenslet = nPx/obj.lenslets.nLenslet;
                mPxLenslet = mPx/obj.lenslets.nLenslet/nLensletArray;
                try
                    buffer     = obj.camera.frame(obj.indexRasterLenslet);
                catch ME
                    fprintf( '@(shackHartmann)> %s\n',ME.identifier)
                    obj.indexRasterLenslet ...
                        = tools.rearrange([nPx,mPx],[nPxLenslet,mPxLenslet]);
                    v = ~obj.validLenslet(:);
                    v = repmat(v,nLensletArray,1);
                    obj.indexRasterLenslet(:,v) = [];
                    buffer     = obj.camera.frame(obj.indexRasterLenslet);
                end
                lensletIntensity = obj.intensityFunction(buffer);
            end
        end
        
        %% Computes the finite difference wavefront
        function out = get.finiteDifferenceWavefront(obj)
            add(obj.log,obj,'Computing the finite differerence wavefront!')
            if isempty(obj.iG)
                add(obj.log,obj,'Computing the finite differerence wavefront!')
                G = sparseGradientMatrix(obj);
                modes = speye((obj.lenslets.nLenslet+1)^2);
                obj.iG = calibrationVault(full(G),...
                    modes(:,obj.validActuator),obj.validActuator,...
                    'noshow',true);
                obj.iG.cond = 100;
            end
            out = obj.iG.M*obj.slopes;
            %             if size(obj.slopes,2)>1
            %                 out = obj.iG.M*obj.slopes;
            %             else
            %                 out = zeros( obj.lenslets.nLenslet+1 );
            %                 out(obj.validActuator) = obj.iG.M*obj.slopes;
            %             end
            out = out*obj.wavefrontUnits;
        end
        
        function setValidLenslet(obj,pupilIntensity)
            %% SETVALIDLENSLET Valid lenslet mask
            %
            % setValidLenslet(obj,pupilIntensity) sets the mask of valid
            % lenslet based on the value of minLightRatio in the lenslets
            % object providing the pupil intensity map
            
            if nargin<2
                pupilIntensity = obj.lensletIntensity./max(obj.lensletIntensity);
            else
                n = length(pupilIntensity);
                nL = n/obj.lenslets.nLenslet;
                pupilIntensity = reshape(pupilIntensity,nL,n*obj.lenslets.nLenslet);
                pupilIntensity = sum(pupilIntensity);
                pupilIntensity = reshape(pupilIntensity,obj.lenslets.nLenslet,obj.lenslets.nLenslet*nL);
                pupilIntensity = reshape(pupilIntensity',nL,obj.lenslets.nLenslet^2);
                pupilIntensity = sum(pupilIntensity);
                surfPx         = nL^2;
                pupilIntensity = pupilIntensity/surfPx;
            end
            obj.validLenslet  = logical( ...
                reshape( pupilIntensity>=obj.lenslets.minLightRatio , ...
                obj.lenslets.nLenslet,obj.lenslets.nLenslet));
            %             obj.referenceSlopes = zeros(2*obj.nValidLenslet,1);
            %             obj.p_referenceSlopes = ...
            %                 repmat(obj.p_referenceSlopes,obj.lenslets.nArray,1);
            %             figure('Name',sprintf('%s valid lenslet',obj.tag)), spy(obj.p_validLenslet)
            dataProcessing(obj)
        end
        
        
        function varargout = dataProcessing(obj)
            %% DATAPROCESSING Processing a SH-WFS detector frame
            %
            % dataProcessing(obj) computes the WFS slopes
            %
            % out = dataProcessing(obj) computes and returns the WFS slopes
            
            [nPx,mPx,nFrame]  = size(obj.camera.frame);
            nLensletArray = obj.lenslets.nArray;
            nPxLenslet = nPx/obj.lenslets.nLenslet;
            mPxLenslet = mPx/obj.lenslets.nLenslet/nLensletArray;
            %             siz(obj.indexRasterLenslet)
            %             obj.nValidLenslet*nLensletArray*nFrame
            %             if numel(obj.indexRasterLenslet)~=(nPxLenslet*mPxLenslet*obj.nValidLenslet*nLensletArray*nFrame)
            if size(obj.indexRasterLenslet,1)~=(nPxLenslet*mPxLenslet) || ...
                    size(obj.indexRasterLenslet,2)~=(obj.nValidLenslet*nLensletArray*nFrame)
                %             try
                % %                 u = obj.indexRasterLenslet;
                % %                 if nFrame>1
                % %                     u = repmat(u,[1,1,nFrame]);
                % %                 end
                %                 buffer     = obj.camera.frame(obj.indexRasterLenslet);
                %             catch ME
                fprintf( '@(shackHartmann)> Setting the raster index \n')
                % get lenslet index
                obj.indexRasterLenslet ...
                    = tools.rearrange([nPx,mPx/nLensletArray,nLensletArray*nFrame],[nPxLenslet,mPxLenslet]);
                % remove index from non-valid lenslets
                v = ~obj.validLenslet(:);
                v = repmat(v,nLensletArray,1);
                v = repmat(v,nFrame,1);
                obj.indexRasterLenslet(:,v) = [];
                %                 u = obj.indexRasterLenslet;
                %                 if nFrame>1
                %                     u = repmat(u,[1,1,nFrame]);
                %                 end
            end
            % Buffer pre-processing
            buffer     = obj.camera.frame(obj.indexRasterLenslet);
            buffer = (buffer - obj.flatField)./obj.pixelGains;
            if isscalar(obj.centroidingMask) || isempty(obj.centroidingMask)
                buffer = bsxfun( @times, obj.centroidingMask(:) , buffer);
            else
                if mPx == nPx
                    buffer = reshape(buffer, [mPxLenslet.^2, obj.nValidLenslet, nFrame]);
                    buffer = bsxfun( @times, obj.centroidingMask(obj.indexRasterLenslet(:,1:obj.nValidLenslet)) , buffer);
                    buffer = reshape(buffer, mPxLenslet.^2, []);
                else
                    buffer = bsxfun( @times, obj.centroidingMask(obj.indexRasterLenslet) , buffer);
                end
            end
            
            %             % Thresholding
            %             if isfinite(obj.framePixelThreshold)
            %                 buffer           = buffer - obj.framePixelThreshold;
            %                 buffer(buffer<0) = 0;
            %             end
            % Thresholding
            if isfinite(obj.framePixelThreshold)
                if numel(obj.framePixelThreshold)>1
                    % intensity based thresholding
                    maxIntensity = max(buffer);
                    threshold    = maxIntensity*obj.framePixelThreshold(2);
                    threshold(threshold<obj.framePixelThreshold(1)) = obj.framePixelThreshold(1);
                    v = obj.validLenslet(:);
                    v = repmat(v,nLensletArray,1);
                    %                     q = zeros(size(v));
                    %                     q(v) = threshold;
                    %                     figure,imagesc(reshape(q,obj.lenslets.nLenslet,[]));set(gca,'clim',[min(threshold),max(threshold)])
                    buffer       = bsxfun( @minus , buffer , threshold);
                else
                    % usual thresholding
                    buffer           = buffer - obj.framePixelThreshold;
                end
                buffer(buffer<0) = 0;
                %                 q = zeros(size(obj.camera.frame));
                %                 q(obj.indexRasterLenslet) = buffer;
                %                 figure,imagesc(q);
            end
            % Centroiding
            if obj.quadCell
                massLenslet ...
                    = sum(buffer)';
                xBuffer = buffer'*obj.quadCellX./massLenslet;
                yBuffer = buffer'*obj.quadCellY./massLenslet;
                xBuffer = reshape(xBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                yBuffer = reshape(yBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                sBuffer ...
                    = bsxfun(@minus,[xBuffer ; yBuffer],obj.referenceSlopes).*obj.slopesUnits;
                index = isnan(sBuffer);
                if any(index(:)) % if all pixels threshold
                    warning('OOMAO:shackHartmann:dataProcessing',...
                        'Threshold (%f) is probably too high or simply there is no light on some of the lenslets',obj.framePixelThreshold)
                    if ~isempty(obj.slopes) && all(size(sBuffer)==size(obj.slopes))
                        sBuffer(index) = obj.slopes(index);
                    end
                end
                %                 obj.slopes = sBuffer;
            elseif obj.centroiding
                massLenslet         = sum(buffer);
                %                 massLenslet(~index) = [];
                %                 buffer(:,~index)    = [];
                %                 size(buffer)
                [x,y]               = ...
                    meshgrid((0:(nPxLenslet-1)),(0:(mPxLenslet-1)));
                if ~obj.rotatedFrame % rotated frame for measurement along and perpendicular to elongated spots
                    %                 xyBuffer  ...
                    %                     = zeros(2*obj.nValidLenslet,1);
                    xBuffer             = bsxfun( @times , buffer , x(:) )  ;
                    xBuffer             = sum( xBuffer ) ./ massLenslet  ;
                    %                 xBuffer             = squeeze((xBuffer));
                    yBuffer             = bsxfun( @times , buffer , y(:) )  ;
                    yBuffer             = sum( yBuffer ) ./ massLenslet  ;
                    %                 yBuffer             = squeeze((yBuffer));
                    %                 xyBuffer = squeeze([xBuffer  yBuffer]);
                    %                 size(xyBuffer)
                    xBuffer = reshape(xBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                    yBuffer = reshape(yBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                    sBuffer = bsxfun(@minus,[xBuffer ; yBuffer],obj.referenceSlopes).*obj.slopesUnits;
                else
                    x = x - (nPxLenslet-1)/2;
                    y = y - (nPxLenslet-1)/2;
                    xPrime = x(:) * cos(obj.theta) + y(:) * sin(obj.theta);
                    %xPrime = xPrime - mean(xPrime(:)) + (nPxLenslet+1)/2;
                    yPrime = -x(:) * sin(obj.theta) + y(:) * cos(obj.theta);
                    %yPrime = yPrime - mean(yPrime(:)) +  (nPxLenslet+1)/2;
                    xPrimeBuffer             = bsxfun( @times , buffer , xPrime )  ;
                    xPrimeBuffer             = sum( xPrimeBuffer ) ./ massLenslet ;
                    %                 xBuffer             = squeeze((xBuffer));
                    yPrimeBuffer             = bsxfun( @times , buffer , yPrime ) ;
                    yPrimeBuffer             = sum( yPrimeBuffer ) ./ massLenslet ;
                    xPrimeBuffer = reshape(xPrimeBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                    yPrimeBuffer = reshape(yPrimeBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                    sBuffer = bsxfun(@minus,[xPrimeBuffer ; yPrimeBuffer],obj.referenceSlopes*0).*obj.slopesUnits;
                end
                index = isnan(sBuffer);
                if any(index(:)) % if all pixels threshold
                    warning('OOMAO:shackHartmann:dataProcessing',...
                        'Threshold (%f) is probably too high or simply there is no light on some of the lenslets',obj.framePixelThreshold)
                    if ~isempty(obj.slopes) && all(size(sBuffer)==size(obj.slopes))
                        sBuffer(index) = obj.slopes(index);
                    end
                end
                %                 obj.slopes = sBuffer;
            elseif obj.correlation
                fprintf('Correlation algorithm\n')
                apod = tukeywin(nPxLenslet, 0.5) * tukeywin(nPxLenslet, 0.5)';
                if isempty(obj.correlationRefFrame) % no offline correlation ref frame
                    obj.correlationRefFrame = reshape(buffer(:, floor(length(obj.validLenslet)/2)+1), nPxLenslet, nPxLenslet);
                    obj.correlationRefFT = conj(ifftshift(fft2(fftshift(obj.correlationRefFrame .* apod))));
                end
                
                if length(obj.correlationRefFT) ~= nPxLenslet% && isempty(obj.correlationRefFT) % multiple correlation ref frames
                    correlationRefFTBuffer = obj.correlationRefFT(obj.indexRasterLenslet);
                    parfor ii = 1:obj.nValidLenslet*nLensletArray*nFrame
                        subap = reshape(buffer(:, ii), nPxLenslet, nPxLenslet);
                        f_mask = reshape(correlationRefFTBuffer(:, ii), nPxLenslet, nPxLenslet);
                        [dX dY] = correlationSH(subap.*apod,f_mask);
                        xBuffer(ii,1) = dX;
                        yBuffer(ii,1) = dY;
                    end
                else
                    correlationRefFTBuffer = conj(ifftshift(fft2(fftshift(obj.correlationRefFrame)))); % single correlation ref frame
                    f_mask = correlationRefFTBuffer;
                    parfor ii = 1:obj.nValidLenslet*nLensletArray*nFrame
                        subap = reshape(buffer(:, ii), nPxLenslet, nPxLenslet);
                        [dX dY] = correlationSH(subap.*apod,f_mask);
                        xBuffer(ii,1) = dX;
                        yBuffer(ii,1) = dY;
                    end
                end
                
                xBuffer = reshape(xBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                yBuffer = reshape(yBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                sBuffer = bsxfun(@minus,[xBuffer ; yBuffer],obj.referenceSlopes).*obj.slopesUnits;
                index = isnan(sBuffer);
                if any(index(:)) % if all pixels threshold
                    warning('OOMAO:shackHartmann:dataProcessing',...
                        'Threshold (%f) is probably too high or simply there is no light on some of the lenslets',obj.framePixelThreshold)
                    if ~isempty(obj.slopes) && all(size(sBuffer)==size(obj.slopes))
                        sBuffer(index) = obj.slopes(index);
                    end
                end
            elseif obj.matchedFilter
                fprintf('Matched filter algorithm\n')
                if isempty(obj.matchedFilterR) % no offline matched filter
                    I0 = reshape(buffer(:, floor(length(obj.validLenslet)/2)+1), nPxLenslet, nPxLenslet);
                    Gx = tools.crop( (tools.shift(I0, 1, 0) - I0) / 1, nPxLenslet);
                    Gy = tools.crop( (tools.shift(I0, 0, 1) - I0) / 1, nPxLenslet);
                    Ixp = tools.shift(I0, 1, 0);
                    Ixm = tools.shift(I0,-1, 0);
                    Iyp = tools.shift(I0, 0, 1);
                    Iym = tools.shift(I0, 0,-1);
                    H = [Gx(:) Gy(:) I0(:) Ixp(:) Ixm(:) Iyp(:) Iym(:)];
                    M = [...
                        1 0 0 1 -1 0  0
                        0 1 0 0  0 1 -1];
                    Cn = diag(I0(:)+obj.camera.readOutNoise.^2);
                    obj.matchedFilterR = sum(I0(:)) * M * pinv(H'*pinv(Cn) * H) * H' * pinv(Cn);
                end
                
                if numel(size(obj.matchedFilterR)) == 3 % && isempty(obj.correlationRefFT) % multiple matchedFilter ref frames
                    matchedFilterRs = obj.matchedFilterR ; % multiple matchedFilter ref frame
                    parfor ii = 1:obj.nValidLenslet*nLensletArray*nFrame
                        subap = reshape(buffer(:, ii), nPxLenslet, nPxLenslet);
                        sMF = matchedFilterRs(:,:,ii)*subap(:) / sum(subap(:));
                        xBuffer(ii,1) = sMF(1);
                        yBuffer(ii,1) = sMF(2);
                    end
                else % single matchedFilter ref frame
                    R = obj.matchedFilterR;
                    parfor ii = 1:obj.nValidLenslet*nLensletArray*nFrame
                        subap = reshape(buffer(:, ii), nPxLenslet, nPxLenslet);
                        sMF = R*subap(:) / sum(subap(:));
                        xBuffer(ii,1) = sMF(1);
                        yBuffer(ii,1) = sMF(2);
                    end
                end
                
                xBuffer = reshape(xBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                yBuffer = reshape(yBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                sBuffer = bsxfun(@minus,[xBuffer ; yBuffer],obj.referenceSlopes).*obj.slopesUnits + (nPxLenslet-1)/2;
                index = isnan(sBuffer);
                if any(index(:)) % if all pixels threshold
                    warning('OOMAO:shackHartmann:dataProcessing',...
                        'Threshold (%f) is probably too high or simply there is no light on some of the lenslets',obj.framePixelThreshold)
                    if ~isempty(obj.slopes) && all(size(sBuffer)==size(obj.slopes))
                        sBuffer(index) = obj.slopes(index);
                    end
                end
                
            elseif obj.brightestPixel
                fprintf('Brightest pixel algorithm\n')
                fprintf('%i%s', obj.nBrightestPixels, ' brightest pixels')
                fprintf('\n')
                nBrightPixels = obj.nBrightestPixels;
                parfor ii = 1:obj.nValidLenslet*nLensletArray*nFrame
                    subap = reshape(buffer(:, ii), nPxLenslet, nPxLenslet);
                    sortSubap = sort(subap(:));
                    brightestPixelThreshold = sortSubap(end-(nBrightPixels));
                    % subap = subap - min(sortSubap(end-(nBrightPixels-1):end));
                    subap(subap<=brightestPixelThreshold) = 0;
                    [dX dY] = cog(subap);
                    xBuffer(ii,1) = dX;
                    yBuffer(ii,1) = dY;
                end
                
                xBuffer = reshape(xBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                yBuffer = reshape(yBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                sBuffer = bsxfun(@minus,[xBuffer ; yBuffer],obj.referenceSlopes).*obj.slopesUnits;
                index = isnan(sBuffer);
                if any(index(:)) % if all pixels threshold
                    warning('OOMAO:shackHartmann:dataProcessing',...
                        'Threshold (%f) is probably too high or simply there is no light on some of the lenslets',obj.framePixelThreshold)
                    if ~isempty(obj.slopes) && all(size(sBuffer)==size(obj.slopes))
                        sBuffer(index) = obj.slopes(index);
                    end
                end
            end
            
            if obj.rmMeanSlopes % remove mean slopes
                obj.meanSlopes = (obj.meanProjection'*sBuffer)/obj.nValidLenslet;
                obj.slopes = sBuffer - obj.meanProjection*obj.meanSlopes;
            else
                obj.slopes = sBuffer;
            end
            
            if nargout>0
                varargout{1} = obj.slopes;
            end
        end
        
        function zern = getZernike(obj,radialOrder)
            zern = zernike(1:zernike.nModeFromRadialOrder(radialOrder),...
                'resolution',obj.lenslets.nLenslet,...
                'pupil',double(obj.validLenslet));
            dzdxy = [zern.xDerivative(obj.validLenslet,:);zern.yDerivative(obj.validLenslet,:)];
            zern.c = dzdxy\obj.slopes;
            %             zern.c = dzdxy'*obj.slopes;
        end
        
        function varargout = grabAndProcess(obj)
            %% GRABANDPROCESS Frame grabbing and processing
            %
            % grabAndProcess(obj) grabs a frame and computes the slopes
            %
            % out = grabAndProcess(obj) grabs a frame, computes and returns
            % the slopes
            
            warning('OOMAO;shackHartmann:grabAndProcess',...
                'DEPRECATED! Instead use the uplus operator (+obj)')
            grab(obj.camera)
            dataProcessing(obj);
            if nargout>0
                varargout{1} = obj.slopes;
            end
        end
        function varargout = uplus(obj)
            %% UPLUS + Update operator
            %
            % +obj grabs a frame and computes the slopes
            %
            % obj = +obj returns the shackHartmann object
            
            grab(obj.camera)
            dataProcessing(obj);
            if nargout>0
                varargout{1} = obj;
            end
        end
        
        function relay(obj,src,~)
            %% RELAY shackhartmann to source relay
            %
            % relay(obj,src) propagates the source through the
            % Shack-Hartmann lenslet array, grab a frame from the detector
            % adding noise if any and process the frame to get the
            % wavefront slopes
            
            if ~isempty(obj.pointingDirection)
                
                nSrc = length(src);
                m_directionVector = obj.directionVector;
                if size(m_directionVector,2)<nSrc
                    m_directionVector = repmat( m_directionVector(:) , 1 , nSrc );
                end
                for kSrc = 1:nSrc
                    delta = m_directionVector(:,kSrc) - src(kSrc).directionVector;
                    if any(delta)
                        warning('oomao:shackHartmann:relay',...
                            'WFS not align on source, this induces an additional tip-tilt error!')
                        D = src(kSrc).opticalPath{1}.D;
                        
                        nOutWavePx    = obj.lenslets.nLensletWavePx*obj.lenslets.fftPad;    % Pixel length of the output wave
                        evenOdd       = rem(obj.lenslets.nLensletWavePx,2);
                        if ~rem(nOutWavePx,2) && evenOdd
                            nOutWavePx = nOutWavePx + evenOdd;
                        end
                        nOutWavePx = max(nOutWavePx,obj.lenslets.nLensletWavePx);
                        %                 fprintf(' - new nOutWavePx = %d\n',nOutWavePx);
                        
                        alpha_p = (src(kSrc).wavelength/D)*obj.lenslets.nLensletWavePx*obj.lenslets.nLenslet/nOutWavePx;
                        delta = delta/alpha_p;
                        
                        [u,v]         = ndgrid((0:(obj.lenslets.nLensletWavePx*obj.lenslets.nLenslet-1)));
                        
                        u = u*delta(2)/nOutWavePx;
                        v = v*delta(1)/nOutWavePx;
                        src(kSrc).phase = 2*pi*(u+v);
                    end
                end
                
            end
            
            %             if isempty(src(1).magnitude)
            %                 obj.camera.photonNoise = false;
            %             else
            %                 obj.camera.photonNoise = true;
            %             end
            %             propagateThrough(obj.lenslets,src)
            if nargin<3
                relay(obj.lenslets,src)
            end
            %             grabAndProcess(obj)
            spotsSrcKernelConvolution(obj,src)
            grab(obj.camera)
            
            if obj.camera.frameCount==0
                dataProcessing(obj);
            else
                obj.slopes = zeros(obj.nSlope,1);
            end
        end
        
        function spotsSrcKernelConvolution(obj,src)
            
            if ~isempty(src(1).extent)% || ~isempty(obj.spotsLgsSrcKernel)
                
                add(obj.log,obj,'Convolution of the spots by source kernel!')
                
                picture   = obj.lenslets.imagelets;
                
                [nPx,mPx,nPicture] = size(picture);
                nPxLenslet = nPx/obj.lenslets.nLenslet;
                mPxLenslet = mPx/obj.lenslets.nLenslet/obj.lenslets.nArray;
                mPxNPicture = mPx*nPicture;
                picture = reshape(picture,nPx,mPxNPicture);
                nLensletArray = obj.lenslets.nArray*nPicture;
                
                indexRasterLenslet_ ...
                    = tools.rearrange(size(picture),[nPxLenslet,mPxLenslet]);
                v = ~obj.validLenslet(:);
                v = repmat(v,nLensletArray,1);
                indexRasterLenslet_(:,v) = [];
                buffer     = picture(indexRasterLenslet_);
                
                buffer     = reshape(buffer,nPxLenslet,nPxLenslet,[]);
                %                 tic
                %                 if isempty(obj.spotsLgsSrcKernel)
                srcExtent = src(1).extent;
                parfor kLenslet=1:size(buffer,3)
                    buffer(:,:,kLenslet) = conv2(buffer(:,:,kLenslet),srcExtent,'same');
                end
                %                 else
                %                     parfor kLenslet=1:size(buffer,3)
                %                         buffer(:,:,kLenslet) = conv2(buffer(:,:,kLenslet),obj.spotsLgsSrcKernel(:,:,kLenslet),'same');
                %                     end
                %                 end
                %                 toc
                
                % cropping in case of elongatedFieldStopSize;
                % L. Blanco 2017/01/02
                %                 if ~isempty(obj.lenslets.elongatedFieldStopSize)
                %                     nCrop = obj.lenslets.fieldStopSize * 2 * obj.lenslets.nyquistSampling;
                %                     buffer = tools.crop(buffer, nCrop);
                %                     nLenslet = size(buffer, 3);
                %                     nRows = sqrt(nLenslet);
                %                     croppedPicture = zeros(nCrop*nRows);
                %                     %conflict with size of picture(indexRasterLenslet_)
                %                     indexRasterLenslet_ = tools.rearrange(size(croppedPicture),[nCrop,nCrop]);
                %                     v = ~obj.validLenslet(:);
                %                     v = repmat(v,nLensletArray,1);
                %                     indexRasterLenslet_(:,v) = [];
                %                     croppedPicture(indexRasterLenslet_) = buffer;
                %                     obj.lenslets.imagelets = reshape( croppedPicture , nCrop*nRows , nCrop*nRows , nPicture);
                %                 else
                %end cropping
                picture(indexRasterLenslet_) = buffer;
                obj.lenslets.imagelets = reshape( picture , nPx , mPx , nPicture);
                %                 end
                
            end
            
        end
        
        
        function out = framelets(obj,lensletI,lensletJ,lensletArrayK)
            %% FRAMELETS Per lenslet detector frame
            %
            % out = framelets(obj,lensletI,lensletJ,lensletArrayK) returns
            % the detector frame restricted to lenslet (I,J) of array # K
            
            nLenslet  = obj.lenslets.nLenslet;
            cameraRes = obj.camera.resolution;
            %             nArray    = obj.lenslets.nArray;
            if nargin<4
                lensletArrayK = 1;
            end
            nPxLenslet = cameraRes/nLenslet;
            u = (1:nPxLenslet(1)) + (lensletI-1)*nPxLenslet(1);
            v = (1:nPxLenslet(2)) + (lensletJ-1)*nPxLenslet(2) + (lensletArrayK-1)*cameraRes(1);
            out = obj.camera.frame(u,v);
        end
        
        function varargout = slopesDisplay(obj,varargin)
            %% SLOPESDISPLAY WFS slopes display
            %
            % slopesDisplay(obj) displays image plot of the slopes
            %
            % slopesDisplay(obj,'PropertyName',PropertyValue) displays
            % image plot of the slopes and set the properties of the
            % graphics object quiver
            %
            % h = slopesDisplay(obj,...) returns the graphics handle
            
            if obj.lenslets.nLenslet>1
                
                nSlopes   = size(obj.slopes,2);
                slopesMap = zeros(2*obj.lenslets.nLenslet^2,nSlopes);
                p         = repmat(obj.validLenslet(:),2,nSlopes);
                slopesMap(p) = obj.slopes;
                slopesMap    = reshape(slopesMap,obj.lenslets.nLenslet,[]);
                
                if nSlopes>3 && rem(nSlopes,2)==0
                    slopesMap = cell2mat(reshape( mat2cell( ...
                        slopesMap,obj.lenslets.nLenslet,2*obj.lenslets.nLenslet*ones(1,nSlopes)) , 2, []));
                end
                
                if ishandle(obj.slopesDisplayHandle)
                    nm = size(slopesMap);
                    if ~all(size(get(obj.slopesDisplayHandle(1),'Cdata'))==nm)
                        set(get(obj.slopesDisplayHandle(1),'parent'),'xlim',0.5+[0 nm(2)],'ylim',0.5+[0 nm(1)])
                    end
                    set(obj.slopesDisplayHandle,'CData',slopesMap)
                else
                    obj.slopesDisplayHandle = imagesc(slopesMap,varargin{:});
                    ax = gca;
                    pos = get(ax,'position');
                    axis xy equal tight
                    ylabel(colorbar('location','EastOutside'),'Pixel')
                    set(ax,'position',pos)
                    
                    hu = findobj(gcf,'Type','uimenu','Label','OOMAO');
                    if isempty(hu)
                        hu = uimenu('Label','OOMAO');
                    end
                    hus  = uimenu(hu,'Label','Slopes Listener Off','Callback',@oomaoMenu);
                    if isvalid(obj.slopesListener) && obj.slopesListener.Enabled
                        set(hus,'Label','Slopes Listener On')
                    end
                    
                    %                     set(gcf,'WindowButtonMotionFcn',@wbmcb)
                    
                end
                
            else
                
                
                if ishandle(obj.slopesDisplayHandle)
                    set(obj.slopesDisplayHandle(1),'XData',obj.slopes(1),'YData',obj.slopes(2))
                    obj.spotTrail = circshift(obj.spotTrail,[0,-1]);
                    obj.spotTrail(:,end) = obj.slopes;
                    set(obj.slopesDisplayHandle(2),'XData',obj.spotTrail(1,:),'YData',obj.spotTrail(2,:))
                else
                    obj.slopesDisplayHandle(1) = plot(obj.slopes(1),obj.slopes(2),'+',...
                        'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],...
                        'MarkerSize',10,'LineWidth',2);
                    obj.spotTrail = zeros(2,10);
                    obj.spotTrail(:,end) = obj.slopes;
                    obj.slopesDisplayHandle(2) = ...
                        line(obj.spotTrail(1,:),obj.spotTrail(2,:),'color','r');
                    %                     set(gca,'xlim',[-1,1],'ylim',[-1,1])
                    grid on
                    axis square
                    
                    hu = findobj(gcf,'Type','uimenu','Label','OOMAO');
                    if isempty(hu)
                        hu = uimenu('Label','OOMAO');
                    end
                    hus  = uimenu(hu,'Label','Slopes Listener Off','Callback',@oomaoMenu);
                    if obj.slopesListener.Enabled
                        set(hus,'Label','Slopes Listener On')
                    end
                end
                
            end
            
            if nargout>0
                varargout{1} = obj.slopesDisplayHandle;
            end
            
            function oomaoMenu(src,~)
                obj.slopesListener.Enabled = ~obj.slopesListener.Enabled;
                if obj.slopesListener.Enabled
                    set(src,'Label','Slopes Listener On')
                else
                    set(src,'Label','Slopes Listener Off')
                end
            end
            
            %         function wbmcb(src,evnt)
            %            cp = get(get(obj.slopesDisplayHandle,'parent'),'CurrentPoint');
            %            disp(round([cp(1,1),cp(1,2)]))
            % %            xdat = [xinit,cp(1,1)];
            % %            ydat = [yinit,cp(1,2)];
            % %            set(hl,'XData',xdat,'YData',ydat);drawnow
            %         end
            
        end
        
        function varargout = intensityDisplay(obj,varargin)
            %% INTENSITYDISPLAY WFS lenslet intensity display
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
            
            intensity = zeros(obj.lenslets.nLenslet,obj.lenslets.nLenslet*obj.lenslets.nArray);
            v = obj.validLenslet(:);
            v = repmat(v,obj.lenslets.nArray,1);
            intensity(v) = obj.lensletIntensity;
            if ishandle(obj.intensityDisplayHandle)
                set(obj.intensityDisplayHandle,...
                    'Cdata',intensity,varargin{:})
            else
                obj.intensityDisplayHandle = imagesc(intensity,varargin{:});
                ax = gca;
                pos = get(ax,'position');
                axis equal tight xy
                colorbar('location','EastOutside')
                set(ax,'position',pos)
            end
            set(get(obj.intensityDisplayHandle,'parent'),'Clim',[floor(min(intensity(v))),ceil(max(intensity(v)))])
            if nargout>0
                varargout{1} = obj.intensityDisplayHandle;
            end
        end
        
        function slopesAndFrameDisplay(obj,varargin)
            imagesc(obj.camera,varargin{:});
            slopesDisplay(obj,'matrix',...
                makehgtform('translate',[1,1,0]),varargin{:});
            %             if isinf(obj.framePixelThreshold)
            %                 clim = get(gca,'clim');
            %                 set(gca,'clim',[obj.framePixelThreshold,clim(2)])
            %             end
            
        end
        
        function slopesAndIntensityDisplay(obj,varargin)
            intensityDisplay(obj,varargin{:});
            n  = obj.lenslets.nLensletImagePx;
            slopesDisplay(obj,'matrix',...
                makehgtform('translate',-[(n-1)/2,(n-1)/2,0]/n,'scale',1/n,'translate',[1,1,0]*2),varargin{:});
        end
        %%
        function subaperturesDisplay(obj,tel,dm)
            nSubap = size(obj.validLenslet,1);
            nPxSubap = tel.resolution/nSubap;
            
            dSubap = tel.D/nSubap;
            xLeftLim = meshgrid(linspace(-tel.D/2, tel.D/2-dSubap,tel.resolution/nPxSubap));
            %xRightLim = xLeftLim + dSubap;
            yLeftLim = xLeftLim';
            %yRightLim = xRightLim';
            figure(345), hold on
            u = linspace(-tel.D/2, tel.D/2,tel.resolution);
            %gridMask = wfs.validLensletSamplingMask(nPxSubap);
            %imagesc(u,u,gridMask(1:end-1,1:end-1))
            imagesc(u,u,tel.pupil)
            pos11 = xLeftLim(obj.validLenslet) + 1i*yLeftLim(obj.validLenslet);
            %scatter(real(pos11), imag(pos11))
            pos21 = xLeftLim(obj.validLenslet)+dSubap + 1i*yLeftLim(obj.validLenslet);
            %scatter(real(pos21), imag(pos21),'r')
            pos12 = xLeftLim(obj.validLenslet) + 1i*(yLeftLim(obj.validLenslet)+dSubap);
            %scatter(real(pos12), imag(pos12),'g')
            pos22 = xLeftLim(obj.validLenslet)+dSubap + 1i*(yLeftLim(obj.validLenslet)+dSubap);
            %scatter(real(pos22), imag(pos22),'k')
            fill(real([pos11 pos12 pos21  pos22])', imag([pos11 pos12 pos22 pos21])','r','FaceColor','none')
            if nargin == 3
                scatter(real(dm.modes.actuatorCoord(dm.validActuator)), imag(dm.modes.actuatorCoord(dm.validActuator)))
            end
            axis square
            title('geometric xLocs and yLocs')
            xlabel('position [m]')
            ylabel('position [m]')
        end
        
        function wavefrontDisplay(obj,varargin)
            wft = zeros( obj.lenslets.nLenslet+1 );
            wft__ = obj.finiteDifferenceWavefront;
            wftRms = std(wft__);
            wft(obj.validActuator) = wft__;
            
            if any(ishandle(obj.finiteDifferenceWavefrontHandle))
                set(obj.finiteDifferenceWavefrontHandle(1),...
                    'CData',wft)
                set(obj.finiteDifferenceWavefrontHandle(2),...
                    'String',sprintf('F.D. Wavefront: %5.2f',wftRms))
            else
                obj.finiteDifferenceWavefrontHandle(1) = ...
                    imagesc(wft,...
                    varargin{:});
                ax = gca;
                pos = get(ax,'position');
                axis equal tight xy
                colorbar('location','EastOutside')
                set(ax,'position',pos)
                obj.finiteDifferenceWavefrontHandle(2) = ...
                    title(sprintf('F.D. Wavefront: %5.2f',wftRms))
            end
        end
        
        function bar(obj,varargin)
            if ishandle(obj.zernCoefsHandle)
                set(obj.zernCoefsHandle,...
                    'YData',obj.zernCoefs)
            else
                obj.zernCoefsHandle = ...
                    bar((1:length(obj.zernCoefs))+1,obj.zernCoefs,...
                    varargin{:});
                grid on
                xlabel('Zernike modes')
                ylabel('Zern. Coefs.')
            end
        end
        
        function [G,mask] = sparseGradientMatrix(obj)
            %% SPARSEGRADIENTMATRIX
            %
            % Gamma = sparseGradientMatrix(obj) computes the sparse
            % gradient such as a wavefront in wavelength units multiply
            % by the gradient matrix gives the centroid in pixel units
            
            nLenslet = obj.lenslets.nLenslet;
            nPxPhase = obj.lenslets.nLenslet + 1;
            nElPhase = nPxPhase^2;
            
            % non zeros row index
            rows = repmat( 1:nLenslet^2 , 4 , 1);
            
            % non zeros column index
            cols = [(1:2) (1:2) + nPxPhase]'; % 1st lenslet stencil
            cols = bsxfun(@plus, cols, 0:nLenslet-1 ); % 1st lenslet column
            p(1,1,1:nLenslet) = (0:nLenslet-1)*nPxPhase;
            cols = bsxfun( @plus, cols , p); % lenslet rows
            
            
            % non zeros values
            nzv = repmat( [-1 1 -1 1]', 1 , nLenslet^2);
            Gy = sparse(rows(:),cols(:),nzv(:),nLenslet^2,nElPhase);
            nzv = repmat( [-1 -1 1 1]', 1 , nLenslet^2);
            Gx = sparse(rows(:),cols(:),nzv(:),nLenslet^2,nElPhase);
            
            mask = ~obj.validActuator;
            Gx(:,mask) = [];
            Gy(:,mask) = [];
            mask = ~obj.validLenslet;
            Gx(mask,:)  = [];
            Gy(mask,:)  = [];
            
            % scaling factor such as a phase in units of wavelength
            % multiply by the gradient matrix will return the centroid in
            % units of pixels
            nOutWavePx    = obj.lenslets.nLensletWavePx*obj.lenslets.fftPad;    % Pixel length of the output wave
            evenOdd       = rem(obj.lenslets.nLensletWavePx,2);
            if ~rem(nOutWavePx,2) && evenOdd
                nOutWavePx = nOutWavePx + evenOdd;
            end
            nOutWavePx = max(nOutWavePx,obj.lenslets.nLensletWavePx);
            a = obj.lenslets.nLensletWavePx/nOutWavePx;
            
            G = 0.5*[Gx;Gy]/a;
            %             figure
            %             spy(G)
            
        end
        
        function varargout = sparseGradientMatrix3x3Stentil(obj)
            %% SPARSEGRADIENTMATRIX
            %
            % Gamma = sparseGradientMatrix(obj)
            %
            % [Gamma,gridMask] = sparseGradientMatrix(obj)
            
            nLenslet = obj.lenslets.nLenslet;
            nMap     = 2*nLenslet+1;
            nValidLenslet ...
                = obj.nValidLenslet;
            
            i0x = [1:3 1:3]; % x stencil row subscript
            j0x = [ones(1,3) ones(1,3)*3]; % x stencil col subscript
            i0y = [1 3 1 3 1 3]; % y stencil row subscript
            j0y = [1 1 2 2 3 3]; % y stencil col subscript
            s0x = [-1 -2 -1  1 2  1]/2; % x stencil weight
            s0y = -[ 1 -1  2 -2 1 -1]/2; % y stencil weight
            %             s0x = [-1 -1 -1  1 1  1]/3; % x stencil weight
            %             s0y = -[ 1 -1  1 -1 1 -1]/3; % y stencil weight
            
            i_x = zeros(1,6*nValidLenslet);
            j_x = zeros(1,6*nValidLenslet);
            s_x = zeros(1,6*nValidLenslet);
            i_y = zeros(1,6*nValidLenslet);
            j_y = zeros(1,6*nValidLenslet);
            s_y = zeros(1,6*nValidLenslet);
            
            [iMap0,jMap0] = ndgrid(1:3);
            gridMask = false(nMap);
            
            u   = 1:6;
            
            % Accumulation of x and y stencil row and col subscript and weight
            for jLenslet = 1:nLenslet
                jOffset = 2*(jLenslet-1);
                for iLenslet = 1:nLenslet
                    
                    if obj.validLenslet(iLenslet,jLenslet)
                        
                        iOffset= 2*(iLenslet-1);
                        i_x(u) = i0x + iOffset;
                        j_x(u) = j0x + jOffset;
                        s_x(u) = s0x;
                        i_y(u) = i0y + iOffset;
                        j_y(u) = j0y + jOffset;
                        s_y(u) = s0y;
                        u = u + 6;
                        
                        gridMask( iMap0 + iOffset , jMap0 + jOffset ) = true;
                        
                    end
                    
                end
            end
            indx = sub2ind([nMap,nMap],i_x,j_x); % mapping the x stencil subscript into location index on the phase map
            indy = sub2ind([nMap,nMap],i_y,j_y); % mapping the y stencil subscript into location index on the phase map
            % row index of non zero values in the gradient matrix
            v = 1:2*nValidLenslet;
            v = v(ones(6,1),:);
            % sparse gradient matrix
            Gamma = sparse(v,[indx,indy],[s_x,s_y],2*nValidLenslet,nMap^2);
            Gamma(:,~gridMask) = [];
            
            varargout{1} = Gamma;
            if nargout>1
                varargout{2} = gridMask;
            end
            
        end
        %%
        function varargout = sparseGradientMatrixAmplitudeWeighted(obj,amplMask)
            %% SPARSEGRADIENTMATRIX
            %
            % Gamma = sparseGradientMatrixAmplitudeWeighted(obj,amplMask)
            %
            % [Gamma,gridMask] = sparseGradientMatrixAmplitudeWeighted(obj,amplMask)
            
            
            nLenslet = obj.lenslets.nLenslet;
            factor = 2;
            
            if nargin == 1
                amplMask = ones(factor*nLenslet+1);
            end
            nMap     = factor*nLenslet+1;
            nValidLenslet ...
                = obj.nValidLenslet;
            dsa = 1; %(meters) Raven subaperture width = 0.8m; NFIRAOS subaperture width = 0.5 m....
            if factor == 2
                i0x = [1:3 1:3 1:3]; % x stencil row subscript
                j0x = [ones(1,3) ones(1,3)*2 ones(1,3)*3]; % x stencil col subscript
                i0y = [1 2 3 1 2 3 1 2 3]; % y stencil row subscript
                j0y = [1 1 1 2 2 2 3 3 3]; % y stencil col subscript
                s0x = [-1/4 -1/2 -1/4 0 0 0 1/4 1/2  1/4]*(1/dsa); % x stencil weight. Raven subaperture width = 0.8m
                % NFIRAOS subaperture width = 0.5 m....
                s0y = -[ 1/4 0 -1/4  1/2 0 -1/2 1/4 0 -1/4]*(1/dsa); % y stencil weight
                %             s0x = [-1 -1 -1  1 1  1]/3; % x stencil weight
                %             s0y = -[ 1 -1  1 -1 1 -1]/3; % y stencil weight
                Gv = [-2 2 -1 1; -2 2 -1 1; -1 1 -2 2; -1 1 -2 2];
                
                i_x = zeros(1,9*nValidLenslet);
                j_x = zeros(1,9*nValidLenslet);
                s_x = zeros(1,9*nValidLenslet);
                i_y = zeros(1,9*nValidLenslet);
                j_y = zeros(1,9*nValidLenslet);
                s_y = zeros(1,9*nValidLenslet);
                
                [iMap0,jMap0] = ndgrid(1:3);
                gridMask = false(nMap);
                
                u   = 1:9;
            elseif factor == 4
                
                i0x = [1:5 1:5]; % x stencil row subscript
                j0x = [ones(1,5) ones(1,5)*5]; % x stencil col subscript
                i0y = [1 5 1 5 1 5 1 5 1 5]; % y stencil row subscript
                j0y = [1 1 2 2 3 3 4 4 5 5]; % y stencil col subscript
                s0x = [-1 -1.5 -2 -1.5 -1  1 1.5 2 1.5 1]/factor; % x stencil weight
                s0y = -[ 1 -1 1.5 -1.5 2 -2 1.5 -1.5 1 -1]/factor; % y stencil weight
                %             s0x = [-1 -1 -1  1 1  1]/3; % x stencil weight
                %             s0y = -[ 1 -1  1 -1 1 -1]/3; % y stencil weight
                
                i_x = zeros(1,10*nValidLenslet);
                j_x = zeros(1,10*nValidLenslet);
                s_x = zeros(1,10*nValidLenslet);
                i_y = zeros(1,10*nValidLenslet);
                j_y = zeros(1,10*nValidLenslet);
                s_y = zeros(1,10*nValidLenslet);
                
                [iMap0,jMap0] = ndgrid(1:5);
                gridMask = false(nMap);
                
                u   = 1:10;
            end
            % Accumulation of x and y stencil row and col subscript and weight
            for jLenslet = 1:nLenslet
                jOffset = factor*(jLenslet-1);
                for iLenslet = 1:nLenslet
                    
                    if obj.validLenslet(iLenslet,jLenslet)
                        
                        I = (iLenslet-1)*factor+1;
                        J = (jLenslet-1)*factor+1;
                        
                        a = amplMask(I:I+factor,J:J+factor);
                        numIllum = sum(a(:));
                        
                        if numIllum == 9
                            
                            iOffset= factor*(iLenslet-1);
                            i_x(u) = i0x + iOffset;
                            j_x(u) = j0x + jOffset;
                            s_x(u) = s0x;
                            i_y(u) = i0y + iOffset;
                            j_y(u) = j0y + jOffset;
                            s_y(u) = s0y;
                            u = u + (factor+1)*3;
                            
                            gridMask( iMap0 + iOffset , jMap0 + jOffset ) = true;
                            
                        elseif numIllum ~= 9
                            a11 = a(1:2,1:2);
                            a21 = a(2:3,1:2);
                            a12 = a(1:2,2:3);
                            a22 = a(2:3,2:3);
                            %
                            
                            du11 = (a11(:))'*Gv/(3*sum(a11(:)));
                            du21 = (a21(:))'*Gv/(3*sum(a21(:)));
                            du12 = (a12(:))'*Gv/(3*sum(a12(:)));
                            du22 = (a22(:))'*Gv/(3*sum(a22(:)));
                            du11 = reshape(du11,2,2);
                            du12 = reshape(du12,2,2);
                            du21 = reshape(du21,2,2);
                            du22 = reshape(du22,2,2);
                            
                            s(1,1) = du11(1,1);
                            s(2,1) = du11(2,1) + du21(1,1);
                            s(3,1) = du21(2,1);
                            s(1,2) = du11(1,2) + du12(1,1);
                            s(2,2) = du11(2,2) + du12(2,1) + du21(1,2) + du22(1,1);
                            s(3,2) = du21(2,2) + du22(2,1);
                            s(1,3) = du12(1,2);
                            s(2,3) = du12(2,2) + du22(1,2);
                            s(3,3) = du22(2,2);
                            
                            mySx = s';
                            mySy = mySx';
                            
                            %                             sPy = [s11 s21 s31; s12 s22 s32; s13 s23 s33]*0.5;
                            %                             sPx = sPy';
                            %                             sPx = sPx(:);
                            %                             sPy = sPy(:);
                            
                            iOffset= factor*(iLenslet-1);
                            i_x(u) = i0x + iOffset;
                            j_x(u) = j0x + jOffset;
                            s_x(u) = mySx(:);
                            i_y(u) = i0y + iOffset;
                            j_y(u) = j0y + jOffset;
                            s_y(u) = mySy(:);
                            u = u + (factor+1)*3;
                            
                            gridMask( iMap0 + iOffset , jMap0 + jOffset ) = true;
                            
                            
                        end
                    end
                    
                end
            end
            indx = sub2ind([nMap,nMap],i_x,j_x); % mapping the x stencil subscript into location index on the phase map
            indy = sub2ind([nMap,nMap],i_y,j_y); % mapping the y stencil subscript into location index on the phase map
            % row index of non zero values in the gradient matrix
            v = 1:2*nValidLenslet;
            v = v(ones((factor+1)*3,1),:);
            % sparse gradient matrix
            Gamma = sparse(v,[indx,indy],[s_x,s_y],2*nValidLenslet,nMap^2);
            Gamma(:,~gridMask) = [];
            
            varargout{1} = Gamma;
            if nargout>1
                varargout{2} = gridMask;
            end
            
        end
        
        function gridMask = validLensletSamplingMask(obj,sample)
            %% VALIDLENSLETSAMPLINGMASK
            %
            % mask = validLensletSamplingMask(obj,n) computes the mask
            % corresponding to the nXn pixel sampling of the lenslets
            
            nLenslet = obj.lenslets.nLenslet;
            nMap     = (sample-1)*nLenslet+1;
            
            [iMap0,jMap0] = ndgrid(1:sample);
            gridMask = false(nMap);
            
            % Accumulation of x and y stencil row and col subscript and weight
            for jLenslet = 1:nLenslet
                jOffset = (sample-1)*(jLenslet-1);
                for iLenslet = 1:nLenslet
                    
                    if obj.validLenslet(iLenslet,jLenslet)
                        
                        iOffset= (sample-1)*(iLenslet-1);
                        
                        gridMask( iMap0 + iOffset , jMap0 + jOffset ) = true;
                        
                    end
                    
                end
            end
            
        end
        
        function out = validLensletArea(obj,R)
            nLenslet ...
                = obj.lenslets.nLenslet;
            validLenslet ...
                = obj.validLenslet;
            nValidLenslet ...
                = obj.nValidLenslet;
            d   = obj.lenslets.pitch;
            
            u   = d*(1-nLenslet:2:nLenslet-1)/2;
            [xLenslet,yLenslet] = meshgrid( u );
            xLenslet = xLenslet(validLenslet);
            yLenslet = yLenslet(validLenslet);
            
            lensletCoordX = ...
                [xLenslet+d/2 xLenslet-d/2 xLenslet-d/2 xLenslet+d/2];
            lensletCoordY = ...
                [yLenslet+d/2 yLenslet+d/2 yLenslet-d/2 yLenslet-d/2];
            rLenslet = hypot(lensletCoordX,lensletCoordY);
            
            out = zeros(nValidLenslet,1);
            for kLenslet = 1:nValidLenslet
                
                plot(lensletCoordX(kLenslet,:),lensletCoordY(kLenslet,:))
                pause
                if any(rLenslet(kLenslet,:)>R)
                    xA = lensletCoordX(kLenslet,2);
                    xB = lensletCoordX(kLenslet,4);
                    yA = lensletCoordY(kLenslet,2);
                    out(kLenslet) = yA*(xB-xA);
                    xA = xA/R;
                    xB = xB/R;
                    out(kLenslet) = out(kLenslet) - ...
                        R.^2.*( asin(xB) - asin(xA) - ...
                        xA*sqrt(1-xA^2) + xB*sqrt(1-xB^2) )/2;
                else
                    out(kLenslet) = obj.lenslets.pitch^2;
                end
                
            end
            
        end
        
        function out = phaseToSlopes(obj,xo,yo,do,sample,alpha,beta,tel)
            %% PHASETOSLOPES
            
            xo = xo(:)';
            yo = yo(:)';
            
            if nargin<6
                alpha = 1;
                beta  = zeros(1,2);
            end
            
            nLenslet ...
                = obj.lenslets.nLenslet;
            validLenslet ...
                = obj.validLenslet;
            nValidLenslet ...
                = obj.nValidLenslet;
            d   = alpha*obj.lenslets.pitch;
            
            u   = d*(1-nLenslet:2:nLenslet-1)/2;
            [xLenslet,yLenslet] = meshgrid( u );
            xLenslet = xLenslet(validLenslet);
            yLenslet = yLenslet(validLenslet);
            
            edge   = linspace(-d/2,d/2,sample)';
            unit   = ones(1,sample-1)*d/2;
            % lenslet contour (x coordinates)
            contourX_ = ...
                [fliplr(edge(2:end)') -unit  edge(1:end-1)'  unit];
            % lenslet contour (y coordinates)
            contourY_ = ...
                [unit   fliplr(edge(2:end)') -unit edge(1:end-1)'];
            % normal to lenslet contour ( X gradient)
            xNormal_ = [ zeros(1,sample-1) -ones(1,sample-1) ...
                zeros(1,sample-1) ones(1,sample-1)];
            % normal to lenslet contour ( Y gradient)
            yNormal_ = [ ones(1,sample-1) zeros(1,sample-1)  ...
                -ones(1,sample-1) zeros(1,sample-1)];
            
            z_  = cell(2*nValidLenslet,1);
            partialLenslet = 0;
            fprintf(' @(shackHartmann.phaseToSlopes)> #%4d/    ', nValidLenslet);
            figure('Name','Truncated lenslet')
            u = linspace(-1,1,obj.lenslets.nLenslet+1).*tel.R;
            axes('xlim',[-1,1]*tel.R,'ylim',[-1,1]*tel.R,...
                'xtick',u,'ytick',u,...
                'xtickLabel',[],'ytickLabel',[],...
                'xgrid','on','ygrid','on')
            axis square
            for kLenslet=1:nValidLenslet
                
                fprintf('\b\b\b\b%4d',kLenslet)
                
                %                 vertex = hypot( ...
                %                     xLenslet(kLenslet) + [ +d -d -d +d]/2 , ...
                %                     yLenslet(kLenslet) + [ +d +d -d -d]/2 );
                
                %                 if any(vertex>tel.R)
                %                     partialLenslet = partialLenslet + 1;
                % lenslet contour (x coordinates)
                contourX = xLenslet(kLenslet) + contourX_;
                % lenslet contour (y coordinates)
                contourY = yLenslet(kLenslet) + contourY_;
                % normal to lenslet contour ( X gradient)
                xNormal = xNormal_;
                % normal to lenslet contour ( Y gradient)
                yNormal = yNormal_;
                % polar coordinate contour
                [contourO,contourR] = cart2pol( contourX , contourY);
                % contour out of pupil
                index   = contourR>tel.R;
                if any(index)
                    contourR(index) = tel.R;
                    xNormal(index) = cos(contourO(index));
                    yNormal(index) = sin(contourO(index));
                    [contourX,contourY] = pol2cart(contourO,contourR);
                    partialLenslet = partialLenslet + 1;
                    line(contourX,contourY,'Marker','.','MarkerEdgeColor','r')
                    drawnow
                end
                % contour steps
                ds = abs( diff( ...
                    complex( ...
                    [contourX contourX(1)] , ...
                    [contourY contourY(1)] ) ... % complex
                    ) ... % diff
                    ); % abs
                A(kLenslet,:) = [ sum(ds.*contourX.*xNormal) , sum(ds.*contourY.*yNormal) ];
                
                
                x_ = bsxfun( @minus , contourX' + beta(1), xo )/do;
                y_ = bsxfun( @minus , contourY' + beta(2), yo )/do;
                zs = bsxfun( @times , ds' , linearSpline(x_).*linearSpline(y_) );
                % contour sum ( X gradient)
                z_{kLenslet} = ...
                    sum( bsxfun( @times , zs , xNormal' ) );
                % contour sum ( Y gradient)
                z_{kLenslet+nValidLenslet} = ...
                    sum( bsxfun( @times , zs , yNormal' ) );
                
                %                 end
                
                %                 x_m = (xLenslet(kLenslet) - d/2 - xo)/do;
                %                 x_p = (xLenslet(kLenslet) + d/2 - xo)/do;
                %                 y_ = bsxfun( @minus , edge + yLenslet(kLenslet), yo )/do;
                %                 z_{kLenslet} = sum(linearSpline(y_)).*...
                %                     ( linearSpline(x_p) - linearSpline(x_m) );
                %
                %                 y_m = (yLenslet(kLenslet) - d/2 - yo)/do;
                %                 y_p = (yLenslet(kLenslet) + d/2 - yo)/do;
                %                 x_ = bsxfun( @minus , edge + xLenslet(kLenslet), xo )/do;
                %                 z_{kLenslet+nValidLenslet} = sum(linearSpline(x_)).*...
                %                     ( linearSpline(y_p) - linearSpline(y_m) );
                
            end
            
            fprintf(' (%d truncated lenslet)\n',partialLenslet)
            out = cell2mat(z_)/d;
            
            %             [n,m] = size(out);
            %             [i,j,s] = find(out);
            %             as = abs(s);
            %             index = as > 1e-12*as;
            %             out = sparse(i(index),j(index),s(index),n,m);
            
        end
        
        %% generate the LGS elongated spot objects as seen by every subaperture
        function  o = generateElongatedKernel(obj,tel,lgs, binningFactor)
            % L.Blanco 02/28/2017 Modified to remove kernels fft
            % computations
            
            %[o fftO] = generateElongatedKernel(obj,tel,lgs, binningFactor)
            nSources =length([lgs(1,:,1).zenith]);
            obj.camera.binningFactor = binningFactor;
            %compute max elongation to create kernel array o
            dH = [lgs(:,1,:).height] - [lgs(:,1,:).objectiveFocalLength];
            d = tel.D/obj.lenslets.nLenslet;
            pixelScale = lgs(1).wavelength/d/2*obj.lenslets.nyquistSampling;
            % maxElong = ceil(max(abs(tel.D*((dH)./([lgs(:,1,:).height].^2+[lgs(:,1,:).height].*dH))/pixelScale)));
            maxElong = ceil(max(abs(tel.D*((dH)./([lgs(:,1,:).height] .* [lgs(1).objectiveFocalLength]))/pixelScale)));
            
            if nSources > 1 % recursive call in case of multiple LGS
                for iSrc = 1:nSources
                    a = tic;
                    o = generateElongatedKernel(obj,tel,lgs(1,iSrc,:), binningFactor);
                    % [o, fftO] = generateElongatedKernel(obj,tel,lgs(1,iSrc,:), binningFactor);
                    obj.spotsLgsSrcKernel(:,:,:,iSrc) = o;
                    % obj.fftSpotsLgsSrcKernel(:,:,:,iSrc) = fftO;
                    obj.camera.spotsLgsSrcKernel(:,:,:,iSrc) = o;
                    % obj.camera.fftSpotsLgsSrcKernel(:,:,:,iSrc) = fftO;
                    toc(a)
                end
            end
            p = [lgs.nPhoton]/sum([lgs.nPhoton]); % flux from Na profile
            xL = lgs(1).viewPoint(1);
            yL = lgs(1).viewPoint(2);
            d = tel.D/obj.lenslets.nLenslet;
            uLenslet = linspace(-1,1,obj.lenslets.nLenslet)*(tel.D/2-d/2);
            [xLenslet,yLenslet] = meshgrid(uLenslet);
            %             maskLenslet = obj.validLenslet;
            %             xLenslet = xLenslet(maskLenslet);
            %             yLenslet = yLenslet(maskLenslet);
            xx = xLenslet(:)-xL;
            yy = yLenslet(:)-yL;
            pixelScale = lgs(1).wavelength/d/2*obj.lenslets.nyquistSampling;
            dH = [lgs.height] - [lgs.objectiveFocalLength];
            %sx = xx*((dH)./([lgs.height].^2+[lgs.height].*dH))/pixelScale; % relative horizontal shifts due to vertical elongation
            %sy = yy*((dH)./([lgs.height].^2+[lgs.height].*dH))/pixelScale; % relative vertical shifts due to vertical elongation
            
            sx = -xx*((dH)./([lgs.height] .* [lgs.objectiveFocalLength]))/pixelScale; % relative horizontal shifts due to vertical elongation
            sy = -yy*((dH)./([lgs.height] .* [lgs.objectiveFocalLength]))/pixelScale; % relative vertical shifts due to vertical elongation
            sx(isnan(sx)) = 0; % in case we have a point source in z (dH=Nan in that case which creates an error)
            sy(isnan(sy)) = 0; % in case we have a point source in z (dH=Nan in that case which creates an error)
            if isempty(lgs(1).extent)
                ps = zeros(4);
                ps(2,2) = 1;
            else
                ps = lgs(1).extent;
                %modif TFu - 2019/04/29
                %%%%%%%%
                nSubap = size(obj.validLenslet,1);
                nPxSubap = tel.resolution/nSubap;
                ps = tools.crop(ps, nPxSubap*binningFactor);
                %%%%%%%%
                % instead of ...
                %ps = tools.crop(ps, 32);
                %%%%%%%
                %                 if mod(size(ps, 1),2) == 1
                %                     ps = ps(2:end, 2:end); % even number of pixels in src.extent
                %                 end
            end
            
            %expand the kernel size to avoid circularization when shifting
            maxElong(isnan(maxElong)) = 0;
            nPxElong = size(ps, 1) + 2 * maxElong;
                
            nPxElong = 4 * binningFactor * ceil(nPxElong / (4*binningFactor));
            %nPxElong = binningFactor * ceil(nPxElong / (binningFactor));

 
            %             twos = 2.^linspace(1,10,10);
            %             nPxElong = twos(find( (abs(nPxElong-twos)) == min(abs(nPxElong-twos))));
           
            ps = tools.crop(ps, nPxElong); %modif ou pas je sais plus 
            
            %binning to match detector sampling
            %binningFactor = obj.lenslets.nLensletsImagePx ./ obj.camera.resolution(1);
          
            nPxElong = nPxElong / binningFactor; % the spot kernels are generated using the WFS binning value (vs lenslet pixels)
           
            %%%%%%%%
           o = zeros([nPxElong nPxElong obj.lenslets.nLenslet^2]);
            % unbinnedO = zeros(length(ps));
            % fftO = zeros([nPxElong nPxElong obj.lenslets.nLenslet^2]);
            nLenslets = obj.lenslets.nLenslet^2;
            lgsHeight = [lgs.height];
            lgsExtent = lgs(1).extent;
            parfor iLenslet = 1:nLenslets
                unbinnedO = zeros(length(ps));
                for iHeight = 1:length(lgsHeight)
                    if isempty(lgsExtent) % no lateral extensioon of the source
                        %o(:,:,iLenslet) = o(:,:,iLenslet) + p(iHeight) .* tools.shift(ps,sx(iLenslet,iHeight), sy(iLenslet,iHeight));
                        defaultShift = 0.; %MODIF Thierry
                        shiftedKernel = tools.shift(ps,sx(iLenslet,iHeight)+defaultShift, sy(iLenslet,iHeight)+defaultShift);
                        shiftedKernel = shiftedKernel .* (shiftedKernel>0);
                        o(:,:,iLenslet) = o(:,:,iLenslet) + tools.binning(p(iHeight) .* shiftedKernel, [nPxElong,nPxElong]);
                    else % non-zero lateral extension of the source
                        defaultShift = 0.;%MODIF Thierry
                        shiftedKernel = tools.shift(ps,sx(iLenslet,iHeight)+defaultShift, sy(iLenslet,iHeight)+ defaultShift);
                        shiftedKernel = shiftedKernel .* (shiftedKernel>=0);
                        %o(:,:,iLenslet) = o(:,:,iLenslet) + p(iHeight) .* tools.shift(ps,sx(iLenslet,iHeight)-0.5, sy(iLenslet,iHeight)-0.5);
                        %o(:,:,iLenslet) = o(:,:,iLenslet) + tools.binning(p(iHeight) .* shiftedKernel, [nPxElong,nPxElong]);
                        unbinnedO = unbinnedO + p(iHeight) .* shiftedKernel;
                        % fftO(:,:,iLenslet) = fftO(:,:,iLenslet) + fftshift(fft2(tools.binning(p(iHeight) .* shiftedKernel, [nPxElong,nPxElong]), nPxElong*2, nPxElong*2));
                    end
                end
                o(:,:,iLenslet) = tools.binning(unbinnedO, [nPxElong,nPxElong]);
                % fftO(:,:,iLenslet) = fftshift(fft2(o(:,:,iLenslet), nPxElong, nPxElong));
            end
            
            if isempty(obj.spotsLgsSrcKernel)
                obj.spotsLgsSrcKernel = o;
                % obj.fftSpotsLgsSrcKernel = fftO;
                obj.camera.spotsLgsSrcKernel = o;
                % obj.camera.fftSpotsLgsSrcKernel = fftO;
                obj.lenslets.convKernel = 'true';
            end
        end
        
        function varargout = theoreticalNoise(obj,tel,atm,gs,ss,varargin)
            %% THEORETICALNOISE WFS theoretical noise
            %
            % noiseVar = theoreticalNoise(obj,tel,atm,gs,ss) computes the
            % theoretical noise variance for a telescope, an atmosphere, a
            % guide star and a science star objects
            %
            % noiseVar = theoreticalNoise(obj,tel,atm,gs,ss,nd) computes the
            % theoretical noise variance for a telescope, an atmosphere, a
            % guide star, a science star objects and the fwhm of a
            % diffraction limited spot in pixel (default: nd=2)
            %
            % noiseVar = theoreticalNoise(obj,...,'skyBackgroundMagnitude',sky)
            % computes the theoretical noise variance including background
            % noise specified with the sky backgroung magnitude at the
            % wavelength of the guide star
            %
            % noiseVar = theoreticalNoise(obj,...,'soao',true) computes the
            % theoretical noise variance for AO corrected WFS
            %
            % noiseVar = theoreticalNoise(obj,...,'naParam',[deltaNa,naAltidude])
            % computes the theoretical noise variance for each leanslet
            % according to the spot elongation derived from the Na layer
            % parameters ; the lgs is launched on-axis
            %
            % noiseVar = theoreticalNoise(obj,...,'naParam',[deltaNa,naAltidude],'lgsLaunchCoord',[xL,yL])
            % computes the theoretical noise variance for Na LGS WFS which
            % the LGS launch telescope location is given by the coordinates
            % [xL,yL]
            %
            % If  lgsLaunchCoord has multiple entries (rows), then the
            % function calls itself recursively, outputing a noise
            % covariance matrix in as many cells, each following the [X;Y]
            % slopes numbering and not the concatenation of 2x2 blocks
            
            
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'shackHartmann') );
            inputs.addRequired('tel',@(x) isa(x,'telescopeAbstract') );
            inputs.addRequired('atm',@(x) isa(x,'atmosphere') );
            inputs.addRequired('gs',@(x) isa(x,'source') );
            inputs.addRequired('ss',@(x) isa(x,'source') );
            inputs.addOptional('ND',2,@isnumeric);
            inputs.addParameter('skyBackground',[],@isnumeric);
            inputs.addParameter('soao',false,@islogical);
            inputs.addParameter('lgsLaunchCoord',[0,0],@isnumeric);
            inputs.addParameter('naParam',[],@isnumeric);
            inputs.addParameter('verbose',true,@islogical);
            inputs.addParameter('NS',[],@isnumeric);
            inputs.addParameter('ensquaredEnergy',1,@isnumeric);
            inputs.addParameter('centroidingAlgorithm','',@ischar);
            inputs.addParameter('emccd',0,@isnumeric);
            inputs.addParameter('rotatedFrame',false,@islogical);
            inputs.addParameter('computeInverseCovMat',false,@islogical);
            inputs.parse(obj,tel,atm,gs,ss,varargin{:});
            
            obj    = inputs.Results.obj;
            tel    = inputs.Results.tel;
            atm    = inputs.Results.atm;
            gs     = inputs.Results.gs;
            ss     = inputs.Results.ss;
            skyBackground ...
                = inputs.Results.skyBackground;
            soao   = inputs.Results.soao;
            ND     = inputs.Results.ND;
            NS     = inputs.Results.NS;
            launchCoord...
                = inputs.Results.lgsLaunchCoord;
            naParam= inputs.Results.naParam;
            verbose= inputs.Results.verbose;
            ensquaredEnergy = inputs.Results.ensquaredEnergy;
            centroidingAlgorithm = inputs.Results.centroidingAlgorithm;
            emccd = inputs.Results.emccd;
            rotatedFrame = inputs.Results.rotatedFrame;
            computeInverseCovMat = inputs.Results.computeInverseCovMat;
            naLgs = false;
            nLgs = size(launchCoord,1);
            nLgs = length(gs);
            noiseVar = cell(nLgs,1);
            % Recursive call to theoreticalNoise if there are more than 1
            % LGS at a time
            if  nLgs > 1
                for iLgs = 1:nLgs
                    noiseVar{iLgs} = obj.theoreticalNoise(tel, atm, gs(iLgs), ss,...
                        'naParam',naParam(iLgs,:),'lgsLaunchCoord',launchCoord(iLgs,:),'centroidingAlgorithm',centroidingAlgorithm,'rotatedFrame',rotatedFrame);
                    if computeInverseCovMat
                        buffer = pinv(full(noiseVar{iLgs}));
                        buffer(abs(buffer)<1e-7) = 0;
                        noiseVar{iLgs} = sparse(buffer);
                    end
                end
                varargout{1} = noiseVar;
                if nargout>1
                    varargout{2} = nph(1);
                end
                return
            end
            % SET WFSs WAVELENGTH
            if gs(1).wavelength ~= atm.wavelength
                originalAtmWavelength = atm.wavelength;
                atm.wavelength = gs(1).wavelength;
                fprintf('ATTENTION, ATM wavelength changed to match guide-star wavelength.\n', obj.camera.exposureTime);
            else
                originalAtmWavelength = atm.wavelength;
            end
            
            nLenslet = obj.lenslets.nLenslet;
            % WFS Pitch
            d = tel.D/nLenslet;
            
            if ~isempty(naParam)
                
                naLgs = true;
                
                deltaNa    = naParam(1);
                naAltitude = naParam(2);
                
                xL = launchCoord(1);
                yL = launchCoord(2);
                
                uLenslet = linspace(-1,1,nLenslet)*(tel.D/2-d/2);
                [xLenslet,yLenslet] = meshgrid(uLenslet);
                maskLenslet = obj.validLenslet;
                xLenslet = xLenslet(maskLenslet);
                yLenslet = yLenslet(maskLenslet);
                
                [oe,re] = cart2pol(xLenslet-xL,yLenslet-yL);
                %                 re = hypot(xLenslet,yLenslet);
                %thetaNa = re*deltaNa/naAltitude^2;
                thetaNa = re*deltaNa/(naAltitude^2 + naAltitude*deltaNa);
            end
            
            if obj.camera.exposureTime ~= tel.samplingTime
                fprintf('ATTENTION, telescope sampling and camera sampling DIFFERENT. Using %1.1f sampling time.\n', obj.camera.exposureTime);
            end
            % Photon #
            if obj.camera.exposureTime ~= 1
                nph = ensquaredEnergy.*obj.lenslets.throughput*obj.camera.quantumEfficiency.*...
                    [gs.nPhoton]*obj.camera.exposureTime*min(tel.area,d^2);
            else
                nph = ensquaredEnergy.*obj.lenslets.throughput*obj.camera.quantumEfficiency.*...
                    [gs.nPhoton]*(obj.camera.exposureTime/obj.camera.clockRate)*tel.samplingTime*min(tel.area,d^2);
            end
            %             nph = obj.lenslets.throughput*obj.camera.quantumEfficiency.*...
            %                 [gs.nPhoton]*obj.camera.exposureTime*obj.lenslets.nLensletImagePx^2*...
            %                 tel.area/tel.pixelArea;
            
            if verbose
                add(obj.log,obj,sprintf('lenslet pitch  : %4.2f cm',d*1e2))
                add(obj.log,obj,sprintf('Fried parameter: %4.2f cm',atm.r0*1e2))
                add(obj.log,obj,sprintf('Number of source photon: %g per subaperture per frame',nph(1)))
            end
            
            % Atmosphere WFS wavelength scaling
            atmWavelength = atm.wavelength;
            atm.wavelength = gs(1).wavelength;
            
            % FWHM in diffraction unit
            if soao
                fwhm = ones(obj.nValidLenslet,1)/d;
            elseif naLgs
                dNa   = gs(1).wavelength./thetaNa;
                if verbose
                    add(obj.log,obj,sprintf('dNa max-min: [%4.2f , %4.2f]cm',max(dNa)*1e2,min(dNa)*1e2))
                end
                index = dNa>min(d,atm.r0);
                dNa(index)...
                    = min(d,atm.r0);
                %                 fwhm  = sqrt(1./atm.r0^2+1./dNa.^2);
                fwhm  = [1./min(d,atm.r0) ; 1./dNa];
                seeingNa = atm.seeingInArcsec*constants.arcsec2radian;
            else
                fwhm = ones(obj.nValidLenslet,1)./min(d,atm.r0);
            end
            
            % Sky backgound photon #
            if isempty(skyBackground)
                nbg = 0;
            else
                skyBackground = source('wavelength',gs(1).photometry,'magnitude',skyBackground);
                nbg = obj.lenslets.throughput*obj.camera.quantumEfficiency.*...
                    skyBackground.nPhoton*obj.camera.exposureTime*tel.area*...
                    obj.camera.pixelScale^2*...
                    prod(obj.camera.resolution/obj.lenslets.nLenslet);
                fprintf(' @(shackHartmann:theoreticalNoise)> Number of background photon %4.2f per frame\n',nbg)
            end
            % WFS phase diff. noise variance
            ron = obj.camera.readOutNoise;
            
            nGs =length(gs);
            noiseVar = zeros(length(fwhm),nGs);
            for kGs = 1:nGs
                
                if nLenslet>1
                    snr = sqrt(2*nph(kGs).^2./( nph(kGs) + ...
                        (2/3)*(gs(kGs).wavelength./ss.wavelength).^2.*(4*ron*d.*fwhm*ND).^2 + ...
                        8*nbg/3) );
                else % quad-cell SNR
                    snr = nph(kGs)./sqrt(nph(kGs) + 4*ron.^2. + nbg);
                end
                noiseVar(:,kGs) = ...
                    (gs(kGs).wavelength./ss.wavelength).^2.*(pi.*d.*fwhm./snr).^2;
                if obj.lenslets.nLenslet==1
                    noiseVar(:,kGs) = (3*pi/16)^2*noiseVar(:,kGs)/4; % To comply with Hardy and Tyler formulaes
                end
                
            end
            
            if naLgs
                %for iPos = 1:size(launchCoord,1)
                if isempty(centroidingAlgorithm)
                    %                 noiseVar = (1/(8*log(2)))*(2*atm.r0.*fwhm).^2/nph + ...
                    %                     (ron/nph).^2.*NS.^2/12;
                    %                 thetaNa
                    %                 seeingNa*constants.radian2arcsec
                    NS = 2*ceil(2*thetaNa/seeingNa);
                    fprintf('NS max-min: [%d,%d]\n',max(NS),min(NS))
                    sigma2X = (1/(8*log(2)))*(2*atm.r0.*fwhm(1)).^2/nph + ...
                        (ron/nph).^2.*NS.^2/12;
                    sigma2Y = (1/(8*log(2)))*(2*atm.r0.*fwhm(2:end)).^2/nph + ...
                        (ron/nph).^2.*NS.^2/12;
                elseif strcmp(centroidingAlgorithm,'cog') % Thomas08, Study of optimal wavefront sensing with elongated laser guide stars, Eq. 3 and following
                    pixelScale = gs(1).wavelength/...
                        (2*d*obj.lenslets.nyquistSampling);
                    NS = max(seeingNa, thetaNa)*seeingNa/pixelScale^2;
                    Nsamp = seeingNa/pixelScale;
                    Ax = pi^2/(2*log(2))*(atm.r0.*fwhm(1)).^2;
                    Bx = pi^2/3*(NS/Nsamp).^2;
                    Ay = pi^2/(2*log(2))*(atm.r0.*fwhm(2:end)).^2;
                    By = pi^2/3*(NS/Nsamp).^2;
                    sigma2X = Ax/nph + Bx*(ron/nph).^2;
                    sigma2Y = Ay/nph + By*(ron/nph).^2;
                end
                
                %                 figure
                %                 map = zeros(nLenslet);
                % %                 size(map(obj.validLenslet))
                % %                 size(sigma2Y)
                %                 map(obj.validLenslet) = sigma2Y + sigma2X;
                %                 imagesc(map)
                
                if ~rotatedFrame
                    B = zeros(obj.nSlope*nGs,3);
                    noiseCovarDiag = [ ...
                        sigma2X.*sin(oe).^2 + sigma2Y.*cos(oe).^2  ...
                        sigma2X.*cos(oe).^2 + sigma2Y.*sin(oe).^2]';
                    noiseCovarDiagP1 = ...
                        -(sigma2X.*ones(obj.nValidLenslet,1) - sigma2Y).*...
                        cos(oe).*sin(oe);
                    B(:,1) = noiseCovarDiag(:);
                    B(1:2:end,2) = noiseCovarDiagP1;
                    B(2:2:end,3) = noiseCovarDiagP1;
                    noiseVar = spdiags(B,[0,-1,1],obj.nSlope*nGs,obj.nSlope*nGs);
                    % noiseVar = bsxfun( @plus , noiseVar(1,:) , noiseVar(2:end,:) );
                    
                    % --- show figure ---
                    if ishandle(obj.noiseDisplayHandle)
                        figure(obj.noiseDisplayHandle.Number)
                        hold on
                    else
                        obj.noiseDisplayHandle = figure;
                    end
                    di = diag(noiseVar);
                    ra = di(1:2:end);
                    rb = di(2:2:end);
                    ra = sqrt(ra.^2 + rb.^2);
                    plotLineColor = {'b', 'k', 'r', 'g', 'm', 'c', 'y'};
                    ellipse(ra/max(ra)/4,ra/max(ra)/20,oe, xLenslet, yLenslet,plotLineColor{randi(7)});
                    hold on
                    scatter(xL, yL, 'ro')
                    box on
                    title('Noise on LGS sub-apertures','fontsize',18)
                    ylabel('distance, [m]')
                    xlabel('distance, [m]')
                    % -----------------------
                    
                    % change output format to comply with [X;Y] slopes
                    % concatenation
                    B = reshape(noiseCovarDiag',2*length(noiseCovarDiag),1);
                    B(1:1:end/2,2) = noiseCovarDiagP1;
                    B(end/2+1:1:end,3) = noiseCovarDiagP1;
                    noiseVar = spdiags(B,[0,-obj.nValidLenslet,obj.nValidLenslet],obj.nSlope*nGs,obj.nSlope*nGs);
                else
                    noiseVar = diag([sigma2X;sigma2Y]);
                end
                
            else
                nGs =length(gs);
                noiseVar = zeros(length(fwhm),nGs);
                for kGs = 1:nGs
                    if isempty(centroidingAlgorithm)
                        if nLenslet>1
                            snr = sqrt(2*nph(kGs).^2./( nph(kGs) + ...
                                (2/3)*(gs(kGs).wavelength./ss.wavelength).^2.*(4*ron*d.*fwhm*ND).^2 + ...
                                8*nbg/3) );
                        else % quad-cell SNR
                            snr = nph(kGs)./sqrt(nph(kGs) + 4*ron.^2. + nbg);
                        end
                        noiseVar(:,kGs) = ...
                            (gs(kGs).wavelength./ss.wavelength).^2.*(pi.*d.*fwhm./snr).^2;
                        if obj.lenslets.nLenslet==1
                            noiseVar(:,kGs) = (3*pi/16)^2*noiseVar(:,kGs)/4; % To comply with Hardy and Tyler formulaes
                        end
                    elseif strcmp(centroidingAlgorithm,'wcog') %thomas06_ComparisonCentroidAlgos_mnras.pdf Eq (23) and (24)
                        Nd = 2*obj.lenslets.nyquistSampling;% #pixels in diffraction-limited spot (may need some fixing for spots sampled below nyquist)
                        %Computation of Nt
                        if d/atm.r0 < 2.
                            Nt = Nd;
                        else
                            Nt = Nd*d/atm.r0*(sqrt(1-(atm.r0/d)^(1./3))) ; % formule de la PSF corrige du tilt
                        end
                        %;Computation of Nw
                        Nw = 2*Nt;
                        if nph(kGs)/ron <= 2
                            Nw           = Nt;
                        end
                        if nph(kGs)/ron > 5,
                            Nw           = 1.5*Nt;
                        end
                        if nph(kGs)/ron > 10,
                            Nw           = 2.*Nt;
                        end
                        if nph(kGs)/ron > 20,
                            Nw           = 3.*Nt;
                        end
                        alphar       =  Nw^2/(Nt^2+Nw^2); %gain lie au WCOG
                        photonNoise  = pi^2/(2*log(2.))*1./nph(kGs)*(Nt/Nd)^2*((Nt^2+Nw^2)/(2*Nt^2+Nw^2))^2/alphar^2;
                        readOutNoise = pi^3/(32*(log(2.))^2)*(ron^2/nph(kGs)^2)*((Nt^2+Nw^2)/Nd)^2/alphar^2;
                        %gain1        = gain*alphar;   %<= prise en compte du gain optique WCOG
                        
                        if emccd
                            noiseVar(:,kGs) = 2*photonNoise + readOutNoise;
                        else
                            noiseVar(:,kGs) = photonNoise + readOutNoise;
                        end
                    elseif strcmp(centroidingAlgorithm,'cog') %thomas06_ComparisonCentroidAlgos_mnras.pdf Eq (23) and (24)
                        Nd = 2*obj.lenslets.nyquistSampling;% #pixels in diffraction-limited spot (may need some fixing for spots sampled below nyquist)
                        %Computation of Nt
                        if d/atm.r0 < 2.
                            Nt = Nd;
                        else
                            Nt = Nd*d/atm.r0*(sqrt(1-(atm.r0/d)^(1./3))) ; % formule de la PSF corrige du tilt
                        end
                        %Nt = 2.5*Nd;
                        %Computation of Ns <= cas d'une formule de CDG classique
                        Ns = 1.5*Nt; % nombre de pixels pour le calcul du CDG
                        photonNoise = pi^2/(2*log(2.))*1./nph(kGs)*(Nt/Nd)^2;
                        readOutNoise = pi^2/3*(ron/nph(kGs))^2*(Ns^2/Nd)^2;
                        %gain1 = gain;
                        if emccd
                            noiseVar(:,kGs) = 2*photonNoise + readOutNoise;
                        else
                            noiseVar(:,kGs) = photonNoise + readOutNoise;
                        end
                    elseif strcmp(centroidingAlgorithm,'quadcell')%thomas06_ComparisonCentroidAlgos_mnras.pdf Eqs (28-29)
                        kappa = 1; %case of a diffraction-limited spot
                        photonNoise = pi^2*kappa*1./nph(kGs);
                        readOutNoise = 4*pi^2*kappa^2*(ron/nph(kGs))^2;
                        if emccd
                            noiseVar(:,kGs) = 2*photonNoise + readOutNoise;
                        else
                            noiseVar(:,kGs) = photonNoise + readOutNoise;
                        end
                    end
                    
                end
                
            end
            
            % Resetting atmosphere wavelength
            atm.wavelength = atmWavelength;
            varargout{1} = noiseVar;
            if nargout>1
                varargout{2} = nph(1);
            end
            atm.wavelength = originalAtmWavelength;
        end
        function gainCalibrationDev(obj,tel,ngs)
            %% GAINCALIBRATION
            % calibrate gain of cetner of gravity
            %
            % wfs.gainCalibration(tel,ngs)
            %
            % the WFS output is in units of \lambda/D/2, i.e. half the
            % diffraction limit for the aperture/sub-aperture
            %
            % When the pixel is \lambda/D (i.e. 2x \lambda/D/2) the gain
            % becomes ~2. For higher binning factors the scaling factor
            % grows in a non-linear fashion
            
            nPx = obj.camera.resolution(1)/size(obj.validLenslet,1);
            d   = tel.D/size(obj.validLenslet,1);
            
            ngs = ngs.*tel*obj;
            obj.pointingDirection = zeros(2,1);
            
            ps = obj.lenslets.pixelScale(ngs,tel);
            fprintf('The current lenslet pixel-scale is %f mas\n',ps.convert('mas')) % sampling of the electic field
            binFactor = 2*obj.lenslets.fieldStopSize/obj.lenslets.nyquistSampling/obj.camera.resolution(1);
            lo2DInMas = ngs.wavelength/(2*d)*constants.radian2mas;
            lo2D = ngs.wavelength/(2*d);
            detectorPixelSizeInMas = lo2DInMas*binFactor;
            detectorPixelSize = lo2D*binFactor;
            fprintf('The current detector pixel-scale is %f mas\n',detectorPixelSizeInMas); % sky-angle of the detector pixels after binning
            
            % the next three lines of code provide better LTAO performance
            % at the centre of the field in the case of KAPA (+7% in H for ZA=50, median MK profile!).
            pixelScale = obj.lenslets.fieldStopSize*ngs.wavelength/d/nPx;
            fprintf('The current delta %f mas\n',pixelScale*constants.radian2mas);
            tipStep = pixelScale/2/2;
            nStep   = floor(nPx/3)*2;
            sx      = zeros(1,nStep+1);
            u       = 0:nStep;
            obj.camera.frameListener.Enabled = false;
            obj.slopesListener.Enabled = false;
            
            warning('off','oomao:shackHartmann:relay')
            for kStep=u
                ngs.zenith = -tipStep*kStep;
                +ngs;
                drawnow
                sx(kStep+1) = median(obj.slopes(1:end/2));
            end
            warning('on','oomao:shackHartmann:relay')
            
            Ox_in  = u*tipStep*constants.radian2arcsec;
            Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
            figure
            plot(Ox_in, Ox_out)
            hold
            plot(Ox_in, Ox_in,'k:')
            xlabel('input')
            ylabel('output')
            %plot(u*tipStep/pixelScale,(u*tipStep/pixelScale)./sx);
            
            slopesLinCoef = polyfit(Ox_in,Ox_out,1);
            title(['gain = ' num2str(1/slopesLinCoef(1))])
            obj.slopesUnits = 1/slopesLinCoef(1);
            ngs.zenith = 0;
            obj.pointingDirection = [];
            
            
        end
        
        function gainCalibration(obj,tel,ngs,unit)
            %% GAINCALIBRATION
            % calibrate gain of cetner of gravity
            %
            % wfs.gainCalibration(tel,ngs)
            %
            % the WFS output is default in units of \lambda/D/2, i.e. half the
            % diffraction limit for the aperture/sub-aperture
            %
            % When the pixel is \lambda/D (i.e. 2x \lambda/D/2) the gain
            % becomes ~2. For higher binning factors the scaling factor
            % grows in a non-linear fashion
            %
            % unit toggles between units of aperture/sub-aperture
            % diffraction or units of pixels such that for a source located
            % one pixelSize off-axis the measurement is exactly one
            if nargin < 4
                unit = 'lo2D';
            end
            
            % reset the calibrated gain
            obj.slopesUnits = 1;
            
            
            nPx = obj.camera.resolution(1)/obj.lenslets.nLenslet;
            d   = tel.D/size(obj.validLenslet,1);
                        ngs = ngs.*tel*obj;
            obj.pointingDirection = zeros(2,1);
            
            ps = obj.lenslets.pixelScale(ngs,tel);
            fprintf('The current lenslet pixel-scale is %f mas\n',ps.convert('mas')) % sampling of the electic field
            binFactor = 2*obj.lenslets.fieldStopSize/obj.lenslets.nyquistSampling/nPx;
            lo2DInMas = ngs.wavelength/(2*d)*constants.radian2mas;
            lo2D = ngs.wavelength/(2*d);
            detectorPixelSizeInMas = lo2DInMas*binFactor;
            detectorPixelSize = lo2D*binFactor;
            fprintf('The current detector pixel-scale is %f mas\n',detectorPixelSizeInMas); % sky-angle of the detector pixels after binning
            
            
            switch unit
                case 'lo2D'
                    
                    % the next three lines of code provide better LTAO performance
                    % at the centre of the field in the case of KAPA (+7% in H for ZA=50, median MK profile!).
                    pixelScale = obj.lenslets.fieldStopSize*ngs.wavelength/d/nPx;
                    fprintf('The current delta %f mas\n',pixelScale*constants.radian2mas);
                    tipStep = pixelScale/2/2;
                    nStep   = floor(nPx/3)*2;
                    sx      = zeros(1,nStep+1);
                    u       = 0:nStep;
                    obj.camera.frameListener.Enabled = false;
                    obj.slopesListener.Enabled = false;
                    
                    warning('off','oomao:shackHartmann:relay')
                    for kStep=u
                        ngs.zenith = -tipStep*kStep;
                        +ngs;
                        drawnow
                        sx(kStep+1) = median(obj.slopes(1:end/2));
                    end
                    warning('on','oomao:shackHartmann:relay')
                    
                    Ox_in  = u*tipStep*constants.radian2arcsec;
                    Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
                    figure
                    plot(Ox_in, Ox_out)
                    hold
                    plot(Ox_in, Ox_in,'k:')
                    xlabel('input')
                    ylabel('output')
                    %plot(u*tipStep/pixelScale,(u*tipStep/pixelScale)./sx);
                    
                    slopesLinCoef = polyfit(Ox_in,Ox_out,1);
                    title(['gain = ' num2str(1/slopesLinCoef(1))])
                    obj.slopesUnits = 1/slopesLinCoef(1);
                    ngs.zenith = 0;
                case 'pixel'
                    % opticalGain fine tuning :: make it exactly 1 at a delta of 1 detector pixel
                    zern = zernike(2,tel.D,'resolution',tel.resolution,'pupil',tel.pupil);
                    zern.c = detectorPixelSizeInMas/constants.radian2mas*tel.D/4;% tel.D/4*ngs.wavelength/(2*tel.D);
                    ngs = ngs.*zern*obj;
                    obj.slopesUnits = 1/median(obj.slopes(1:end/2));% wfs.slopesUnits/wfs.slopes(1)*pixelSize/lo2DInMas;
            end
            
            % reset the wfs pointing
            obj.pointingDirection = [];
            
        end
        
        %         function gainCalibration(obj,tel,ngs)
        %             %% GAINCALIBRATION
        %             % calibrate gain of cetner of gravity
        %             %
        %             % wfs.gainCalibration(tel,ngs)
        %             %
        %             % the WFS output is in units of \lambda/D/2, i.e. half the
        %             % diffraction limit for the aperture/sub-aperture
        %             %
        %             % When the pixel is \lambda/D (i.e. 2x \lambda/D/2) the gain
        %             % becomes ~2. For higher binning factors the scaling factor
        %             % grows in a non-linear fashion
        %
        %
        %             nPx = obj.camera.resolution(1)/size(obj.validLenslet,1);
        %             d   = tel.D/size(obj.validLenslet,1);
        %
        %             ngs = ngs.*tel*obj;
        %             obj.pointingDirection = zeros(2,1);
        %
        %
        %             binFactor = max(obj.lenslets.nLensletsImagePx/obj.camera.resolution(1),1);
        %             pixelScale = ngs.wavelength/d/2/obj.lenslets.nyquistSampling*binFactor;
        %             fprintf('Pixel Scale: %12.2f mas\n', pixelScale*constants.radian2mas)
        %             nStep = 10;
        %             tipStep = pixelScale/nStep/2/binFactor; % on-sky off-set up to 1/2 pixel/binFactor (division by bin factor is to "avoid" non-linear operation)
        %
        %             % the next three lines of code provide better LTAO performance
        %             % at the centre of the field in the case of KAPA (+7% in H for ZA=50, median MK profile!).
        %             pixelScale = obj.lenslets.fieldStopSize*ngs.wavelength/nPx;
        %             tipStep = pixelScale/2/2;
        %             nStep   = floor(nPx/3)*2;
        %
        %             sx      = zeros(1,nStep+1);
        %             u       = 0:nStep;
        %             obj.camera.frameListener.Enabled = false;
        %             obj.slopesListener.Enabled = false;
        %             u*tipStep*constants.radian2mas
        %             warning('off','oomao:shackHartmann:relay')
        %             for kStep=u
        %                 ngs.zenith = -tipStep*kStep;
        %                 +ngs;
        %                 drawnow
        %                 sx(kStep+1) = median(obj.slopes(1:end/2));
        %             end
        %             warning('on','oomao:shackHartmann:relay')
        %
        %             Ox_in  = u*tipStep*constants.radian2arcsec;
        %             Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
        %             figure
        %             plot(Ox_in, Ox_out)
        %             hold
        %             plot(Ox_in, Ox_in,'k:')
        %             xlabel('input')
        %             ylabel('output')
        %             %plot(u*tipStep/pixelScale,(u*tipStep/pixelScale)./sx);
        %
        %             slopesLinCoef = polyfit(Ox_in,Ox_out,1);
        %             title(['gain = ' num2str(1/slopesLinCoef(1))])
        %             obj.slopesUnits = 1/slopesLinCoef(1);
        %             ngs.zenith = 0;
        %             obj.pointingDirection = [];
        %
        %
        %         end
    end
    
    methods (Static)
        
        function obj = loadobj(obj)
            %% LOADOBJ
            add(obj.log,obj,'Load!')
            setSlopesListener(obj)
            obj.log = logBook.checkIn(obj);
        end
        
                %% get reconstruction grid - used in caes where a phase is reconstructed where the actuators are located, or a oversampling is applied
                function val = reconstructionGrid(obj,os)
                % os stands for overSampling. Can be 1 or 2. 
                % If os=1, reconstructionGrid pitch = subaperturePitch, 
                % if os=2, reconstructionGrid pitch = subaperturePitch/2
            if nargin == 0
                val = get.validActuator(obj);
            elseif os ==2
                nElements            = os*obj.lenslets.nLenslet+1; % Linear number of lenslet+actuator
                validLensletActuator = zeros(nElements);
                index                = 2:2:nElements; % Lenslet index
                validLensletActuator(index,index) = obj.validLenslet;
                for xLenslet = index
                    for yLenslet = index
                        if validLensletActuator(xLenslet,yLenslet)==1
                            xActuatorIndice = ones(3,1)*(xLenslet-1:xLenslet+1);
                            yActuatorIndice =ones(3,1)*(yLenslet-1:yLenslet+1) ;
                            validLensletActuator(xActuatorIndice,yActuatorIndice) = 1;
                        end
                    end
                end
                val = logical(validLensletActuator);
            end
        end
    end
    
    methods (Access=private)
        
        function setSlopesListener(obj)
            %% SETSLOPESLISTENER Slopes listener
            obj.slopesListener = addlistener(obj,'slopes','PostSet',...
                @(src,evnt) obj.slopesDisplay );
            obj.slopesListener.Enabled = false;
        end
        
    end
    
end

function y = linearSpline(x)
%% LINEARSPLINE Linear spline function
%
% y = linearSpline(x) computes the function y = 1 - |x| for |x|<1 and y = 0
% elsewhere

[m,n] = size(x);
x     = abs(x);
index = x < 1;
[i,j] = find(index);
s     = 1 - x(index);
y = sparse(i,j,s,m,n);
% y = zeros(size(x));
% y(index) = 1 - x(index);

end


function y = linearSplineInt(x)
%% LINEARSPLINEINT Linear spline integral
%
% y = linearSplineInt(x) computes the function y = -(x-sign(x))^2/(2sign(x))

y = -(x-sign(x)).^2./(2.*sign(x));

end