classdef detector < handle
    % DETECTOR Create a detector object
    %
    % obj = detector(resolution) creates a detector object from the detector
    % resolution
    
    properties
        % Add or not photon noise to the frame        
        photonNoise = false;
        % Photon background noise (# of photon per frame)
        nPhotonBackground = 0;
        % # of photo-electron per pixel rms
        readOutNoise = 0;
        % quantum efficiency
        quantumEfficiency = 1;
        % dark background (scalar or 2D matrix)
        darkBackground = 0;
        %pixelGains
        pixelGains = 1;
        % units of one pixel
        pixelScale;
        % detector resolution
        resolution;
        % detector region of interest
        regionOfInterest;
        roiSouthWestCorner;
        % frame rate [Hz]
        frameRate = 1;
        % clock rate [Hz]
        clockRate = 1;
        %  Delay after which the camera start integrating
        startDelay = 0;
        % detector timer
        paceMaker;
        % frame update listener
        frameListener;
        % frame grabber callback function
        frameGrabber;
        % detector tag
        tag= 'DETECTOR';
        frameBuffer = 0;
        frameCount  = 0;
        % LGS elongated spot kernel
        spotsLgsSrcKernel;
        % LGS elongated spot kernel in Fourier domain
        fftSpotsLgsSrcKernel;
        % detector binningFactor
        binningFactor;
        %spatial offset the pixel intensity output 
        offSetSpotsLgsSrcKernel = [];
    end
    
    properties (Dependent)
        % Exposure time [second]
        exposureTime;
    end
    
    properties (SetObservable=true)
        % detector frame
        frame;
    end
    
    properties (Access=protected)
        frameHandle;
        log;
    end
    
    properties (Access=private)
        p_exposureTime;
    end
    
    methods
        
        %% Constructor
        function obj = detector(resolution,pixelScale)
            narginchk(1, 2)
            if numel(resolution)==1
                obj.resolution = resolution*ones(1,2);
            else
                obj.resolution = resolution;
            end
            obj.regionOfInterest   = obj.resolution;
            obj.roiSouthWestCorner = [1,1];
            if nargin>1
                obj.pixelScale  = pixelScale;
            end
            obj.exposureTime = 1;

            setFrameListener(obj)
            
            % Timer settings
            obj.paceMaker = timer;
            obj.paceMaker.name = 'Detector';
            obj.paceMaker.TimerFcn = @(src,evnt) obj.grab;% {@timerCallBack, obj};
            obj.paceMaker.ExecutionMode = 'FixedSpacing';
            obj.paceMaker.BusyMode = 'error';
            obj.paceMaker.Period = 1;
            obj.paceMaker.ErrorFcn = 'disp('' @detector: frame rate too high!'')';
%             function timerCallBack( timerObj, event, a)
%                 %                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
%                 a.grab;
%             end
            %             obj.frameRate = 1;
            obj.log = logBook.checkIn(obj);
        end
        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.frameHandle) && ishandle(obj.frameHandle(1))
                delete(get(obj.frameHandle(1),'Parent'));
            end
            if isvalid(obj.paceMaker)
                if strcmp(obj.paceMaker.Running,'on')
                    stop(obj.paceMaker)
                end
                delete(obj.paceMaker)
            end
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % disp(obj) prints information about the detector object
          
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' %dx%d pixels camera \n',...
                obj.resolution)
            if ~isempty(obj.pixelScale)
            fprintf('  . pixel scale: %4.2f milli-arcsec \n',...
                obj.pixelScale*constants.radian2arcsec*1000)                
            end            
            fprintf('  . quantum efficiency: %3.1f \n',...
                obj.quantumEfficiency)
            if obj.photonNoise
                fprintf('  . photon noise enabled\n')
            else
                fprintf('  . photon noise disabled\n')
            end
            fprintf('  . %.1f photo-events rms read-out noise \n',...
                obj.readOutNoise)
            fprintf('  . %3.1fms exposure time and %3.1fHz frame rate \n',...
                obj.exposureTime*1e3,obj.frameRate)
            fprintf('----------------------------------------------------\n')
            
        end
        
        function obj = saveobj(obj)
            %% SAVEOBJ
            delete(obj.frameListener)
            add(obj.log,obj,'Save!')
        end                
        
        %% Set/Get spotsLgsSrcKernel
        function set.spotsLgsSrcKernel(obj,val)
            obj.spotsLgsSrcKernel = val;
%             if obj.clockRate==1
%                 obj.clockRate = 1/obj.p_exposureTime;
%                 fprintf('Clock rate is: %.2fHz\n',obj.clockRate)
%             end
        end
        function val = get.spotsLgsSrcKernel(obj)
            val = obj.spotsLgsSrcKernel;
        end
        
        %% Set/Get extendedResolution
        function set.binningFactor(obj,val)
            obj.binningFactor = val;
%             if obj.clockRate==1
%                 obj.clockRate = 1/obj.p_exposureTime;
%                 fprintf('Clock rate is: %.2fHz\n',obj.clockRate)
%             end
        end
        function val = get.binningFactor(obj)
            val = obj.binningFactor;
        end
        
        %% Set/Get exposureTime
        function set.exposureTime(obj,val)
            obj.p_exposureTime = val;
%             if obj.clockRate==1
%                 obj.clockRate = 1/obj.p_exposureTime;
%                 fprintf('Clock rate is: %.2fHz\n',obj.clockRate)
%             end
        end
        function val = get.exposureTime(obj)
            val = obj.p_exposureTime;
        end
        
        %% Displays detector frame
        function imagesc(obj,varargin)
            % IMAGESC Display the detector frame
            %
            % imagesc(obj) displays the frame of the detector object
            %
            % imagesc(obj,'PropertyName',PropertyValue) displays the frame of
            % the detector object and set the properties of the graphics object
            % imagesc
            %
            % h = imagesc(obj,...) returns the graphics handle
            %
            % See also: imagesc
            
            if ismatrix(obj.frame) % frame must be a 2-D array for plotting
                
                if isempty(obj.frame)
                    m_frame = obj.frameBuffer;
                    titleColor = 'r';
                else
                    m_frame = obj.frame;
                    [n,m] = size(m_frame);
                    try
                    if m>2*n
                        m_frame = cell2mat(reshape( mat2cell( ...
                            m_frame,n,n*ones(1,m/n) ) , 2, []));
                    end
                    end
                    titleColor = 'k';
                end
                
                if all(ishandle(obj.frameHandle)) && ~isempty(obj.frameHandle)
                    nm = size(m_frame);
                    if ~all(size(get(obj.frameHandle(1),'Cdata'))==nm)
                        set(get(obj.frameHandle(1),'parent'),'xlim',0.5+[0 nm(2)],'ylim',0.5+[0 nm(1)])
                    end
                    set(obj.frameHandle(1),'Cdata',m_frame,varargin{:});
                    set(obj.frameHandle(2),'String',...
                        sprintf('Frame #%d - Exposure: %gs',obj.frameCount,obj.exposureTime),...
                        'color',titleColor)
                else
                    if isempty(obj.pixelScale)
                        obj.frameHandle(1) = image(m_frame,...
                            'CDataMApping','Scaled',...
                            varargin{:});
                    else
                        [n,m] = size(m_frame);
                        x = 0.5*linspace(-1,1,n)*...
                            n*obj.pixelScale*constants.radian2arcsec;
                        y = 0.5*linspace(-1,1,m)*...
                            m*obj.pixelScale*constants.radian2arcsec;
                        obj.frameHandle(1) = image(x,y,m_frame,...
                            'CDataMApping','Scaled',...
                            varargin{:});
                    end
                    obj.frameHandle(2) = title(sprintf('Frame #%d - Exposure: %gs',obj.frameCount,obj.exposureTime),...
                        'color',titleColor);
                    colormap(pink)
                    axis xy equal tight
                    colorbar('location','EastOutside')
                    hu = findobj(gcf,'Type','uimenu','Label','OOMAO');
                    if isempty(hu)
                        hu = uimenu('Label','OOMAO');
                    end
                    hus  = uimenu(hu,'Label','Frame Listener Off','Callback',@oomaoMenu);
                    if obj.frameListener.Enabled
                        set(hus,'Label','Frame Listener On')
                    end
                end
            
            end
            function oomaoMenu(src,~)
                obj.frameListener.Enabled = ~obj.frameListener.Enabled;
                if obj.frameListener.Enabled
                    set(src,'Label','Frame Listener On')
                else
                    set(src,'Label','Frame Listener Off')
                end
            end
        end
        
        function varargout = grab(obj,~)
            %% GRAB Frame grabber
            %
            % grab(obj) grabs a frame
            %
            % out = grab(obj) grabs a frame and returns it
            
            switch class(obj.frameGrabber)
                case {'lensletArray','gpuLensletArray'}
                    readOut(obj,obj.frameGrabber.imagelets)
                case {'pyramid','metaPyramid'}
                    readOut(obj,obj.frameGrabber.lightMap)
                case 'function_handle'
                    buffer = obj.frameGrabber();
                    [n,m] = size(buffer);
                    u = obj.roiSouthWestCorner(1):...
                        min(obj.roiSouthWestCorner(1)+obj.regionOfInterest(1)-1,n);
                    v = obj.roiSouthWestCorner(2):...
                        min(obj.roiSouthWestCorner(2)+obj.regionOfInterest(2)-1,m);
                    obj.frame = buffer(u,v);
                otherwise
            end
            if nargout>0
                varargout{1} = obj.frame;
            end
        end
        
        function relay(obj,src)
            
            % Here we check the last object the source went through before
            % the detector
            srcLastPath = src.opticalPath{end-1};
            switch class(srcLastPath) 
                case 'telescope'
                    f = tools.cartAndPol(obj.resolution(1),...
                        'output','radius');
                    % pixel scale in radian
                    f = obj.pixelScale*f.*(obj.resolution(1)-1)./src.wavelength/2;
                    obj.frame = psf(srcLastPath,f);
                otherwise
                    if src.timeStamp>=obj.startDelay
                        obj.startDelay = -Inf;
                        obj.frameBuffer = obj.frameBuffer + src.intensity;
                        if src.timeStamp>=obj.exposureTime
                            src.timeStamp = 0;
                            %disp(' @(detector:relay)> reading out and emptying buffer!')
                            readOut(obj,obj.frameBuffer)
                            obj.frameBuffer = 0*obj.frameBuffer;
                        end
                    end
            end
            
        end
        
    end
    
    methods (Access=protected)
        
        function readOut(obj,image)
            %% READOUT Detector readout
            %
            % readOut(obj,image) adds noise to the image: photon noise if
            % photonNoise property is true and readout noise if
            % readOutNoise property is greater than 0
            
            %             image = image;%This is now done in telescope.relay (.*obj.exposureTime;) % flux integration
            
            % convolution by extended LGS kernel
            if ~isempty(obj.spotsLgsSrcKernel) || ~isempty(obj.fftSpotsLgsSrcKernel)
                [n,m,p] = size(image);
                nLenslet = sqrt(size(obj.spotsLgsSrcKernel, 3));
                nPxDetector = obj.resolution(1) / nLenslet;
                if mod(nPxDetector,2)
                    compShift = 0;
                else
                    compShift = 0;
                end
                if length(size(image)) == 2 % regular case
                    nArray = m/n;
                elseif length(size(image)) == 3 % calibration matrix case
                    nArray = p;
                end
                
                for iArray = 1:nArray
                    if length(size(image)) == 2 % regular case
                        tmp = image(:, (iArray-1)*n + [1:n]);
                    elseif length(size(image)) == 3 % DM calibration matrix case
                        tmp = image(:,:, iArray);
                    end
                    
                    [nPx,mPx,nPicture] = size(tmp);
                    nPxLenslet = nPx/nLenslet;
                    mPxLenslet = mPx/nLenslet;
                    mPxNPicture = mPx*nPicture;
                    tmp = reshape(tmp,nPx,mPxNPicture);
                    nLensletArray = nArray*nPicture;
                    
                    indexRasterLenslet_ ...
                        = tools.rearrange(size(tmp),[nPxLenslet,mPxLenslet]);
                    %                 v = ~obj.validLenslet(:);
                    %                 v = repmat(v,nLensletArray,1);
                    %                 indexRasterLenslet_(:,v) = [];
                    buffer     = tmp(indexRasterLenslet_);
                    buffer     = reshape(buffer,nPxLenslet,nPxLenslet,[]);
                    %buffer = tools.crop(buffer, nPxDetector);
                    bufferDetect = zeros(nPxDetector, nPxDetector, nLenslet^2);
                    
                    % create local variables to avoid object deletion in
                    % parfor loop
                    binningFactor = obj.binningFactor;
                    offSetSpotsLgsSrcKernel = obj.offSetSpotsLgsSrcKernel;
                    spotsLgsSrcKernel = obj.spotsLgsSrcKernel;
                    fftSpotsLgsSrcKernel = obj.fftSpotsLgsSrcKernel;
                    
                    parfor ii = 1:nLenslet^2
                        if ~isempty(offSetSpotsLgsSrcKernel)
                            % shifts from geometric SH
                            subap = buffer(:,:,ii);
                            nPhotons = sum(subap(:));
                            if length(size(image)) == 2 % regular case
                                bufferDetect(:,:,ii) = tools.crop(tools.shift(spotsLgsSrcKernel(:,:,ii,iArray),...
                                    offSetSpotsLgsSrcKernel(ii,iArray)+compShift / binningFactor,...
                                    offSetSpotsLgsSrcKernel(ii+nLenslet^2,iArray)+compShift / binningFactor), nPxDetector);
                                % check compShift / binningFactor
                                bufferDetect(:,:,ii) = bufferDetect(:,:,ii) .* (bufferDetect(:,:,ii)>=0);
                            elseif length(size(image)) == 3 % for calibration matrix
                                bufferDetect(:,:,ii) = tools.crop(tools.shift(spotsLgsSrcKernel(:,:,ii),...
                                    offSetSpotsLgsSrcKernel(ii)+compShift / binningFactor,...
                                    offSetSpotsLgsSrcKernel(ii+nLenslet^2)+compShift / binningFactor), nPxDetector);
                                bufferDetect(:,:,ii) = bufferDetect(:,:,ii) .* (bufferDetect(:,:,ii)>=0);
                            end
                            
                            bufferDetect(:,:,ii) = bufferDetect(:,:,ii) / sum(sum(bufferDetect(:,:,ii))) * nPhotons;
                        else
                            convMethod = 'conv';
                            if length(size(image)) == 2 % regular case   
                                switch convMethod
                                    case 'truncatedFFT'
%                                         convTime = ones(1000, 1);
%                                         for iTime = 1:1000
%                                         a = tic;
                                        kTilda = fftSpotsLgsSrcKernel(:,:,ii, iArray);
                                        h = buffer(:,:,ii);
                                        fftSize = length(kTilda)*binningFactor;
                                        hTilda = tools.crop(fftshift(fft2(tools.crop(h, fftSize))), length(kTilda));
                                        % hTilda = tools.crop(fftshift(fft2(tools.crop(h, length(kTilda)*binningFactor))), length(kTilda));
                                        subap = fftshift(abs(ifft2(hTilda .* kTilda)));%, fftSize, fftSize));
                                        bufferDetect(:,:,ii) = tools.crop(subap, nPxDetector);
                                        % subapTrunc = bufferDetect(:,:,ii);
%                                         convTime(iTime) = toc(a);
%                                         end
%                                         mean(convTime)*1000
                                        
                                    case 'FFT'                                        
                                        % TO BE DEBUGGED
                                        binnedKernel = padArray2(spotsLgsSrcKernel(:,:,ii, iArray), binningFactor, binningFactor) / binningFactor^2; % kernel is expanded to lenslets pixelScale
                                        convKernel = abs(fftshift(ifft2(fft2(buffer(:,:,ii)) .* fft2(binnedKernel, nPxLenslet, nPxLenslet))));% to be implemented
                                        nCrop = nPxDetector * binningFactor;
                                        convKernel = tools.crop(convKernel, nCrop(1)); % the convolved spot is binned to the WFS pixelScale
                                        bufferDetect(:,:,ii) = tools.crop(tools.binning(convKernel, size(convKernel) / binningFactor), nPxDetector);
                                        %subapFFT = bufferDetect(:,:,ii);
                                    case 'conv'
%                                          convTime = ones(1000, 1);
%                                         for iTime = 1:1000
%                                             a = tic;
                                        binnedKernel = padArray2(spotsLgsSrcKernel(:,:,ii, iArray), binningFactor, binningFactor) / binningFactor^2; % kernel is expanded to lenslets pixelScale
                                        if length(buffer(:,:,ii)) >= length(binnedKernel)
                                            convKernel = conv2(buffer(:,:,ii), binnedKernel, 'same'); % convolution is performed at the lenslet pixelScale to preserve the lenslet PSF resolution;
                                        else
                                            convKernel = conv2(binnedKernel, buffer(:,:,ii), 'same');
                                        end
                                        % convKernel = tools.shift(convKernel, (0.5-compShift) * binningFactor, (0.5-compShift) * binningFactor);
                                        %convKernel = tools.shift(convKernel, -4/binningFactor, -4/binningFactor);
                                        nCrop = nPxDetector * binningFactor;
                                        convKernel = tools.crop(convKernel, nCrop(1)); % the convolved spot is binned to the WFS pixelScale
                                        bufferDetect(:,:,ii) = tools.crop(tools.binning(convKernel, size(convKernel) / binningFactor), nPxDetector);
                                        %subapConv = bufferDetect(:,:,ii);
%                                           convTime(iTime) = toc(a);
%                                         end
%                                         mean(convTime)*1000
                                end
                                
                            elseif length(size(image)) == 3 % for calibration matrix
                                
                                switch convMethod
                                    case 'truncatedFFT'
                                        kTilda = fftSpotsLgsSrcKernel(:,:,ii);
                                        h = buffer(:,:,ii);
                                        fftSize = length(kTilda)*binningFactor;
                                        hTilda = tools.crop(fftshift(fft2(tools.crop(h, fftSize))), length(kTilda));
                                        % hTilda = tools.crop(fftshift(fft2(tools.crop(h, length(kTilda)*binningFactor))), length(kTilda));
                                        subap = fftshift(abs(ifft2(hTilda .* kTilda)));%, fftSize, fftSize));
                                        bufferDetect(:,:,ii) = tools.crop(subap, nPxDetector);
                                        
                                    case 'FFT'
                                        % TO BE DEBUGGED
                                        binnedKernel = padArray2(spotsLgsSrcKernel(:,:,ii), binningFactor, binningFactor) / binningFactor^2; % kernel is expanded to lenslets pixelScale
                                        convKernel = abs(fftshift(ifft2(fft2(buffer(:,:,ii)) .* fft2(binnedKernel, nPxLenslet, nPxLenslet))));% to be implemented
                                        nCrop = nPxDetector * binningFactor;
                                        convKernel = tools.crop(convKernel, nCrop(1)); % the convolved spot is binned to the WFS pixelScale
                                        bufferDetect(:,:,ii) = tools.crop(tools.binning(convKernel, size(convKernel) / binningFactor), nPxDetector);
                                        % subapFFT = bufferDetect(:,:,ii);
                                    case 'conv'
                                        binnedKernel = padArray2(spotsLgsSrcKernel(:,:,ii), binningFactor, binningFactor) / binningFactor^2; % kernel is expanded to lenslets pixelScale
                                        convKernel = conv2(buffer(:,:,ii), binnedKernel); % convolution is performed at the lenslet pixelScale to preserve the lenslet PSF resolution;
                                        % convKernel = tools.shift(convKernel, (0.5-compShift) * binningFactor, (0.5-compShift) * binningFactor);
                                        nCrop = nPxDetector * binningFactor;
                                        convKernel = tools.crop(convKernel, nCrop(1)); % the convolved spot is binned to the WFS pixelScale
                                        bufferDetect(:,:,ii) = tools.crop(tools.binning(convKernel, size(convKernel) / binningFactor), nPxDetector);
                                        % subapConv = bufferDetect(:,:,ii);
                                end
                                
%                                 binnedKernel = padArray2(spotsLgsSrcKernel(:,:,ii), binningFactor, binningFactor) / binningFactor^2;
%                                 convKernel = conv2(buffer(:,:,ii), binnedKernel);
%                                 %convKernel = tools.shift(convKernel, (0.5-compShift) * binningFactor, (0.5-compShift) * binningFactor);
%                                 nCrop = nPxDetector * binningFactor;
%                                 convKernel = tools.crop(convKernel, nCrop(1));
%                                 bufferDetect(:,:,ii) = tools.crop(tools.binning(convKernel, size(convKernel) / binningFactor), nPxDetector);
                            end
                        end
                    end
                    
                    
                    indexRasterDetector = tools.rearrange(obj.resolution,[nPxDetector,nPxDetector]);
                    tmpDetector(indexRasterDetector) = bufferDetect;
                    tmpDetector = reshape(tmpDetector , obj.resolution(1) , obj.resolution(2) , nPicture);
                    if length(size(image)) == 2
                        output(:, (iArray-1)*obj.resolution(1) + [1:obj.resolution(1)]) = tmpDetector;
                    elseif length(size(image)) == 3
                        output(:,:, iArray) = tmpDetector;
                    end
                end
                image = output;
            end
            % end of convolution
            
            [n,m,~] = size(image);
            if any(n>obj.resolution(1))
                %                 disp('Binning')
                image = tools.binning(image,obj.resolution.*[n,m]/n);
            end
            
            
            %             % convolution by extended LGS kernel
            %             if ~isempty(obj.spotsLgsSrcKernel)
            %                 [n,m,~] = size(image);
            %                 nLenslet = sqrt(size(obj.spotsLgsSrcKernel, 3));
            %                 nPxDetector = obj.resolution(1) / nLenslet;
            %                 nArray = m/n;
            %                 for iArray = 1:nArray
            %                     tmp = image(:, (iArray-1)*n + [1:n]);
            %                     [nPx,mPx,nPicture] = size(tmp);
            %                     nPxLenslet = nPx/nLenslet;
            %                     mPxLenslet = mPx/nLenslet;
            %                     mPxNPicture = mPx*nPicture;
            %                     tmp = reshape(tmp,nPx,mPxNPicture);
            %                     nLensletArray = nArray*nPicture;
            %
            %                     indexRasterLenslet_ ...
            %                         = tools.rearrange(size(tmp),[nPxLenslet,mPxLenslet]);
            %                     %                 v = ~obj.validLenslet(:);
            %                     %                 v = repmat(v,nLensletArray,1);
            %                     %                 indexRasterLenslet_(:,v) = [];
            %                     buffer     = tmp(indexRasterLenslet_);
            %                     buffer     = reshape(buffer,nPxLenslet,nPxLenslet,[]);
            %                     buffer = tools.crop(buffer, nPxDetector);
            %                     for ii = 1:nLenslet^2
            %                         buffer(:,:,ii) = tools.crop(conv2(buffer(:,:,ii), obj.spotsLgsSrcKernel(:,:,ii, iArray), 'same'), nPxDetector);
            %                     end
            %                     indexRasterDetector = tools.rearrange(obj.resolution,[nPxDetector,nPxDetector]);
            % %                     tmp(indexRasterLenslet_) = buffer;
            % %                     tmp = reshape(tmp , nPx , mPx , nPicture);
            %                     tmpDetector(indexRasterDetector) = buffer;
            %                     tmpDetector = reshape(tmpDetector , obj.resolution(1) , obj.resolution(2) , nPicture);
            %                     output(:, (iArray-1)*obj.resolution(1) + [1:obj.resolution(1)]) = tmpDetector;
            %                     %end
            %                 end
            %                 image = output;
            %             end
            %             % end of convolution
            
            obj.frameCount = obj.frameCount + 1;
            if obj.frameCount<obj.p_exposureTime*obj.clockRate && obj.frameCount>obj.startDelay*obj.clockRate
                obj.frameCount = obj.frameCount - obj.startDelay*obj.clockRate;
                obj.startDelay = 0;
                obj.frameBuffer = obj.frameBuffer + image;
                obj.frame = [];
            else
                image = obj.frameBuffer + image;
                
                if license('checkout','statistics_toolbox')
                    if obj.photonNoise
                        image = poissrnd(image + obj.nPhotonBackground) - obj.nPhotonBackground;
                    end
                    image = obj.quantumEfficiency*image;
                    if obj.readOutNoise>0
                        image = image + randn(size(image)).*obj.readOutNoise;
                    end
                else
                    if obj.photonNoise
                        buffer    = image + obj.nPhotonBackground;
                        image = image + randn(size(image)).*(image + obj.nPhotonBackground);
                        index     = image<0;
                        image(index) = buffer(index);
                        image = obj.quantumEfficiency*image;
                    end
                    image = obj.quantumEfficiency*image;
                    if obj.readOutNoise>0
                        image = image + randn(size(image)).*obj.readOutNoise;
                    end
                end
                if any(obj.darkBackground > 0)
                    image = image + obj.darkBackground;
                end
                if any(obj.pixelGains > 0)
                    image = image.*obj.pixelGains;
                end
                
                if obj.frameCount>obj.startDelay*obj.clockRate
                    obj.frameCount = 0;
                end
                obj.frame = image;
                obj.frameBuffer = 0;
            end
        end
        
    end
    
    methods (Static)
            
        function obj = loadobj(obj)
            %% LOADOBJ
            add(obj.log,obj,'Load!')
            setFrameListener(obj)
            obj.log = logBook.checkIn(obj);
        end
        
    end
    
    methods (Access=private)
        
        function setFrameListener(obj)
            %% SETFRAMELISTENER % Frame listener
            obj.frameListener = addlistener(obj,'frame','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.frameListener.Enabled = false;
        end
        
    end
    
end