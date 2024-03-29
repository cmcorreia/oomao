classdef lensletArray < handle
% LENSLETARRAY Create a lenslet array object
%
% obj = lensletArray(nLenslet) creates a lenslet array object from the
% number of lenslet on one side of the array

    properties
        % the # of lenslet on one side of the lenslet array
        nLenslet;
        % lenslet size
        pitch;
        % the minimum amount of light per lenslet
        minLightRatio;
        % lenslets conjugation altitude
        conjugationAltitude = Inf;
        % the input wave
        wave;
        % imagelets listener
        imageletsListener
        % stacked imagelets sum
        sumStack = false;
        % the number of lenslet arrays
        nArray = 1;
        % optical throughput
        throughput=1;
        % lenslet array tag
        tag='LENSLET ARRAY';
        
        % extra field of view given in diffraction fwhm units for elongated
        % spots convolution. The sky-psf on the lenlslet focal plane is
        % extended using zero-padding, which is an approximation
        elongatedFieldStopSize = [];
        convKernel = false;

        % amount of opticalAberration on each lenslet to produce an aberrated spot
        opticalAberration;
        
    end
    
    properties (SetObservable=true)
        % the image at each lenslet focus
        imagelets;
    end
    
    properties (SetAccess=private)
        fftPad;        
    end

    properties (Dependent)
        % the # of pixel per lenslet
        nLensletWavePx;
        % the Nyquist sampling factor
        nyquistSampling;
        % the lenslet field of view given in diffraction fwhm units
        fieldStopSize;
        % offsets in fraction of the telescope diameter (i.e. units of 1/D)
        offset;
        % rotation in radians
        rotation
    end
    
    properties (Dependent, SetAccess=private)
        % the # of pixel per imagelet
        nLensletImagePx;
        % the total # of pixel per imagelet
        nLensletsImagePx;

    end

    properties (Access=private)
        p_nLensletWavePx;
        p_nyquistSampling;
        p_fieldStopSize;
        fftPhasor;
        imageHandle;
        zoomTransmittance;
        log;
        % vector of offsets
        p_offset = [0;0;0];
    end
    methods

        %% Constructor
        function obj = lensletArray(nLenslet)
            narginchk(1, 3)
            obj.nLenslet        = nLenslet;
            obj.minLightRatio   = 0;
            obj.nyquistSampling = 1;
            
            setImageletsListener(obj)
            
            obj.log = logBook.checkIn(obj);
        end

        %% Destructor
        function delete(obj)
            if ishandle(obj.imageHandle)
                delete(get(obj.imageHandle,'parent'));
            end
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % disp(obj) prints information about the lenslet array object
          
            fprintf('___ %s ___\n',obj.tag)
            if obj.nArray>1
                fprintf(' %d %dx%d lenslet array: \n',...
                    obj.nArray,obj.nLenslet*ones(1,2))
            else
                fprintf(' %dx%d lenslet array: \n',...
                    obj.nLenslet*ones(1,2))
            end
            fprintf('  . %3.1f pixels across the diffraction limited spot fwhm\n',...
                obj.nyquistSampling*2);
            fprintf('  . %d pixels across the square lenslet field stop size\n',...
                obj.fieldStopSize*obj.nyquistSampling*2);            
            fprintf('  . optical throughput coefficient: %3.1f\n',...
                obj.throughput);            
            fprintf('----------------------------------------------------\n')
            
        end
        
        function obj = saveobj(obj)
            %% SAVEOBJ
            delete(obj.imageletsListener)
            add(obj.log,obj,'Save!')
        end        

        %% Lenslet image sampling
        function nLensletImagePx = get.nLensletImagePx(obj)
            nLensletImagePx = ceil(obj.fieldStopSize.*obj.nyquistSampling*2);
%             if ~isempty(obj.elongatedFieldStopSize)
%                 nLensletImagePx = ceil((obj.fieldStopSize + obj.elongatedFieldStopSize)  .*obj.nyquistSampling*2);
%             end
        end

        %% Whole lenslet array image sampling
        function nLensletsImagePx = get.nLensletsImagePx(obj)
            nLensletsImagePx = obj.nLenslet.*obj.nLensletImagePx;
        end

        %% Get and Set rotation
        function set.rotation(obj,val)
            if isempty(val)
                obj.p_offset(3,:) = 0;
            else
                obj.p_offset(3,length(val)) = 0;
                obj.p_offset(3,:) = val;
            end
            %if length(obj.rotation) == 1 && length(obj.rotation) < obj.nSrc
            %    obj.rotation = repmat(obj.rotation,1,obj.nSrc);
            %end
        end
        
        function out = get.rotation(obj)
            out = obj.p_offset(3,:);
        end
        %% Get and Set Offsets
        function set.offset(obj,val)
            % check offset's magnitude
%             if any(abs(val)) > 1
%                 warning('Offsets larger than a sub-aperture not supported. Truncating the offsets to 1')
%                 val(val > 1) = 1;
%                 val(val < -1) = -1;
%             end
            % check offset's size consistency
            if any(size(val)) == 1
                if size(val,1) >  size(val,2)
                    val = val';
                end
            end
            if isscalar(val) % create a 2x1 vector
                offset_ = [1;1]*val;
            elseif size(val,2) >  size(val,1)
                if size(val,2) == 2 % if a 1x2 vector, transpose
                    offset_ =val';
                elseif size(val,1) == 1 % if a 1 x N vector, replicate along the other direction
                    offset_ = repmat(val, 2,1);
                else
                    offset_ = val;
                end
            else
                offset_ = val;
            end
            obj.p_offset(1:2,:) = offset_;
            rotation_ =  obj.p_offset(3,:);
            if any(rotation_) %rotation already set
                if length(rotation_) < size(obj.p_offset,2)
                     obj.p_offset(3,size(obj.p_offset,2)) = 0;
                end
            else
                 obj.p_offset(3,size(obj.p_offset,2)) = 0;
            end               
            %if size(obj.offset,2) == 1 && size(obj.offset,2) < obj.nSrc
            %    obj.offset = repmat(obj.offset,1,obj.nSrc);
            %end
        end
        function out = get.offset(obj)
            out = obj.p_offset(1:2,:);
        end
        %% Nyquist sampling
        function out = get.nyquistSampling(obj)
            out = obj.p_nyquistSampling;
        end
        function set.nyquistSampling(obj,val)
            obj.p_nyquistSampling = val;
            obj.fftPad = ceil(obj.p_nyquistSampling*2);
%             if ~isempty(obj.nLensletWavePx) % refreshing phasor
%                 setPhasor(obj);
%             end
        end
        
        %% Field stop
        function out = get.fieldStopSize(obj)
            out = obj.p_fieldStopSize;
        end
        function set.fieldStopSize(obj,val)
            obj.p_fieldStopSize = val;
%             if ~isempty(obj.nLensletWavePx) % refreshing phasor
%                 setPhasor(obj);
%             end
        end
        function out = get.elongatedFieldStopSize(obj)
            out = obj.elongatedFieldStopSize;
        end
        function set.elongatedFieldStopSize(obj,val)
            obj.elongatedFieldStopSize = val;
        end
        %% Number of wave pixel per lenslet
        function out = get.nLensletWavePx(obj)
            out = obj.p_nLensletWavePx;
        end
        function set.nLensletWavePx(obj,val)
            obj.p_nLensletWavePx = val;
            if isempty(obj.p_fieldStopSize)
                fprintf(' @(lensletArray)> Setting the lenslet field stop size!\n')
                obj.fieldStopSize  = obj.p_nLensletWavePx/obj.p_nyquistSampling/2;
            end
            setPhasor(obj);
            %             % set the wave reshaping index (2D to 3D)
            %             logBook.add(obj,'Set the wave reshaping index (2D to 3D)')
        end
        
        function out = pixelScale(obj,src,tel)
            %% PIXELSCALE Sky pixel scale
            %
            % out = pixelScale(obj,src,tel) returns the pixel scale
            % projected in the sky for the specified source and telescope
            % object
            %
            % ATTENTION: This is the "numerical" pixel scale, not the
            % detector pixel scale which after binning may be different
            % from this one
            
            nOutWavePx    = obj.nLensletWavePx*obj.fftPad;    % Pixel length of the output wave
            evenOdd       = rem(obj.nLensletWavePx,2);
            if ~rem(nOutWavePx,2) && evenOdd
                nOutWavePx = nOutWavePx + evenOdd;
            end
            nOutWavePx = max(nOutWavePx,obj.nLensletWavePx);
            out = skyAngle( (src.wavelength/tel.D)*obj.nLensletWavePx*obj.nLenslet/nOutWavePx);
            % show also in mas 
            out.unit = 'mas';
        end

        function wavePrgted = propagateThrough(obj,src_in)
            %% PROPAGATETHROUGH Fraunhoffer wave propagation to the lenslets focal plane
            %
            % propagateThrough(obj) progates the object wave throught the
            % lenslet array and computes the imagelets
            
            % for a given wave of size [nxn], 
            % nyquistSampling=1 means 2 pixels per fhwm obtained by setting
            % fftPad=2, fieldStopSize should be set to default n/nLenslet
            % fwhm
            % nyquistSampling=0.5 means 1 pixels per fhwm obtained by
            % setting fftPad=1, fieldStopSize should be set to default
            % n/nLenslet/nyquistSampling/2
            
            % Check for LGS asterisms
            [n1,n2,n3] = size(src_in);
            n12 = n1*n2;
            if ndims(src_in)==3 && n12>1
                m_imagelets = zeros(obj.nLensletsImagePx,obj.nLensletsImagePx,n12);
                count       = 0;
                for k1 = 1:n1
                    for k2 = 1:n2
                        count = count + 1;
                        add(obj.log,obj,sprintf('Processing LGS #%d/%d',count,n12))
                        propagateThrough(obj,src_in(k1,k2,:));
                        m_imagelets(:,:,count) = obj.imagelets;
                    end
                end
                obj.imagelets = reshape(m_imagelets,size(m_imagelets,1),[]);
                obj.nArray = n12;
                return
            else
                src = src_in;
            end
            
            % apply any offsets
            if ~isempty(obj.offset) && any(obj.offset(:) ~= 0) || ~isempty(obj.rotation) && any(obj.rotation(:) ~= 0)
                buf = src.phaseUnmasked;
                [nx, ny, nz] = size(src(1).phaseUnmasked);
                R = 0.5;
                [x0,y0] = meshgrid(linspace(-R,R,nx), linspace(-R,R,ny));
                nSrc_ = numel(src);
                phasei = zeros(nx, ny, nz);
                val = zeros(nx, ny, nSrc_);
                for iSrc = 1:nSrc_
                    if size(obj.offset,2) >= iSrc
                        xi = x0 - obj.offset(1,iSrc);
                        yi = y0 - obj.offset(2,iSrc);
                        ri = obj.rotation(iSrc);
                        for kWave = 1:nz
                            if ri ~= 0
                                if license('checkout','signal_toolbox')
                                    phasei(:,:,kWave) = imrotate(src(iSrc).phaseUnmasked(:,:,kWave), ri*180/pi,'crop');
                                else
                                    phasei(:,:,kWave) = rotate(src(iSrc).phaseUnmasked(:,:,kWave),ri*180/pi);
                                end
                                phasei(:,:,kWave) = interp2(x0,y0,phasei(:,:,kWave),xi,yi,'cubic',0);
                            else
                                phasei(:,:,kWave) = interp2(x0,y0,src(iSrc).phaseUnmasked(:,:,kWave),xi,yi,'cubic',0);
                            end
                        end
                        src(iSrc).resetPhase(phasei);
                        
                        % offset amplitude
                        if ri ~= 0
                            if license('checkout','signal_toolbox')
                                amp = imrotate(src(iSrc).amplitudeUnmasked, ri*180/pi,'crop');
                            else
                                amp = rotate(src(iSrc).amplitudeUnmasked,ri*180/pi);
                            end
                            amp = interp2(x0,y0,amp,xi,yi,'cubic',0);
                        else
                            amp = interp2(x0,y0,src(iSrc).amplitude,xi,yi,'cubic',0);
                        end
                    else
                        src(iSrc).resetPhase(buf);
                        amp = src(iSrc).amplitude;
                    end
                    val(:,:,iSrc) = amp.*exp(1i*src(iSrc).phaseUnmasked);
                end
                val = reshape(val, nx, ny*nSrc_);
                
            else % no offset case
                val = src.catWave; % get complex amplitudes
            end
            
            %             % apply any offsets
%             if ~isempty(obj.offset) && any(obj.offset(:) ~= 0) || ~isempty(obj.rotation) && any(obj.rotation(:) ~= 0)
%                 buf = src.phaseUnmasked;
%                 [nx, ny, nz] = size(src(1).phaseUnmasked);
%                 R = 0.5;
%                 [x0,y0] = meshgrid(linspace(-R,R,nx), linspace(-R,R,ny));
%                 nSrc_ = numel(src);
%                 phasei = zeros(nx, ny, nz);
%                 val = zeros(nx, ny, nSrc_);
%                 for iSrc = 1:nSrc_
%                     if size(obj.offset,2) >= iSrc                        
%                         xi = x0 - obj.offset(1,iSrc);
%                         yi = y0 - obj.offset(2,iSrc);
%                         ri = obj.rotation(iSrc);
%                         for kWave = 1:nz
%                             if ri ~= 0
%                                 if license('checkout','signal_toolbox')
%                                     phasei(:,:,kWave) = imrotate(src(iSrc).phaseUnmasked(:,:,kWave), ri*180/pi,'crop');
%                                 else
%                                     phasei(:,:,kWave) = rotate(src(iSrc).phaseUnmasked(:,:,kWave),ri*180/pi);
%                                 end
%                                 phasei(:,:,kWave) = interp2(x0,y0,phasei(:,:,kWave),xi,yi,'cubic',0);
%                             else
%                                 phasei(:,:,kWave) = interp2(x0,y0,src(iSrc).phaseUnmasked(:,:,kWave),xi,yi,'cubic',0);
%                             end
%                         end
%                         src(iSrc).resetPhase(phasei);
%                     else
%                         src(iSrc).resetPhase(buf);
%                     end
%                     val(:,:,iSrc) = src(iSrc).amplitudeUnmasked.*exp(1i*src(iSrc).phaseUnmasked);
%                 end
%                 val = reshape(val, nx, ny*nSrc_);
%                 
%             else
%                 val = src.catWave; % get complex amplitudes
%             end
            

            
            if ndims(src)==3 % if src object array is 3D then it is an LGS hence collapse the 3D imagelets 
                obj.sumStack = true;
            end
            if isscalar(val) % if source amplitude and phase not set, set a default one
                n = obj.nLensletWavePx*obj.nLenslet;
                set(src,...
                    'mask',true(n),...
                    'amplitude',ones(n),...
                    'phase',zeros(n));
                if ~isempty(src(1).opticalPath)
                    cellfun(@(x)relay(x,src),src(1).opticalPath,'uniformOutput',false)
                end
                val = src.catWave;
            end
            [nLensletsWavePx,nLensletsWavePxNGuideStar,nWave] = size(val);
            if ndims(src)==3 || nWave>1 
                obj.nArray = 1;
            else
                obj.nArray = numel(src);
            end
            % Resize the 3D input into a 2D input
            nLensletsWavePxNGuideStar = nLensletsWavePxNGuideStar*nWave;
            val = reshape(val,nLensletsWavePx,nLensletsWavePxNGuideStar);
            nLensletWavePx = nLensletsWavePx./obj.nLenslet;
            if isempty(obj.nLensletWavePx) || all(obj.nLensletWavePx~=nLensletWavePx)
                obj.nLensletWavePx = nLensletWavePx;
            end
            nLensletArray = nLensletsWavePxNGuideStar/nLensletsWavePx;
            
            if ~isempty(obj.opticalAberration)
                val = obj.opticalAberration.*val;
            end
%             obj.nArray = nLensletArray;
            % Invocation of the zoom optics for conjugation to finite
            % distance
%             if isfinite(obj.conjugationAltitude)
% %                 val                 = obj.zoomTransmittance.*val;
%                 val = repmat(obj.zoomTransmittance,1,nLensletArray).*val;
%             end
%             nOutWavePx    = obj.nLensletImagePx*obj.fftPad;    % Pixel length of the output wave
%             evenOdd       = rem(obj.nLensletImagePx,2);
            nOutWavePx    = obj.nLensletWavePx*obj.fftPad;    % Pixel length of the output wave
            evenOdd       = rem(obj.nLensletWavePx,2);
            if ~rem(nOutWavePx,2) && evenOdd
                nOutWavePx = nOutWavePx + evenOdd;
            end
            nOutWavePx = max(nOutWavePx,obj.nLensletWavePx);
            nLensletSquareWavePx    = obj.nLensletWavePx*obj.nLenslet^2*nLensletArray;
            wavePrgted = zeros(nOutWavePx,nLensletSquareWavePx);
            val        = val./nOutWavePx;
%             nLensletWavePx   = obj.nLensletWavePx;
            if obj.convKernel == false
                nLensletImagePx  = obj.nLensletImagePx;
                nLensletsImagePx = obj.nLensletsImagePx;
            else
                nLensletImagePx  = obj.nLensletWavePx*obj.fftPad;
                nLensletsImagePx = obj.nLensletWavePx*obj.fftPad*obj.nLenslet;
            end
            
            %%% ODD # OF PIXELS PER LENSLET
            if isempty(obj.fftPhasor)
%                 fprintf('ODD # OF PIXELS PER LENSLET (%d) Phasor empty!\n',obj.nLensletImagePx)
                % Shape the wave per columns of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletSquareWavePx);
                u         = any(val); % Index of non-zeros columns
                wavePrgted(:,u) = fftshift(fft(val(:,u),nOutWavePx),1);
                % Select the field of view
                v = [];
                if nOutWavePx>nLensletImagePx
%                     disp('Cropping!')
%                     centerIndex = (nOutWavePx+1)/2;
%                     halfLength  = (nLensletImagePx-1)/2;
                    centerIndex = ceil((nOutWavePx+1)/2);
                    halfLength  = floor(nLensletImagePx/2);
                    v           = true(nOutWavePx,1);
%                     v((-halfLength:halfLength)+centerIndex) ...
%                                 = false;
                    v((0:nLensletImagePx-1)-halfLength+centerIndex) ...
                                = false;
                elseif nOutWavePx<nLensletImagePx
                    if obj.convKernel == false
                    error('lensletArray:propagateThrough:size','The computed image is smaller than the expected image!')
                    end
                end
                wavePrgted(v,:) = [];
                % Back to transpose 2D
                val       = reshape( wavePrgted ,...
                    nLensletsImagePx,obj.nLensletWavePx*obj.nLenslet*nLensletArray).';
                % Shape the wave per rows of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                u         = any(val); % Index of non-zeros columns
                wavePrgted = zeros(nOutWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                wavePrgted(:,u)  = fftshift(fft(val(:,u),nOutWavePx),1);
                wavePrgted(v,:) = [];
            else
            %%% EVEN # OF PIXELS PER LENSLET
%                 fprintf('EVEN # OF PIXELS PER LENSLET (%d) Phasor exist!\n',obj.nLensletImagePx)
                val       = val.*repmat(obj.fftPhasor,1,nLensletArray);
                % Shape the wave per columns of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletSquareWavePx);
                u         = any(val); % Index of non-zeros columns
                wavePrgted(:,u) = fft(val(:,u),nOutWavePx);
                v = [];
                if nOutWavePx>nLensletImagePx
%                     disp('Cropping!')
%                     centerIndex = nOutWavePx/2+1;
%                     halfLength  = nLensletImagePx/2;
                    centerIndex = ceil((nOutWavePx+1)/2) + rem(nLensletWavePx,2);
                    halfLength  = floor(nLensletImagePx/2);
                    v           = true(nOutWavePx,1);
                    v((0:nLensletImagePx-1)+centerIndex-halfLength) ...
                                = false;
                elseif nOutWavePx<nLensletImagePx;
                    if obj.convKernel == false
                    error('lensletArray:propagateThrough:size','The computed image is smaller than the expected image!')
                    end
                end
                wavePrgted(v,:) = [];
                % Back to transpose 2D
                val       = reshape( wavePrgted ,...
                    nLensletsImagePx,obj.nLensletWavePx*obj.nLenslet*nLensletArray).';
                % Shape the wave per rows of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                u         = any(val); % Index of non-zeros columns
                wavePrgted = zeros(nOutWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                wavePrgted(:,u)  = fft(val(:,u),nOutWavePx);
                wavePrgted(v,:) = [];
            end
            
                % Back to transpose 2D
                wavePrgted  = reshape(wavePrgted,nLensletsImagePx*nLensletArray,nLensletsImagePx).';

            % and back to input wave array shape
            [n,m] = size(wavePrgted);
            wavePrgted = reshape(wavePrgted,[n,m/nWave,nWave]);

                        
            % extend "artificially" the lenslet fieldStopSize
            tmp = wavePrgted;
            if ~isempty(obj.elongatedFieldStopSize)
                nPadPix = obj.elongatedFieldStopSize/2*obj.nyquistSampling*2;
                [nLensletsWavePx,nLensletsWavePxNGuideStar,nWave] = size(wavePrgted);
                if nLensletsWavePx ~= nLensletsWavePxNGuideStar
                    nFrames = nLensletsWavePxNGuideStar / nLensletsWavePx;
                    wavePrgted = zeros(obj.nLenslet * (obj.nLensletImagePx),obj.nLenslet * (obj.nLensletImagePx),nFrames);
                    for iFrame = 1:nFrames
                        wavePrgted(:,:,iFrame) = tmp(:,(iFrame-1)*nLensletsWavePx+(1:nLensletsWavePx));
                    end
                else
                    nFrames = nWave;
                end
                
                tmp = wavePrgted;
                extendedFramePx = obj.nLenslet * (obj.nLensletImagePx+2*nPadPix);
                wavePrgted = zeros(extendedFramePx,extendedFramePx,nFrames);
                for iFrame = 1:nFrames
                    %reshape wavePrgted in 'line'
                    currentFrame = tmp(:,:,iFrame);
                    currentFrame =  reshape(currentFrame,obj.nLensletImagePx,obj.nLensletImagePx*obj.nLenslet^2);
                    %tic
                    nPx = obj.nLensletImagePx;
                    nL = obj.nLenslet;
                    indices1 = zeros(nPx*nL,1);
                    indices2 = zeros(nPx*nL^2,1);
                    for ii = 1:nL
                        indices1((ii-1)*nPx+1:ii*nPx) = linspace(ii,(nPx-1)*nL+ii, nPx);
                    end
                    for ii = 1:nL
                        indices2((ii-1)*nL*nPx+1:ii*nL*nPx) = indices1 + (ii-1)*nL*nPx;
                    end
                    currentFrame = currentFrame(:, indices2);
                    % 2D->3D
                    currentFrame = reshape(currentFrame, nPx, nPx, nL^2);
                    
                    currentFrame = padarray(currentFrame,[nPadPix,nPadPix,0],0,'both');
                    nPx = nPx + 2*nPadPix;
                    %back to 2d
                    currentFrame = reshape(currentFrame,nPx, nPx*nL^2);
                    
                    indices1 = zeros(nPx*nL,1);
                    indices2 = zeros(nPx*nL^2,1);
                    for ii = 1:nPx
                        indices1((ii-1)*nL+1:ii*nL) = linspace(ii,(nL-1)*nPx+ii, nL);
                    end
                    
                    for ii = 1:nL
                        indices2((ii-1)*nL*nPx+1:ii*nL*nPx) = indices1 + (ii-1)*nL*nPx;
                    end
                    currentFrame = currentFrame(:,indices2);
                    %back to square frame
                    currentFrame = reshape(currentFrame, nPx*nL, nPx*nL);
                    %toc
                    wavePrgted(:,:,iFrame) = currentFrame;
                end
                
                
                if nLensletsWavePx ~= nLensletsWavePxNGuideStar
                    tmp = wavePrgted;
                    wavePrgted = zeros(extendedFramePx, extendedFramePx * nFrames);
                    for iFrame = 1:nFrames
                        wavePrgted(:,(iFrame-1) * extendedFramePx + (1:extendedFramePx)) = tmp(:,:, iFrame);
                    end
                    
                    
                end
                
            end
            % end of artificially extending the lenslet fieldStopSize
            
            if obj.sumStack
                wavePrgted = mean(wavePrgted,3);
            end
            obj.sumStack = false;
            
            % reset phase from any offsets 
            if ~isempty(obj.offset) && any(obj.offset(:) ~= 0) || ~isempty(obj.rotation) && any(obj.rotation(:) ~= 0)
                src(iSrc).resetPhase(buf);
            end
        end
        
        function relay(obj,src)
            wavePrgted = propagateThrough(obj,src);
            obj.imagelets = wavePrgted.*conj(wavePrgted)*obj.throughput;            
        end

        function varargout = imagesc(obj,varargin)
        %% IMAGESC Display the lenslet imagelets
        %
        % imagesc(obj) displays the imagelets of the lenslet array object 
        %
        % imagesc(obj,'PropertyName',PropertyValue) displays the imagelets
        % of the lenslet array object and set the properties of the
        % graphics object imagesc
        %
        % h = imagesc(obj,...) returns the graphics handle
        %
        % See also: imagesc
        
 %             if ishandle(obj.imageHandle)
%                 axisLim = [0,obj.nLensletsImagePx]+0.5;
%                 set( get(obj.imageHandle,'parent') , ...
%                     'xlim',axisLim,'ylim',axisLim);
%             end
       if ishandle(obj.imageHandle)
            set(obj.imageHandle,'Cdata',obj.imagelets,varargin{:});
        else
            %                 figure
            obj.imageHandle = image(obj.imagelets,...
                'CDataMApping','Scaled',varargin{:});
            colormap(pink)
            axis xy square
            axisXLim = [0,size(obj.imagelets,2)]+0.5;
            axisYLim = [0,size(obj.imagelets,1)]+0.5;
            set( get(obj.imageHandle,'parent') , ...
                'xlim',axisXLim,'ylim',axisYLim);
            axis equal tight
            colorbar
        end
        if nargout>0
            varargout{1} = obj.imageHandle;
        end
        end


    end
    
    methods (Access=private)     
        
        function setPhasor(obj)
        % Phasor settings for Fraunhoffer propagation in order to have the
        % spots centered in between pixels for even pixel sampling
%             nOutWavePx    = obj.nLensletImagePx*obj.fftPad;    % Pixel length of the output wave
%             evenOdd       = rem(obj.nLensletImagePx,2);
%             if ~rem(nOutWavePx,2) && evenOdd
%                 nOutWavePx = nOutWavePx + evenOdd;
%             end
%             fprintf(' 0ld nOutWavePx = %d',nOutWavePx);
            nOutWavePx    = obj.nLensletWavePx*obj.fftPad;    % Pixel length of the output wave
            evenOdd       = rem(obj.nLensletWavePx,2);
            if ~rem(nOutWavePx,2) && evenOdd
                nOutWavePx = nOutWavePx + evenOdd;
            end
            nOutWavePx = max(nOutWavePx,obj.nLensletWavePx);
%             fprintf(' - new nOutWavePx = %d\n',nOutWavePx);

            if ~rem(obj.nLensletImagePx,2)
                % shift the intensity of half a pixel for even sampling
                fprintf(' @(lensletArray)> Set phasor (shift the intensity of half a pixel\n for even intensity sampling)\n')
                [u,v]         = ndgrid((0:(obj.nLensletWavePx-1)).*(~rem(obj.nLensletWavePx,2)-nOutWavePx)./nOutWavePx);
                obj.fftPhasor = repmat( exp(-1i.*pi.*(u+v)) , obj.nLenslet, obj.nLenslet );
            else
                fprintf(' @(lensletArray)> Reset phasor\n')
                obj.fftPhasor = [];
            end            
        end
        
        function setImageletsListener(obj)
            %% SETIMAGELETSLISTENERS Imagelets listener
            obj.imageletsListener = addlistener(obj,'imagelets','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.imageletsListener.Enabled = false;
        end

    end
    
    methods (Static)
            
        function obj = loadobj(obj)
            %% LOADOBJ
            add(obj.log,obj,'Load!')
            setImageletsListener(obj)
            obj.log = logBook.checkIn(obj);
        end
        
    end
    
end