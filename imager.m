classdef imager < detector
    %% Imaging camera
    %
    % imgr = imager(resolution) creates an imaging camera object from the
    % detector resolution
    %
    % Example:
    % tel = telescope(1,'resolution',21,'samplingTime',1);
    % imgr = imager(21);
    % src = source.*tel*imgr;
    % figure
    % imagesc(imgr.frame)

    properties
        % reference frame
        referenceFrame;
        % imaging lens
        imgLens;
        % Strehl ratio
        strehl;
        % entrapped energy
        ee;
        % entrapped energy slit width in # pixels (conversion to physical
        % units done by multiplying by the pixel scale)  -
        % ccorreia ~May2019
        eeWidth;
        % telescope diameter
        diameter
    end
        
    properties (Access=private)
%         % integration count;
%         frameCount=0;
        % telescope
        tel;
        %diameter
    end
    
    methods
        
        %% Constructor
        function obj = imager(varargin)
            
              p = inputParser;
              addParameter(p,'nyquistSampling',4,@isnumeric)
              addParameter(p,'fieldStopSize',10,@isnumeric)
              addParameter(p,'exposureTime',1,@isnumeric)
              addParameter(p,'clockRate',1,@isnumeric)
              addParameter(p,'diameter',1,@isnumeric)
              %addOptional(p,'tel', telescope(-1), @(x) isa(x,'telescopeAbstract') );
              
              
              parse(p,varargin{:})
              
              %             if isa(in,'telescopeAbstract')
              %                 resolution = in.resolution;
              %             elseif isnumeric(in)
              %                 resolution = in;
              %             else
              %                 error('oomao:imager','Inputer is either numeric or a telescope class')
              %             end
              
              resolution = 2*p.Results.nyquistSampling*p.Results.fieldStopSize;
              obj = obj@detector(resolution);
              obj.imgLens = lensletArray(1);
              obj.imgLens.nyquistSampling = p.Results.nyquistSampling;
              obj.imgLens.fieldStopSize = p.Results.fieldStopSize;
              %             if isa(in,'telescopeAbstract')
              %                 obj.tel = in;
              %                 obj.exposureTime = in.samplingTime;
              %             end
              obj.exposureTime = p.Results.exposureTime;
              obj.clockRate = p.Results.clockRate;
              if ~isempty(p.Results.diameter)
                  obj.diameter = p.Results.diameter;
              else
                  obj.diameter = 1;
                  fprintf('Caution, a default 1m telescope is associated to the class\n')
              end
              
               % Frame listener
%             obj.frameListener = addlistener(obj,'frameBuffer','PostSet',...
%                 @(src,evnt) obj.imagesc );
%             obj.frameListener.Enabled = false;
        end
        
        function relay(obj,src)
            %% RELAY source propagation
            
            relay(obj.imgLens,src)
%             if all([src.timeStamp]>=obj.startDelay)
%                 if obj.frameCount==0
%                     disp(' @(detector:relay)> starting imager integration!')
%                 end
%                 obj.startDelay = -Inf;
%                 obj.frameCount = obj.frameCount + 1;
%                 obj.frameBuffer = obj.frameBuffer + cat(2,src.intensity);
%                 if src.timeStamp>=obj.exposureTime
%                     src.timeStamp = 0;
%                     flush(obj,src);
%                 end
%             end
              readOut(obj,obj.imgLens.imagelets)
              if obj.frameCount==0 && obj.startDelay==0
                  flush(obj,length(src))
                  %flush(obj)
              end
        end
        
        function flush(obj,nSrc)
            %fprintf(' @(detector:relay)> reading out and emptying buffer (%d frames)!\n',obj.frameCount)
            if nargin<2
                nSrc = 1;
            end
            if ~isempty(obj.referenceFrame) && ~isempty(obj.frame)
                obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize*2;
                [n1,n2] = size(obj.referenceFrame);%n = length(obj.referenceFrame);
                nSrc = n2/n1; % if referenceFrame is square, nSrc=1 otherwise, nSrc = n2/n1 
                src_ = source.*obj.referenceFrame;
                wavePrgted = propagateThrough(obj.imgLens,src_);
                otf =  abs(wavePrgted);
                otf = mat2cell(otf,size(otf,1),size(otf,2)/nSrc*ones(1,nSrc));
                nFrame = obj.exposureTime*obj.clockRate;
                
                
                m_frame = obj.frame/nFrame;
                m_frame = reshape(m_frame, size(m_frame,1), size(m_frame,2)*size(m_frame,3));
                nf = [nSrc size(m_frame,2)/n2]; %nf = size(m_frame)/n;
                m_frame = mat2cell( m_frame , n1, n2*ones(1,nf(2))*size(m_frame,3)); % m_frame = mat2cell( m_frame , n*ones(1,nf(1)), n*ones(1,nf(2)));
                
                obj.strehl = zeros(nf);% obj.strehl = zeros(1,length(m_frame));
                if ~isempty(obj.eeWidth)
                    obj.ee     = zeros([nf length(obj.eeWidth)]);%zeros(1,length(m_frame));
                    obj.ee = squeeze(obj.ee); % remove non-singleton dimensions
                end
                
                if ~isempty(obj.tel)
                    D = obj.tel.D;
                else
                    D = obj.diameter;
                end
                for kFrame = 1:length(m_frame)
                    src_ = src_.*m_frame{kFrame};
                    wavePrgted = propagateThrough(obj.imgLens,src_);
                    otfAO =  abs(wavePrgted);
                    otfAO = mat2cell(otfAO,size(otfAO,1),size(otfAO,2)/nSrc*ones(1,nSrc));
                    
                    for kobj=1:nSrc
                        %a = [a otfAO];
                        % strehl ratio
                        obj.strehl(kobj,kFrame) = sum(otfAO{kobj}(:))/sum(otf{kobj}(:));
                        % entrapped energy
                        
                        if ~isempty(obj.eeWidth) % EE in eeWidth # of pixels
                            for kIntegBoxSize = 1:length(obj.eeWidth)
                                %a      = (obj.eeWidth(kIntegBoxSize)/(src_.wavelength/D*constants.radian2arcsec))/D;
                                a = obj.eeWidth(kIntegBoxSize);
                                nOtf   = length(otfAO{kobj});
                                %u      = linspace(-1,1,nOtf).*D;
                                u      = linspace(-1,1,nOtf); % removed normalisation by D. The integration box is assumed in # of pixels. Conversion to physical units taken care outside of this function
                                [x,y]  = meshgrid(u);
                                eeFilter ...
                                    = a^2*(sin(pi.*x.*a)./(pi.*x.*a)).*...
                                    (sin(pi.*y.*a)./(pi.*y.*a));
                                otfAO{kobj} = otfAO{kobj}/max(otfAO{kobj}(:));
                                obj.ee(kFrame,kIntegBoxSize) = real(trapz(u,trapz(u,otfAO{kobj}.*eeFilter)));
                                % go to
                                % http://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/PSFtheory.pdf
                                % for checking the accuracy of this function using the DL PSF. When
                                % done delete this comment. ccorreia 7/12/2016
                            end
                        end
                    end
                end
                obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize/2;
            end
            obj.frameCount = 0;
        end
        
%         function imagesc(obj,varargin)
%             %% IMAGESC Display the detector frame
%             %
%             % imagesc(obj) displays the frame of the detector object
%             %
%             % imagesc(obj,'PropertyName',PropertyValue) displays the frame of
%             % the detector object and set the properties of the graphics object
%             % imagesc
%             %
%             % h = imagesc(obj,...) returns the graphics handle
%             %
%             % See also: imagesc
%             
%             disp('Got it')
%             if ishandle(obj.frameHandle)
%                 set(obj.frameHandle,'Cdata',obj.frameBuffer,varargin{:});
%                 %                 xAxisLim = [0,size(obj.frame,2)]+0.5;
%                 %                 yAxisLim = [0,size(obj.frame,1)]+0.5;
%                 %                 set( get(obj.frameHandle,'parent') , ...
%                 %                     'xlim',xAxisLim,'ylim',yAxisLim);
%             else
%                 obj.frameHandle = image(obj.frameBuffer,...
%                     'CDataMApping','Scaled',...
%                     varargin{:});
%                 colormap(pink)
%                 axis xy equal tight
%                 colorbar('location','SouthOutside')
%             end
%         end

   
    end

end