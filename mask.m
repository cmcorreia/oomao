classdef mask < handle
    %% mask
    % Creates object containg transparency maks for fourier-masked
    % WFSensors
    %
    %  m = mask(n,nPx); creates a 4-faced pyramid (centred, no rotation).
    %
    %  m = mask(n,nPx,'alpha',alpha); creates a 4-faced pyramid (centred,
    %  no rotation) with a angle given by alpha (default=pi/2). This allows
    %  to simulate the flattened pyramid
    %
    %  m =
    %  mask(n,nPx,'alpha',pi,'type','pyramidSLM','rotation',-pi/5,'centrePos',55);
    %  creates a pyramid mask n-faced pyramid mask with rotation and centrePosition different from zero (user-defined) 
    
    
    properties
        
        % # of points
        resolution;
        
        % centre postion
        centrePos
        
        % type :: pyramid
        type
        
        % # faces
        nFaces
        
        % angle of the faces
        alpha
        
        % rotation with respect to centre (in radians)
        rotation
        
        % phase map
        phaseMask
        
        % transparency phase map
        theMask
        
        % tag
        tag = 'FOURIER TRANSPARENCY MASK';
        
    end
    
    
    properties (Dependent)%, SetAccess=private)
        
        
        
    end
    
    properties (Access=private)
        %p_p;
    end
    
    methods
        
        %% Constructor
        function obj = mask(nFaces,resolution,varargin)
            p = inputParser;
            addRequired(p,'nFaces', @isnumeric);
            addRequired(p,'resolution', @isnumeric);
            addOptional(p,'alpha',pi/2,@isnumeric);
            
            addOptional(p,'centrePos',0,@isnumeric);
            addOptional(p,'rotation',0,@isnumeric);
            
            paramName = 'type';
            validationFcn = @(x) validateattributes(x,{'char'},{'nonempty'});
            addOptional(p,paramName,'pyramid',validationFcn);
            p.parse(nFaces,resolution,varargin{:})
            
            
            obj.nFaces      = p.Results.nFaces;
            obj.resolution  = p.Results.resolution;
            if length(obj.resolution) == 1
                obj.resolution = [1;1]*obj.resolution;
            end
            obj.type        = p.Results.type;
            
            obj.centrePos   = p.Results.centrePos;
            if obj.centrePos == 0
                obj.centrePos = obj.resolution/2+0.5;
            end
            if length(obj.centrePos) == 1
                obj.centrePos = [1;1]*obj.centrePos; 
            end
            obj.alpha = p.Results.alpha;
            obj.rotation = p.Results.rotation;
            switch obj.type
                case 'pyramid'
                    pyramid(obj, obj.alpha)
                case 'pyramidSLM'
                    pyramidSLM(obj, obj.nFaces, obj.alpha, obj.centrePos, obj.rotation)
                case 'zelda'
                case 'FQMP'
                otherwise
                    'NO OPTION CHOSE. PLEASE RESTART'
            end
            
            
        end
        
        %% Destructor
        function delete(obj)
            
        end
        %% Display
        function display(obj)
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' %dX%d pixels transparency mask: \n ',...
                obj.resolution(1),obj.resolution(2))
            fprintf('----------------------------------------------------\n')
        end
        %% Displays object phaseMask
        function imagesc(obj,varargin)
            % IMAGESC Display the object phase map
            %
            % imagesc(obj) displays the frame of the mask object
            %
            % imagesc(obj,'PropertyName',PropertyValue) displays the frame of
            % the mask object and set the properties of the graphics object
            % imagesc
            %
            % h = imagesc(obj,...) returns the graphics handle
            %
            % See also: imagesc
            imagesc(obj.phaseMask)
        end
        %% Set the DEFAULT 4-sided pyramid mask
        function obj = pyramid(obj, alpha,rooftop)
            if ~exist('rooftop','var') || isempty(rooftop)
                rooftop = [0 0];
            end
            if isscalar(alpha)
                realAlpha = ones(1,4)*alpha;
            else
                realAlpha = alpha;
            end
            % Complex alpha for PYR surface slope error in 2 directions...
            imagAlpha = imag(realAlpha);
            realAlpha = real(realAlpha);
            
            nx = rooftop(1);
            ny = rooftop(2);
            
            [fx,fy] = freqspace(obj.resolution,'meshgrid');
            fx = fx.*floor(obj.resolution(1)/2);
            fy = fy.*floor(obj.resolution(2)/2);
            %pym = zeros(pxSide);
            
            % pyramid face transmitance and phase for fx>=0 & fy>=0
            mask  = graduatedHeaviside(fx,nx).*graduatedHeaviside(fy,nx);
            phase = -realAlpha(1).*(fx+fy) + -imagAlpha(1).*(-fx+fy);
            pym   = mask.*exp(1i.*phase);
            pyp   = mask.*phase;
            % pyramid face transmitance and phase for fx>=0 & fy<=0
            mask  = graduatedHeaviside(fx,ny).*graduatedHeaviside(-fy,-ny);
            phase = -realAlpha(2).*(fx-fy) + -imagAlpha(2).*(fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            pyp   = pyp + mask.*phase;
            
            % pyramid face transmitance and phase for fx<=0 & fy<=0
            mask  = graduatedHeaviside(-fx,-nx).*graduatedHeaviside(-fy,-nx);
            phase = realAlpha(3).*(fx+fy) + -imagAlpha(3).*(fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            pyp   = pyp + mask.*phase;
            % pyramid face transmitance and phase for fx<=0 & fy>=0
            mask  = graduatedHeaviside(-fx,-ny).*graduatedHeaviside(fy,ny);
            phase = -realAlpha(4).*(-fx+fy) + -imagAlpha(4).*(-fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            pyp   = pyp + mask.*phase;
            obj.phaseMask = pyp;
            obj.theMask = fftshift(pym)./sum(abs(pym(:)));
        end
        %% Set pyramid phase map
        function obj = pyramidSLM(obj, nFaces, angle, centrePos, rotation, rooftop)
            
            if ~exist('rooftop','var') || isempty(rooftop)
                rooftop = [0 0];
            end
            nx = obj.resolution(1);
            ny = obj.resolution(2);
            
            cx = centrePos(1);
            cy = centrePos(2);
            
            x = linspace(-cx+1, nx-cx, nx);
            y = linspace(-cy+1, ny-cy, ny);
            [xGrid, yGrid] = meshgrid(x, y);
            xGrid = xGrid';
            yGrid = yGrid';
            
            % DEFINE THE SLOPE GRID
            angleGrid = atan2(yGrid*sin(-rotation) + ...
                xGrid*cos(-rotation), ...
                yGrid*cos(-rotation) - ...
                xGrid*sin(-rotation));
            
            % INITIALIZE PYRAMID MASK
            pyp = zeros(nx, ny);
            for kFaces=0:nFaces-1
                theta = pi*(1/nFaces - 1) + kFaces*2*pi/nFaces + rotation;
                slope = sin(theta)*xGrid + cos(theta)*yGrid;

                %slope((-pi+kFaces*2*pi/nFaces >= angleGrid) |...
                %    (angleGrid >= (-pi+(kFaces + 1)*2*pi/nFaces))) = 0;

                % Take into account the last tile of the pyramid mask
                if kFaces == nFaces-1
                    slope((-pi+kFaces*2*pi/nFaces <= angleGrid) &...
                        (angleGrid <= (-pi+(kFaces + 1)*2*pi/nFaces))) = 0;
                else
                    slope((-pi+kFaces*2*pi/nFaces <= angleGrid) &...
                        (angleGrid < (-pi+(kFaces + 1)*2*pi/nFaces))) = 0;
                end
                
                pyp = pyp + angle*slope;
            end
            
            obj.phaseMask = pyp;
            pym = exp(1i*pyp);
            obj.theMask = fftshift(pym)./sum(abs(pym(:)));
        end
        
        %% Set pyramid phase map
        function obj = zelda(obj, center, radius, phase_shift)
            
            Nx = obj.size(1);
            Ny = obj.size(2);
            
            Cx = center(1);
            Cy = center(2);
            
            x = linspace(-Cx+1, Nx-Cx, Nx);
            y = linspace(-Cy+1, Ny-Cy, Ny);
            [x_grid, y_grid] = meshgrid(x, y);
            x_grid = x_grid';
            y_grid = y_grid';
            
            % DEFINE THE POSITION GRID
            positionGrid = sqrt(x_grid.^2 + y_grid.^2);
            
            % INITIALIZE ZELDA MASK
            Zelda = zeros(Nx, Ny);
            Zelda(positionGrid <= radius) = phase_shift;
            obj.phasemap = Zelda;
            
        end
        %% Set FQPM phase map
        function obj = FQPM(obj, phase_shift, center)
            
            Cx = center(1);
            Cy = center(2);
            Nx = obj.size(1);
            Ny = obj.size(2);
            
            FQPM = phase_shift*[zeros(Cx, Cy), ...
                ones(Cx, Ny-Cy); ...
                ones(Nx-Cx, Cy), ...
                zeros(Nx-Cx, Ny-Cy)];
            
            obj.phasemap = FQPM;
        end
    end
end

