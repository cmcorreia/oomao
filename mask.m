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
    %
    % m = mask(4,256,'rotation',[0.25 0]*pi);
    % creates a mask with the vertical edge rotated by pi/4 degrees and the
    % horizontal axis remaining unaffected
    %{
    %% individual defects
    figure (1)
    % reference pyramid mask
    m = mask(4,256);
    subplot(2,3,1)
    imagesc(m)
    title('Defect-free phase mask')
    axis off
    % decentered mask
    m = mask(4,256,'centrePos',[-10 20]);
    subplot(2,3,2)
    imagesc(m)
    title('Decentred phase mask')
    axis off
    % rotated edges mask
    m = mask(4,256,'rotation',[0.3 1]*pi/4);
    subplot(2,3,3)
    imagesc(m)
    title('Rotated edge phase mask')
    axis off
    % rooftop mask
    m = mask(4,256,'rooftop',[50 0]);
    subplot(2,3,4)
    imagesc(m)
    title('Rooftopped phase mask')
    axis off
    % facet output angle
    alpha = [1 1 1 0.5]*pi/2 + 1i*[0 0.0 0.0 0.0]*pi/2;
    m = mask(4,256,'alpha',alpha);
    subplot(2,3,5)
    imagesc(m)
    title('Radial Output-angle P.M.')
    axis off
    
    % facet output angle
    alpha = [1 1 1 1]*pi/2 + 1i*[0 0.0 0.0 0.9]*pi/2;
    m = mask(4,256,'alpha',alpha);
    subplot(2,3,6)
    imagesc(m)
    title('Azimuthal Output-angle P.M.')
    axis off
%     % edgewidth
%     m = mask(4,256,'edgeWidth',[10 0]);
%     subplot(2,3,6)
%     imagesc(m)
%     title('Edge-width phase mask')
%     axis off

    figs(gca,'pwfsPhaseMasIndividualDefects')
    %% all effects at once
    figure(2)
    alpha = [1 1 1 0.8]*pi/2 + 1i*[0 0.0 0.0 0.2]*pi/2;
    m = mask(4,256,'edgeWidth',[10 0], 'alpha',alpha, 'rooftop',[50 0], 'rotation',[0.3 1]*pi/4,'centrePos',[-10 20]);
    imagesc(m)
    axis square
    axis off
    title('All-defects phase mask')
    set(gca,'fontSize',16)
    figs(gca,'pwfsPhaseMaskAllDefects')
        %}
    properties
        
        % # of points
        resolution;
        
        % centre postion
        centrePos
        
        % type :: pyramid
        type
        
        % # faces
        nFaces
        
        % angle of the faces. Default pi/2 separates re-imaged pupils by 2x
        % the diameter (centre to centre)
        alpha
        
        % rotation with respect to centre (in radians) \in [-pi/4:pi/4]. If
        % two-dimensional, each value rotates one of the pyramid axis
        % individualy
        rotation
        
        % edgeWidth in pixels
        edgeWidth
        
        %rooftop in pixels
        rooftop
        
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
            addOptional(p,'edgeWidth',[0,0],@isnumeric);
            addOptional(p,'rooftop',[0,0],@isnumeric);
            
            paramName = 'type';
            validationFcn = @(x) validateattributes(x,{'char'},{'nonempty'});
            addOptional(p,paramName,'pyramidIntersectPlanes',validationFcn);
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
            else
                if size(obj.centrePos,1) > size(obj.centrePos,2)
                    obj.centrePos = obj.centrePos + (obj.resolution/2+0.5);
                else
                    obj.centrePos = obj.centrePos' + (obj.resolution/2+0.5);
                end
            end
            if length(obj.centrePos) == 1
                obj.centrePos = [1;1]*obj.centrePos;
            end
            obj.alpha       = p.Results.alpha;
            obj.rotation    = p.Results.rotation;
            obj.edgeWidth   = p.Results.edgeWidth;
            obj.rooftop     = p.Results.rooftop;
            switch obj.type
                case 'pyramidIntersectPlanes'
                    pyramidIntersectPlanes(obj, obj.centrePos, obj.alpha, obj.rotation, obj.rooftop)
                case 'pyramid'
                    pyramid(obj, obj.centrePos, obj.alpha, obj.edgeWidth, obj.rotation, obj.rooftop)
                case 'pyramidSLM'
                    pyramidSLM(obj, obj.nFaces, obj.alpha, obj.centrePos, obj.rotation,obj.edgeWidth)
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
        function obj = pyramidIntersectPlanes(obj, centrePos, alpha, rotation, rooftop)
            if ~exist('centrePos','var') || isempty(centrePos)
                cx = 0;
                cy = 0;
            else
                cx = centrePos(1) - (obj.resolution(1)/2+0.5);
                cy = centrePos(2) - (obj.resolution(2)/2+0.5);
            end
            
            
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
            
            if ~exist('rotation','var') || isempty(rotation)
                rotation = [0 0];
            end
            if isscalar(rotation)
                r1 = rotation/(pi/4);
                r2 = rotation/(pi/4);
            else
                r1 = rotation(1)/(pi/4);
                r2 = rotation(2)/(pi/4);
            end
            
            [fx,fy] = freqspace(obj.resolution,'meshgrid');
            fx0 = fx.*floor(obj.resolution(1)/2)-cx;
            fy0 = fy.*floor(obj.resolution(2)/2)-cy;
            
            fx = fx0 +r1*fy0;
            fy = fy0 -r2*fx0;
            
            planes(:,:,1) = -realAlpha(1).*(fx+fy) + -imagAlpha(1).*(-fx+fy);
            planes(:,:,2) = -realAlpha(2).*(fx-fy) + -imagAlpha(2).*(fx+fy);
            planes(:,:,3) = -realAlpha(3).*(-fx-fy) + -imagAlpha(3).*(fx-fy);
            planes(:,:,4) = -realAlpha(4).*(-fx+fy) + -imagAlpha(4).*(-fx-fy);
            pyp = min(planes,[],3);
            
            % apply rooftop
            pyp = pyp + rooftop(1);
            pyp(pyp >0) = 0;
            pym = exp(1i.*pyp);
            
            obj.phaseMask = pyp;
            obj.theMask = fftshift(pym)./sum(abs(pym(:)));
        end
        
        %% Set the DEFAULT 4-sided pyramid mask
        function obj = pyramid(obj, centrePos, alpha, edgeWidth, rotation, rooftop)
            if ~exist('centrePos','var') || isempty(centrePos)
                cx = 0;
                cy = 0;
            else
                cx = centrePos(1) - obj.resolution(1)/2+0.5;
                cy = centrePos(2) - obj.resolution(2)/2+0.5;
            end
            
            if ~exist('edgeWidth','var') || isempty(edgeWidth)
                edgeWidth = [0 0];
            end
            
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
            
            nx = edgeWidth(1);
            ny = edgeWidth(2);
            
            
            if ~exist('rotation','var') || isempty(rotation)
                rotation = [0 0];
            end
            if isscalar(rotation)
                r1 = rotation/(pi/4);
                r2 = rotation/(pi/4);
            else
                r1 = rotation(1)/(pi/4);
                r2 = rotation(2)/(pi/4);
            end
            
            [fx,fy] = freqspace(obj.resolution,'meshgrid');
            fx0 = fx.*floor(obj.resolution(1)/2)-cx;
            fy0 = fy.*floor(obj.resolution(2)/2)-cy;
            
            fx = fx0 +r1*fy0;
            fy = fy0 -r2*fx0;
            %ny = -ny;
            %pym = zeros(pxSide);
            Hx = graduatedHeaviside(fx,nx);
            Hy = graduatedHeaviside(fy,ny);
            % pyramid face transmitance and phase for fx>=0 & fy>=0
            %mask  = graduatedHeaviside(fx,nx).*graduatedHeaviside(fy,ny);
            mask = Hx.*Hy;
            phase = -realAlpha(1).*(fx+fy) + -imagAlpha(1).*(-fx+fy);
            pym   = mask.*exp(1i.*phase);
            pyp   = mask.*phase;
            % pyramid face transmitance and phase for fx>=0 & fy<=0
            %mask  = graduatedHeaviside(fx,nx).*graduatedHeaviside(-fy,-ny);
            mask = Hx.*(1-Hy);
            phase = -realAlpha(2).*(fx-fy) + -imagAlpha(2).*(fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            pyp   = pyp + mask.*phase;
            % pyramid face transmitance and phase for fx<=0 & fy<=0
            %mask  = graduatedHeaviside(-fx,-nx).*graduatedHeaviside(-fy,-ny);
            mask = (1-Hx).*(1-Hy);
            phase = -realAlpha(3).*(-fx-fy) + -imagAlpha(3).*(fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            pyp   = pyp + mask.*phase;
            % pyramid face transmitance and phase for fx<=0 & fy>=0
            %mask  = graduatedHeaviside(-fx,-nx).*graduatedHeaviside(fy,ny);
            mask = (1-Hx).*Hy;
            phase = -realAlpha(4).*(-fx+fy) + -imagAlpha(4).*(-fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            pyp   = pyp + mask.*phase;
            % apply rooftop
            pyp = pyp + rooftop(1);
            pyp(pyp >0) = 0;
            pym = exp(1i.*pyp);
            
            obj.phaseMask = pyp;
            obj.theMask = fftshift(pym)./sum(abs(pym(:)));
        end
        %% Set pyramid phase map
        function obj = pyramidSLM(obj, nFaces, angle, centrePos, rotation, edgeWidth)
            
            if ~exist('edgeWidth','var') || isempty(edgeWidth)
                edgeWidth = [0 0];
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

