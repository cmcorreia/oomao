classdef segment
    % Create a telescope segment giving number of sides, size in px, etc...
    
    properties
        % Segment Properties
        nSides; % number of sides
        radius; % m
        nPxInterp;% matrix size in px
        nPxBase;
        keepBase;
        % user can make a segment in nPxBase resolution. 
        % if nPxInterp is not equal to nPxBase, segment will be interpolate
        matrix; % double matrix containing shape of segment after interpolation
        matrixBase; % matrix containing shape of segment before interpolation
        posX;   % position X in a pupil
        posY;   % position Y in a pupil
        angleSegment; % angle of the segment
        pixelRatio;   % in px/m
        
        % Errors properties
        posXError;
        posYError;
        angleError;
        phaseError; % matrix
        coeffPhaseModes; % array 1 x nModes
        reflexivity;
        sizeError;
    end
    
    methods
        function obj = segment(nsides,R,nPxBase,varargin)
            % CHECKING INPUTS
            inputs = inputParser;
            addRequired(inputs,'nsides',@isnumeric);
            addRequired(inputs,'radius',@isnumeric);
            addRequired(inputs,'nPxBase',@isnumeric);
            addParameter(inputs,'nPxInterp',[],@isnumeric);
            addParameter(inputs,'keepBase',0,@isnumeric);
            addParameter(inputs,'angle',0, @isnumeric);
            addParameter(inputs,'posX',0, @isnumeric);
            addParameter(inputs,'posY',0, @isnumeric);
            addParameter(inputs,'angleError',0, @isnumeric);
            addParameter(inputs,'posXError',0, @isnumeric);
            addParameter(inputs,'posYError',0, @isnumeric);
            addParameter(inputs,'phaseError',0, @isnumeric);
            addParameter(inputs,'phaseCoeff',0, @isnumeric);
            addParameter(inputs,'reflexivity',1, @isnumeric);
            addParameter(inputs,'sizeError',1, @isnumeric);
            inputs.parse(nsides,R,nPxBase,varargin{:});
            
            % INSTANTIATION
            obj.nSides          = inputs.Results.nsides;
            obj.radius          = inputs.Results.radius;
            obj.nPxInterp       = inputs.Results.nPxInterp;
            obj.nPxBase         = inputs.Results.nPxBase;
            obj.keepBase        = inputs.Results.keepBase;
            obj.angleSegment    = inputs.Results.angle;
            obj.posX            = inputs.Results.posX;
            obj.posY            = inputs.Results.posY;
            obj.angleError      = inputs.Results.angleError;
            obj.posXError       = inputs.Results.posXError;
            obj.posYError       = inputs.Results.posYError;
            obj.phaseError      = inputs.Results.phaseError;
            obj.coeffPhaseModes = inputs.Results.phaseCoeff;
            obj.reflexivity     = inputs.Results.reflexivity;
            obj.sizeError       = inputs.Results.sizeError;
            
            % Construction
            if isempty(obj.nPxInterp)
                obj.nPxInterp=obj.nPxBase;
            end
            [obj.matrix, obj.matrixBase]=mkSegment(obj);
            if (~obj.keepBase) || (obj.nPxInterp==obj.nPxBase)
                obj.matrixBase=[];
            end
            obj.pixelRatio=(obj.nPxInterp)/(2*obj.radius);
        end

        function [outInterp,outBase] = mkSegment(obj)
            nP=obj.nPxBase;
            if (obj.nSides>25)
                [rr ,cc] = meshgrid(linspace(-nP,nP,obj.nPxBase));
                tmp=(sqrt((rr).^2+(cc).^2)<=obj.sizeError*(nP+0.1));
            else
                theta=linspace(0,2*pi,obj.nSides+1);
                xv=(nP+1)*obj.sizeError*cos(theta+obj.angleSegment+obj.angleError);
                yv=(nP+1)*obj.sizeError*sin(theta+obj.angleSegment+obj.angleError);
                [X,Y]=meshgrid(linspace(-nP,nP,nP));
                [in,on]=inpolygon(X,Y,xv,yv);
                tmp=in|on;
            end
            if obj.nPxBase~=obj.nPxInterp
                xi       = linspace(-1,1,length(tmp));
                [Xi,Yi]  = meshgrid(xi);
                xo       = linspace(-1,1,obj.nPxInterp);
                [Xo,Yo]  = meshgrid(xo);
                tmp = interp2(Xi,Yi,double(tmp),Xo,Yo,'linear');
                outInterp=tmp.*obj.reflexivity.*exp(1i*obj.phaseError);
            else
                outInterp=tmp.*obj.reflexivity.*exp(1i*obj.phaseError);
            end
            outBase=tmp.*obj.reflexivity.*exp(1i*obj.phaseError);
        end
        
        function [out, phase]=applyPhase(obj,modesCube)
            % WARNING : the "phase" is in modesCube and coeffPhaseModes unit. 
            % 
            nM=length(obj.coeffPhaseModes);
            if size(modesCube,3)>=nM
                phase=sum(bsxfun(@times,modesCube(:,:,1:nM),reshape(obj.coeffPhaseModes,[1 1 nM])),3);
                out=logical(abs(obj.matrix)).*obj.reflexivity.*exp(1i*phase);
            else
                error('Number of coefficients doesn''t match the size of the cube of Modes');
            end
        end
        
        function [out,reflexivity]=applyReflexivity(obj,R)
            out=logical(abs(obj.matrix)).*R.*exp(1i*obj.phaseError);
            reflexivity=R;
        end
        %% Overloading disp, plus and minus operators
        
        function M=plus(M,obj)
            inT=tic;
            [LM,CM]=size(M);
            if mod(LM,2) %odd
                cxloc=(LM-1)/2+1;
            else % even
                cxloc=LM/2;
            end
            if mod(CM,2) %odd
                cyloc=(CM-1)/2+1;
            else % even
                cyloc=CM/2;
            end
            
            [sx,sy]=size(obj.matrix);
            if mod(sx,2) %impair
                Ax=cxloc+(-(sx-1)/2:(sx-1)/2);
            else %pair
                Ax=cxloc+(-sx/2+1:sx/2);
            end
            if mod(sy,2) %impair
                Ay=cyloc+(-(sy-1)/2:(sy-1)/2);
            else %pair
                Ay=cyloc+(-sy/2+1:sy/2);
            end
            
            posx=Ax+round(obj.pixelRatio*obj.posX);
            posy=Ay+round(obj.pixelRatio*obj.posY);
            M(posy,posx)=M(posy,posx)+obj.matrix;
            inT=toc(inT);
            disp(['Time inside : ' num2str(inT)]);
        end
        
        function M=minus(M,obj)
            [LM,CM]=size(M);
            if mod(LM,2) %odd
                cxloc=(LM-1)/2+1;
            else % even
                cxloc=LM/2;
            end
            if mod(CM,2) %odd
                cyloc=(CM-1)/2+1;
            else % even
                cyloc=CM/2;
            end
            
            [sx,sy]=size(obj.matrix);
            if mod(sx,2) %impair
                Ax=cxloc+(-(sx-1)/2:(sx-1)/2);
            else %pair
                Ax=cxloc+(-sx/2+1:sx/2);
            end
            if mod(sy,2) %impair
                Ay=cyloc+(-(sy-1)/2:(sy-1)/2);
            else %pair
                Ay=cyloc+(-sy/2+1:sy/2);
            end
            
            posx=Ax+round(obj.pixelRatio*obj.posX);
            posy=Ay+round(obj.pixelRatio*obj.posY);
            M(posy,posx)=M(posy,posx).*(~logical(abs(obj.matrix)));
        end
        
        function disp(obj,what)
            
            if nargin==1
                figure
                colormap(parula(1024));
                imagesc(abs(obj.matrix));
                title('Reflexivity');
                axis equal
                colorbar
                
                figure
                colormap(jet(1024));
                tmp=obj.matrix;
                tmp((tmp==0))=NaN;
                imagesc(angle(tmp));
                title('Phase');
                axis equal
                colorbar
                
            elseif nargin==2
                switch what
                    case 'r'
                        figure
                        colormap(parula(1024));
                        imagesc(abs(obj.matrix));
                        title('Reflexivity');
                        axis equal
                        colorbar
                        
                    case 'p'
                        figure
                        colormap(jet(1024));
                        tmp=obj.matrix;
                        tmp((tmp==0))=NaN;
                        imagesc(angle(tmp));
                        title('Phase');
                        axis equal
                        colorbar
                    otherwise
                        disp('not valid argument');
                end
            end
        end
        
        %         function display(obj)
        %             disp('___SEGMENT___');
        %             disp(['Number of sides : ', num2str(obj.nSides)]);
        %             disp(['Size of the matrix : ', num2str(obj.nPixels)]);
        %
        %         end
        
    end
    
end

