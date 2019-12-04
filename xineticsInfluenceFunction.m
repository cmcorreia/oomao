classdef xineticsInfluenceFunction  < handle
    %    Returns an accurate model of Xinetics Inc. deformable mirror
    %   influence function shape. The model is empiric and based on high
    %   spatial resolution measurements of the average influence function of
    %   the deformable mirror of ALTAIR, Gemini North telescope adaptive
    %   optics system and on the very same DM on the Keck AO system.
    %
    %   Mechanical coupling for the Xinetics influence function is 11% and
    %   cannot be modified in the model.
    %
    % INPUTS name | type | unit
    %
    %   RESOLUTION | SCALAR REAL | PIXELS
    %   This is the actuator pitch length expressed in pixels. In other
    %   words, RESOLUTION=pitch/dxp where dxp would be the pixel size in
    %   the DM plane. Does not have to be an integer, can be anything > 0.
    %
    %   DIMMAT | SCALAR INTEGER | PIXELS
    %   Matrix size. Must be at least 5 times the RESOLUTION, otherwise the influence
    %   function is truncated on the edges and is not complete.
    %
    % OPTIONAL INPUTS name | type | unit | default
    %
    %   influenceCentre | SCALAR REAL | PIXELS | 0,0
    %   Indicates where must be located the maximum of the influence function,
    %   relative to the pixel coordinate influenceCentre/2+1, which is the
    %   default location of the maximum (=1).
    %
    % OUTPUT type | unit
    %
    %   REAL ARRAY(DIMMAT,DIMMAT) | max=1
    %   The influence function model. Its maximum is located on the pixel
    %   [DIMMAT/2,DIMMAT/2], unless DX and DY are set to a specific value.
  
    properties       
        % mechanicalCoupling
        mechCoupling;
        % spline polynomials coefs
        splineP;
        % path listener
        xineticsListener;
        % modes
        modes;
        % actuators coordinates
        actuatorCoord;
        % condition for influence function centring
        influenceCentre;
        % influence function tag
        tag = 'XINETICS INFLUENCE FUN';
    end
    
    properties (SetObservable=true)
        % path
        func;
        pitch; % actuator pitch in physical units
        p; % function parameter
    end
    
    
    properties (Access=private)
        % influence function display handle;
        displayHandle
        log
    end
    
    methods
        
        %% Constructor
        function obj = xineticsInfluenceFunction(pitch,varargin)
            inputs = inputParser;
            inputs.addRequired('pitch',@isnumeric);
            inputs.addParameter('influenceCentre',0,@isnumeric);
            inputs.parse(pitch,varargin{:});

            
            obj.pitch           = inputs.Results.pitch;
            obj.influenceCentre = inputs.Results.influenceCentre;   
            obj.mechCoupling    = 0.11;
            
            obj.xineticsListener = addlistener(obj,'func','PostSet',...
                @(src,evnt) obj.show );
            obj.xineticsListener.Enabled = false;
            obj.log = logBook.checkIn(obj);
            
            
            % List of parameters: see KAON 1151
            m = 0.180267421;
            obj.p = [2.24506, 6.28464*m^2,-18.1956*m^4,31.2025*m^6,76.9336,...
                -39.7956,m];
            obj.splineP         = 0.11;
        end
        
        %% Destructor
        function delete(obj)
            if ishandle(obj.displayHandle)
                delete(get(obj.displayHandle,'parent'));
            end
            checkOut(obj.log,obj)
        end
        function display(obj)
            %% DISP Display object information
            %
            % disp(obj) prints information about the influenceFunction
            % object
            fprintf('___ %s ___\n',obj.tag)
            fprintf('  . mechanical coupling: %3.1f\n',...
                obj.mechCoupling);  
            fprintf('  . pitch in meters: %3.1f\n',...
                obj.pitch);  
            fprintf('----------------------------------------------------\n')
        end
        
        function set.splineP(obj,val)
                                                
            % Define the geometry
            c = 1/sqrt(log(1/val));
            df = 1e-10;
            mx = sqrt(-log(df)*c^2);
            r = linspace(-mx,mx,1001);          
            %[X,Y] = meshgrid(x);
            %r     = hypot(X,Y);
            
            % Define the M(r) function
            tmp     = -150*obj.p(7)^(8)*r.^8;
            w       = find(r.^(8) < 1/(3*obj.p(7)^8));
            mask    = 0*tmp;
            mask(w) = exp(tmp(w));
            % Define sub function
            e  = (obj.p(1) +obj. p(2)*r.^2 + obj.p(3)*r.^4 + obj.p(4)*r.^6).*mask;
            re = (abs(r).^e).*obj.p(7).^e;
            % Get the influence function model
            f  = exp(-obj.p(5)*re).*(1+obj.p(6)*re).*mask;     
            
            % Define the interpolation class
            obj.func(:,1) = r*obj.pitch;
            obj.func(:,2) = f;
            obj.splineP   = spline(r*obj.pitch,f);
            
        end
        
        
        function out = mtimes(obj,c)
            %% MTIMES Influence function multiplication
            %
            % v = obj*c multiply the influence functions in obj by the
            % vector c
            
            out = obj.modes*c;
        end
        function out = mldivide(obj,c)
            %% MLDIVIDE Influence function projection
            %
            % v = obj\c project the influence functions onto the vector c
            
            out = obj.modes\c;
        end
              
        function show(obj,varargin)
            %% SHOW 
            
            %figure
                plot(obj.func(:,1),obj.func(:,2))
                hold on
                plot(obj.pitch*ones(100,1), linspace(0,1,100),'k--')
                plot(-obj.pitch*ones(100,1), linspace(0,1,100),'k--')
                plot(obj.func(:,1), obj.mechCoupling*ones(size(obj.func(:,1))),'k--')
                axis tight
                xlabel('position [m]')
                ylabel('normalised influence function')
            
        end
        
        function setInfluenceFunction(obj,nIF,resolution,validActuator,ratioTelDm,offset, diameter)
            %% SETINFLUENCEFUNCTION
                         
            
            if isempty(obj.actuatorCoord)
                if obj.influenceCentre==0
                    xIF         = linspace(-1,1,nIF)*(nIF-1)/2*obj.pitch - offset(1);
                    yIF         = linspace(-1,1,nIF)*(nIF-1)/2*obj.pitch - offset(2);
                    [xIF2,yIF2] = ndgrid(xIF,yIF);
                    obj.actuatorCoord = yIF2 + 1i*flip(xIF2);                                        
                    u0 = ratioTelDm.*linspace(-1,1,resolution)*(nIF-1)/2*obj.pitch; % scaled by telescope diamter
                else
                    xIF               = 1:nIF;
                    yIF               = 1:nIF;
                    [xIF2,yIF2]       = ndgrid(xIF,yIF);
                    obj.actuatorCoord = xIF2 + 1i*yIF2;                    
                    u0                = 1:nIF;
                end
                
                nValid = sum(validActuator(:));
                kIF = 0;
                
                u           = bsxfun( @minus , u0' , xIF );
                wu          = zeros(resolution,nIF);
                index_u     = u >= -obj.func(end,1) & u <= obj.func(end,1);
                nu          = sum(index_u(:));
                wu(index_u) = ppval(obj.splineP,u(index_u));
                
                v           = bsxfun( @minus , u0' , yIF);
                wv          = zeros(resolution,nIF);
                index_v     = v >= -obj.func(end,1) & v <= obj.func(end,1);
                nv          = sum(index_v(:));
                wv(index_v) = ppval(obj.splineP,v(index_v));
                
                
                % Sparse allocation
                m_modes = spalloc(resolution^2,nValid,nu*nv);
                indIF                   = 1:nIF^2;
                indIF(~validActuator)   = [];
                [iIF,jIF]               = ind2sub([nIF,nIF],indIF);
                kIF                     = 1:nValid;
                wv                      = sparse(wv(:,iIF(kIF)));
                wu                      = sparse(wu(:,jIF(kIF)));
                
                fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    \n',nValid)
                for kIF = 1:nValid % parfor doesn't work with sparse matrix!
                    fprintf('\b\b\b\b%4d',kIF)
                    buffer = wv(:,kIF)*wu(:,kIF)';
                    m_modes(:,kIF) = buffer(:);
                end
                fprintf('\n')
                obj.modes = m_modes;
                
            else
                
                xIF = real(obj.actuatorCoord(:)') - offset(1);
                yIF = imag(obj.actuatorCoord(:)') - offset(2);
                               
                if ~isempty(diameter)
                    u0 = linspace(-1,1,resolution)*diameter/2; % normalised by the equivalent optical diameter of the DM
                else
                    u0 = linspace(-1,1,resolution)*max(xIF); % normalised by max of the x coordinate of the DM
                end
                nValid = sum(validActuator(:));
                kIF = 0;
                
                u       = bsxfun( @minus , u0' , xIF );
                wu      = zeros(resolution,nIF);
                index_u = u >= -obj.func(end,1) & u <= obj.func(end,1);                
                wu(index_u) = ppval(obj.splineP,u(index_u));
               
                v       = bsxfun( @minus , u0' , yIF);
                wv      = zeros(resolution,nIF);
                index_v = v >= -obj.func(end,1) & v <= obj.func(end,1);               
                wv(index_v) = ppval(obj.splineP,v(index_v));
                
                expectedNnz = max(sum(index_u))*max(sum(index_v))*nIF;
                add(obj.log,obj,sprintf('Expected non-zeros: %d',expectedNnz))
                add(obj.log,obj,sprintf('Computing the %d 2D DM zonal modes...',nValid))
               
                fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    ',nValid)
                s_i = zeros(expectedNnz,1);
                s_j = zeros(expectedNnz,1);
                s_s = zeros(expectedNnz,1);
                index = 0;
                for kIF = 1:nIF
                    fprintf('\b\b\b\b%4d',kIF)
                    buffer      = wv(:,kIF)*wu(:,kIF)';
                    [i_,~,s_]   = find(buffer(:));
                    n           = length(i_);
                    index       = (1:n) + index(end);
                    s_i(index)  = i_;
                    s_s(index)  = s_;
                    s_j(index)  = ones(n,1)*kIF;
                    
                end
                fprintf('\n')
                                
                index       = 1:index(end);
                obj.modes = sparse(s_i(index),s_j(index),s_s(index),resolution^2,nValid);
                add(obj.log,obj,sprintf('Actual non-zeros: %d',nnz(obj.modes)))
            end                
        end    
    end
end