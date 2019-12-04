classdef gaussianInfluenceFunction < handle
    % INFLUENCEFUNCTION Create an influence function object
    %
    % obj = influenceFunction('monotonic',mech) creates a cubic Bezier
    % influence function monotically decreasing from 1 to 0 with the
    % mechanical coupling mech
    %
    % obj = influenceFunction('overshoot',mech) creates a cubic Bezier
    % influence function with a negative overshoot and the mechanical
    % coupling mech
    %
    % Try show(obj) to view a cut of the influence function
    %
    % Edit 19/11/2018 C. T. Heritier
    % It is now possible to apply the following transformation to the DM by
    % using the parameter 'misReg' considering that the Inflence Functions 
    % coordinates are input. Otherwise, misReg has no effect.
    % misReg must be a structure and can contain the followgin fields:
    %
    %   misReg.rotAngle     Global rotation of the coordinates in degree
    %   misReg.shiftX       Shift in X in m
    %   misReg.shiftY       Shift in Y in m
    %   misReg.magnRadial   Radial Magnification in percentage of diameter
    %   misReg.magnNormal   Normal Magnification in percentage of diameter
    %   misReg.anamAngle    Anamorphosis angle in degrees
    %
    %   In a case of a magnification, the shape of the influence function
    %   is updated to ensure a constant mechanical coupling of the IF.
    properties
        % mechanicalCoupling
        mechCoupling;
        % spline polynomials coefs
        splineP;
        % path listener
        gaussianListener;
        % modes
        modes
        % actuators coordinates
        actuatorCoord;
        % condition for influence function centring
        influenceCentre;
        % influence function tag
        tag = 'GAUSSIAN INFLUENCE FUN';
    end
    
    properties (SetObservable=true)
        % path
        gaussian;
        pitch; % actuator pitch in physical units
    end
    
    
    properties (Access=private)
        % influence function display handle;
        displayHandle
        log
    end
    
    methods
        
        %% Constructor
        function obj = gaussianInfluenceFunction(mech,pitch, influenceCentre)
            if ~exist('pitch','var')
                pitch = 1; % default pitch value
            end
            if ~exist('influenceCentre','var')
                influenceCentre = 0;
            end
            
            obj.mechCoupling = mech;
            obj.pitch = pitch;
            obj.influenceCentre = influenceCentre;
            obj.gaussianListener = addlistener(obj,'gaussian','PostSet',...
                @(src,evnt) obj.show );
            obj.gaussianListener.Enabled = false;
            obj.log = logBook.checkIn(obj);
            obj.splineP = mech;
            
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
        
        %% Set spline
        function set.splineP(obj,val)
            
            % Notes: only using simple setup: x axis is linearly defined,
            % unlike in bezier (influenceFunction) class.  May change this
            % later.
            
            c = 1/sqrt(log(1/val));
            df = 1e-10;
            mx = sqrt(-log(df)*c^2);
            x = linspace(-mx,mx,1001);
            f = exp(-x.^2/c^2);
            
            obj.gaussian(:,1) = x*obj.pitch;
            obj.gaussian(:,2) = f;
            obj.splineP = spline(x*obj.pitch,f);
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
            plot(obj.gaussian(:,1),obj.gaussian(:,2))
            hold on
            plot(obj.pitch*ones(100,1), linspace(0,1,100),'k--')
            plot(-obj.pitch*ones(100,1), linspace(0,1,100),'k--')
            plot(obj.gaussian(:,1), obj.mechCoupling*ones(size(obj.gaussian(:,1))),'k--')
            axis tight
            xlabel('position [m]')
            ylabel('normalised influence function')
            
        end
        
        function setInfluenceFunction(obj,nIF,resolution,validActuator,ratioTelDm,offset, diameter,misReg)
            %% SETINFLUENCEFUNCTION
            
            %             if nargin<5
            %                 ratioTelDm = 1;
            %             end
            %             z = linspace(-1,1,nIF)*(nIF-1)/2;
            
            if ~isempty(misReg) && nnz(offset)~=0
                error('You cannot use both parameters offset and misReg to apply a shift, use only one of them')
            end
            
            if ~isempty(misReg)
                % default values
                if ~isfield(misReg,'shiftX')
                    misReg.shiftX=0;
                end
                if ~isfield(misReg,'shiftY')
                    misReg.shiftY=0;
                end
                
                if ~isfield(misReg,'rotAngle')
                    misReg.rotAngle=0;
                end
                if ~isfield(misReg,'magnNorm')
                    misReg.magnNorm=100;
                end
                if ~isfield(misReg,'magnRadial')
                    misReg.magnRadial=100;
                end
                if ~isfield(misReg,'anamAngle')
                    misReg.anamAngle=0;
                end
            end
            
            if isempty(misReg)
                % backward compatibility with 'offset' parameter
                if nnz(offset)~=0
                    misReg.shiftX=offset(1);        %Shift X
                    misReg.shiftY=offset(2);        %Shift Y
                    misReg.rotAngle=0;              %Rotation angle
                    misReg.magnNorm=100;              %Radial magnification
                    misReg.magnRadial=100;            %Normal magnification
                    misReg.anamAngle=0;             %Anamorphose Angle
                    
                %default case
                else
                    misReg.shiftX=0;                %Shift X
                    misReg.shiftY=0;                %Shift Y
                    misReg.rotAngle=0;              %Rotation angle
                    misReg.magnNorm=100;              %Radial magnification
                    misReg.magnRadial=100;            %Normal magnification
                    misReg.anamAngle=0;             %Anamorphose Angle
                end
            end
            
            if isempty(obj.actuatorCoord)
                
                if obj.influenceCentre==0
                    xIF = linspace(-1,1,nIF)*(nIF-1)/2*obj.pitch - offset(1);
                    yIF = linspace(-1,1,nIF)*(nIF-1)/2*obj.pitch - offset(2);
                    [xIF2,yIF2] = ndgrid(xIF,yIF);
                    obj.actuatorCoord = yIF2 + 1i*flip(xIF2);
                    
                    %u0 = ratioTelDm.*linspace(-1,1,resolution)*(nIF-1)/2*(resolution-1)/resolution*obj.pitch;
                    u0 = ratioTelDm.*linspace(-1,1,resolution)*(nIF-1)/2*obj.pitch; % scaled by telescope diamter
                else
                    xIF = 1:nIF;
                    yIF = 1:nIF;
                    [xIF2,yIF2] = ndgrid(xIF,yIF);
                    obj.actuatorCoord = xIF2 + 1i*yIF2;
                    
                    u0 = 1:nIF;
                end
                
                nValid = sum(validActuator(:));
                kIF = 0;
                
                u = bsxfun( @minus , u0' , xIF );
                wu = zeros(resolution,nIF);
                index_v = u >= -obj.gaussian(end,1) & u <= obj.gaussian(end,1);
                nu = sum(index_v(:));
                wu(index_v) = ppval(obj.splineP,u(index_v));
                
                v = bsxfun( @minus , u0' , yIF);
                wv = zeros(resolution,nIF);
                index_v = v >= -obj.gaussian(end,1) & v <= obj.gaussian(end,1);
                nv = sum(index_v(:));
                wv(index_v) = ppval(obj.splineP,v(index_v));
                
                m_modes = spalloc(resolution^2,nValid,nu*nv);
                
                %             fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    \n',nValid)
                %             for jIF = 1:nIF
                %
                %                 for iIF = 1:nIF
                %                     if validActuator(iIF,jIF)
                %                         buffer = wv(:,iIF)*wu(:,jIF)';
                %                         kIF = kIF + 1;
                %                         obj.modes(:,kIF) = buffer(:);
                %                                                 fprintf('\b\b\b\b%4d',kIF)
                %                     end
                %                 end
                %
                %             end
                %             fprintf('\n')
                
                indIF = 1:nIF^2;
                indIF(~validActuator) = [];
                [iIF,jIF] = ind2sub([nIF,nIF],indIF);
                kIF = 1:nValid;
                wv = sparse(wv(:,iIF(kIF)));
                wu = sparse(wu(:,jIF(kIF)));
                fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    \n',nValid)
                for kIF = 1:nValid % parfor doesn't work with sparse matrix!
                    fprintf('\b\b\b\b%4d',kIF)
                    buffer = wv(:,kIF)*wu(:,kIF)';
                    m_modes(:,kIF) = buffer(:);
                end
                fprintf('\n')
                obj.modes = m_modes;
                
                
            else                    
                % initial coordinates
                xIF0 = real(obj.actuatorCoord(:)');
                yIF0 = imag(obj.actuatorCoord(:)');
                % rotation for anamorphosis
                [xIF1,yIF1]=lamTools.rotateDM(xIF0,yIF0,misReg.anamAngle*pi/180);
                
                % Non-Symmetric Magnification Normal and radial
                xIF2=xIF1*misReg.magnNorm/100;
                yIF2=yIF1*misReg.magnRadial/100;
                
                % de-rotation for anamorphosis
                [xIF3,yIF3]=lamTools.rotateDM(xIF2,yIF2,-misReg.anamAngle*pi/180);
                
                % rotation of the coordinates
                [xIF4,yIF4]=lamTools.rotateDM(xIF3,yIF3,misReg.rotAngle*pi/180);
                
                % shift of the coordinates
                xIF = xIF4 - misReg.shiftX;
                yIF = yIF4 - misReg.shiftY;
                obj.actuatorCoord = xIF + 1i*(yIF);
                fprintf(' -- Rotation %.2f °-- Shift X %.0f %% -- Shift Y -- %.0f %%',misReg.rotAngle,100*misReg.shiftX/obj.pitch,100*misReg.shiftY/obj.pitch)
                fprintf(' -- Radial Magnificaton %.0f %% -- Normal Magnification %.0f %% -- Anamorposis angle %.2f \n',misReg.magnRadial,misReg.magnNorm,misReg.anamAngle)
                
                % case with no magnification
                if misReg.magnNorm== 100 && misReg.magnRadial==100;
                    %u0 = linspace(-1,1,resolution)*max(abs(obj.actuatorCoord(:)));
                    if ~isempty(diameter)
                        u0 = linspace(-1,1,resolution)*diameter/2; % normalised by the equivalent optical diameter of the DM
                    else
                        u0 = linspace(-1,1,resolution)*max(xIF); % normalised by max of the x coordinate of the DM
                    end
                    nValid = sum(validActuator(:));
                    kIF = 0;
                    
                    u = bsxfun( @minus , u0' , xIF );
                    wu = zeros(resolution,nIF);
                    index_u = u >= -obj.gaussian(end,1) & u <= obj.gaussian(end,1);
                    %                 nu = sum(index_u(:));
                    wu(index_u) = ppval(obj.splineP,u(index_u));
                    
                    if ~isempty(diameter)
                        u0 = linspace(-1,1,resolution)*diameter/2; % normalised by the equivalent optical diameter of the DM
                    else
                        u0 = linspace(-1,1,resolution)*max(yIF); % normalised by max of the x coordinate of the DM
                    end
                    v = bsxfun( @minus , u0' , yIF);
                    wv = zeros(resolution,nIF);
                    index_v = v >= -obj.gaussian(end,1) & v <= obj.gaussian(end,1);
                    %                 nv = sum(index_v(:));
                    wv(index_v) = ppval(obj.splineP,v(index_v));
                    expectedNnz = max(sum(index_u))*max(sum(index_v))*nIF;
                    add(obj.log,obj,sprintf('Expected non-zeros: %d',expectedNnz))
                    add(obj.log,obj,sprintf('Computing the %d 2D DM zonal modes...',nValid))
                    %                 m_modes = spalloc(resolution^2,nValid,expectedNnz);
                    fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    ',nValid)
                    s_i = zeros(expectedNnz,1);
                    s_j = zeros(expectedNnz,1);
                    s_s = zeros(expectedNnz,1);
                    index = 0;
                    for kIF = 1:nIF
                        fprintf('\b\b\b\b%4d',kIF)
                        buffer = wv(:,kIF)*wu(:,kIF)';
                        [i_,~,s_] = find(buffer(:));
                        n = length(i_);
                        index = (1:n) + index(end);
                        s_i(index) = i_;
                        s_s(index) = s_;
                        s_j(index) = ones(n,1)*kIF;
                        %                     m_modes(:,kIF) = buffer(:);
                    end
                    fprintf('\n')
                    %                 add(obj.log,obj,sprintf('Actual non-zeros: %d',nnz(m_modes)))
                    %                 obj.modes = m_modes;
                    index = 1:index(end);
                    obj.modes = sparse(s_i(index),s_j(index),s_s(index),resolution^2,nValid);
                    add(obj.log,obj,sprintf('Actual non-zeros: %d',nnz(obj.modes)))
                % case with magnification
                else                    
                    if ~isempty(diameter)
                        %recentered
                        u0x = (resolution)/2+xIF*(resolution)/diameter;
                        u0y = (resolution)/2+yIF*(resolution)/diameter;
                    else
                        %recenterered
                        u0x = resolution/2+xIF*resolution/max(xIF);
                        u0y = resolution/2+yIF*resolution/max(yIF);
                        
                    end
                    %
                    nValid = sum(validActuator(:));
                    fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    ',nValid)
                    modes2D=zeros(resolution*resolution,nIF);
                    for kIF=1:nIF
                        
                        fprintf('\b\b\b\b%4d',kIF)
                        
                        x0=u0x(kIF);
                        y0=u0y(kIF);
                        % get the number of actuator along the diameter
                        nActAlongDiameter=(diameter)/obj.pitch;
                        % compute FWHM of the 2D gaussian in each direction
                        cx=(misReg.magnRadial/100)*(resolution/nActAlongDiameter)/sqrt(2*log(1/obj.mechCoupling));
                        cy=(misReg.magnNorm/100)*(resolution/nActAlongDiameter)/sqrt(2*log(1/obj.mechCoupling));
                        % Anamorphosis Angle
                        theta=misReg.anamAngle*pi/180+pi/2;
                        
                        X=linspace(0,1,resolution)*resolution;
                        Y=linspace(0,1,resolution)*resolution;
                        [X,Y]=meshgrid(X,Y);
                        
                        % compute 2D gaussian
                        a = cos(theta)^2/(2*cx^2) + sin(theta)^2/(2*cy^2);
                        b = -sin(2*theta)/(4*cx^2) + sin(2*theta)/(4*cy^2);
                        c = sin(theta)^2/(2*cx^2) + cos(theta)^2/(2*cy^2);
                        
                        Z = exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2));
                        
                        % create 2D mode
                        modes2D(:,kIF)=reshape(Z,1,resolution*resolution);
                        
                    end
                    fprintf(')\n')
                    obj.modes=modes2D;
                end
            end
        end
    end
    
end