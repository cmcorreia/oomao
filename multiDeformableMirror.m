classdef multiDeformableMirror < handle
    % DEFORMABLEMIRROR Create a stack of deformableMirror object at
    % different altitude
    
    properties
        % number of DMs
        nDm;
        % Array of DMs
        dms;
        % DMs Altitude
        zLocations;
        % # of total valid actuators (full dms stack)
        nValidActuators;
        % overSampling factor for the DM Grids
        overSampling;
        % telescope
        tel;
        % nDeformableMirror tag
        tag = 'N DEFORMABLE MIRROR';
        % The Poke Matrix
        thePokeMatrix;
        % commands
        coefs
        %valid actuators
        validActuators
    end
    
    properties (Access=private)
        log;
    end
    
    methods
        
        %% Constructor
        function obj = multiDeformableMirror(zLocations,nActuators,tel,varargin)
            p = inputParser;
            addRequired(p,'zLocations', @isnumeric);
            addRequired(p,'nActuators', @isnumeric);
            addOptional(p,'tel', telescope(-1), @(x) isa(x,'telescopeAbstract') );
            addParameter(p,'overSampling',[], @isnumeric);
            addParameter(p,'modes', {}, @(x) isnumeric(x) || ...
                (isa(x,'influenceFunction') || isa(x,'gaussianInfluenceFunction') ...
                || isa(x,'splineInfluenceFunction') || isa(x,'zernike') ...
                || isa(x,'hexagonalPistonTipTilt')) || isa(x,'differentialGaussianInfluenceFunction'));
            addParameter(p,'resolution', [], @isnumeric);
            addParameter(p,'validActuators', {}, @iscell);
            addParameter(p,'offsets', [] , @isnumeric);
            addParameter(p,'diameters', [], @isnumeric);
            parse(p,zLocations,nActuators, tel, varargin{:});
            
            obj.nDm = length(p.Results.nActuators);
            obj.zLocations =  p.Results.zLocations;
            obj.tel = p.Results.tel;
            
            obj.overSampling = p.Results.overSampling;
            if isempty(obj.overSampling)
                obj.overSampling = 2*ones(obj.nDm,1);
            end
            if isempty(p.Results.validActuators)
                validActuators = {};
                for kDm = 1:obj.nDm
                    validActuators{kDm} = true(nActuators(kDm));
                end
            else
                for kDm = 1:obj.nDm
                    validActuators{kDm} = p.Results.validActuators{kDm};
                end
            end
            
            for kDM = 1:obj.nDm
                D_m = obj.tel.D + 2*obj.zLocations.*tan(0.5*obj.tel.fieldOfView);
                if isempty(p.Results.resolution)
                    do = obj.tel.D/(obj.tel.resolution-1);
                    layersNPixel = 1 + round(D_m./do);
                else
                    layersNPixel = p.Results.resolution;
                end
                obj.dms{kDM} = deformableMirror(p.Results.nActuators(kDM),...
                    'modes',p.Results.modes(kDM),...
                    'validActuator',validActuators{kDM},...
                    'zLocation',obj.zLocations(kDM), ...
                    'resolution',layersNPixel(kDM));
                
                
                %  'resolution',p.Results.resolutions(kDM),...
                % 'validActuator',p.Results.validActuators(kDM),...
                %,...
                %  'offset',p.Results.offsets(kDM),...
                %  'diameter',p.Results.diameters(kDM));
            end
            if isempty(p.Results.validActuators)
                obj.setValidActuators
            end
            nValidActuators = 0;
            for kDm = 1:obj.nDm
                nValidActuators = nValidActuators + obj.dms{kDm}.nValidActuator;
            end
            obj.nValidActuators = nValidActuators;
            obj.log = logBook.checkIn(obj);
        end
        %
        %% Destructor
        function delete(obj)
            %             if isa(obj.modes,'influenceFunction')
            %                 delete(obj.modes)
            %             end
            checkOut(obj.log,obj)
        end
        
        %% SETS AND GETS
        function set.coefs(obj,val)
            obj.coefs = [];
            for kDm = 1:obj.nDm
                if kDm == 1
                    init = 1;
                end
                range = init:init-1+obj.dms{kDm}.nValidActuator;
                if ~isscalar(val);
                    obj.dms{kDm}.coefs = val(range ,1);
                else
                    obj.dms{kDm}.coefs = ones(length(range),1)*val;
                end
                obj.coefs = [obj.coefs; obj.dms{kDm}.coefs];
                init = init + obj.dms{kDm}.nValidActuator;
            end
            
        end
        function relay(obj,srcs)
            %% RELAY the stack of deformable mirror to source relay
            %
            % relay(obj,srcs) crops the phase of the deformableMirror stack
            % in the src directions and relay the corresponding phase into
            % the properties of the source object(s)
            
            nSrc       = numel(srcs);
            wavenumber = 2*pi/srcs(1).wavelength;
            dmPhase = {};
            for kDm = 1:obj.nDm
                dmPhase{kDm} = -2*obj.dms{kDm}.surface*wavenumber;
            end
            for kSrc=1:nSrc % Browse the srcs array
                src = srcs(kSrc);
                nDm_            = obj.nDm;
                altitude_m      = obj.zLocations;
                tel             = src.opticalPath{1};
                sampler_m       = linspace(-1,1,tel.resolution);
                phase_m         = dmPhase;
                R_              = tel.D/2;
                D_m = tel.D + 2*altitude_m.*tan(0.5*tel.fieldOfView);
                do = tel.D/(tel.resolution-1);
                layersNPixel = 1 + round(D_m./do);
                srcDirectionVector1 = src.directionVector(1);
                srcDirectionVector2 = src.directionVector(2);
                srcHeight = src.height;
                m_origin = tel.origin;
                out = zeros(size(src.amplitude,1),size(src.amplitude,2),obj.nDm);
                nOut = size(out,1);
                try % check whether a parpool exist or can be created
                    p = gcp('nocreate');
                catch
                    p = [];
                end
                if isempty(p)
                    for kDm = 1:nDm_
                        layerSampling_m = D_m(kDm)*0.5*linspace(-1,1,layersNPixel(kDm));
                        %                             disp(kLayer)
                        height = altitude_m(kDm);
                        %                             sampling = { layerSampling_m{kLayer} , layerSampling_m{kLayer} };
                        [xs,ys] = meshgrid(layerSampling_m);
                        layerR = R_*(1-height./srcHeight);
                        u = sampler_m*layerR;
                        xc = height.*srcDirectionVector1;
                        yc = height.*srcDirectionVector2;
                        [xi,yi] = meshgrid(u+xc+m_origin(1),u+yc+m_origin(2));
                        %                                 disp( [ xc yc ]+layerR )
                        % out(:,:,kLayer) = spline2(sampling,phase_m{kLayer},{u+yc,u+xc});
                        % size(linear(xs,ys,phase_m{kLayer},xi,yi))
                        % size(out(:,:,kLayer))
                        out(:,:,kDm) = utilities.linear(xs,ys,phase_m{kDm},xi,yi);
                        
                        %                                 F = TriScatteredInterp(xs(:),ys(:),phase_m{kLayer}(:));
                        % out(:,:,kLayer) = F(xi,yi);
                    end
                else
                    parfor kDm = 1:nDm_
                        layerSampling_m = D_m(kDm)*0.5*linspace(-1,1,layersNPixel(kDm));
                        %                             disp(kLayer)
                        height = altitude_m(kDm);
                        %                             sampling = { layerSampling_m{kLayer} , layerSampling_m{kLayer} };
                        [xs,ys] = meshgrid(layerSampling_m);
                        layerR = R_*(1-height./srcHeight);
                        u = sampler_m*layerR;
                        xc = height.*srcDirectionVector1;
                        yc = height.*srcDirectionVector2;
                        [xi,yi] = meshgrid(u+xc+m_origin(1),u+yc+m_origin(2));
                        %                                 disp( [ xc yc ]+layerR )
                        % out(:,:,kLayer) = spline2(sampling,phase_m{kLayer},{u+yc,u+xc});
                        % size(linear(xs,ys,phase_m{kLayer},xi,yi))
                        % size(out(:,:,kLayer))
                        
                        out(:,:,kDm) = utilities.linear(xs,ys,phase_m{kDm},xi,yi);
                    end
                end
                out = sum(out,3);
                %out = (tel.phaseScreenWavelength/src.wavelength)*out;
                
                % relay integrated phase to the source
                srcs(kSrc).phase = out;
                srcs(kSrc).amplitude = 1;
            end
        end
        
        function displayActuatorLayout(obj,src)
            figure
            color = {'r','k','g','b','c'};
            for kDm = obj.nDm:-1:1
                c = obj.dms{kDm}.modes.actuatorCoord;
                scatter(real(c(:)), imag(c(:)),color{kDm})
                hold on
                polar(linspace(0,2*pi,100),0.5*obj.tel.diameterAt(obj.zLocations(kDm)*ones(1,100)),color{kDm});
                
                if nargin > 1 & kDm > 1
                    o = linspace(0,2*pi,101)';
                    xP = obj.tel.D*cos(o)/2;
                    yP = obj.tel.D*sin(o)/2;
                    plot(xP,yP,'color',ones(1,3)*0.8)
                    if ~isempty(src)
                        kLayer = 1;
                        for kSrc=1:numel(src)
                            q = 1 - obj.zLocations(kDm)/src(kSrc).height;
                            xSrc = src(kSrc).directionVector(1).*...
                                obj.zLocations(kDm);
                            ySrc = src(kSrc).directionVector(2).*...
                                obj.zLocations(kDm);
                            plot(xSrc+xP*q,ySrc+yP*q,'color',ones(1,3)*0.8)
                        end
                    end
                end
            end
        end
        
        %% 
        function validActuators = get.validActuators(obj)
            for kDm = 1:obj.nDm
                validActuators{kDm} = obj.dms{kDm}.validActuator;
            end
        end
        %%
        function setValidActuators(obj,val)
            if nargin == 2
                for kDm = 1:obj.nDm
                    obj.dms{kDm}.validActuator = val{kDm};
                end
            else
                validActuators = logical(abs(obj.dms{kDm}.modes.actuatorCoord)<(obj.tel.diameterAt(obj.dms{kDm}.zLocation)/2)*1.2);
                obj.dms{kDm}.validActuator = validActuators;
            end
        end
        function calibMultiDmCell = calibrationMultiDm(obj,sensor,srcs,calibDmStroke,steps,varargin)
            calibMultiDmCell = calibration(obj,sensor,srcs,calibDmStroke,steps,varargin{:});
        end
        
        function calibMultiDmCell = calibration(obj,sensor,srcs,calibDmStroke,steps,varargin)
            nSrcs = numel(srcs);
            calibMultiDmCell = cell(obj.nDm,nSrcs);
            if nargin<5 || isempty(steps)
                steps = 1;
                if obj.dms{1}.nValidActuator>1000
                    steps  = 25;
                end
            end
            for kDm = 1:obj.nDm
                obj.dms{kDm}.coefs = 0;
            end
            %% Calibration
            for kDm = 1:obj.nDm
                if ~(obj.dms{kDm}.zLocation == 0)
                    nSrc = numel(srcs);
                else
                    nSrc = 1;
                end
                for kSrc = 1:nSrc
                    src = srcs(kSrc);
                    if isa(sensor,'pyramid') || sensor.lenslets.nLenslet>1
                        src = src*obj.dms{kDm}*sensor;
                        calibDmCommands = speye(obj.dms{kDm}.nValidActuator)*calibDmStroke;
                        tId = tic;
                        if steps==1
                            obj.dms{kDm}.coefs = calibDmCommands;
                            +src;
                            sp = sensor.slopes;
                            obj.dms{kDm}.coefs = -calibDmCommands;
                            +src;
                            sm = sensor.slopes;
                            pokeMatrix = 0.5*(sp-sm);
                            
                            nMode = size(calibDmCommands,2);
                            slopesStats = zeros(nMode,3);
                            slopesStats(:,1) = max(sp);
                            slopesStats(:,2) = min(sp);
                            slopesStats(:,3) = std(sp);
                        else
                            nMode = size(calibDmCommands,2);
                            nC              = floor(nMode/steps);
                            u               = 0;
                            pokeMatrix  = zeros(sensor.nSlope,nMode);
                            slopesStats = zeros(3,nMode);
                            %                     fprintf(' . actuators range:          ')
                            h = waitbar(0,'DM/WFS calibration ...');
                            while u(end)<nMode
                                u = u(end)+1:min(u(end)+nC,nMode);
                                waitbar(u(end)/nMode)
                                fprintf('\b\b\b\b\b\b\b\b\b%4d:%4d',u(1),u(end))
                                obj.dms{kDm}.coefs = calibDmCommands(:,u);
                                +src;
                                sp = sensor.slopes;
                                obj.dms{kDm}.coefs = -calibDmCommands(:,u);
                                +src;
                                sm = sensor.slopes;
                                pokeMatrix(:,u) = 0.5*(sp-sm);
                                slopesStats(u,1) = max(sp);
                                slopesStats(u,2) = min(sp);
                                slopesStats(u,3) = std(sp);
                            end
                            close(h)
                        end
                        elt = toc(tId);
                        
                        %                 pokeMatrix = src.wavelength*pokeMatrix./calibDmStroke;
                        fprintf('__ Poke Matrix Stats ___\n')
                        fprintf(' . computing time: %5.2fs\n',elt)
                        fprintf(' . size: %dx%d\n',size(pokeMatrix))
                        nonZerosIndex = pokeMatrix(:)~=0;
                        nonZeros = sum(nonZerosIndex);
                        fprintf(' . non zeros values: %d i.e. %4.2f%%\n',nonZeros,100*nonZeros/numel(pokeMatrix))
                        pokeMatrixNZ = pokeMatrix(nonZerosIndex);
                        fprintf(' . min. and max. values: [%5.2f,%5.2f]\n',max(pokeMatrixNZ),min(pokeMatrixNZ))
                        fprintf(' . mean and median of absolute values: [%5.2f,%5.2f]\n',mean(abs(pokeMatrixNZ)),median(abs(pokeMatrixNZ)))
                        fprintf('________________________\n')
                        
                        pokeMatrix = pokeMatrix*diag((1./diag(calibDmStroke)));
                        
                        if isnumeric(obj.dms{kDm}.modes)
                            calibMultiDmCell{kDm,kSrc} = calibrationVault(pokeMatrix,obj.dms{kDm}.modes,src.mask,slopesStats, varargin{:});
                        else
                            calibMultiDmCell{kDm,kSrc} = calibrationVault(pokeMatrix,obj.dms{kDm}.modes.modes,src.mask,slopesStats, varargin{:});
                        end
                    else
                        %%% Tip-Tilt sensor calibration
                        
                        tel = src.opticalPath{1};
                        zern = zernike(tel,2:3);
                        zern.c = eye(2)*calibDmStroke;%src.wavelength/4;
                        src = src.*zern;
                        %                 buf = reshape(src.phase,tel.resolution,2*tel.resolution);
                        
                        % zernike projection onto DM influence functions
                        obj.dms{kDm} = obj.dms{kDm}\src;
                        %                 src = src.*tel*obj;
                        dmTtCoefs = obj.dms{kDm}.coefs;
                        %
                        %                 buf = [buf;reshape(src.phase,tel.resolution,2*tel.resolution)];
                        %                 figure
                        %                 imagesc(buf)
                        %                 axis square
                        %                 colorbar
                        
                        %                 obj.coefs = dmTtCoefs;
                        src = src.*tel*obj.dms{kDm}*sensor;
                        pokeTipTilt = sensor.slopes/calibDmStroke;
                        calib = calibrationVault(pokeTipTilt);
                        calib.spaceJump = dmTtCoefs/calibDmStroke;
                        calibMultiDmCell{kDm,kSrc} = calib;
                    end
                    obj.dms{kDm}.coefs = 0;
                end
            end
        end
    end
end
