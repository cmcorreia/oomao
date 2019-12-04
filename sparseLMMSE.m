classdef sparseLMMSE < handle
    %% SPARSELMMSE MMSE wavefront estimation from SH-WFS slopes
    %
    % sparseLMMSE(wfs,tel,atmModel,guideStar,'mmseStar',science)
    % computes the tomographic reconstructor to estimate the directional
    % phase in the science direction(s) from SH-WFSs meurements based on
    % turbulence model given by atmModel.
    %
    % The tomography process is done iteratively using (minres/cgs)  in
    % MATLAB and you can change method and parameters for the itelative
    % computation by options
    %   'solver' : Iterative method (minres / cgs), default = minres
    %   'MAXIT' : The maximum number of Iteration
    %   'RTOL'  : The tolerance
    % The defaults are MAXIT=50 and RTOL=1e-3.
    %
    % Also, you can get the tomographic reconstructor explicitly by a
    % function of 'getReconstructionMatrix' as
    %    mmse = slopesLinearMMSE(wfs,tel,atmModel,guideStar,'mmseStar',science, other options...);
    %    R = mmse.getReconstructionMatrix();
    %    phi = R*s; % reconstruction 
    % where the input s should be given by in [pixel] and the output phi is
    % in [nm]. Don't use this function for large-scale AO systems!!!!
    %
    % If you want to estimate layered phase at each altitudes (not
    % directional phase for science directions), you change a flag
    %    'layerEstimation' : true/false, default is false
    %
    % In case with LGSs, you have to remove tip/tilt/(focus) from
    % measurements and the reconstructor should be created with taking
    % into account of this removement. This can be done by the option
    %    'isTTRM' : removing tip/tilt (true/false)
    %    'isFocusRM : removing focus (true/false)
    % When these parameter are true, tip/til/focus are removed from
    % measurements, recontructor and reconstructed phase.
    % The default is no removing i.e. isTTRM=false and isFocusRM=false
    %
    % For spatial sampling of esitimated layered phase at different
    % altitudes can be changed with option of 'overSampling'.
    %    spatial sampling = sub-aperture size / overSampling
    %    Default value is overSampling=2 for all altitudes i.e. the spatial
    %    sampling is half of sub-aperture diameter.
    %    overSampling = 1 means the sub-aperture sampling for all altitudes and
    %    overSampling = [1 2 2] means the sub-aperture sampling for the
    %    first layer and half of sub-aperutre sampling for the rest.
    %
    %  Gradient operator assumes 3x3 stensil
    %     (see  shackHartmann.sparseGradientMatrix3x3Stentil)
    %
    % For the regularization tarm, you can choose from
    %    'phaseCovariance' = 'L2' : sparse approximated matrix
    %    'phaseCovariance' = 'real' : actual phase-covariance (slow and
    %    heavey, Don't use this option for large-scale AO systems)
    %
    % Demoscripts for test
    %     sparseLMMSE.demo_OL_LTAO; % Open-loop LTAO script
    %     sparseLMMSE.demo_LayerEstimation; % Open-loop tomography for atmopheric volume
    
    properties
        tag = 'sparseLMMSE';
        % diameter of sub-aperture [m]
        dSub
        % number of sub-aperture along diameter
        nSub
        % diameter of telescope
        D
        % slopes pupil mask
        slopesMask
        % output wavefront pupil mask
        outputWavefrontMask
        % atmosphere mode
        atmModel
        % wave-front sensor object
        wfs;
        % guide star direction
        guideStar
        % MINRES relative tolerance
        RTOL;
        % MINRES maximum number of iteration
        MAXIT;
        % MINRES initial guess
        x0
        % warm start flag
        warmStart = false;
        % flag for TT removal
        isTTRM
        % flag for Focus removal
        isFocusRM
        %
        tPredict
        tMulti
        maskAmplitudeWighted;
        
        % flag for computing the rotation matrix
        rotatedFrame;
        % rotator: Matrix that rotates the slopes from the original x-y
        % referential to a x'-y' per sub-aperture (depends on the LLT)
        rotationMatrix;
        
        % discrete gradient operator (for the SH WFS)
        Gamma;
        GammaT;
        % propagator from layers to gridMask pupil locations in guide-star
        % direction
        H;
        HT;
        % propagator from layers to gridMask pupil locations in science
        % directions
        Hss;
        % phase covariance matrix
        L2;
        % mask for atmosphere grid
        atmGrid;
        % overSampling factor
        overSampling
        % type of phase covariance
        phaseCovariance
        % output phase grid
        outputPhaseGrid
        % solver
        solver
        phaseMask
        % flag for estiamtion of layered phase
        layerEstimation
    end
    
    properties (Dependent)
        % mmse estimates directions
        mmseStar;
        % inverse of noise variance
        iNoiseVar;
    end
    
    properties (Dependent,SetAccess=private)
        % number of guide stars
        nGuideStar
        % number of science stars
        nMmseStar
    end
    
    properties (Access=private)
        log
        p_mmseStar
        p_iNoiseVar = 0;
        % handle to private function A*x
        AmultFun
        wavefrontToMeter;
        c;
        ttrmfunSlope; % function for the TT or/and focus removal
        ttrmfunPhase; % function for the TT or/and focus removal
        ttrmfunSlopeTT; % function for the TT or/and focus removal
    end
    
        
    methods
            %% Constructor
        function obj = sparseLMMSE(wfs,tel,atmModel,guideStar,varargin)
            
            inputs = inputParser;
            inputs.addRequired('wfs', @(x) isa(x,'shackHartmann') );
            inputs.addRequired('tel', @(x) isa(x,'telescope'));
            inputs.addRequired('atmModel', @(x) isa(x,'atmosphere'));
            inputs.addRequired('guideStar', @(x) isa(x,'source'));
            inputs.addOptional('mmseStar', [],@(x) isa(x,'source'));
            inputs.addOptional('solver','minres',@ischar);
            inputs.addOptional('slopesMask', [], @islogical );
            inputs.addOptional('iNoiseVar', 0, @isnumeric );
            inputs.addOptional('isTTRM', false, @islogical );
            inputs.addOptional('isFocusRM', false, @islogical );
            inputs.addOptional('RTOL', 1e-3, @isnumeric );
            inputs.addOptional('MAXIT', 50, @isnumeric );
            inputs.addOptional('isWarmStart', true, @islogical );
            inputs.addOptional('tPredict', 0, @isnumeric );
            inputs.addOptional('tMulti', 0, @isnumeric );
            inputs.addOptional('maskAmplitudeWighted', [], @islogical); % need to be tested....
            inputs.addOptional('overSampling',[], @isnumeric);
            inputs.addOptional('phaseCovariance','L2',@ischar);
            inputs.addOptional('outputWavefrontMask',[],@isnumeric);
            inputs.addOptional('layerEstimation', false, @islogical);
            inputs.addOptional('rotatedFrame', false, @islogical );
            
            inputs.parse(wfs,tel,atmModel,guideStar,varargin{:});

            obj.nSub          = inputs.Results.wfs.lenslets.nLenslet;
            obj.D             = inputs.Results.tel.D;
            obj.dSub          = obj.D/obj.nSub;
            obj.wfs           = inputs.Results.wfs;
            obj.atmModel      = inputs.Results.atmModel;
            obj.guideStar     = inputs.Results.guideStar;
            obj.iNoiseVar     = inputs.Results.iNoiseVar;
            obj.isTTRM        = inputs.Results.isTTRM;
            obj.isFocusRM     = inputs.Results.isFocusRM;
            obj.tPredict      = inputs.Results.tPredict;
            obj.tMulti        = inputs.Results.tMulti(:);
            obj.RTOL          = inputs.Results.RTOL;
            obj.MAXIT         = inputs.Results.MAXIT;
            obj.warmStart     = inputs.Results.isWarmStart;
            obj.maskAmplitudeWighted = inputs.Results.maskAmplitudeWighted;
            obj.phaseCovariance = inputs.Results.phaseCovariance;
            obj.solver        = inputs.Results.solver;
            obj.layerEstimation= inputs.Results.layerEstimation;
            obj.rotatedFrame  = inputs.Results.rotatedFrame;
            
            if obj.tMulti(1) ~= 0
                obj.tMulti(1) = [0; obj.tMulti];
            end
            if ~isempty(inputs.Results.slopesMask)
                obj.slopesMask = inputs.Results.slopesMask;
            else
                obj.slopesMask = repmat( repmat( inputs.Results.wfs.validLenslet(:) , 2, 1), 1, obj.nGuideStar );
            end;
            makePhaseMask(obj);
            obj.p_mmseStar    = inputs.Results.mmseStar;
            if isempty(obj.p_mmseStar)
                obj.p_mmseStar = obj.guideStar;
            end
            
            obj.overSampling = inputs.Results.overSampling;
            if isempty(obj.overSampling)
                obj.overSampling = ones(obj.atmModel.nLayer,1)*2;
            end
            
            obj.log = logBook.checkIn(obj);
            obj.log.verbose = false;
            
            
            if isempty(inputs.Results.outputWavefrontMask)
                obj.outputWavefrontMask = inputs.Results.wfs.validActuator;
            end
            
            [x,y] = meshgrid(linspace(-1,1,obj.nSub+1)*obj.D/2);
            obj.outputPhaseGrid = complex(x(obj.outputWavefrontMask),y(obj.outputWavefrontMask));
            
            %Compute individual compouds
            defineAtmGrid(obj);
            if inputs.Results.rotatedFrame
                setRotationMatrix(obj);
            else
                obj.rotationMatrix = [];
            end
            setGamma(obj);
            setH(obj);
            setHss(obj);
            setL2(obj);
            obj.createTTRMfunc;

            obj.wavefrontToMeter = wfs.lenslets.fieldStopSize*obj.guideStar(1).wavelength...
                /(tel.D/wfs.lenslets.nLenslet)/(wfs.camera.resolution(1)/wfs.lenslets.nLenslet);
        end
        
        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.log)
                checkOut(obj.log,obj)
            end
        end
        
        %% Set/Get mmseStar
        function set.mmseStar(obj,val)
            obj.p_mmseStar = val;
            %add(obj.log,obj,'Computing the mmse/guide stars covariance matrix')
            %phaseToSlopesCovariance(obj);
        end
        function val = get.mmseStar(obj)
            val = obj.p_mmseStar;
        end
        
        %% Set/Get noiseVar
        function set.iNoiseVar(obj,val)
            if isscalar(val) % for scalar
                obj.p_iNoiseVar = val;
            elseif isvector(val) % for vector
                obj.p_iNoiseVar = val(:);
            else % for matrix
                obj.p_iNoiseVar = val;
            end
        end
        function val = get.iNoiseVar(obj)
            val = obj.p_iNoiseVar;
        end
        
        %% Set nGuideStar
        function val = get.nGuideStar(obj)
            val = length(obj.guideStar(1,:,1));
        end
        
        %% Set nMmseStar
        function val = get.nMmseStar(obj)
            val = length(obj.mmseStar(1,:,1));
        end
        
        function makePhaseMask(obj)
            %%  define mask for wavefront of WFS
            nMap     = 2*obj.nSub+1;
            
            [iMap0,jMap0] = ndgrid(1:3);
            gridMask = false(nMap^2,obj.nGuideStar);
            for jLenslet = 1:obj.nSub
                jOffset = 2*(jLenslet-1);
                for iLenslet = 1:obj.nSub
                    index1 = jLenslet+obj.nSub*(iLenslet-1);
                    iOffset = 2*(iLenslet-1);
                    for kGs = 1:obj.nGuideStar
                        if obj.slopesMask(index1,kGs)
                            index2 = jMap0+jOffset+nMap*(iMap0+iOffset-1);
                            gridMask(index2(:),kGs) = true;
                        end
                    end
                end
            end
            obj.phaseMask = gridMask;
        end
        
        %% set a rotation matrix
        function setRotationMatrix(obj)
            rTheta = cell(obj.nGuideStar,1);
            for kGs = 1:obj.nGuideStar
                xL = obj.guideStar(kGs).viewPoint(1);
                yL = obj.guideStar(kGs).viewPoint(2);
                d = obj.D/obj.wfs.lenslets.nLenslet;
                uLenslet = linspace(-1,1,obj.wfs.lenslets.nLenslet)*(obj.D/2-d/2);
                [xLenslet,yLenslet] = meshgrid(uLenslet);
                maskLenslet = obj.wfs.validLenslet;
                xLenslet = xLenslet(maskLenslet);
                yLenslet = yLenslet(maskLenslet);
                
                [oe,~] = cart2pol(xLenslet-xL,yLenslet-yL);
                
                rTheta{kGs} = [diag(cos(oe)) diag(sin(oe)); -diag(sin(oe)) diag(cos(oe)) ];
                %rTheta{kGs} = rTheta{kGs}
            end
            obj.rotationMatrix = blkdiag(rTheta{:});
        end
        
        
        %% set the sparse gradient matrix using a 3x3 stencil
        function setGamma(obj)
            [p_Gamma,phaseMask3x3] = sparseGradientMatrix3x3Stentil(obj.wfs);
            p_Gamma = p_Gamma/2/obj.dSub;
            if ~isempty(obj.maskAmplitudeWighted)
                [p_Gamma,phaseMask3x3] = sparseGradientMatrixAmplitudeWeighted(obj.wfs, obj.maskAmplitudeWighted);
                p_Gamma = p_Gamma/obj.dSub;
            end
            obj.Gamma = cell(obj.nGuideStar,1);
            vL = repmat(obj.wfs.validLenslet(:),2,1);
            for kGs = 1:obj.nGuideStar
                if ~isempty(obj.rotationMatrix)
                    obj.Gamma{kGs} = p_Gamma(:,...
                    obj.phaseMask(phaseMask3x3(:),kGs));
                else
                obj.Gamma{kGs} = p_Gamma(obj.slopesMask(vL,kGs),...
                    obj.phaseMask(phaseMask3x3(:),kGs));
                end
            end
            obj.Gamma = blkdiag(obj.Gamma{:});
            % Model the rotated slopes
            if ~isempty(obj.rotationMatrix)
                obj.Gamma = obj.rotationMatrix*obj.Gamma;
                obj.Gamma = obj.Gamma(obj.slopesMask(obj.slopesMask(:)),:);
            end
            obj.GammaT = obj.Gamma';
        end

        %% set propagator H from GS to WFS
        function setH(obj)
            [x,y] = meshgrid(linspace(-.5,.5,2*obj.nSub+1)*obj.D);
            grid = x+y*1i;
            p_H = cell(obj.nGuideStar,obj.atmModel.nLayer);
            for kGs = 1:obj.nGuideStar
                a = p_bilinearSplineInterpMat(obj,...
                    obj.guideStar(1,kGs,1),...
                    grid(obj.phaseMask(:,kGs)),...
                    0);
                p_H(kGs,:) = a;
            end
            obj.H = cell2mat(p_H);
            obj.HT = obj.H';
        end
        
        %% set propagator Hss from science star to pupil
        function setHss(obj)
            p_Hss = p_bilinearSplineInterpMat(obj,...
                obj.mmseStar,...
                obj.outputPhaseGrid,...
                0);
            obj.Hss = cell2mat(p_Hss);
        end
        
        % set bi-harmonic operator (appox to inverse phase covariance
        % matrix)
        function setL2(obj)
            switch obj.phaseCovariance
                case 'real'
                    p_L2 = cell(obj.atmModel.nLayer,1);
                    for kLayer=1:obj.atmModel.nLayer
                        [x,y] = meshgrid(obj.atmGrid{kLayer,1},obj.atmGrid{kLayer,2});
                        rho = complex(x,y);
                        tmp = phaseStats.covarianceMatrix(rho,slab(obj.atmModel, kLayer));
                        tmp = tmp*(obj.atmModel.wavelength/(2*pi))^2;
                        p_L2{kLayer} = inv(tmp);
                    end
                case 'L2'
                    p_L2 = p_sparseInverseCovarianceMatrix(obj);
            end
            obj.L2 = blkdiag(p_L2{:});
        end
        
        %% Wavefront reconstruction (new version by Y.O 18/3/2017)
        function out = mtimes(obj,wfs)
            if isa(wfs,'shackHartmann')
                obj.c = wfs.slopes(:);
            else
                obj.c = wfs(:);
            end
            
            % remove tip/tilt/focus from slopes if needed
            v = 0;
            for k1=1:length(obj.tMulti)
                for k2=1:obj.nGuideStar
                    v = (1:sum(obj.slopesMask(:,k2)))+v(end);
                    obj.c(v) = obj.ttrmfunSlope{k2}(obj.c(v));
                end
            end
            
            % estimation
            switch obj.solver
                case 'cg'
                    [yy] = cgs(@(x) obj.leftSide(x),obj.rightSide(obj.c),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart); % conjugate-gradient method
                case 'minres'
                    [yy] = minres(@(x) obj.leftSide(x),obj.rightSide(obj.c),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart); % minres
                case 'bicg'
                    [yy] = bicgstab(@(x) obj.leftSide(x),obj.rightSide(obj.c),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart); % minres
            end
            obj.x0 = yy;
            
            % remove tip/tilt/focus from the estimated phase if needed
            if ~obj.layerEstimation
                out = obj.ttrmfunPhase(obj.Hss*yy);
            else
                out = yy;
            end
            
            out = out*obj.wavefrontToMeter;            
        end
        
        %% Compute reconstruction matrix explicitly
        function [R,Left,Right]=getReconstructionMatrix(obj)
                
            % mask
            mask = cell(obj.nGuideStar,1);
            for k=1:obj.nGuideStar
                mask{k} = obj.slopesMask(1:end/2,k);
            end
            
            % filtering
            filter = cell(obj.nGuideStar,2);
            if obj.isTTRM 
                mode = [2 3];
                if obj.isFocusRM
                    mode = [2 3 4];
                end
                for k=1:obj.nGuideStar
                    zer = zernike(1:max(mode),'resolution',obj.nSub,'pupil',mask{k});
                    zx = zer.xDerivative(mask{k},mode);
                    zy = zer.yDerivative(mask{k},mode);
                    
                    filter{k,1} = eye(sum(mask{k}(:)))-zx*pinv(zx);
                    filter{k,2} = eye(sum(mask{k}(:)))-zy*pinv(zy);
                end
            else
                for k=1:obj.nGuideStar
                    filter{k,1} = eye(sum(mask{k}(:)));
                    filter{k,2} = eye(sum(mask{k}(:)));
                end
            end
            filter = blkdiag(filter{:});
            
            if isvector(obj.iNoiseVar)
                Cnn = diag(obj.iNoiseVar);
            elseif isscalar(obj.iNoiseVar)
                Cnn = obj.iNoiseVar*eye(sum(obj.slopesMask(:)));
            else
                Cnn = obj.iNoiseVar;
            end
            
            Right = obj.HT*obj.GammaT*filter'*Cnn;
            
            Left=obj.HT*obj.GammaT*filter'*Cnn*filter*obj.Gamma*obj.H + obj.L2;
            
            if ~obj.layerEstimation
                R = obj.Hss*pinv(Left)*Right*obj.wavefrontToMeter;
            else
                R = pinv(Left)*Right*obj.wavefrontToMeter;
            end
        end
    end
    
    
    %methods (Access=private)
    methods
        
        function createTTRMfunc(obj)
            %% create functins to remove TT or/and Focus

            obj.ttrmfunSlope = cell(length(obj.guideStar),1);
            if ~obj.isTTRM && ~obj.isFocusRM
                for iT=1:length(obj.tMulti)
                    for iGs=1:obj.nGuideStar
                        iN = iGs + obj.nGuideStar*(iT-1);
                        obj.ttrmfunSlope{iN} = @(x) x;
                        obj.ttrmfunSlopeTT{iN} = @(x) x;
                    end
                end
                obj.ttrmfunPhase = @(x) x;
            else
                mode = [];
                if obj.isTTRM
                    mode = [mode 2 3];
                end
                if obj.isFocusRM
                    mode = [mode 4];
                end
                for iT=1:length(obj.tMulti)
                    for iGs=1:obj.nGuideStar
                        mask = reshape(obj.slopesMask(1:end/2,iGs),obj.nSub,obj.nSub);
                        zers = zernike(1:4,...
                            obj.D,...
                            'resolution',obj.nSub,...
                            'pupil',double(mask));
                        RMs = [zers.xDerivative(:,mode); zers.yDerivative(:,mode);];
                        RMs(~obj.slopesMask(:,iGs),:) = [];
                        iRMs = pinv(RMs);
                        iN = iGs + obj.nGuideStar*(iT-1);
                        obj.ttrmfunSlope{iN} = @(x) x-RMs*(iRMs*x);
                        obj.ttrmfunSlopeTT{iN} = @(x) x-iRMs'*(RMs'*x);
                    end
                end
                zerp = zernike(1:4,...
                    obj.D,...
                    'resolution',obj.nSub+1,...
                    'pupil',double(obj.outputWavefrontMask));
                RMp = zerp.modes(:,mode);
                RMp(~obj.outputWavefrontMask,:) = [];
                iRMp = pinv(RMp);
                obj.ttrmfunPhase = @(x) x-RMp*(iRMp*x);
            end
            
        end
        
        function out = leftSide(obj,x)
            %% Computes y=(G'*H'*Cnn^-1*G*H+L'*L) *x
            out = obj.H*x;
            out = obj.Gamma*out;
            v = 0;
            for k=1:obj.nGuideStar
                v = (1:sum(obj.slopesMask(:,k)))+v(end);
                out(v) = obj.ttrmfunSlope{k}(out(v));
            end
            if isvector(obj.iNoiseVar)
                out = obj.iNoiseVar.*out;
            else
                out = obj.iNoiseVar*out;
            end
            v = 0;
            for k=1:obj.nGuideStar
                v = (1:sum(obj.slopesMask(:,k)))+v(end);
                out(v) = obj.ttrmfunSlopeTT{k}(out(v));
            end
            out = obj.GammaT*out;
            out = obj.HT*out;
            out = out + obj.L2*x;
        end
 
        function out = rightSide(obj,x)
            %% Computes  G'*H'*Cnn^-1*s
            if isvector(obj.iNoiseVar)
                out = obj.iNoiseVar.*x;
            else
                out = obj.iNoiseVar*x;
            end
            v = 0;
            for k=1:obj.nGuideStar
                v = (1:sum(obj.slopesMask(:,k)))+v(end);
                out(v) = obj.ttrmfunSlopeTT{k}(out(v));
            end
            out = obj.GammaT*out;
            out = obj.HT*out;
        end
        
        
        function defineAtmGrid(obj)
            %% Define spatial grid for layered phase at each altitude

            vg = [obj.guideStar(1,:,1).directionVector];
            vs = [obj.mmseStar(1,:,1).directionVector];
            vg = vg(1:2,:);
            vs = vs(1:2,:);
            
            obj.atmGrid = cell(obj.atmModel.nLayer,2);
            
            for kLayer = 1:obj.atmModel.nLayer
                
                pitchLayer = obj.dSub/obj.overSampling(kLayer);
                height = obj.atmModel.layer(kLayer).altitude;
                m = 1-height/obj.guideStar(1).height;

                dDirecG = vg*height;
                dDirecS = vs*height;
                dmin = min([dDirecG-obj.D/2*m dDirecS-obj.D/2],[],2);
                dmax = max([dDirecG+obj.D/2*m dDirecS+obj.D/2],[],2);
                
                nPxLayerX = floor((dmax(1)-dmin(1))/pitchLayer)+2;
                nPxLayerY = floor((dmax(2)-dmin(2))/pitchLayer)+2;
                
                Dx = (nPxLayerX-1)*pitchLayer;
                Dy = (nPxLayerY-1)*pitchLayer;
                
                sx = dmin(1)-(Dx-(dmax(1)-dmin(1)))/2;
                sy = dmin(2)-(Dy-(dmax(2)-dmin(2)))/2;
                
                obj.atmGrid{kLayer,1} = linspace(0,1,nPxLayerX)*Dx+sx;
                obj.atmGrid{kLayer,2} = linspace(0,1,nPxLayerY)*Dy+sy;
                
            end
            
        end
        
        
        function H = p_bilinearSplineInterpMat(obj,star,outputCoord,dt)
            %% Bi-linear interpolation to compute propagation matrix
            
            nStar = length(star);
            H     = cell(nStar, obj.atmModel.nLayer);
            
            fprintf('___ BI-LINEAR INTERPOLATION OPERATOR ___\n')

            for kLayer = 1:obj.atmModel.nLayer;
                
                pitchLayer = obj.dSub/obj.overSampling(kLayer);
                height = obj.atmModel.layer(kLayer).altitude;
                
                for kGs=1:nStar
                    
                    fprintf(' [%d,%d]',kGs,kLayer)
                    
                    % pupil center in layer
                    beta  = star(kGs).directionVector*height;
                    scale = 1-height/star(kGs).height;
                    H{kGs,kLayer} = p_bilinearSplineInterp(...
                        obj.atmGrid{kLayer,1}(:),...
                        obj.atmGrid{kLayer,2}(:),...
                        pitchLayer,...
                        real(outputCoord)*scale+beta(1),...
                        imag(outputCoord)*scale+beta(2));
                end
                
                fprintf('\n');

            end
            
            % local function
            function P = p_bilinearSplineInterp(xo,yo,do,xi,yi)
                
                ni = length(xi);
                
                nxo = length(xo);
                nyo = length(yo);
                no = nxo*nyo;
                
                % remove the interaporated points out of the original grid
                mask = xi>=xo(1) & yi>=yo(1) & xi<=xo(end) & yi<=yo(end);
                
                % index for the inteporated grid
                index = (1:ni)';
                index = index(mask);
                
                % x & y index for the original grid
                ox = floor((xi-xo(1))/do)+1;
                ox = ox(mask);
                oy = floor((yi-yo(1))/do)+1;
                oy = oy(mask);

                % bi-linear inteporation value
                fxo = abs(xi(mask)-(xo(1)+do*(ox-1)))/do;
                fyo = abs(yi(mask)-(yo(1)+do*(oy-1)))/do;
                s1 = (1-fxo).*(1-fyo);
                s2 = fxo.*(1-fyo);
                s3 = (1-fxo).*fyo;
                s4 = fxo.*fyo;
                
                % vectoraized index for the original grid
                o1 = oy+nyo*(ox-1);
                o2 = oy+nyo*ox;
                o3 = oy+1+nyo*(ox-1);
                o4 = oy+1+nyo*ox;
                
                % masking
                o1(s1==0)=[];
                i1 = index(s1~=0);
                s1(s1==0)=[];
                o2(s2==0)=[];
                i2 = index(s2~=0);
                s2(s2==0)=[];
                o3(s3==0)=[];
                i3 = index(s3~=0);
                s3(s3==0)=[];
                o4(s4==0)=[];
                i4 = index(s4~=0);
                s4(s4==0)=[];
                
                % intepolation matrix
                P1 = sparse(i1,o1,s1,ni,no);
                P2 = sparse(i2,o2,s2,ni,no);
                P3 = sparse(i3,o3,s3,ni,no);
                P4 = sparse(i4,o4,s4,ni,no);
                
                P = P1+P2+P3+P4;
            end
        end
        
        function p_L2 = p_sparseInverseCovarianceMatrix(obj)
            %% Computes sparse-approximated phase covariance
            
            % 5x5 stencil
            a = [0 sqrt(2) 1 sqrt(2) 2;...
                sqrt(2) 0 1 2 sqrt(2);...
                1 1 0 1 1;...
                sqrt(2) 2 1 0 sqrt(2);...
                2 sqrt(2) 1 sqrt(2) 0;];

            p_L2 = cell(obj.atmModel.nLayer,1);
            for kLayer = 1:obj.atmModel.nLayer;
                
                m = length(obj.atmGrid{kLayer,2});
                n = length(obj.atmGrid{kLayer,1});
                N = m*n;
                
                e = ones(N,1);
                
                ex1 = e;
                ex1(end-m+1:end) = 2;
                ex1 = circshift(ex1,-m);
                
                ex2 = e;
                ex2(1:m) = 2;
                ex2 = circshift(ex2,m);
               
                ey1 = e;
                ey1(mod((1:N),m)==1) = 0;
                ey1(mod((1:N),m)==0) = 2;
                ey1 = circshift(ey1,-1);
                
                ey2 = e;
                ey2(mod((1:N),m)==1) = 2;
                ey2(mod((1:N),m)==0) = 0;
                ey2 = circshift(ey2,1);
                
                p_L = spdiags([ex1 ey1 -4*e ey2 ex2],[-m -1 0 1 m],N,N);
                p_L2{kLayer} = p_L'*p_L;
                

                ex1 = e;
                ex1(end-m+1:end) = 2;
                ex1(1:m) = 0;
                
                ex2 = e;
                ex2(1:m) = 2;
                ex2(end-m+1:end) = 0;
               
                ey1 = e;
                ey1(mod((1:N),m)==1) = 0;
                ey1(mod((1:N),m)==0) = 2;
                
                ey2 = e;
                ey2(mod((1:N),m)==1) = 2;
                ey2(mod((1:N),m)==0) = 0;

                E = [ex1 ey1 -4*e ey2 ex2];
                pitchLayer = obj.dSub/obj.overSampling(kLayer);
                b = phaseStats.structureFunction(a*pitchLayer,slab(obj.atmModel,kLayer));
                b = b*(obj.atmModel.wavelength/2/pi)^2;
                
                E = mat2cell(E,ones(size(E,1),1),5);
                c=sum(cellfun(@(x) x*b*x',E));
                
                p_L2{kLayer} = N/(-1/2*c)*p_L2{kLayer};
            end
            
        end
        
        
    end
    
    methods (Static)
        
        function demo_OL_LTAO()
            %% Demo scripts for open-loop LTAO
            
            % Telescope
            nL                  = 16;               % number of lenslets
            nPx                = 10;               % number of pixels per lenslet (backwards compatibility with SH WFS)
            nRes              = nL*nPx;           % resolution on the pupil plane (no of pixels)
            D                   = 8;               % this needs be done to ensure influence functions are at the corners of the sub-apertures
            d                  = D/nL;             % lenslet pitch
            samplingFreq       = 100;              % WFS sampling time
            obstructionRatio   = 0.3;              % central obscuration ratio
            fieldOfViewInArcsec= 90;               % fieldOfViewInArcsec
            tel = telescope(D,'resolution',nRes,...
                'obstructionRatio',obstructionRatio,...
                'fieldOfViewInArcsec',fieldOfViewInArcsec,...
                'samplingTime',1/samplingFreq);
            
            % Calibration Source
            ngs = source('wavelength',photometry.R);
            
            % Shack Hartmann WFS
            wfs = shackHartmann(nL,2*nPx*nL,0.6);
            ngs = ngs.*tel*wfs;
            wfs.INIT
            +wfs;
            wfs.camera.frameListener.Enabled = false;
            wfs.slopesListener.Enabled = false;
            wfs.gainCalibration(tel,ngs);
    
            % DM
            dmCrossCouplingCoeff = 0.4;
            bifa = gaussianInfluenceFunction(dmCrossCouplingCoeff,0.5);
            dm = deformableMirror(nL+1,'modes',bifa,...
                'resolution',tel.resolution,...
                'validActuator',wfs.validActuator);
            
            % Low-resolution DM
            bifaLowRes = gaussianInfluenceFunction(dmCrossCouplingCoeff,0.5);
            dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,...
                'resolution',nL+1,...
                'validActuator',wfs.validActuator);
            F = 2*bifaLowRes.modes(wfs.validActuator,:);
            iF = pinv(full(F));
    
            % DM calibration
            ngs = ngs.*tel;
            calibDm = calibration(dm,wfs,ngs,ngs.wavelength/100,nL+1,'cond',1e2);
            
            % Guide stars
            ast = source('asterism',{[3,arcsec(30),0]},...
                'wavelength', photometry.R);

            % On-axis science source
            science = source('zenith',zeros(1,1),...
                'azimuth',zeros(1,1),...
                'wavelength',photometry.H);

            % Atmosphere
            r0 = 0.15;
            L0 = 30;
            fractionalR0    = [0.5,0.3,0.2];
            altitude        = [0e3,5e3,12e3];
            windSpeed       = [10,5,20];
            windDirection   = [0,pi/2,pi];
            atm = atmosphere(photometry.V0,r0,L0,...
                'fractionnalR0',fractionalR0,'altitude',altitude,...
                'windSpeed',windSpeed,'windDirection',windDirection);
            tel = tel + atm;
            
            % Tomoraphic Reconstructor
            mmse =sparseLMMSE(wfs,tel,atm,ast,'mmseStar',science,...
                'iNoiseVar', 1/1e-14);
            %R = mmse.getReconstructionMatrix; % MVM case

%             % MMSE with explicit noise covariance matrix
%             iNoiseCovMat = ...
%                 wfs.theoreticalNoise(tel, atm, ast, ngs,'naParam',repmat([10000,90000],3,1),'lgsLaunchCoord',zeros(3,2),...
%                 'computeInverseCovMat',true);
%             mmse =sparseLMMSE(wfs,tel,atm,ast,'mmseStar',science,...
%                 'iNoiseVar', blkdiag(iNoiseCovMat{:}));
% 
%             
%             % MMSE with rotated slopes and explicit noise covariance matrix
%             iNoiseCovMat = ...
%                 wfs.theoreticalNoise(tel, atm, ast, ngs,'naParam',repmat([10000,90000],3,1),'lgsLaunchCoord',zeros(3,2),...
%                 'rotatedFrame',true,'computeInverseCovMat',true);
%                 
%             mmse =sparseLMMSE(wfs,tel,atm,ast,'mmseStar',science,...
%                 'iNoiseVar', blkdiag(iNoiseCovMat{:}),'rotatedFrame',true);
            
            % Science camera
            cam = imager;
            cam.frameListener.Enabled = false;
            tel = tel - atm;
            science = science.*tel*cam;
            cam.referenceFrame = cam.frame;
            tel = tel + atm;
            
            % Loop initialization
            reset(tel);
            ast = ast.*tel*wfs;
            science = science.*tel*dm*cam;
            flush(cam)
            cam.frame = cam.frame*0;
            cam.clockRate    = 1;
            cam.exposureTime = 50;
            cam.startDelay   = 10;
            nIt = cam.startDelay+cam.exposureTime;
            dm.coefs = 0;

            % Display initialization
            figure(4410);
            subplot(2,2,1);
            h1 = imagesc(zeros(nRes),[-1 1]*4e-6);
            axis xy equal tight
            title('Input WF');
            h1.Parent.XTickLabel = '';
            h1.Parent.YTickLabel = '';

            subplot(2,2,2);
            h2 = imagesc(zeros(nL+1),[-1 1]*4e-6);
            axis xy equal tight
            title('Reconstructed WF');
            h2.Parent.XTickLabel = '';
            h2.Parent.YTickLabel = '';

            subplot(2,2,3);
            h3 = imagesc(zeros(nRes),[-1 1]*4e-6);
            axis xy equal tight
            title('Residual WF');
            h3.Parent.XTickLabel = '';
            h3.Parent.YTickLabel = '';

            subplot(2,2,4);
            h4 = imagesc(zeros(size(cam.referenceFrame)));
            title('PSF');
            axis equal tight xy
            h4.Parent.XTickLabel = '';
            h4.Parent.YTickLabel = '';
            
            % LTAO loop (open loop)
            rec = zeros(nL+1);
            for k=1:nIt
                tic
                fprintf('Frame %d/%d\n',k,nIt);
                
                +tel;
                
                science=science.*tel;
                set(h1,'CData',science.catMeanRmPhase/science.waveNumber);
                
                +ast;
                rec(wfs.validActuator) = mmse*wfs; % iterative case
                %rec(wfs.validActuator) = R*wfs.slopes(:); % MVM case
                
                dm.coefs= iF*rec(wfs.validActuator);
                set(h2,'CData',rec);
                
                science=science*dm*cam;
                set(h3,'CData',science.catMeanRmPhase/science.waveNumber);
                
                if cam.frameCount == 0 || isscalar(cam.frameBuffer)
                    set(h4,'CData',cam.frame);
                else
                    set(h4,'CData',cam.frameBuffer);
                end
                
                drawnow;
            end
            
            fprintf('\n\nStrehl Ratio = %.2f\n\n',cam.strehl*100);
            % Strehl Ratio = 44.71
        end
        
        
        function demo_LayerEstimation()
            %% Demo scripts for tomographic estimation of layered phase 
            
            % Telescope
            nL                  = 16;               % number of lenslets
            nPx                = 10;               % number of pixels per lenslet (backwards compatibility with SH WFS)
            nRes              = nL*nPx;           % resolution on the pupil plane (no of pixels)
            D                   = 8;               % this needs be done to ensure influence functions are at the corners of the sub-apertures
            d                  = D/nL;             % lenslet pitch
            samplingFreq       = 100;              % WFS sampling time
            obstructionRatio   = 0.3;              % central obscuration ratio
            fieldOfViewInArcsec= 90;               % fieldOfViewInArcsec
            tel = telescope(D,'resolution',nRes,...
                'obstructionRatio',obstructionRatio,...
                'fieldOfViewInArcsec',fieldOfViewInArcsec,...
                'samplingTime',1/samplingFreq);
            
            % Calibration Source
            ngs = source('wavelength',photometry.R);
            
            % Shack Hartmann WFS
            wfs = shackHartmann(nL,2*nPx*nL,0.6);
            ngs = ngs.*tel*wfs;
            wfs.INIT
            +wfs;
            wfs.camera.frameListener.Enabled = false;
            wfs.slopesListener.Enabled = false;
            wfs.gainCalibration(tel,ngs);
    
            % Guide stars
            ast = source('asterism',{[3,arcsec(30),0]},...
                'wavelength', photometry.R);

            % On-axis science source
            science = source('zenith',zeros(1,1),...
                'azimuth',zeros(1,1),...
                'wavelength',photometry.H);

            % Atmosphere
            r0 = 0.15;
            L0 = 30;
            fractionalR0    = [0.5,0.3,0.2];
            altitude        = [0e3,5e3,12e3];
            windSpeed       = [10,5,20];
            windDirection   = [0,pi/2,pi];
            atm = atmosphere(photometry.V0,r0,L0,...
                'fractionnalR0',fractionalR0,'altitude',altitude,...
                'windSpeed',windSpeed,'windDirection',windDirection);
            tel = tel + atm;
            
            % Tompraphic Reconstructor
            mmse = sparseLMMSE(wfs,tel,atm,ast,'mmseStar',science,...
                'iNoiseVar', 1/1e-14,...
                'layerEstimation',true);
            %R = mmse.getReconstructionMatrix; % MVM case
            
            % Loop initialization
            reset(tel);
            ast = ast.*tel*wfs;
            nIt = 60;

            % Compute reconstructed area at each altitude
            ex1 = [mmse.atmGrid{1,1}(1), mmse.atmGrid{1,1}(end)]; % xlim for reconstructed phase
            ey1 = [mmse.atmGrid{1,2}(1), mmse.atmGrid{1,2}(end)]; % ylim for reconstructed phase
            ex2 = [mmse.atmGrid{2,1}(1), mmse.atmGrid{2,1}(end)]; % xlim for reconstructed phase
            ey2 = [mmse.atmGrid{2,2}(1), mmse.atmGrid{2,2}(end)]; % ylim for reconstructed phase
            ex3 = [mmse.atmGrid{3,1}(1), mmse.atmGrid{3,1}(end)]; % xlim for reconstructed phase
            ey3 = [mmse.atmGrid{3,2}(1), mmse.atmGrid{3,2}(end)]; % ylim for reconstructed phase
            eDx1 = ex1(2)-ex1(1);
            eDy1 = ey1(2)-ey1(1);
            eDx2 = ex2(2)-ex2(1);
            eDy2 = ey2(2)-ey2(1);
            eDx3 = ex3(2)-ex3(1);
            eDy3 = ey3(2)-ey3(1);
            
            aD1 = (size(atm.layer(1).phase,1)-1)*D/(nRes-1); % diameter of meta-pupil
            aD2 = (size(atm.layer(2).phase,1)-1)*D/(nRes-1); % diameter of meta-pupil
            aD3 = (size(atm.layer(3).phase,1)-1)*D/(nRes-1); % diameter of meta-pupil
            [aGx1,aGy1] = meshgrid(linspace(-1,1,size(atm.layer(1).phase,1))*aD1/2); % grid
            [aGx2,aGy2] = meshgrid(linspace(-1,1,size(atm.layer(2).phase,1))*aD2/2); % grid
            [aGx3,aGy3] = meshgrid(linspace(-1,1,size(atm.layer(3).phase,1))*aD3/2); % grid

             % Display initialization
            figure(4410);
            subplot(2,3,1);
            h1 = imagesc(aGx1(:),aGy1(:),atm.layer(1).phase);
            axis xy equal tight
            title('Input 1');
            h1.Parent.XTickLabel = '';
            h1.Parent.YTickLabel = '';
            xlim(ex1);
            ylim(ey1);
            pbaspect([eDx1/eDx1 eDy1/eDy1 1]);

            subplot(2,3,2);
            h2 = imagesc(aGx2(:),aGy2(:),atm.layer(2).phase);
            axis xy equal tight
            title('Input 2');
            h2.Parent.XTickLabel = '';
            h2.Parent.YTickLabel = '';
            xlim(ex2);
            ylim(ey2);
            pbaspect([eDx2/eDx1 eDy2/eDy1 1]);
            
            subplot(2,3,3);
            h3 = imagesc(aGx3(:),aGy3(:),atm.layer(3).phase);
            axis xy equal tight
            title('Input 3');
            h3.Parent.XTickLabel = '';
            h3.Parent.YTickLabel = '';
            xlim(ex3);
            ylim(ey3);
            pbaspect([eDx3/eDx1 eDy3/eDy1 1]);
            
            subplot(2,3,4);
            nx1 = length(mmse.atmGrid{1,1});
            ny1 = length(mmse.atmGrid{1,2});
            h4 = imagesc(zeros(nx1,ny1));
            axis xy equal tight
            title('Reconstructed  1');
            h4.Parent.XTickLabel = '';
            h4.Parent.YTickLabel = '';
            pbaspect([eDx1/eDx1 eDy1/eDy1 1]);
            
            subplot(2,3,5);
            nx2 = length(mmse.atmGrid{2,1});
            ny2 = length(mmse.atmGrid{2,2});
            h5 = imagesc(zeros(nx2,ny2));
            axis xy equal tight
            title('Reconstructed 2');
            h5.Parent.XTickLabel = '';
            h5.Parent.YTickLabel = '';
            pbaspect([eDx2/eDx1 eDy2/eDy1 1]);

            subplot(2,3,6);
            nx3 = length(mmse.atmGrid{3,1});
            ny3 = length(mmse.atmGrid{3,2});
            h6 = imagesc(zeros(nx3,ny3));
            axis xy equal tight
            title('Reconstructed 3');
            h6.Parent.XTickLabel = '';
            h6.Parent.YTickLabel = '';
            pbaspect([eDx3/eDx1 eDy3/eDy1 1]);
            
            % LTAO loop (open loop)
            est1 = zeros(nx1,ny1);
            est2 = zeros(nx2,ny2);
            est3 = zeros(nx3,ny3);
            for k=1:nIt
                tic
                fprintf('Frame %d/%d\n',k,nIt);
                
                +tel;
                set(h1,'CData',atm.layer(1).phase);
                set(h2,'CData',atm.layer(2).phase);
                set(h3,'CData',atm.layer(3).phase);
                
                +ast;
                est = mmse*wfs; % iterative case
                %est = R*wfs.slopes(:); % MVM case
                
                est1 = reshape(est(1:nx1*ny1),nx1,ny1);
                set(h4,'CData',est1);
                est2 = reshape(est(nx1*ny1+(1:nx2*ny2)),ny2,nx2);
                set(h5,'CData',est2);
                est3 = reshape(est(nx1*ny1+nx2*ny2+(1:nx3*ny3)),ny3,nx3);
                set(h6,'CData',est3);

                drawnow;
            end
            
        end
    end
end