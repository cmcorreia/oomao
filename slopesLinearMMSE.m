classdef slopesLinearMMSE < handle
    %% SLOPESLINEARMMSE MMSE wavefront estimation from SH-WFS slopes
    %
    % slopesLinearMMSE(wfs,tel,atmModel,guideStar,'mmseStar',science)
    % computes the tomographic reconstructor to estimate the directional
    % phase in the science direction(s) from SH-WFSs meurements based on
    % turbulence model given by atmModel.
    %
    % The tomography process is done iteratively using minres function in
    % MATLAB and you can change parameters for the itelative computation
    % by option
    %   'MAXIT' : The maximum number of Iteration
    %   'RTOL'  : The tolerance
    % The defaults are MAXIT=20 and RTOL=1e-2.
    %
    % Also, you can get the tomographic reconstructor explicitly by a
    % function of 'getReconstructionMatrix' as
    %    mmse = slopesLinearMMSE(wfs,tel,atmModel,guideStar,'mmseStar',science, other options...);
    %    R = mmse.getReconstructionMatrix();
    %    phi = R*s; % reconstruction 
    % where the input s should be given by in [pixel] and the output phi is
    % in [nm]. Don't use this function for large-scale AO systems!!!!
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
    % You can choose slopes Model to derive theoretical covariance matrix.
    %   'covModel' : FFT or Hudgin or Fried, default id FFT model
    %
    % Demoscripts for test
    %     slopeLinearMMSE.demo_OL_LTAO; % Open-loop LTAO script
    
    properties
        tag = 'slopesLinearMMSE';
        % diameter of sub-aperture [m]
        dSub
        % number of sub-aperture along diameter
        nSub
        % diameter of telescope
        D
        % slopes pupil mask
        slopesMask
        % output wavefront pupil mask
        wavefrontMask
        % wavefront array size
        wavefrontSize;
        % atmosphere mode
        atmModel
        % guide star covariance matrix
        Cxx
        % mmse star covariance matrix
        Cox
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
        % bilinear sparse operator
        B = [];
        % covariance Model
        covModel = 'FFT';
        % flag for TT removal
        isTTRM
        % flag for Focus removal
        isFocusRM
        % Time lag for linear prediction control
        tPredict
        % Time difference for multi-time-step reconstruction
        tMulti
        % converson factor
        wavefrontToMeter;
        % Parameter for spot elongation (by Leo)
        binningFactor = 1;
        convKernel = false;
        % solver
        solver
        NP
    end
    
    properties (Dependent)
        % mmse estimates directions
        mmseStar;
        % noise variance
        noiseVar;
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
        p_noiseVar = 0;
        fun;
        funp;
        c;
        %wavefrontToMeter;
        ttrmfunSlope; % function for the TT or/and focus removal
        ttrmfunPhase; % function for the TT or/and focus removal
        
        p_NF; % Fourier spectrum resolution
        p_sf; % Fourier spectrum oversampling factor
    end
    
    methods
        
        %% Constructor
        function obj = slopesLinearMMSE(wfs,tel,atmModel,guideStar,varargin)
            
            inputs = inputParser;
            inputs.addRequired('wfs', @(x) isa(x,'shackHartmann') ); % shack hartmann object
            inputs.addRequired('tel', @(x) isa(x,'telescope')); % telescope object
            inputs.addRequired('atmModel', @(x) isa(x,'atmosphere')); % atmosphere object
            inputs.addRequired('guideStar', @(x) isa(x,'source')); % source object for guide stars
            inputs.addOptional('mmseStar', [],@(x) isa(x,'source')); % source object for science stars
            inputs.addOptional('solver','minres',@ischar);
            inputs.addOptional('slopesMask', [], @islogical ); % mask for WFSs
            inputs.addOptional('noiseVar', 0, @isnumeric ); % Noise covariance matrix in unit of [radian^2]
            inputs.addOptional('covModel', 'FFT', @ischar ); % Slope Model for theoretical covariance
            inputs.addOptional('isTTRM', false, @islogical ); % flag for removing tip/tilt
            inputs.addOptional('isFocusRM', false, @islogical ); % flag for removing focus
            inputs.addOptional('RTOL', 5e-2, @isnumeric ); % Tolerance for iterative method (minres)
            inputs.addOptional('MAXIT', 200, @isnumeric ); % Muximum number of iteration  (minres)
            inputs.addOptional('isWarmStart', true, @islogical ); % flag for warm start
            inputs.addOptional('tPredict', 0, @isnumeric ); % Time lag for linear prediction control
            inputs.addOptional('tMulti', 0, @isnumeric ); % Time difference for mulit-time-step reconstruction
            inputs.addOptional('NF', 1024, @isnumeric ); % Tuning parameter for FFT model
            inputs.addOptional('sf', 4, @isnumeric ); % Tuning parameter for FFT model
            inputs.addOptional('wavefrontMask', [], @islogical ); % phse mask for reconstructed phase
            
            inputs.parse(wfs,tel,atmModel,guideStar,varargin{:});
            
            obj.nSub    = inputs.Results.wfs.lenslets.nLenslet;
            obj.D      = inputs.Results.tel.D;
            obj.dSub      = obj.D/obj.nSub;
            if isempty(inputs.Results.wavefrontMask)
                obj.wavefrontMask = inputs.Results.wfs.validActuator;
            else
                obj.wavefrontMask = inputs.Results.wavefrontMask;
            end
            
            obj.atmModel      = inputs.Results.atmModel;
            obj.guideStar     = inputs.Results.guideStar;
            obj.covModel      = inputs.Results.covModel;
            obj.isTTRM        = inputs.Results.isTTRM;
            obj.isFocusRM     = inputs.Results.isFocusRM;
            obj.tPredict      = inputs.Results.tPredict;
            obj.tMulti        = inputs.Results.tMulti(:);
            obj.RTOL          = inputs.Results.RTOL;
            obj.MAXIT         = inputs.Results.MAXIT;
            obj.warmStart     = inputs.Results.isWarmStart;
            obj.solver        = inputs.Results.solver;
            obj.binningFactor = inputs.Results.wfs.camera.binningFactor;
            obj.convKernel    = inputs.Results.wfs.lenslets.convKernel;
            
            if obj.tMulti(1) ~= 0
                obj.tMulti(1) = [0; obj.tMulti];
            end
            if ~isempty(inputs.Results.slopesMask)
                obj.slopesMask = inputs.Results.slopesMask;
            else
                obj.slopesMask = repmat( repmat( inputs.Results.wfs.validLenslet(:) , 2, 1), 1, obj.nGuideStar );
            end
            obj.noiseVar      = inputs.Results.noiseVar;
            obj.p_mmseStar    = inputs.Results.mmseStar;
            if isempty(obj.p_mmseStar)
                obj.p_mmseStar = obj.guideStar;
            end
            obj.p_NF            = inputs.Results.NF;
            obj.p_sf            = inputs.Results.sf;
            obj.createTTRMfunc;
            
            obj.log = logBook.checkIn(obj);
            obj.log.verbose = false;
            %add(obj.log,obj,'atmosphere wavelength set to mmse star wavelength!')
            %atmWavelength = obj.atmModel.wavelength;
            %obj.atmModel.wavelength = obj.p_mmseStar(1).wavelength;
            
            
            %------ Computing the covariance matrix ------%
            if isinf(obj.guideStar(1).height)
                add(obj.log,obj,sprintf('Computing the covariance matrices for NGS with %s Model',obj.covModel))
            else
                add(obj.log,obj,sprintf('Computing the covariance matrices for LGS with %s Model',obj.covModel))
            end
            
            % slopes-slopes covariance
            dv = [obj.guideStar.directionVector];
            switch (obj.covModel)
                case {'HUDGIN','hudgin','Hudgin'}
                    scovfun = @(dx,dy,lag) slopestoSlopesCovarianceMapHudgin(obj, dx, dy, lag);
                case {'FRIED','fried','Fried'}
                    scovfun = @(dx,dy,lag) slopestoSlopesCovarianceMapFried(obj, dx, dy, lag);
                otherwise % FFT
                    scovfun = @(dx,dy,lag) slopestoSlopesCovarianceMapFFT(obj, dx, dy, lag);
            end
            nm = [obj.nSub obj.nSub];
            m_Cxx = cell(obj.nGuideStar*length(obj.tMulti));
            cmap = scovfun(0, 0, 0);
            aCxx{1,1} = toeplitzBlockToeplitz( nm, nm, cmap{1,1} );
            aCxx{2,2} = toeplitzBlockToeplitz( nm, nm, cmap{2,2} );
            aCxx{1,2} = toeplitzBlockToeplitz( nm, nm, cmap{1,2} );
            aCxx{2,1} = toeplitzBlockToeplitz( nm, nm, cmap{2,1} );

            nss = sum(1:obj.nGuideStar-1);
            t = 0;
            fprintf('Computing slope-slope covariance : 000%%...');
            for iGs=1:obj.nGuideStar
                m_Cxx{iGs,iGs} = aCxx;
                for jGs=iGs+1:obj.nGuideStar
                    
                    deltaSrc_x = dv(1,jGs)-dv(1,iGs);
                    deltaSrc_y = dv(2,jGs)-dv(2,iGs);
                    
                    cmap = scovfun(deltaSrc_x, deltaSrc_y, 0);
                    
                    m_Cxx{iGs,jGs}{1,1} = toeplitzBlockToeplitz( nm, nm, cmap{1,1} );
                    m_Cxx{iGs,jGs}{2,2} = toeplitzBlockToeplitz( nm, nm, cmap{2,2} );
                    m_Cxx{iGs,jGs}{1,2} = toeplitzBlockToeplitz( nm, nm, cmap{1,2} );
                    m_Cxx{iGs,jGs}{2,1} = toeplitzBlockToeplitz( nm, nm, cmap{2,1} );
                    
                    m_Cxx{jGs,iGs}{1,1} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{1,1},2) );
                    m_Cxx{jGs,iGs}{2,2} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{2,2},2) );
                    m_Cxx{jGs,iGs}{1,2} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{1,2},2) );
                    m_Cxx{jGs,iGs}{2,1} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{2,1},2) );
                    
                    t=t+1;
                    fprintf('\b\b\b\b\b\b\b%03.0f%%...',t/nss*100);
                end
            end
            fprintf('done\n');
            
            % for multi-time-step
            for iT=1:length(obj.tMulti)
                if iT>1
                    for iGs=1:obj.nGuideStar
                        for jGs=1:obj.nGuideStar
                            iN = iGs + obj.nGuideStar*(iT-1);
                            jN = jGs + obj.nGuideStar*(iT-1);
                            m_Cxx{iN,jN} = m_Cxx{iGs,jGs};
                        end
                    end
                end
                for jT=iT+1:length(obj.tMulti)
                    for iGs=1:obj.nGuideStar
                        for jGs=1:obj.nGuideStar
                            
                            iN = iGs + obj.nGuideStar*(iT-1);
                            jN = jGs + obj.nGuideStar*(jT-1);

                            deltaSrc_x = dv(1,jGs)-dv(1,iGs);
                            deltaSrc_y = dv(2,jGs)-dv(2,iGs);
                            
                            cmap = scovfun(deltaSrc_x, deltaSrc_y, -obj.tMulti(iT)+obj.tMulti(jT));
                            
                            m_Cxx{iN,jN}{1,1} = toeplitzBlockToeplitz( nm, nm, cmap{1,1} );
                            m_Cxx{iN,jN}{2,2} = toeplitzBlockToeplitz( nm, nm, cmap{2,2} );
                            m_Cxx{iN,jN}{1,2} = toeplitzBlockToeplitz( nm, nm, cmap{1,2} );
                            m_Cxx{iN,jN}{2,1} = toeplitzBlockToeplitz( nm, nm, cmap{2,1} );
                            
                            m_Cxx{jN,iN}{1,1} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{1,1},2) );
                            m_Cxx{jN,iN}{2,2} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{2,2},2) );
                            m_Cxx{jN,iN}{1,2} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{1,2},2) );
                            m_Cxx{jN,iN}{2,1} = toeplitzBlockToeplitz( nm, nm, rot90(cmap{2,1},2) );
                        end
                    end
                end
            end
            obj.Cxx = m_Cxx;
            
            % phase-slopes covariance
            fprintf('Computing phase-slope covariance : 000%%...');
            if isinf(obj.guideStar(1).height)
                switch (obj.covModel)
                    case {'HUDGIN','hudgin','Hudgin'}
                        pcovfun = @(dx,dy,lag) phaseToSlopesCovarianceMapHudginAllLayer(obj,dx,dy, lag);
                    case {'FRIED','fried','Fried'}
                        pcovfun = @(dx,dy,lag) phaseToSlopesCovarianceMapFriedAllLayer(obj,dx,dy, lag);
                    otherwise % FFT
                        pcovfun = @(dx,dy,lag) phaseToSlopesCovarianceMapFFTAllLayer(obj,dx,dy, lag);
                end
            else
                switch (obj.covModel)
                    case {'HUDGIN','hudgin','Hudgin'}
                        pcovfun = @(dx,dy,pad,atmLayer,gl,lag) phaseToSlopesCovarianceMapHudginSingleLayer(obj,pad,atmLayer,dx,dy,gl, lag);
                    case {'FRIED','fried','Fried'}
                        pcovfun = @(dx,dy,pad,atmLayer,gl,lag) phaseToSlopesCovarianceMapFriedSingleLayer(obj,pad,atmLayer,dx,dy,gl, lag);
                    otherwise % FFT
                        pcovfun = @(dx,dy,pad,atmLayer,gl,lag) phaseToSlopesCovarianceMapFFTSingleLayer(obj,pad,atmLayer,dx,dy,gl, lag);
                end
            end
            if isinf(obj.guideStar(1).height) % for NGS
                nps = obj.nMmseStar*length(obj.tMulti)*obj.nGuideStar;
                t=0;
                
                dvm = [obj.mmseStar.directionVector];
                m_Cox = cell(obj.nMmseStar,1);
                nm = [obj.nSub+1 obj.nSub];
                for iSci = 1:obj.nMmseStar
                    m_Cox{iSci} = cell(1,1);
                    m_Cox{iSci}{1} = cell(1,obj.nGuideStar*length(obj.tMulti));
                    for iT=1:length(obj.tMulti)
                        for iGs = 1:obj.nGuideStar
                            iN = iGs + obj.nGuideStar*(iT-1);
                            deltaSrc_x = dv(1,iGs)-dvm(1,iSci);
                            deltaSrc_y = dv(2,iGs)-dvm(2,iSci);
                            cmap = pcovfun(deltaSrc_x,deltaSrc_y,obj.tPredict+obj.tMulti(iT));
                            m_Cox{iSci}{1}{iN}{1} = toeplitzBlockToeplitz( nm, nm, cmap{1} );
                            m_Cox{iSci}{1}{iN}{2} = toeplitzBlockToeplitz( nm, nm, cmap{2} );
                            
                            t=t+1;
                            fprintf('\b\b\b\b\b\b\b%03.0f%%...',t/nps*100);
                            
                        end
                    end
                end
                obj.Cox = m_Cox;
                
            else  % for LGS
                nps = obj.nMmseStar*length(obj.tMulti)*obj.nGuideStar*obj.atmModel.nLayer;
                t=0;
                
                height = obj.guideStar(1).height;
                dvm = [obj.mmseStar.directionVector];
                m_Cox = cell(obj.nMmseStar,1);
                for iSci = 1:obj.nMmseStar
                    m_Cox{iSci} = cell(obj.atmModel.nLayer,1);
                    for kLayer = 1:obj.atmModel.nLayer
                        m_Cox{iSci}{kLayer} = cell(1,obj.nGuideStar*length(obj.tMulti));
                        atmLayer = slab(obj.atmModel,kLayer);
                        g = 1 - atmLayer.layer.altitude/height;
                        pad = ceil( 0.5*obj.nSub*(1-g)/g );
                        nm = [obj.nSub+1+2*pad obj.nSub];
                        for iT=1:length(obj.tMulti)
                            for iGs = 1:obj.nGuideStar
                                iN = iGs + obj.nGuideStar*(iT-1);
                                deltaSrc_x = dv(1,iGs)-dvm(1,iSci);
                                deltaSrc_y = dv(2,iGs)-dvm(2,iSci);
                                cmap = pcovfun(deltaSrc_x,deltaSrc_y,pad,atmLayer,g,obj.tPredict+obj.tMulti(iT));
                                m_Cox{iSci}{kLayer}{iN}{1} = toeplitzBlockToeplitz( nm, nm, cmap{1} );
                                m_Cox{iSci}{kLayer}{iN}{2} = toeplitzBlockToeplitz( nm, nm, cmap{2} );
                                
                                t=t+1;
                                fprintf('\b\b\b\b\b\b\b%03.0f%%...',t/nps*100);
                                
                            end
                        end
                        if iSci == 1
                            NI = obj.nSub+1;
                            NP = NI+2*pad;
                            obj.NP(kLayer) = NP;
                            obj.B{kLayer} = tools.bilinearSparseInterpolator(NI,NP,1,g);
                        end
                    end
                end
                obj.Cox = m_Cox;
                obj.B = cell2mat( obj.B );
                
            end
            fprintf('done\n');

            obj.fun = @(x) mtimes4squareBlocksVaryingSlopesMaskYoshito(obj.Cxx,x,repmat(obj.slopesMask(:),length(obj.tMulti),1),obj.p_noiseVar,obj.ttrmfunSlope);
            %{
            if isempty(inputs.Results.slopesMask)
                obj.fun = @(x) mtimes4squareBlocks(obj.Cxx,x,wfs.validLenslet(:));
            else
                obj.fun = @(x) mtimes4squareBlocksVaryingSlopesMask(obj.Cxx,x,obj.slopesMask);
            end
            %}
                            
            obj.funp = @(x) mtimes4precond(obj.Cxx,x,obj.slopesMask);
            %obj.funp = @(x) mtimes4precond(obj.Cxx,x,wfs.validLenslet(:));
            obj.c   = zeros(obj.nSub^2*2*obj.nGuideStar*length(obj.tMulti),1 );
            if ~obj.convKernel
                obj.wavefrontToMeter = ...
                    obj.guideStar(1).wavelength*obj.atmModel.wavelength/2/pi/(obj.dSub*wfs.lenslets.fftPad)/wfs.lenslets.nyquistSampling;
            else
                obj.wavefrontToMeter = ...
                    obj.guideStar(1).wavelength*obj.atmModel.wavelength/2/pi/(obj.dSub*wfs.lenslets.fftPad)/wfs.lenslets.nyquistSampling*obj.binningFactor;
                fprintf('\n Spots created with convolution kernels \n')
            end
            %obj.wavefrontToMeter = wfs.lenslets.fieldStopSize*obj.guideStar(1).wavelength/obj.dSub...
            %    /(wfs.camera.resolution(1)/obj.nSub)*obj.atmModel.wavelength/(2*pi);
            obj.wavefrontSize    = sum(obj.wavefrontMask(:));
            obj.x0               = zeros(size(obj.c));
            obj.log.verbose = true;
            %obj.atmModel.wavelength = atmWavelength;
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
        % size of noise covairance should be Ns where Ns is the total
        % number of sub-aperture. If the input is masked noise covariance
        % matrix (or vector), this function reassigns the input to the
        % non-masked covariance matrix (or vector).
        function set.noiseVar(obj,val)
            ns = length(obj.slopesMask(:));
            if isscalar(val) % for scalar
                obj.p_noiseVar = val;
            elseif isvector(val) % for vector
                if length(val) == ns
                    obj.p_noiseVar = val(:);
                else % reassigning
                    obj.p_noiseVar = zeros(ns,1);
                    obj.p_noiseVar(obj.slopesMask(:)==1) = val(:);
                end
            else % for matrix
                if size(val,1) == ns
                    obj.p_noiseVar = val;
                else  % reassigning
                    nsEach = sum(obj.slopesMask,1)/2;
                    nsEach = [nsEach; nsEach];
                    nsEach = nsEach(:);
                    val=mat2cell(val,nsEach,nsEach);
                    obj.p_noiseVar = cell(obj.nGuideStar,1);
                    for k=1:obj.nGuideStar
                        tmp = cell(2,2);
                        ind = 1:2*obj.nSub^2;
                        ind(~obj.slopesMask(:,k)) = [];
                        % x-x
                        [c,r,v]=find(val{2*k-1,2*k-1});
                        tmp{1,1} = sparse(ind(c),ind(r),v,obj.nSub^2,obj.nSub^2);
                        % x-y
                        [c,r,v]=find(val{2*k-1,2*k});
                        tmp{1,2} = sparse(ind(c),ind(r),v,obj.nSub^2,obj.nSub^2);
                        % y-x
                        [c,r,v]=find(val{2*k,2*k-1});
                        tmp{2,1} = sparse(ind(c),ind(r),v,obj.nSub^2,obj.nSub^2);
                        % y-y
                        [c,r,v]=find(val{2*k,2*k});
                        tmp{2,2} = sparse(ind(c),ind(r),v,obj.nSub^2,obj.nSub^2);
                        obj.p_noiseVar{k} = cell2mat(tmp);
                    end
                    obj.p_noiseVar = blkdiag(obj.p_noiseVar{:});
                end
            end
            %slopestoSlopesCovariance(obj);
        end
        function val = get.noiseVar(obj)
            val = obj.p_noiseVar;
        end
        
        %% Set nGuideStar
        function val = get.nGuideStar(obj)
            val = length(obj.guideStar(1,:,1));
        end
        
        %% Set nMmseStar
        function val = get.nMmseStar(obj)
            val = length(obj.mmseStar(1,:,1));
        end
        
        %{
        %% Wavefront reconstruction (Old version)
        function out = mtimes(obj,wfs)
            if isa(wfs,'shackHartmann')
                obj.c( obj.slopesMask ) = wfs.slopes(:);
            else
                obj.c( obj.slopesMask ) = wfs;
            end
            [yy] = minres(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart);
            %[yy,~] = minres(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,obj.fun,[],obj.x0*obj.warmStart);
            obj.x0 = yy;
            m_Cox = obj.Cox{1};
            n = length(m_Cox);
            m = length(m_Cox{1});
            out = cell( n , 1 );
            ns = obj.nSub^2*2;
            %             nFig = 100 + randi(100);
            for kn=1:n
                out{kn} = zeros( m_Cox{kn}{1}{1}.nRow^2 , 1 );
                u  = 1:ns;
                for km=1:m
                    yyy = yy(u);
                    out{kn} = out{kn} + ...
                        m_Cox{kn}{km}{1}*yyy(1:end/2) + m_Cox{kn}{km}{2}*yyy(1+end/2:end);
                    u = u + ns;
                end
                
                %                 figure(nFig)
                %                 subplot(1,n,kn)
                %                 imagesc( reshape( out{kn} , ones(1,2)*m_Cox{kn}{1}{1}.nRow ) )
                %                 axis square
                %                 colorbar
            end
            out = cell2mat( out );
            if ~isempty(obj.B)
                out = obj.B*out;
            end

            if length(out)>prod(obj.wavefrontSize)
                out(obj.wavefrontMask) = [];
            else
                if numel(out)==numel(obj.wavefrontMask)
                    out(obj.wavefrontMask) = 0;
                end
                out = reshape(out,obj.wavefrontSize);
            end
            out = out*obj.wavefrontToMeter;
            
        end
        %}
        
        %% Wavefront reconstruction (new version by Y.O 18/3/2017)
        function out = mtimes(obj,wfs)
            if isa(wfs,'shackHartmann')
                obj.c( obj.slopesMask ) = wfs.slopes(:);
            else
                obj.c( repmat(obj.slopesMask(:),length(obj.tMulti),1) ) = wfs(:);
            end
            
            % remove tip/tilt/focus from slopes if needed
            v=0;
            ns = 2*obj.nSub^2;
            for ii=1:(obj.nGuideStar*length(obj.tMulti))
                v = (1:ns)+v(end);
                obj.c(v) = obj.ttrmfunSlope{ii}(obj.c(v));
            end
            
            % iterative computation
            switch obj.solver
                case 'cg'
                    [yy] = cgs(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart);
                case 'minres'
                    [yy] = minres(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart);
                case 'bicg'
                    [yy] = bicgstab(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart);
            end
            
            %[yy,~] = minres(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,obj.fun,[],obj.x0*obj.warmStart);
            obj.x0 = yy;
            
            nSci = length(obj.Cox);
            nLayer = length(obj.Cox{1});
            nGs = length(obj.Cox{1}{1});
            ns = obj.nSub^2;
            fU = cell(nLayer,1);
            for kl=1:nLayer
                % Topelitz prameters
                nU = obj.Cox{1}{kl}{1}{1}.nU;
                mu = obj.Cox{1}{kl}{1}{1}.mu(:);
                na = obj.Cox{1}{kl}{1}{1}.na;
                
                % Shuffling
                fU{kl} = zeros(nU,2*nGs);
                vx = 1:ns;
                vy = vx + ns;
                for ii=1:nGs
                    fU{kl}(mu,2*ii-1) = yy(vx);
                    fU{kl}(mu,2*ii) = yy(vy);
                    vx = vx + 2*ns;
                    vy = vy + 2*ns;
                end
                
                % pre forward 1-D fft
                fU{kl} = fft(fU{kl},na);
            end

            % vector-vecotr element-wise product, sum and backward fft
            out = cell(nLayer,1);
            for kl=1:nLayer
                xi = obj.Cox{1}{kl}{1}{1}.xi(:);
                na = obj.Cox{1}{kl}{1}{1}.na;
                out{kl} = zeros(na,nSci);
                for ks = 1:nSci
                    for kg=1:nGs
                        out{kl}(:,ks) = out{kl}(:,ks) + ...
                            mtimesFast(obj.Cox{ks}{kl}{kg}{1},fU{kl}(:,2*kg-1))...
                            + mtimesFast(obj.Cox{ks}{kl}{kg}{2},fU{kl}(:,2*kg));
                    end
                    
                end % CC change 27/04/2018 cropping to output size done outside the nSci and nGs loops
                out{kl} = ifft(out{kl});
                out{kl} = out{kl}(xi,:);
            end

            % interpolation 
            out = cell2mat( out );
            if ~isempty(obj.B)
                out = obj.B*out;
            end
            
            % remove tip/tilt/focus from the estimated phase if needed
            out = obj.ttrmfunPhase(out);

            % masking
            out(~obj.wavefrontMask(:),:) = [];
            out = out*obj.wavefrontToMeter;            
        end
        
        %%
        function out = estimation(obj,wfs)
            % [yy,flag,relres,iter,resvec] = minres(fun,c,1e-3,50);
            if isa(wfs,'shackHartmann')
                obj.c( obj.slopesMask ) = wfs.slopes;
            else
                obj.c( obj.slopesMask ) = wfs;
            end
            yy = minres(obj.fun,obj.c,obj.RTOL,obj.MAXIT,[],[],obj.x0);
            if obj.warmStart
                obj.x0 = yy;
            end
            out = obj.Cox{1}*yy(1:end/2) + obj.Cox{2}*yy(1+end/2:end);
            if length(out)>prod(obj.wavefrontSize)
                out(obj.wavefrontMask) = [];
            else
                out(obj.wavefrontMask) = 0;
                out = reshape(out,obj.wavefrontSize);
            end
            out = out*obj.wavefrontToMeter;
        end
        %%
        function varargout = imagesc(obj)
            %             n = length(obj.Cxx);
            m_Cxx = cellfun( ...
                @(x) [x{1,1}.elements, x{1,2}.elements x{2,1}.elements, x{2,2}.elements] ,...
                obj.Cxx,'UniformOutput',false);
            m_Cxx = cell2mat( m_Cxx );
            m_Cox = cellfun( ...
                @(x) [x{1,1}.elements, x{1,2}.elements] ,...
                obj.Cox{1},'UniformOutput',false);
            %             m_Cox = cell2mat( m_Cox );
            %             subplot(2*n+1,2*n,[1,4*n^2])
            %             imagesc(m_Cxx)
            %             axis square
            %             colorbar
            %             subplot(2*n+1,2*n,[1+4*n^2,(2*n+1)*(2*n)])
            %             imagesc(m_Cox)
            %             axis equal tight
            %             colorbar
            if nargout>0
                varargout{1} = m_Cxx;
                varargout{2} = m_Cox;
            end
        end
        
        %% Compute reconstruction matrix explicitly
        function [R,exCxx,exCox,Cnn] = getReconstructionMatrix(obj)
            
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
                    filter{k,1} = 1;
                    filter{k,2} = 1;
                end
            end
            
            % slopes-to-slopes covariance
            exCxx = cell(obj.nGuideStar*2); 
            for k1=1:obj.nGuideStar
                for k2=1:obj.nGuideStar
                    for t1=1:2
                        for t2=1:2
                            exCxx{2*k1-t1+1,2*k2-t2+1} = filter{k1,t1}*tools.covMap2Matrix(obj.Cxx{k1,k2}{-t1+3,-t2+3}.elements,...
                                obj.nSub,obj.nSub,'mask1',mask{k1},'mask2',mask{k2})*filter{k2,t2}';
                        end
                    end
                end
            end
            exCxx = cell2mat(exCxx);

            % phase-to-slope covariance
            if isinf(obj.guideStar(1).height)
                exCox = cell(1,obj.nGuideStar);
                for k2=1:obj.nGuideStar
                    for t2=1:2
                        exCox{1,2*k2-t2+1} = tools.covMap2Matrix(obj.Cox{1}{1}{k2}{-t2+3}.elements,...
                            obj.nSub,obj.nSub+1,'mask1',mask{k2},'mask2',obj.wavefrontMask)*filter{k2,t2}';
                    end
                end
                exCox = cell2mat(exCox);
            else
                exCox = cell(obj.atmModel.nLayer,2*obj.nGuideStar);
                for k1=1:obj.atmModel.nLayer
                    for k2=1:obj.nGuideStar
                        for t2=1:2
                            n = size(obj.Cox{1}{k1}{k2}{-t2+3}.elements,1)-obj.nSub+1;
                            exCox{k1,2*k2-t2+1} = tools.covMap2Matrix(obj.Cox{1}{k1}{k2}{-t2+3}.elements,...
                                obj.nSub,n,'mask1',mask{k2})*filter{k2,t2}';
                        end
                    end
                end
                exCox = obj.B*cell2mat(exCox);
                exCox = exCox(obj.wavefrontMask(:),:);
            end

            % noise
            if isscalar(obj.noiseVar)
                Cnn = speye(size(exCxx,1))*obj.noiseVar;
            elseif isvector(obj.noiseVar)
                Cnn = diag(obj.noiseVar(obj.slopesMask(:)));
            else
                Cnn = obj.noiseVar;
            end
            
            R = exCox*pinv(exCxx+Cnn);
            R = R*obj.wavefrontToMeter;
        end
        
    end
    
    methods (Access=private)

        function createTTRMfunc(obj)
            %% create functins to remove TT or/and Focus

            obj.ttrmfunSlope = cell(length(obj.guideStar),1);
            if ~obj.isTTRM && ~obj.isFocusRM
                for iT=1:length(obj.tMulti)
                    for iGs=1:obj.nGuideStar
                        iN = iGs + obj.nGuideStar*(iT-1);
                        obj.ttrmfunSlope{iN} = @(x) x;
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
                        zers = zernike(1:max(mode),...
                            obj.D,...
                            'resolution',obj.nSub,...
                            'pupil',double(mask));
                        RMs = [zers.xDerivative(:,mode); zers.yDerivative(:,mode);];
                        iRMs = pinv(RMs);
                        iN = iGs + obj.nGuideStar*(iT-1);
                        obj.ttrmfunSlope{iN} = @(x) x-RMs*(iRMs*x);
                    end
                end
                zerp = zernike(1:max(mode),...
                    obj.D,...
                    'resolution',obj.nSub+1,...
                    'pupil',double(obj.wavefrontMask));
                RMp = zerp.modes(:,mode);
                iRMp = pinv(RMp);
                obj.ttrmfunPhase = @(x) x-RMp*(iRMp*x);
            end
            
        end
        
        function  map = slopestoSlopesCovarianceMapFFT(obj,deltaSrc_x,deltaSrc_y,lag)
            %% slopestoSlopesCovarianceMapFFT
            %
            % map =
            % slopestoSlopesCovarianceMapFFT(obj,deltaSrc_x,deltaSrc_y)
            % computes slope-slope covariance model with the FFT model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : x-x covariance
            %     map{1,2} : x-y covariance
            %     map{2,1} : y-x covariance
            %     map{2,2} : y-y covariance

            covxx = zeros( obj.p_NF);
            covyy = zeros( obj.p_NF);
            covxy = zeros( obj.p_NF);
            b0 = obj.p_NF/2+1;
            b  = (1-obj.nSub:obj.nSub-1)*obj.p_sf + b0;
            [fx0,fy0] = freqspace(obj.p_NF,'meshgrid');
            
            for kLayer = 1:obj.atmModel.nLayer
                
                r0 = (obj.atmModel.r0^(-5/3)*obj.atmModel.layer(kLayer).fractionnalR0)^(-3/5);
                L0 = obj.atmModel.layer(kLayer).layeredL0;
                h = obj.atmModel.layer(kLayer).altitude;
                ws = obj.atmModel.layer(kLayer).windSpeed;
                theta = obj.atmModel.layer(kLayer).windDirection;
                
                g = 1-h/obj.guideStar(1).height;
                dg = obj.dSub*g;
                lf = obj.p_sf/(2*dg);
                fx = lf*fx0;
                fy = lf*fy0;
                delta = 2*lf/obj.p_NF;
                dw = ws*[cos(theta);sin(theta)]*lag;
                dx = h*deltaSrc_x + dw(1);
                dy = h*deltaSrc_y + dw(2);
                
                phasor = exp( 2*1i*pi.*( dx*fx + dy*fy ) );
                
                spec = @(fx,fy,u,v) obj.atmModel.wavelength.^2*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
                    delta.^2.*spectrum(hypot(fx,fy),r0,L0).*...
                    (tools.sinc(dg*fx).*tools.sinc(dg*fy)).^2.*phasor;
                
                covxx = covxx + spec(fx,fy,[1,0],[1,0]);
                covyy = covyy + spec(fx,fy,[0,1],[0,1]);
                covxy = covxy + spec(fx,fy,[0,1],[1,0]);
            end
            
            covxx = real( fftshift( fft2( fftshift( covxx ) ) ) );
            map{1,1} = covxx(b,b);
            covyy = real( fftshift( fft2( fftshift( covyy ) ) ) );
            map{2,2} = covyy(b,b);
            covxy = real( fftshift( fft2( fftshift( covxy ) ) ) );
            map{1,2} = covxy(b,b);
            map{2,1} = map{1,2};
            
        end
        
        function  map = slopestoSlopesCovarianceMapHudgin(obj,deltaSrc_x,deltaSrc_y,lag)
            %% slopestoSlopesCovarianceMapHudgin
            %
            % map =
            % slopestoSlopesCovarianceMapFFT(obj,deltaSrc_x,deltaSrc_y)
            % computes slope-slope covariance model with the Hudgin-like
            % 2-points model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : x-x covariance
            %     map{1,2} : x-y covariance
            %     map{2,1} : y-x covariance
            %     map{2,2} : y-y covariance

            map{1,1} = zeros(2*obj.nSub-1);
            map{1,2} = zeros(2*obj.nSub-1);
            map{2,1} = zeros(2*obj.nSub-1);
            map{2,2} = zeros(2*obj.nSub-1);
            z = meshgrid((1-obj.nSub:obj.nSub-1));
            z = z+z'*1i;
            
            height = obj.guideStar(1).height;
            
            for kLayer = 1:obj.atmModel.nLayer
                
                r0 = (obj.atmModel.r0^(-5/3)*obj.atmModel.layer(kLayer).fractionnalR0)^(-3/5);
                L0 = obj.atmModel.layer(kLayer).layeredL0;
                h = obj.atmModel.layer(kLayer).altitude;
                ws = obj.atmModel.layer(kLayer).windSpeed;
                theta = obj.atmModel.layer(kLayer).windDirection;
                
                g = 1-h/height;
                dg = obj.dSub*g;
                dw = ws*[cos(theta);sin(theta)]*lag;
                dx = h*deltaSrc_x + dw(1);
                dy = h*deltaSrc_y + dw(2);
                zh = z*dg - dx - dy*1i;
                
                Dphi = @(z_) structureFunction(abs(z_),r0,L0);
                
                A = Dphi(zh);
                covxx = -A+Dphi(zh+dg)+Dphi(zh-dg)-A;
                covyy = -A+Dphi(zh+dg*1i)+Dphi(zh-dg*1i)-A;
                covxy = Dphi(zh+dg/2*1i+dg/2)-Dphi(zh-dg/2*1i+dg/2)-...
                    Dphi(zh+dg/2*1i-dg/2)+Dphi(zh-dg/2*1i-dg/2);
                
                map{1,1} = map{1,1} + covxx;
                map{2,2} = map{2,2} + covyy;
                map{1,2} = map{1,2} + covxy;
            end
            
            rad2m = obj.atmModel.wavelength/(2*pi);
            map{1,1} = rad2m*map{1,1}*rad2m/2/obj.dSub/obj.dSub;
            map{2,2} = rad2m*map{2,2}*rad2m/2/obj.dSub/obj.dSub;
            map{1,2} = rad2m*map{1,2}*rad2m/2/obj.dSub/obj.dSub;
            map{2,1} = map{1,2};
            
        end
        
        function  map = slopestoSlopesCovarianceMapFried(obj,deltaSrc_x,deltaSrc_y,lag)
            %% slopestoSlopesCovarianceMapFried
            %
            % map =
            % slopestoSlopesCovarianceMapFFT(obj,deltaSrc_x,deltaSrc_y)
            % computes slope-slope covariance model with the Fried 
            % 4-points model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : x-x covariance
            %     map{1,2} : x-y covariance
            %     map{2,1} : y-x covariance
            %     map{2,2} : y-y covariance

            map{1,1} = zeros(2*obj.nSub-1);
            map{1,2} = zeros(2*obj.nSub-1);
            map{2,1} = zeros(2*obj.nSub-1);
            map{2,2} = zeros(2*obj.nSub-1);
            z = meshgrid((1-obj.nSub:obj.nSub-1));
            z = z+z'*1i;
            
            height = obj.guideStar(1).height;
            
            for kLayer = 1:obj.atmModel.nLayer
                
                r0 = (obj.atmModel.r0^(-5/3)*obj.atmModel.layer(kLayer).fractionnalR0)^(-3/5);
                L0 = obj.atmModel.layer(kLayer).layeredL0;
                h = obj.atmModel.layer(kLayer).altitude;
                ws = obj.atmModel.layer(kLayer).windSpeed;
                theta = obj.atmModel.layer(kLayer).windDirection;
                
                g = 1-h/height;
                dg = obj.dSub*g;
                dw = ws*[cos(theta);sin(theta)]*lag;
                dx = h*deltaSrc_x + dw(1);
                dy = h*deltaSrc_y + dw(2); 
                zh = z*dg - dx - dy*1i;
                
                Dphi = @(z_) structureFunction(abs(z_),r0,L0);

                A = Dphi(zh);
                B1 = Dphi(zh+dg)+Dphi(zh-dg);
                B2 = Dphi(zh+dg*1i)+Dphi(zh-dg*1i);
                C1 = Dphi(zh+dg+dg*1i)+Dphi(zh-dg-dg*1i);
                C2 = Dphi(zh+dg-dg*1i)+Dphi(zh-dg+dg*1i);
                covxx = -4*A+2*B1-2*B2+C1+C2;
                covyy = -4*A-2*B1+2*B2+C1+C2;
                covxy = C1-C2;

                map{1,1} = map{1,1} + covxx;
                map{2,2} = map{2,2} + covyy;
                map{1,2} = map{1,2} + covxy;
            end

            rad2m = obj.atmModel.wavelength/(2*pi);
            map{1,1} = rad2m*map{1,1}*rad2m/8/obj.dSub/obj.dSub;
            map{2,2} = rad2m*map{2,2}*rad2m/8/obj.dSub/obj.dSub;
            map{1,2} = rad2m*map{1,2}*rad2m/8/obj.dSub/obj.dSub;
            map{2,1} = map{1,2};
            
        end
        
        function map = phaseToSlopesCovarianceMapFFTSingleLayer(obj,pad,atmLayer,deltaSrc_x,deltaSrc_y,gl,lag)
            %% phaseToSlopesCovarianceMapFFTSingleLayer
            %
            % map =
            % phaseToSlopesCovarianceMapFFTSingleLayer(obj,pad,atmLayer,deltaSrc_x,deltaSrc_y,gl)
            % computes phase-slope covariance model with the FFT model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : phase-x covariance
            %     map{1,2} : phase-y covariance
            
            alpha = 1;
            
            r0 = (atmLayer.r0^(-5/3)*atmLayer.layer.fractionnalR0)^(-3/5);
            L0 = atmLayer.layer.layeredL0;
            h = atmLayer.layer.altitude;
            ws = atmLayer.layer.windSpeed;
            theta = atmLayer.layer.windDirection;
            
            [fx,fy] = freqspace(obj.p_NF,'meshgrid');
            dg = obj.dSub*gl;
            lf = obj.p_sf/(obj.dSub*2*gl);
            fx = lf*alpha*fx;
            fy = lf*alpha*fy;
            delta = 2*lf*alpha/obj.p_NF;
            
            dw = ws*[cos(theta);sin(theta)]*lag;
            dx = h*deltaSrc_x + dw(1);
            dy = h*deltaSrc_y + dw(2);
            
            phasor = exp( 2*1i*pi.*( (dx+0.5*obj.dSub)*fx + (dy+0.5*obj.dSub)*fy ) );
            
            spec = @(u) atmLayer.wavelength.*1i*(fx.*u(1) + fy.*u(2)).*...
                delta.^2.*spectrum(hypot(fx,fy),r0,L0).*...
                tools.sinc(dg*fx).*tools.sinc(dg*fy).*phasor./gl;
            
            b0 = obj.p_NF/2+1;
            b = (1-obj.nSub-pad:obj.nSub+pad)*obj.p_sf + b0;
            
            covx  = fftshift(real( fft2( fftshift( spec([1,0]) ) ) ) );
            covy  = fftshift(real( fft2( fftshift( spec([0,1]) ) ) ) );
            map{1} = covx(b,b);
            map{2} = covy(b,b);
            
        end
           
        function map = phaseToSlopesCovarianceMapHudginSingleLayer(obj,pad,atmLayer,deltaSrc_x,deltaSrc_y,gl,lag)
            %% phaseToSlopesCovarianceMapHudginSingleLayer
            %
            % map =
            % phaseToSlopesCovarianceMapHudginSingleLayer(obj,pad,atmLayer,deltaSrc_x,deltaSrc_y,gl)
            % computes phase-slope covariance model with the FFT model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : phase-x covariance
            %     map{1,2} : phase-y covariance
            
            z = meshgrid(1-obj.nSub-pad:obj.nSub+pad);
            z = z - 0.5;
            z = z+z'*1i;
            
            dg = obj.dSub*gl;
            r0 = (atmLayer.r0^(-5/3)*atmLayer.layer.fractionnalR0)^(-3/5);
            L0 = atmLayer.layer.layeredL0;
            h = atmLayer.layer.altitude;
            ws = atmLayer.layer.windSpeed;
            theta = atmLayer.layer.windDirection;
            dw = ws*[cos(theta);sin(theta)]*lag;
            dx = h*deltaSrc_x + dw(1);
            dy = h*deltaSrc_y + dw(2);
            zh = z*dg - dx - dy*1i;
            
            Dphi = @(z_) covariance(abs(z_),r0,L0);
            
            covx = Dphi(zh+dg/2)-Dphi(zh-dg/2);
            covy = Dphi(zh+dg/2*1i)-Dphi(zh-dg/2*1i);
            
            rad2m = obj.atmModel.wavelength/(2*pi);
            map{1} = -covx*rad2m/obj.dSub;
            map{2} = -covy*rad2m/obj.dSub;
            
        end

        function map = phaseToSlopesCovarianceMapFriedSingleLayer(obj,pad,atmLayer,deltaSrc_x,deltaSrc_y,gl,lag)
            %% phaseToSlopesCovarianceMapFriedSingleLayer
            %
            % map =
            % phaseToSlopesCovarianceMapFriedSingleLayer(obj,pad,atmLayer,deltaSrc_x,deltaSrc_y,gl)
            % computes phase-slope covariance model with the FFT model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : phase-x covariance
            %     map{1,2} : phase-y covariance
            
            z = meshgrid(1-obj.nSub-pad:obj.nSub+pad);
            z = z - 0.5;
            z = z+z'*1i;
            
            dg = obj.dSub*gl;
            r0 = (atmLayer.r0^(-5/3)*atmLayer.layer.fractionnalR0)^(-3/5);
            L0 = atmLayer.layer.layeredL0;
            h = atmLayer.layer.altitude;
            ws = atmLayer.layer.windSpeed;
            theta = atmLayer.layer.windDirection;
            dw = ws*[cos(theta);sin(theta)]*lag;
            dx = h*deltaSrc_x + dw(1);
            dy = h*deltaSrc_y + dw(2);
            zh = z*dg - dx - dy*1i;
            
            Dphi = @(z_) covariance(abs(z_),r0,L0);
            
            A = Dphi(zh+dg/2+dg/2*1i);
            B = Dphi(zh+dg/2-dg/2*1i);
            C = Dphi(zh-dg/2+dg/2*1i);
            D = Dphi(zh-dg/2-dg/2*1i);
            covx = 0.5*(A+B-C-D);
            covy = 0.5*(A-B+C-D);
                        
            rad2m = obj.atmModel.wavelength/(2*pi);
            map{1} = -covx*rad2m/obj.dSub;
            map{2} = -covy*rad2m/obj.dSub;
            
        end
        
        function map = phaseToSlopesCovarianceMapFFTAllLayer(obj,deltaSrc_x,deltaSrc_y,lag)
            %% phaseToSlopesCovarianceMapFFTAllLayer
            %
            % map =
            % phaseToSlopesCovarianceMapFFTAllLayer(obj,deltaSrc_x,deltaSrc_y)
            % computes phase-slope covariance model with the FFT model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : phase-x covariance
            %     map{1,2} : phase-y covariance
            
            alpha = 1;
            [fx,fy] = freqspace(obj.p_NF,'meshgrid');
            lf = obj.p_sf/(obj.dSub*2);
            fx = lf*alpha*fx;
            fy = lf*alpha*fy;
            delta = 2*lf*alpha/obj.p_NF;
            covx = zeros( obj.p_NF );
            covy = zeros( obj.p_NF );
            for kLayer=1:obj.atmModel.nLayer
                r0 = (obj.atmModel.r0^(-5/3)*obj.atmModel.layer(kLayer).fractionnalR0)^(-3/5);
                L0 = obj.atmModel.layer(kLayer).layeredL0;
                h = obj.atmModel.layer(kLayer).altitude;
                ws = obj.atmModel.layer(kLayer).windSpeed;
                theta = obj.atmModel.layer(kLayer).windDirection;
                
                dw = ws*[cos(theta);sin(theta)]*lag;
                dx = h*deltaSrc_x + dw(1);
                dy = h*deltaSrc_y + dw(2);
                
                phasor = exp( 2*1i*pi.*( (dx+0.5*obj.dSub)*fx + (dy+0.5*obj.dSub)*fy ) );
                
                spec = @(u) obj.atmModel.wavelength.*1i*(fx.*u(1) + fy.*u(2)).*...
                    delta.^2.*spectrum(hypot(fx,fy),r0,L0).*...
                    tools.sinc(obj.dSub*fx).*tools.sinc(obj.dSub*fy).*phasor;
                
                b0 = obj.p_NF/2+1;
                b = (1-obj.nSub:obj.nSub)*obj.p_sf + b0;
                covx  = covx + fftshift(real( fft2( fftshift( spec([1,0]) ) ) ) );
                covy  = covy + fftshift(real( fft2( fftshift( spec([0,1]) ) ) ) );
            end
            map{1} = covx(b,b);
            map{2} = covy(b,b);
            
        end
        
        function map = phaseToSlopesCovarianceMapHudginAllLayer(obj,deltaSrc_x,deltaSrc_y,lag)
            %% phaseToSlopesCovarianceMapHudginAllLayer
            %
            % map =
            % phaseToSlopesCovarianceMapHudginAllLayer(obj,deltaSrc_x,deltaSrc_y)
            % computes phase-slope covariance model with the Hudgin-like 2 points model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : phase-x covariance
            %     map{1,2} : phase-y covariance

            z = meshgrid(1-obj.nSub:obj.nSub);
            z = z - 0.5;
            z = z+z'*1i;
            covx = zeros(2*obj.nSub);
            covy = zeros(2*obj.nSub);
            
            for kLayer=1:obj.atmModel.nLayer
                r0 = (obj.atmModel.r0^(-5/3)*obj.atmModel.layer(kLayer).fractionnalR0)^(-3/5);
                L0 = obj.atmModel.layer(kLayer).layeredL0;
                h = obj.atmModel.layer(kLayer).altitude;
                ws = obj.atmModel.layer(kLayer).windSpeed;
                theta = obj.atmModel.layer(kLayer).windDirection;
                
                dw = ws*[cos(theta);sin(theta)]*lag;
                dx = h*deltaSrc_x + dw(1);
                dy = h*deltaSrc_y + dw(2);
                                
                zh = z*obj.dSub - dx - dy*1i;
                
                Dphi = @(z_) covariance(abs(z_),r0,L0);
                
                covx = covx + Dphi(zh+obj.dSub/2)-Dphi(zh-obj.dSub/2);
                covy = covy + Dphi(zh+obj.dSub/2*1i)-Dphi(zh-obj.dSub/2*1i);
            end
            
            rad2m = obj.atmModel.wavelength/(2*pi);
            map{1} = -covx*rad2m/obj.dSub;
            map{2} = -covy*rad2m/obj.dSub;
        end
        
        function map = phaseToSlopesCovarianceMapFriedAllLayer(obj,deltaSrc_x,deltaSrc_y,lag)
            %% phaseToSlopesCovarianceMapFriedAllLayer
            %
            % map =
            % phaseToSlopesCovarianceMapHudginAllLayer(obj,deltaSrc_x,deltaSrc_y)
            % computes phase-slope covariance model with the Fried 4 points model.
            % deltaSrc_x and deltaSrc_y shoud be in angle of radian.
            % Output is covariance maps in squre of angle of radian.
            %     map{1,1} : phase-x covariance
            %     map{1,2} : phase-y covariance

            z = meshgrid(1-obj.nSub:obj.nSub);
            z = z - 0.5;
            z = z+z'*1i;
            covx = zeros(2*obj.nSub);
            covy = zeros(2*obj.nSub);
            
            for kLayer=1:obj.atmModel.nLayer
                r0 = (obj.atmModel.r0^(-5/3)*obj.atmModel.layer(kLayer).fractionnalR0)^(-3/5);
                L0 = obj.atmModel.layer(kLayer).layeredL0;
                h = obj.atmModel.layer(kLayer).altitude;
                ws = obj.atmModel.layer(kLayer).windSpeed;
                theta = obj.atmModel.layer(kLayer).windDirection;
                
                dw = ws*[cos(theta);sin(theta)]*lag;
                dx = h*deltaSrc_x + dw(1);
                dy = h*deltaSrc_y + dw(2);
                                
                zh = z*obj.dSub - dx - dy*1i;
                
                Dphi = @(z_) covariance(abs(z_),r0,L0);
                
                A = Dphi(zh+obj.dSub/2+obj.dSub/2*1i);
                B = Dphi(zh+obj.dSub/2-obj.dSub/2*1i);
                C = Dphi(zh-obj.dSub/2+obj.dSub/2*1i);
                D = Dphi(zh-obj.dSub/2-obj.dSub/2*1i);
                covx = covx + 0.5*(A+B-C-D);
                covy = covy + 0.5*(A-B+C-D);
            end
            
            rad2m = obj.atmModel.wavelength/(2*pi);
            map{1} = -covx*rad2m/obj.dSub;
            map{2} = -covy*rad2m/obj.dSub;
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
            
            % Tompraphic Reconstructor
            mmse =slopesLinearMMSE(wfs,tel,atm,ast,'mmseStar',science,...
                'noiseVar', 1e-14);
            %R = mmse.getReconstructionMatrix; % MVM case
            
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
            % Strehl Ratio = 46.79
        end
        
    end
    
end


function out = mtimes4squareBlocksVaryingSlopesMaskYoshito( OO, w, mask, noiseVar, ttrmfun)

n = length(OO);
nw = length(w);
ns = nw/n/2;

% Topelitz prameters
nU = OO{1,1}{1,1}.nU;
mu = OO{1,1}{1,1}.mu(:);
xi = OO{1,1}{1,1}.xi(:);
na = OO{1,1}{1,1}.na;

% Shuffling
fU = zeros(nU,2*n);
vx = 1:ns;
vy = vx + ns;
for ii=1:n
    fU(mu,2*ii-1) = w(vx);
    fU(mu,2*ii) = w(vy);
    vx = vx + 2*ns;
    vy = vy + 2*ns;
end

% pre forward 1-D fft
fU = fft(fU,na);

%
tmp = zeros(na,2*n);
for ii=1:n
    for jj=1:n
        tmp(:,2*ii-1) = tmp(:,2*ii-1)...
            + mtimesFast(OO{ii,jj}{1,1},fU(:,2*jj-1))...
            + mtimesFast(OO{ii,jj}{1,2},fU(:,2*jj));
        tmp(:,2*ii) = tmp(:,2*ii)...
            + mtimesFast(OO{ii,jj}{2,1},fU(:,2*jj-1))...
            + mtimesFast(OO{ii,jj}{2,2},fU(:,2*jj));
    end
end

tmp = ifft(tmp);
tmp = tmp(xi,:);

% noise
if isvector(noiseVar)
    Cnw = noiseVar.*w;
else
    Cnw = noiseVar*w;
end

% mask
out = mask(:).*(tmp(:)+Cnw(:));


% low-order removal
v=0;
for ii=1:n
    v = (1:nw/n)+v(end);
    out(v) = ttrmfun{ii}(out(v));
end

end


function out = mtimes4squareBlocksVaryingSlopesMask( OO, w, mask )

n = length(OO);
nw = length(w);
u = 1:nw/n;
out = zeros( size(w) );
if nargin>2
    for ii=1:n
        v = 1:nw/n;
        mask_ = mask(1:end/2,ii);
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                maskedMtimes(O{1,1},x,mask_) + maskedMtimes(O{1,2},y,mask_)
                maskedMtimes(O{2,1},x,mask_) + maskedMtimes(O{2,2},y,mask_)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
else
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                mtimes(O{1,1},x) + mtimes(O{1,2},y)
                mtimes(O{2,1},x) + mtimes(O{2,2},y)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
end

end


function out = mtimes4squareBlocks( OO, w, mask )

n = length(OO);
nw = length(w);
u = 1:nw/n;
out = zeros( size(w) );
if nargin>2
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                maskedMtimes(O{1,1},x,mask) + maskedMtimes(O{1,2},y,mask)
                maskedMtimes(O{2,1},x,mask) + maskedMtimes(O{2,2},y,mask)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
else
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                mtimes(O{1,1},x) + mtimes(O{1,2},y)
                mtimes(O{2,1},x) + mtimes(O{2,2},y)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
end

end


function out = mtimes4precondNewCcorreia( OO, w, mask)

n = length(OO);
nw = length(w);
u = 1:nw/n;
out = zeros( size(w) );
if nargin>2
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                invMaskedMtimes(O{1,1},x,mask) + invMaskedMtimes(O{1,2},y,mask)
                invMaskedMtimes(O{2,1},x,mask) + invMaskedMtimes(O{2,2},y,mask)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
else
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                invMtimes(O{1,1},x) + invMtimes(O{1,2},y)
                invMtimes(O{2,1},x) + invMtimes(O{2,2},y)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
end

end


function out = mtimes4precond( OO, w, mask)

x = w(1:end/2);
y = w(1+end/2:end);
O = OO{1};
A = O{1,1};
B = O{2,2};
out = [medfilt1( A\x ) ; medfilt1( B\y )];

end

function out = structureFunction(rho,r0,L0)
%% STRUCTUREFUNCTION Phase structure function
%
% out = phaseStats.structureFunction(rho,atm) computes the
% phase structure function from the baseline rho and an
% atmosphere object
%
% See also atmosphere

if isinf(L0)
    out   = zeros(size(rho));
    index = rho~=0;
    out(index) = 2.*(24.*gamma(6./5)./5).^(5./6).*(rho(index)./r0).^(5./3);
else
    out = 2.*(variance(r0,L0)-covariance(rho,r0,L0));
end
end

function out = spectrum(f,r0,L0)
%% SPECTRUM Phase power spectrum density
%
% out = phaseStats.spectrum(f,atm) computes the phase power
% spectrum density from the spatial frequency f and an
% atmosphere object
%
% See also atmosphere

out = (24.*gamma(6./5)./5).^(5./6).*...
    (gamma(11./6).^2./(2.*pi.^(11./3))).*...
    r0.^(-5./3);
out = out.*(f.^2 + 1./L0.^2).^(-11./6);
end

function out = covariance(rho,r0,L0)
%% COVARIANCE Phase covariance
%
% out = phaseStats.covariance(rho,atm) computes the phase covariance from
% the baseline rho and an atmosphere object
%
% See also atmosphere

L0r0ratio= (L0./r0).^(5./3);
cst      = (24.*gamma(6./5)./5).^(5./6).*...
    (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
    L0r0ratio;
out   = ones(size(rho)).*(24.*gamma(6./5)./5).^(5./6).*...
    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
index         = rho~=0;
u             = 2.*pi.*rho(index)./L0;
out(index) = cst.*u.^(5./6).*besselk(5./6,u);
end

function  out = variance(r0,L0)
%% VARIANCE Phase variance
%
% out = phaseStats.variance(atm) computes the phase variance
% from an atmosphere object
%
% See also atmosphere
L0r0ratio= (L0./r0).^(5./3);
out   = (24.*gamma(6./5)./5).^(5./6).*...
    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
end