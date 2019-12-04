classdef sparseLMMSE_old < handle
    %% SPARSELMMSE_OLD Create a sparseLMMSE object
    % 
    % 
    
    properties
        tag = 'linearMMSE_old';
        % telescope object
        tel;
        % wave-front sensor object
        wfs;
        % deformable mirror object
        dm;
        % atmospheric model for the reconstruction (may be different from
        % the system atmosphere)
        atmModel;
        % guide star asterism (source class)
        guideStar;
        
        % science optimisation directions (source class)
        scienceStar;
        % discrete gradient operator (for the SH WFS)
        Gamma;
        %
        gridMask;
        % propagator from layers to gridMask pupil locations in guide-star
        % direction
        H;
        % propagator from layers to gridMask pupil locations in science
        % directions
        Hss;
        % propagator from reconstructed layers to DM actuator locations
        Ha;
        % inverse of the noise covariance matrix
        iCn;
        % noise covariance matrix
        Cn;
        % user-defined method to compute the noise covariance matrix: theoretical or empirical
        noiseCovarComputeMethod;
        % Na-layer width
        deltaNa;
        % bi-harmonic approximation to iCphi
        L2;
        % inverse phase covariance matrix
        iCphi;
        %
        mask;
        % forward model for linear system of equation Ax=b
        A;
        % left-hand side of Ax=b
        b;
        % Cholesky upper/lower diagonal matrix
        L;
        % string for the algorithm/solver to find solution x in Ax=b
        solver;
        % the wfs input data
        data;
        % used to predict through bilinear interpolation the layered phase
        % lag seconds ahead
        lag;
        % normalisation variable
        q;
        % slopes location mask
        slopesMask;
        % history of estimates for iterative warm-restart
        psEstHist;
        % fitting matrix
        fitMatrix;
        %fitting-step solver
        fitSolver;
    end
    properties (Access=private)%(Dependent)
        %number of guide-star directions
        nGuideStar;
        % number of science optimisation directions
        nScienceStar;
        % handle to private function A*x
        AmultFun
    end
    
    
    methods
        %% Constructor
        function obj = sparseLMMSE_old(tel,wfs,dm,atmModel,guideStar,scienceStar, varargin)
            inputs = inputParser;
            inputs.addRequired('tel',@(x) isa(x,'telescopeAbstract') );
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann') );
            inputs.addRequired('dm',@(x) isa(x,'deformableMirror') );
            inputs.addRequired('atmModel',@(x) isa(x,'atmosphere') );
            inputs.addRequired('guideStar',@(x) isa(x,'source') );
            inputs.addRequired('scienceStar',@(x) isa(x,'source') );
            inputs.addParameter('solver','Cholesky',@ischar);
            inputs.addParameter('sparseMtimes','handle',@ischar);
            inputs.addParameter('lag',0,@isnumeric);
            inputs.addOptional('slopesMask',[], @isnumeric );
            inputs.addOptional('iCn',[],@isnumeric);
            inputs.addOptional('noiseCovarComputeMethod','theoretical',@ischar);
            inputs.addOptional('deltaNa',0,@isnumeric); % default point-source LGS
            inputs.addOptional('fitMatrix',[],@isnumeric);
            inputs.addOptional('fitSolver','backSolve',@ischar);
            inputs.parse(tel,wfs,dm,atmModel,guideStar,scienceStar,varargin{:});
            obj.tel         = inputs.Results.tel;
            obj.wfs         = inputs.Results.wfs;
            obj.dm          = inputs.Results.dm;
            obj.atmModel    = inputs.Results.atmModel;
            obj.guideStar   = inputs.Results.guideStar;
            obj.noiseCovarComputeMethod = inputs.Results.noiseCovarComputeMethod;
            obj.deltaNa     = inputs.Results.deltaNa;
            obj.scienceStar = inputs.Results.scienceStar;
            %obj.nScienceStar= length(obj.scienceStar);
            obj.solver      = inputs.Results.solver;
            if ~isempty(inputs.Results.lag)
                obj.lag     = inputs.Results.lag;
            else
                obj.lag     = 0;
            end
            if ~isempty(inputs.Results.slopesMask)
                obj.slopesMask = inputs.Results.slopesMask;
                % convert to rasterised index for vectorised indexing
                for iGs = 1:obj.nGuideStar
                    obj.slopesMask(:,iGs) = obj.slopesMask(:,iGs) + (iGs-1)*obj.wfs.nValidLenslet*2;
                end
                obj.slopesMask = obj.slopesMask(:);
            else
                obj.slopesMask    = find(repmat( repmat( ones(nnz(inputs.Results.wfs(1).validLenslet(:)),1) , 2, 1), 1, obj.nGuideStar ));
            end
            
            %Compute individual compouds
            setGamma(obj);
            setH(obj);
            setL2(obj);
            
            % inverse noise covariance matrix
            if isempty(inputs.Results.iCn)
                if strcmp(obj.noiseCovarComputeMethod,'theoretical')
                    fprintf('Computing theoretical noise covariance matrix\n');
                    calcTheoreticaliCn(obj);
                else
                    fprintf('Computing empirical noise covariance matrix\n');
                    calcEmpiricaliCn(obj);
                end
            else
                obj.iCn =  inputs.Results.iCn;
                d = obj.tel.D/obj.wfs.lenslets.nLenslet;
                obj.q = 2*pi/d/(obj.wfs.lenslets.nyquistSampling*2);
            end
            

            % Fitting step
            if isempty( inputs.Results.fitMatrix)
                setHa(obj);
                obj.fitSolver = 'backSolve';
            else
                obj.fitMatrix = inputs.Results.fitMatrix;
                obj.fitSolver = 'vmm';
            end
            
            
            if ~isempty(inputs.Results.sparseMtimes)
                % use separate comopunds of A
                obj.AmultFun = @(x) AmultSplit(x);
                obj.b = (obj.Gamma*obj.H)'*obj.iCn;
            else % assemble A and b
                obj.AmultFun = @(x) Amult(x);
                buildSparseModel(obj);
            end
            
            
            if strcmp(obj.solver,'Cholesky')
                obj.L = chol(obj.A,'lower');
            end
            
            obj.data = zeros(size(obj.Gamma,1),1);
            obj.psEstHist = zeros(size(obj.A,1),1);
        end
        
        %% Gets and Sets
%          function set.Cn(obj,Cn)
%              obj.Cn = Cn;%lamTools.invertMeasurementNoiseCovarianceMatrix(Cn);
%              %buildSparseModel(obj);
%          end
        function set.iCn(obj,iCn)
            obj.iCn = iCn;
            buildSparseModel(obj);
        end
        
        %% Get nGuideStar
        function out = get.nGuideStar(obj)
            out  = length(obj.guideStar);
        end
        %% Get nScienceStar
        function out = get.nScienceStar(obj)
            out  = length(obj.scienceStar);
        end
        %% set the sparse gradient matrix using a 3x3 stencil
        function setGamma(obj)
            [p_Gamma,p_gridMask] = sparseGradientMatrix3x3Stentil(obj.wfs);
            d = obj.tel.D/obj.wfs.lenslets.nLenslet;
            p_Gamma = p_Gamma/d;
            p_Gamma = repmat({p_Gamma},1,obj.nGuideStar);
            obj.Gamma = blkdiag(p_Gamma{:});
            obj.Gamma = obj.Gamma(obj.slopesMask,:); % select only the valid sub-apertures from current listing
            obj.gridMask = p_gridMask;
        end
        
        %% set propagator H from GS to WFS and from science star to pupil
        function setH(obj)
            [p_H,obj.mask] = bilinearSplineInterpMat([obj.guideStar,obj.scienceStar],obj.atmModel,obj.tel,obj.gridMask,obj.lag);
            p_Hss = p_H(obj.nGuideStar+1:end,:);
            p_Hss = cell2mat(p_Hss);
            p_Hss(:,~cell2mat(obj.mask)) = [];
            p_H(obj.nGuideStar+1:end,:) = [];
            p_H = cell2mat(p_H);
            p_H(:,~cell2mat(obj.mask)) = [];
            obj.H = p_H;
            obj.Hss = p_Hss;
        end
        
        % set bi-harmonic operator (appox to inverse phase covariance
        % matrix)
        function setL2(obj)
            d = obj.tel.D/obj.wfs.lenslets.nLenslet;
            p_L2 = phaseStats.sparseInverseCovarianceMatrix(obj.mask,d, obj.atmModel);
            obj.L2 = blkdiag(p_L2{:});
        end
        
        %% set inverse of noise covariance matrix
        function calcTheoreticaliCn(obj)
            nGs = length(obj.guideStar);
            lgsLaunchCoord = reshape([obj.guideStar.viewPoint],2,nGs)';
            heightNa = [obj.guideStar.height]';
            buffer = obj.wfs.theoreticalNoise(obj.tel, obj.atmModel, obj.guideStar, obj.scienceStar,...
                'lgsLaunchCoord',lgsLaunchCoord,...
                'naParam',[ones(nGs,1)*obj.deltaNa heightNa],'centroidingAlgorithm','cog');
            % convert cell to 3-D mat
            obj.Cn = zeros([size(buffer{1}), nGs]);
            for iGs = 1:nGs
                obj.Cn(:,:,iGs) = buffer{iGs}/pi^2*obj.wfs.lenslets.nyquistSampling; % convert to pixel displacement at the focal plane of the WFS lenslets;;
            end
            obj.iCn = lamTools.invertMeasurementNoiseCovarianceMatrix(obj.Cn);
           
            d = obj.tel.D/obj.wfs.lenslets.nLenslet; % subaperture pitch (OOMAO will need it) lblanco 12/14/2016
            obj.q = 2*pi/d/(obj.wfs.lenslets.nyquistSampling*2);%lblanco 12/14/2016
        end
        %% set inverse of noise covariance matrix
        function calcEmpiricaliCn(obj,nMeas)
            if ~isempty(obj.tel.opticalAberration)
                opticalAberration = obj.tel.opticalAberration;
                obj.tel.opticalAberration = [];
            else
                opticalAberration = [];
            end
            % WFS noise covariance matrix
            fprintf(' Computing the WFS noise covariance matrix ...\n');
            if nargin == 1
                nMeas = 250;
            end
            d = obj.tel.D/obj.wfs.lenslets.nLenslet;
            slopes = zeros(obj.wfs.nSlope,nMeas);
            %ngs = obj.guideStar(1);
            p_wfs = obj.wfs;
            %ngs = ngs.*obj.tel*obj.wfs;
            for kMeas=1:nMeas
                +p_wfs
                slopes(:,kMeas) = p_wfs.slopes;
            end
            obj.q = 2*pi/d/(p_wfs.lenslets.nyquistSampling*2);
            slopes = slopes.*obj.q;
            Cn = slopes*slopes'/nMeas;
            Cn = diag(diag(Cn));
            % iCn = diag(1./diag(Cn));
            p_iCn = sparse(1:p_wfs.nSlope,1:p_wfs.nSlope, 1./diag(Cn) );
            p_iCn = repmat({p_iCn},1,obj.nGuideStar);
            p_iCn = blkdiag(p_iCn{:});
            obj.iCn = p_iCn(obj.slopesMask, obj.slopesMask); % collect only those sub-apertures in the current listing
            obj.tel.opticalAberration = opticalAberration; % put it back
        end
        %% set fit matrix Ha
        function setHa(obj)
            obj.Ha = obj.dm.modes.modes(obj.gridMask(:),:);
        end
        
        %% build the sparse linear system of equations
        function buildSparseModel(obj)
            G = obj.Gamma*obj.H;
            % M = (G'*iCn*G+L2)\(G'*iCn);
            obj.A = G'*obj.iCn*G+obj.L2;
            obj.b = G'*obj.iCn;
        end
        %%
        function reset(obj)
            obj.psEstHist = zeros(size(obj.A,1),1);
        end
        
        %% Wavefront reconstruction
        function out = mtimes(obj,input)
            if isa(input,'shackHartmann')
                obj.data = input.slopes(obj.slopesMask(:))*obj.q;
            else
                obj.data = input(:)*obj.q;
            end
            %% estimation
            switch obj.solver
                case 'Cholesky'
                    psEst = obj.L'\(obj.L\( obj.b * (obj.data))); % Cholesky back-solve
                case 'cg'
                    psEst = cgs(obj.AmultFun,obj.b*(obj.data),[],50,[],[],obj.psEstHist); % conjugate-gradient method
                case 'GaussianElimination'
                    psEst = obj.A\( obj.b * obj.data); % Gauss elimination
            end
            obj.psEstHist = psEst;
            %% forward-projection to aperture in science directions
            psEst = obj.Hss*psEst;
            
            %% fitting
            switch obj.fitSolver
                case 'backSolve'
                    psEst = reshape(psEst,size(obj.Ha,1),obj.nScienceStar);
                    out = (obj.Ha\psEst)/obj.guideStar(1).waveNumber;
                case 'vmm'
                    psEst = reshape(psEst,size(obj.fitMatrix,2),obj.nScienceStar);
                    out = (obj.fitMatrix*psEst)/obj.guideStar(1).waveNumber;
            end
        end
    end
end

function out = AmultSplit(x)
out = obj.G*x;
out = obj.iCn*out;
out = obj.G'*out;
out = out + obj.L2*x;
end

function out = Amult(x)
out = obj.A*x;end
