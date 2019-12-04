classdef turbulenceProfiler < handle
   
    
    properties (SetObservable = true)
        %Sub-classes
        tel;
        atm;
        atmMod;
        atmFit;
        wfs;
        gs;
        mod;
        %Model
        FUN;
        %Measurements
        sGen;           % Generator class
        generatorFlag;  % Either 'simulation' or 'model'
        nFrames;        %# of frames to be simulated
        Css             %Matrix/Map of measurements
        %Fitting stuffs
        OPT;            % OPTIONS class
        maxIter;        %#iterations max
        tolFun;         %tolerance and diff value on FUN between two steps
        indexFit;       % Index vector of parameters to fit
        UB;             %Upper bounds
        LB;             %Lower bounds
        fitr0;          %Flag to fit r0 values for each layer
        fith;           %Flag to fit altitude values for each layer
        fitL0;          %Flag to fit L0 values for each layer
        fitMis;         %Flag to fit mis-registration for each wfs
        flagTT;
        %Outputs
        X0;             %init guess
        X;              % Results from the fit
        Xf;             % What we should retrieve using a model matrix as input
        resFit;
        resnorm;
        outFit;
        flag;
        
    end
    
    methods
        
        %% %%%%%%%%%%%
        % CONSTRUCTOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        function obj = turbulenceProfiler(tel,atmMod,wfs,gs,varargin)            
                            
            %Checking inputs
            inputs = inputParser;
            inputs.addRequired('tel',@(x) isa(x,'telescope') );
            inputs.addRequired('atmMod',@(x) isa(x,'atmosphere'));
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann'));
            inputs.addRequired('gs',@(x) isa(x,'source'));
            inputs.addParameter('atmosphere',atmMod, @(x) isa(x,'atmosphere') );
            inputs.addParameter('generatorFlag','model', @(x) isa(x,'char') );
            inputs.addParameter('nFrames',2048, @isnumeric );
            inputs.addParameter('maxIter',15, @isnumeric );
            inputs.addParameter('tolFun',1e-10, @isnumeric );
            inputs.addParameter('fitr0',0, @isnumeric );
            inputs.addParameter('fith',0, @isnumeric );
            inputs.addParameter('fitL0',0, @isnumeric );
            inputs.addParameter('fitMis',0, @isnumeric );
            inputs.addParameter('flagTT',false, @islogical);
            inputs.parse(tel,atmMod,wfs,gs,varargin{:});
            
            %Instantiation
            obj.tel      = tel;
            obj.atmMod   = atmMod;
            obj.wfs      = wfs;
            obj.gs       = gs;
            %Results from parser
            obj.atm           = inputs.Results.atmosphere;
            obj.generatorFlag = inputs.Results.generatorFlag;
            obj.nFrames       = inputs.Results.nFrames;
            obj.maxIter       = inputs.Results.maxIter;
            obj.tolFun        = inputs.Results.tolFun;
            obj.fitr0         = inputs.Results.fitr0;
            obj.fith          = inputs.Results.fith;
            obj.fitL0         = inputs.Results.fitL0;
            obj.fitMis        = inputs.Results.fitMis;
            obj.flagTT        = inputs.Results.flagTT;
            
            % ------------------ COVARIANCE MODEL ------------------------%
            obj.mod  = slopesCovarianceModel(tel,atmMod,wfs,gs,'flagTT',obj.flagTT);
                                                                                  
            
            % --------------------- MEASUREMENTS -------------------------%
            obj      = obj.covarianceMatrixGenerator();          
            
            % ------------------- FITTING OPTIONS ------------------------%
            obj      = obj.getOptions();
            obj      = obj.getIndexFit(obj.fitr0,obj.fith,obj.fitL0,obj.fitMis);
            obj      = obj.getBounds();
        end
  
        %% %%%%%%%%%%%
        % Turbulence Generator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = covarianceMatrixGenerator(obj)
            if strcmp(obj.generatorFlag,'simulation')                
                obj.sGen = slopesGenerator(obj.tel,obj.atm,obj.wfs,obj.gs,'nFrames',obj.nFrames);
                obj.sGen = obj.sGen.generateOpenLoopSlopes();
                obj.Css  = obj.sGen.slopes * obj.sGen.slopes'./obj.nFrames;
                
            elseif strcmp(obj.generatorFlag,'model')
                obj.Xf  = [35.0,10.0];
                obj.Css = obj.mod.covarianceModel(obj.Xf);
            end
            
            %Filtering tip-tilt if required
            if obj.mod.flagTT
                obj.Css =  obj.mod.manageTipTilt(obj.Css);
            end            
        end
        
        
        %% %%%%%%%%%%%
        %Fitting procedure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function obj = getOptions(obj)
            %obj.OPT         = optimoptions('lsqcurvefit');
            obj.OPT         = optimoptions('lsqnonlin');
            obj.OPT.MaxIter = obj.maxIter;
            obj.OPT.TolFun  = obj.tolFun;
            obj.OPT.TolX    = 1e-10;
        end
        function obj = getBounds(obj)
            %Grabbing the fitting indexes vector
            idx    = obj.indexFit;
            %Tomographic resolution
            dH     = obj.mod.getTomographicResolution();
            %Defining lower bounds for each parameters
            LBr0   = obj.mod.fr0.*0.1;          
            LBh    = obj.mod.h - 2*dH;
            LBL0   = zeros(1,obj.atmMod.nLayer);
            LBX    = 0.85.*obj.mod.xSrc;
            LBY    = 0.85.*obj.mod.ySrc;
            LBG    = 0.85.*ones(1,obj.mod.nwfs);
            LBth   = -180.*ones(1,obj.mod.nwfs);
            LBdX   = 0.85*obj.mod.pitch.*ones(1,obj.mod.nwfs);
            LBdY   = 0.85*obj.mod.pitch.*ones(1,obj.mod.nwfs);
            LBs    = 0.85.*obj.mod.s;
            LBt    = -ones(1,3);
            %Defining upper bounds
            UBr0   = obj.mod.fr0.*10;
            UBh    = obj.mod.h + 2*dH;
            UBL0   = 100.*ones(1,obj.atmMod.nLayer);
            UBX    = 1.15.*obj.mod.xSrc;
            UBY    = 1.15.*obj.mod.ySrc;
            UBG    = 1.15.*ones(1,obj.mod.nwfs);
            UBth   = 180.*ones(1,obj.mod.nwfs);
            UBdX   = 1.15*obj.mod.pitch.*ones(1,obj.mod.nwfs);
            UBdY   = 1.15*obj.mod.pitch.*ones(1,obj.mod.nwfs);
            UBs    = 1.15.*obj.mod.s;
            UBt    = -ones(1,3);
            %Concatenating useful bounds
            obj.LB = [LBr0,LBh,LBL0,LBX,LBY,LBG,LBth,LBdX,LBdY,LBs,LBt];
            obj.UB = [UBr0,UBh,UBL0,UBX,UBY,UBG,UBth,UBdX,UBdY,UBs,UBt];
            obj.LB = obj.LB(find(idx));
            obj.UB = obj.UB(find(idx)); 
        end
        function obj = getIndexFit(obj,fitr0,fith,fitL0,fitMis)
            %Index for atmospheric parameters
            idxr0    = zeros(1,obj.atmMod.nLayer) + fitr0;
            idxh     = zeros(1,obj.atmMod.nLayer) + fith;
            idxL0    = zeros(1,obj.atmMod.nLayer) + fitL0;
            obj.mod.idxFitAtm = [idxr0,idxh,idxL0];
            %Index for system parameters
            idxxSrc  = zeros(1,obj.mod.nwfs);
            idxySrc  = zeros(1,obj.mod.nwfs);
            idxG     = zeros(1,obj.mod.nwfs) + fitMis;
            idxth    = zeros(1,obj.mod.nwfs) + fitMis;
            idxdX    = zeros(1,obj.mod.nwfs) + fitMis;
            idxdY    = zeros(1,obj.mod.nwfs) + fitMis;
            idxs     = zeros(1,obj.mod.nwfs) + fitMis;
            idxtrack = zeros(1,3);
            obj.mod.idxFitSys = [idxxSrc,idxySrc,idxG,idxth,idxdX,idxdY,idxs,idxtrack];
            %Defining a fake intial guess
            obj.indexFit = [obj.mod.idxFitAtm,obj.mod.idxFitSys];            
        end
        function obj = doLsq1Step(obj)                                   
            % Performing the fitting
            obj         = obj.getIndexFit(obj.fitr0,obj.fith,obj.fitL0,obj.fitMis);
            obj.mod.atm = obj.profiler2atm(obj.X0);
            obj.mod     = obj.mod.atmToInit();
            obj         = obj.getBounds();
            obj.FUN     = @(x0) obj.mod.covarianceModel(x0) - obj.Css;
           
            [tmpX,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]= lsqnonlin(...
                obj.FUN,obj.X0);%,obj.LB,obj.UB,obj.OPT);
                        
            obj.X     = obj.mod.parameters;
            obj.X(obj.indexFit==1) = tmpX;
            % Merging results into an atmosphere class
            obj.atmFit = obj.profiler2atm(obj.X);
        end
        function obj = doLsq2Steps(obj)            
            %1. Get the mean slopes removal matrix
            P                   = obj.mod.groundLayerFilter();
            %2. Apply to measurements
            Css                 = P*obj.Css*P';
            %3. Redefine the residual error
            obj.mod.twoSteps    = 1;
            FUN1                = @(x0) obj.mod.covarianceModel(x0) - Css;            
            %4. Managing the fitting procedure for the 1st step
            obj                 = obj.getIndexFit(obj.fitr0,obj.fith,obj.fitL0,obj.fitMis);
            idxFit              = obj.indexFit;
            idxGL               = 1:obj.atmMod.nLayer:3*obj.atmMod.nLayer;
            obj.indexFit(idxGL) = 0;
            obj.mod.idxFitAtm   = obj.indexFit(1:3*obj.atmMod.nLayer);
            tmpX0               = obj.X0(obj.indexFit==1);
            obj                 = obj.getBounds();
            %5. Getting the altitude parameters
            X1                  = lsqnonlin(FUN1,tmpX0);%,obj.LB,obj.UB,obj.OPT);
            
            %6. Managing the fitting procedure for the 2nd step
            obj.indexFit        = idxFit;
            idxALT              = [2:obj.atm.nLayer ...
                2+obj.atm.nLayer:2*(obj.atm.nLayer)...
                2+2*obj.atm.nLayer:3*obj.atm.nLayer];
            
            obj.indexFit(idxALT) = 0;
            obj.mod.idxFitAtm    = obj.indexFit(1:3*obj.atm.nLayer);
            tmpX0                = obj.X0(obj.indexFit==1);
            obj                  = obj.getBounds();
            %7. Getting the ground layer parameters
            obj.mod.twoSteps = 0;
            FUN2             = @(x0) obj.mod.covarianceModel(x0) - obj.Css;                       
            X2               = lsqnonlin(FUN2,tmpX0);%,obj.LB,obj.UB,obj.OPT);
            
            %8. Merging results
            obj.X         = obj.mod.parameters;
            obj.X(idxGL)  = X2;
            obj.X(idxALT) = X1;
            
            % Merging results into an atmosphere class
            obj.atmFit = obj.profiler2atm(obj.X);
        end
        
        %% %%%%%%%%%%%
        % Display procedures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function atm = profiler2atm(obj,X)
                     
            % Grabbing the fitted atmosphere
            nl  = obj.atmMod.nLayer;
            pf  = X;
            fr0 = pf(1:nl);
            r0  = sum(fr0)^(-3/5.);
            fr0 = fr0/sum(fr0);
            h   = pf(1+nl:2*nl);
            L0  = pf(1+2*nl:3*nl);
            wS  = zeros(nl,1);
            wD  = zeros(nl,1);
            
            % Defining the atmosphere
            atm = atmosphere(photometry.V0,r0,'fractionnalR0',fr0,'altitude',h,...
                'layeredL0',L0,'windSpeed',wS,'windDirection',wD);
        end        
    end   
    
    methods (Static)
        function obj = profilerDemo
            
            clear all;
            close all;
            
            %Atmosphere parameters
            r0            = 0.15;
            L0            = 30;
            fractionnalR0 = [.75,.25];
            layerAlt      = [0,10].*1e3;
            windSpeed     = [10,25];
            windDirection = [0.,0.];
            %System parameters
            D             = 10;
            nL            = 14;
            d             = D/nL;
            obs           = 0;
            foV           = 90;
            Fe            = 500;
            nPxGs         = 10;
            nRes          = nPxGs*nL;
            
            %% --------------------- CREATING SUB-CLASSES ---------------------------%%
            
            atmMod = atmosphere(photometry.V0,r0,L0,'fractionnalR0',fractionnalR0,'altitude',layerAlt,...
                'windSpeed',windSpeed,'windDirection',windDirection);
            
            tel = telescope(D,'resolution',nRes,'obstructionRatio',obs ,...
                'fieldOfViewInArcsec',foV ,'samplingTime',1/Fe);
            
            %On-axis target
            sci = source();
            
            %Guide strs asterism
            ast = source('asterism',{[3,1*constants.arcmin2radian,pi/3]});
            
            %WFS type, defined only one is enough
            wfs = shackHartmann(nL,nRes,0.75);
            sci = sci.*tel*wfs;
            wfs.INIT;
            +wfs;
            wfs.slopesUnits = calibrateWfsPixelScale(wfs,sci,tel.D,nPxGs,0);
            
            %Creating a command matrix
            bif          = influenceFunction('monotonic',0.35);
            dm           = deformableMirror(nL+1,'modes',bif,'resolution',nRes,'validActuator',wfs.validActuator);
            dmCalib      = calibration(dm,wfs,sci,sci.wavelength/40);
            dmCalib.cond = 30;
            MC           = dmCalib.M;
            ps           = constants.radian2arcsec*wfs.slopesUnits*sci.wavelength/(2*d*wfs.lenslets.nyquistSampling);
            unit         = 1e9/ps/sqrt(dm.nValidActuator);
            
            %% --------------- TURBULENCE PROFILER ---------------------%%
            obj     = turbulenceProfiler(tel,atmMod,wfs,ast,'fitr0',1);
            obj.Xf  = 0.25^(-5/3.).*[50,50]./100;
            obj     = obj.covarianceMatrixGenerator();
            obj.FUN = @(x0) obj.mod.covarianceModel(x0) - obj.Css;
            obj     = obj.doLsq1Step();  
            
            
            %% --------------- TOMOGRAPHIC ERRORS -----------------------%%
            %Spatial covariance model of slopes
            mod     = slopesCovarianceModel(tel,atm,wfs,[ast,sci],'Projector',MC,'ProjectorUnits',unit);
            %covariance matrix of real measurements
            Css     = mod.covarianceModel(obj.Xf);
            %Fitted covariance matrix
            Clearn  = mod.covarianceModel(obj.X);
            %Reconstructor
            R       = mod.getMMMSEreconstructor(Clearn);
            %Wavefront error
            Cee     = mod.getErrorCovarianceMatrix(Css,R);
            wfeTomo = mod.getWaveFrontError();
            obj.mod = mod;

        end
    end
end


%% ------------------------ OLD STUFFS --------------------------------
 function  obj = getInitialGuess(obj)
            %Defining the initial guess
            obj.mod.fr0   = 0.2^(-5/3.).*[70,(30/(obj.atmMod.nLayer-1))...
                .*ones(1,obj.atmMod.nLayer - 1)]./100;
            obj.mod.h     = linspace(0,16,obj.atmMod.nLayer).*1e3;
            obj.mod.L0    = 30.*ones(1,obj.atmMod.nLayer);
            obj.mod.G     = ones(1,obj.mod.nwfs);           
            obj.mod.th    = zeros(1,obj.mod.nwfs);
            obj.mod.dX    = zeros(1,obj.mod.nwfs);
            obj.mod.dY    = zeros(1,obj.mod.nwfs);
            obj.mod.s     = ones(1,obj.mod.nwfs);           
            obj.mod.track = [0.,0.,0.]; 
            obj.xdata     = obj.mod.objToData();
 end