classdef newTurbulenceProfiler < handle
%% NEWTURBULENCEPROFILER developed by Y.O.
% This class computes the turbulence profiling with SLODAR-like method. 
% This class is still under development.
%
    
    properties (SetObservable = true)
        %Sub-classes
        covModel;
        wfs;
                
        isDomeSeeing;
        isUnsensed;
        OPT;
        mask;
        isMask4Map;
        boxWidth4mask;
        hmax;
        isMyLMA;
        
        fitOutput;
        fitResnorm;
        fitResidual;
        isFlatL0;
        usedGsIdx;
    end
    
    properties (Dependent)
        % fit parameters and flags (atm)
        nLayer;
        
        h;
        fith;
        LBh;
        UBh;
       
        fr0;
        fitfr0;
        LBfr0;
        UBfr0;
        
        r0;
        fitr0;
        
        L0;
        fitL0;
        LBL0;
        UBL0;
        
        wSpeed;
        fitwSpeed;
        LBwSpeed;
        UBwSpeed;
        
        wDirection;
        fitwDirection;
        LBwDirection;
        UBwDirection;
        
        dpeak;
        fitdpeak;
        LBdpeak;
        UBdpeak;
        
        fr0dome;
        fitfr0dome;
        LBfr0dome;
        UBfr0dome;
        
        r0dome;
        fitr0dome;
        
        L0dome;
        fitL0dome;
        LBL0dome;
        UBL0dome;
        
        wSpeeddome;
        fitwSpeeddome;
        LBwSpeeddome;
        UBwSpeeddome;
        
        wDirectiondome;    % wind direction for dome seeing
        fitwDirectiondome; % fit flag of wind direction for dome seeing
        LBwDirectiondome;  % lower bound of wind direction for dome seeing
        UBwDirectiondome;  % upper bound of wind direction for dome seeing
        
        fr0unsensed;
        fitfr0unsensed;
        LBfr0unsensed;
        UBfr0unsensed;
        
        r0unsensed;
        fitr0unsensed;

        L0unsensed;
        fitL0unsensed;
        LBL0unsensed;
        UBL0unsensed;
        
        noise;
        fitnoise;
        LBnoise;
        UBnoise;
        
        maxIter;        %#iterations max
        tolFun;         %tolerance and diff value on FUN between two steps
        
        isMap;
        isXY;
        isAuto;
        isFFT;
        flagTT;
        flagFocus;
        rmDiagonal;
        twoSteps;
    end
    
    properties (Dependent,SetAccess=private)
        D;
        gs;
        nGs;
        paramAtm;
        fitFlagAtm;
        LBAtm;
        UBAtm;
        nParamAtm;
        nFitAtm;
        fitParam;
        fitLB;
        fitUB;
    end
    
    properties (Access=private)
        % fit parameter (atm)
        p_nLayer;
        p_idxh;
        p_idxfr0;
        p_idxL0;
        p_idxwSpeed;
        p_idxwDirection;
        p_idxdpeak;
        p_idxfr0dome;
        p_idxL0dome;
        p_idxwSpeeddome;
        p_idxwDirectiondome;
        p_idxfr0unsensed;
        p_idxL0unsensed;
        p_idxnoise;

        p_paramAtm;
        p_fitFlagAtm;
        p_LBAtm;
        p_UBAtm;
        
    end
    
    methods
        
        %% %%%%%%%%%%%
        % CONSTRUCTOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        function obj = newTurbulenceProfiler(D,wfs,gs,varargin)            
                            
            %Checking inputs
            inputs = inputParser;
            inputs.addRequired('D', @isnumeric);
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann'));
            inputs.addRequired('gs',@(x) isa(x,'source'));
            
            inputs.addParameter('atm', [], @(x) isa(x,'atmosphere'));
            inputs.addParameter('nLayer', 5, @isnumeric);
            inputs.addParameter('isMap', false, @islogical);
            inputs.addParameter('isFFT', false, @islogical);
            inputs.addParameter('twoSteps', false, @islogical);
            inputs.addParameter('isXY', false, @islogical);
            inputs.addParameter('isAuto', false, @islogical);
            inputs.addParameter('rmDiagonal', false, @islogical);
            inputs.addParameter('flagTT', true, @islogical);
            inputs.addParameter('flagFocus', true, @islogical);
            inputs.addParameter('isMask4Map', false, @islogical);
            inputs.addParameter('boxWidth4mask', 2, @isnumeric);
            inputs.addParameter('resolution', 0, @isnumeric);
            inputs.addParameter('hmax', 0, @isnumeric);
            inputs.addParameter('isMyLMA', false, @islogical);
            inputs.addParameter('isFlatL0', false, @islogical);
            inputs.addParameter('usedGsIdx', [], @isnumeric);

            inputs.addParameter('maxIter',15, @isnumeric );
            inputs.addParameter('tolFun',1e-10, @isnumeric );
            
            inputs.parse(D,wfs,gs,varargin{:});
            
            if isempty(inputs.Results.usedGsIdx)
                obj.usedGsIdx = 1:length(gs);
            else
                obj.usedGsIdx = inputs.Results.usedGsIdx;
            end
            
            %Instantiation
            obj.covModel = slopesCovarianceModel(telescope(D),...
                atmosphere(photometry.V0,0.1),wfs(1),gs(obj.usedGsIdx),...
                'isMap',inputs.Results.isMap,...
                'isFFT',inputs.Results.isFFT,...
                'isXY',inputs.Results.isXY,...
                'flagTT',inputs.Results.flagTT,...
                'flagFocus',inputs.Results.flagFocus,...
                'isAuto',inputs.Results.isAuto,...
                'twoSteps',inputs.Results.twoSteps);
            obj.rmDiagonal = inputs.Results.rmDiagonal;
            obj.isMyLMA = inputs.Results.isMyLMA;
            obj.wfs = wfs;
            obj.isFlatL0 = inputs.Results.isFlatL0;
            obj.isMask4Map = inputs.Results.isMask4Map;
            obj.boxWidth4mask = inputs.Results.boxWidth4mask;

            % if an input is an atmosphere object...
            if ~isempty(inputs.Results.atm)
                atm = inputs.Results.atm;
                obj.nLayer   = atm.nLayer;
                obj.h        = [atm.layer.altitude];
                obj.fr0      = [atm.layer.fractionnalR0].*(atm.r0)^(-5/3.);
                obj.L0       = [atm.layer.layeredL0];
                obj.dpeak    = 0;
                if ~isempty([atm.layer.windSpeed])
                    obj.wSpeed     = [atm.layer.windSpeed];
                    obj.wDirection = [atm.layer.windDirection];
                else
                    obj.wSpeed     = 0;
                    obj.wDirection = 0;
                end
                
            % if an input is the number of layer...
            elseif inputs.Results.nLayer>0
                obj.nLayer        = inputs.Results.nLayer;
                obj.h             = linspace(0,15000,obj.nLayer);
                obj.r0            = 0.2;
                obj.L0            = 30;
                obj.dpeak         = 0;
                obj.wSpeed        = 0;
                obj.wDirection    = 0;
                
            % SLODAR like method
            else
                obj.getHeight(inputs.Results.resolution,inputs.Results.hmax);
            end
            
            % Initialize
            obj.isDomeSeeing     = false;
            obj.isUnsensed       = false;
            obj.p_fitFlagAtm(:)  = false;
            obj.p_LBAtm(:)       = 0;
            obj.p_UBAtm(:)       = Inf;
            obj.UBwDirection     = 2*pi;
            obj.UBwDirectiondome = 2*pi;

            % fitting parameters
            if ~obj.isMyLMA
                obj.OPT           = optimoptions('lsqnonlin','Display','iter');
                %obj.OPT.Algorithm = 'levenberg-marquardt';
                obj.OPT.MaxIter   = inputs.Results.maxIter;
                obj.OPT.TolFun    = inputs.Results.tolFun;
                obj.OPT.TolX      = 1e-10;
            end

        end

        %% Get D
        function out=get.D(obj)
            out = obj.covModel.tel.D;
        end

        %% Get gs
        function out=get.gs(obj)
            out = obj.covModel.gs;
        end
        
        %% Get nGs
        function out=get.nGs(obj)
            out = length(obj.gs);
        end

        %% Get nParamAtm
        function out=get.nParamAtm(obj)
            out = length(obj.p_paramAtm);
        end

        %% Get nFitAtm
        function out=get.nFitAtm(obj)
            out = sum(obj.p_fitFlagAtm);
        end

        %% Get paramAtm
        function out=get.paramAtm(obj)
            out = obj.p_paramAtm;
        end
        
        %% Get fitFlagAtm
        function out=get.fitFlagAtm(obj)
            out = obj.p_fitFlagAtm;
        end

        %% Get LBAtm
        function out=get.LBAtm(obj)
            out = obj.p_LBAtm;
        end
        
        %% Get UBAtm
        function out=get.UBAtm(obj)
            out = obj.p_UBAtm;
        end
        
        %% Get fitInit
        function out=get.fitParam(obj)
            out = obj.p_paramAtm(obj.p_fitFlagAtm);
        end
        
        %% Get fitLB
        function out=get.fitLB(obj)
            out = obj.p_LBAtm(obj.p_fitFlagAtm);
        end
        
        %% Get fitInit
        function out=get.fitUB(obj)
            out = obj.p_UBAtm(obj.p_fitFlagAtm);
        end
        
        %% Get/Set maxIter
        function out = get.maxIter(obj)
            out = obj.OPT.MaxIter;
        end
        function set.maxIter(obj,val)
            obj.OPT.MaxIter = val;
        end
        
        %% Get/Set tolFun
        function out = get.tolFun(obj)
            out = obj.OPT.TolFun;
        end
        function set.tolFun(obj,val)
            obj.OPT.TolFun = val;
        end

        %% Get/Set isMap
        function out = get.isMap(obj)
            out = obj.covModel.isMap;
        end
        function set.isMap(obj,val)
            obj.covModel.isMap = logical(val);
        end
        
        %% Get/Set isXY
        function out = get.isXY(obj)
            out = obj.covModel.isXY;
        end
        function set.isXY(obj,val)
            obj.covModel.isXY = logical(val);
        end

        %% Get/Set isAuto
        function out = get.isAuto(obj)
            out = obj.covModel.isAuto;
        end
        function set.isAuto(obj,val)
            obj.covModel.isAuto = logical(val);
        end
        
        %% Get/Set isMap
        function out = get.isFFT(obj)
            out = obj.covModel.isFFT;
        end
        function set.isFFT(obj,val)
            obj.covModel.isFFT = logical(val);
        end   
        
        %% Get/Set flagTT
        function out = get.flagTT(obj)
            out = obj.covModel.flagTT;
        end
        function set.flagTT(obj,val)
            obj.covModel.flagTT = logical(val);
        end
        
        %% Get/Set flagTT
        function out = get.flagFocus(obj)
            out = obj.covModel.flagFocus;
        end
        function set.flagFocus(obj,val)
            obj.covModel.flagFocus = logical(val);
        end
        
        %% Get/Set rmDiagonal
        function out = get.rmDiagonal(obj)
            out = obj.covModel.rmDiagonal;
        end
        function set.rmDiagonal(obj,val)
            obj.covModel.rmDiagonal = logical(val);
            if ~obj.rmDiagonal
               obj.fitnoise = false;
               obj.noise = 0;
            end
        end
        
        %% Get/Set twoSteps
        function out = get.twoSteps(obj)
            out = obj.covModel.twoSteps;
        end
        function set.twoSteps(obj,val)
            obj.covModel.twoSteps = logical(val);
        end
        
        %% Get/Set nLayer
        function out = get.nLayer(obj)
            out = obj.p_nLayer;
        end
        function set.nLayer(obj,val)
            obj.p_nLayer     = val;
            obj.p_paramAtm   = zeros(obj.nLayer*6+6+2*obj.nGs,1);
            obj.p_fitFlagAtm = true(obj.nLayer*6+6+2*obj.nGs,1);
            obj.p_LBAtm      = zeros(obj.nLayer*6+6+2*obj.nGs,1);
            obj.p_UBAtm      = zeros(obj.nLayer*6+6+2*obj.nGs,1);
            
            % set index for all parameters
            obj.p_idxh              = (1:obj.nLayer);
            obj.p_idxfr0            = (1:obj.nLayer)+obj.nLayer;
            obj.p_idxL0             = (1:obj.nLayer)+obj.nLayer*2;
            obj.p_idxwSpeed         = (1:obj.nLayer)+obj.nLayer*3;
            obj.p_idxwDirection     = (1:obj.nLayer)+obj.nLayer*4;
            obj.p_idxdpeak          = (1:obj.nLayer)+obj.nLayer*5;
            obj.p_idxfr0dome        = 1+obj.nLayer*6;
            obj.p_idxL0dome         = 2+obj.nLayer*6;
            obj.p_idxwSpeeddome     = 3+obj.nLayer*6;
            obj.p_idxwDirectiondome = 4+obj.nLayer*6;
            obj.p_idxfr0unsensed    = 5+obj.nLayer*6;
            obj.p_idxL0unsensed     = 6+obj.nLayer*6;
            obj.p_idxnoise          = (1:2*obj.nGs)+6+obj.nLayer*6;

        end
        
        %% Get/Set for h
        function out = get.h(obj)
            out = obj.p_paramAtm(obj.p_idxh);
        end
        function set.h(obj,val)
            obj.p_paramAtm(obj.p_idxh) = val(:);
        end
        function out = get.fith(obj)
            out = obj.p_fitFlagAtm(obj.p_idxh);
        end
        function set.fith(obj,val)
            obj.p_fitFlagAtm(obj.p_idxh) = logical(val(:));
        end
        function out = get.LBh(obj)
            out = obj.p_LBAtm(obj.p_idxh);
        end
        function set.LBh(obj,val)
            obj.p_LBAtm(obj.p_idxh) = val(:);
        end
        function out = get.UBh(obj)
            out = obj.p_UBAtm(obj.p_idxh);
        end
        function set.UBh(obj,val)
            obj.p_UBAtm(obj.p_idxh) = val(:);
        end
        
        %% Get/Set for fr0
        function out = get.fr0(obj)
            out = obj.p_paramAtm(obj.p_idxfr0);
        end
        function set.fr0(obj,val)
            obj.p_paramAtm(obj.p_idxfr0) = val(:);
        end
        function out = get.fitfr0(obj)
            out = obj.p_fitFlagAtm(obj.p_idxfr0);
        end
        function set.fitfr0(obj,val)
            obj.p_fitFlagAtm(obj.p_idxfr0) = logical(val(:));
        end
        function out = get.LBfr0(obj)
            out = obj.p_LBAtm(obj.p_idxfr0);
        end
        function set.LBfr0(obj,val)
            obj.p_LBAtm(obj.p_idxfr0) = val(:);
        end
        function out = get.UBfr0(obj)
            out = obj.p_UBAtm(obj.p_idxfr0);
        end
        function set.UBfr0(obj,val)
            obj.p_UBAtm(obj.p_idxfr0) = val(:);
        end
        
        %% Get/Set r0
        function out = get.r0(obj)
            out = obj.p_paramAtm(obj.p_idxfr0).^(-3/5);
        end
        function set.r0(obj,val)
            obj.p_paramAtm(obj.p_idxfr0) = (val(:)).^(-5/3);
        end
        function out = get.fitr0(obj)
            out = obj.p_fitFlagAtm(obj.p_idxfr0);
        end
        function set.fitr0(obj,val)
            obj.p_fitFlagAtm(obj.p_idxfr0) = logical(val(:));
        end
        
        %% Get/Set L0
        function out = get.L0(obj)
            out = obj.p_paramAtm(obj.p_idxL0);
        end
        function set.L0(obj,val)
            if ~obj.isFlatL0
                obj.p_paramAtm(obj.p_idxL0) = val(:);
            else
                obj.p_paramAtm(obj.p_idxL0) = val(1);
            end
        end
        function out = get.fitL0(obj)
            if ~obj.isFlatL0
                out = obj.p_fitFlagAtm(obj.p_idxL0);
            else
                out = obj.p_fitFlagAtm(obj.p_idxL0(1));
            end
        end
        function set.fitL0(obj,val)
            if ~obj.isFlatL0
                obj.p_fitFlagAtm(obj.p_idxL0) = logical(val(:));
            else
                obj.p_fitFlagAtm(obj.p_idxL0(1)) = logical(val(1));
                obj.p_fitFlagAtm(obj.p_idxL0(2:end)) = false;
            end
        end
        function out = get.LBL0(obj)
            out = obj.p_LBAtm(obj.p_idxL0);
        end
        function set.LBL0(obj,val)
            obj.p_LBAtm(obj.p_idxL0) = val(:);
        end
        function out = get.UBL0(obj)
            out = obj.p_UBAtm(obj.p_idxL0);
        end
        function set.UBL0(obj,val)
            obj.p_UBAtm(obj.p_idxL0) = val(:);
        end
        
        %% Get/Set wSpeed
        function out = get.wSpeed(obj)
            out = obj.p_paramAtm(obj.p_idxwSpeed);
        end
        function set.wSpeed(obj,val)
            obj.p_paramAtm(obj.p_idxwSpeed) = val(:);
        end
        function out = get.fitwSpeed(obj)
            out = obj.p_fitFlagAtm(obj.p_idxwSpeed);
        end
        function set.fitwSpeed(obj,val)
            obj.p_fitFlagAtm(obj.p_idxwSpeed) = logical(val(:));
        end
        function out = get.LBwSpeed(obj)
            out = obj.p_LBAtm(obj.p_idxwSpeed);
        end
        function set.LBwSpeed(obj,val)
            obj.p_LBAtm(obj.p_idxwSpeed) = val(:);
        end
        function out = get.UBwSpeed(obj)
            out = obj.p_UBAtm(obj.p_idxwSpeed);
        end
        function set.UBwSpeed(obj,val)
            obj.p_UBAtm(obj.p_idxwSpeed) = val(:);
        end
        
        %% Get/Set wDirection
        function out = get.wDirection(obj)
            out = obj.p_paramAtm(obj.p_idxwDirection);
        end
        function set.wDirection(obj,val)
            obj.p_paramAtm(obj.p_idxwDirection) = val(:);
        end
        function out = get.fitwDirection(obj)
            out = obj.p_fitFlagAtm(obj.p_idxwDirection);
        end
        function set.fitwDirection(obj,val)
            obj.p_fitFlagAtm(obj.p_idxwDirection) = logical(val(:));
        end
        function out = get.LBwDirection(obj)
            out = obj.p_LBAtm(obj.p_idxwDirection);
        end
        function set.LBwDirection(obj,val)
            obj.p_LBAtm(obj.p_idxwDirection) = val(:);
        end
        function out = get.UBwDirection(obj)
            out = obj.p_UBAtm(obj.p_idxwDirection);
        end
        function set.UBwDirection(obj,val)
            obj.p_UBAtm(obj.p_idxwDirection) = val(:);
        end
        
        %% Get/Set dpeak
        function out = get.dpeak(obj)
            out = obj.p_paramAtm(obj.p_idxdpeak);
        end
        function set.dpeak(obj,val)
            obj.p_paramAtm(obj.p_idxdpeak) = val(:);
        end
        function out = get.fitdpeak(obj)
            out = obj.p_fitFlagAtm(obj.p_idxdpeak);
        end
        function set.fitdpeak(obj,val)
            obj.p_fitFlagAtm(obj.p_idxdpeak) = logical(val(:));
        end
        function out = get.LBdpeak(obj)
            out = obj.p_LBAtm(obj.p_idxdpeak);
        end
        function set.LBdpeak(obj,val)
            obj.p_LBAtm(obj.p_idxdpeak) = val(:);
        end
        function out = get.UBdpeak(obj)
            out = obj.p_UBAtm(obj.p_idxdpeak);
        end
        function set.UBdpeak(obj,val)
            obj.p_UBAtm(obj.p_idxdpeak) = val(:);
        end
        
        %% Get/Set fr0dome
        function out = get.fr0dome(obj)
            out = obj.p_paramAtm(obj.p_idxfr0dome);
        end
        function set.fr0dome(obj,val)
            obj.p_paramAtm(obj.p_idxfr0dome) = val(:);
        end
        function out = get.fitfr0dome(obj)
            out = obj.p_fitFlagAtm(obj.p_idxfr0dome);
        end
        function set.fitfr0dome(obj,val)
            obj.p_fitFlagAtm(obj.p_idxfr0dome) = logical(val(:));
        end
        function out = get.LBfr0dome(obj)
            out = obj.p_LBAtm(obj.p_idxfr0dome);
        end
        function set.LBfr0dome(obj,val)
            obj.p_LBAtm(obj.p_idxfr0dome) = val(:);
        end
        function out = get.UBfr0dome(obj)
            out = obj.p_UBAtm(obj.p_idxfr0dome);
        end
        function set.UBfr0dome(obj,val)
            obj.p_UBAtm(obj.p_idxfr0dome) = val(:);
        end
        
        %% Get/Set r0dome
        function out = get.r0dome(obj)
            out = obj.p_paramAtm(obj.p_idxfr0dome).^(-3/5);
        end
        function set.r0dome(obj,val)
            obj.p_paramAtm(obj.p_idxfr0dome) = (val(:)).^(-5/3);
        end
        function out = get.fitr0dome(obj)
            out = obj.p_fitFlagAtm(obj.p_idxfr0dome);
        end
        function set.fitr0dome(obj,val)
            obj.p_fitFlagAtm(obj.p_idxfr0dome) = logical(val(:));
        end
        
        %% Get/Set L0dome
        function out = get.L0dome(obj)
            out = obj.p_paramAtm(obj.p_idxL0dome);
        end
        function set.L0dome(obj,val)
            obj.p_paramAtm(obj.p_idxL0dome) = val(:);
        end
        function out = get.fitL0dome(obj)
            out = obj.p_fitFlagAtm(obj.p_idxL0dome);
        end
        function set.fitL0dome(obj,val)
            obj.p_fitFlagAtm(obj.p_idxL0dome) = logical(val(:));
        end
        function out = get.LBL0dome(obj)
            out = obj.p_LBAtm(obj.p_idxL0dome);
        end
        function set.LBL0dome(obj,val)
            obj.p_LBAtm(obj.p_idxL0dome) = val(:);
        end
        function out = get.UBL0dome(obj)
            out = obj.p_UBAtm(obj.p_idxL0dome);
        end
        function set.UBL0dome(obj,val)
            obj.p_UBAtm(obj.p_idxL0dome) = val(:);
        end
        
        %% Get/Set wSpeeddome
        function out = get.wSpeeddome(obj)
            out = obj.p_paramAtm(obj.p_idxwSpeeddome);
        end
        function set.wSpeeddome(obj,val)
            obj.p_paramAtm(obj.p_idxwSpeeddome) = val(:);
        end
        function out = get.fitwSpeeddome(obj)
            out = obj.p_fitFlagAtm(obj.p_idxwSpeeddome);
        end
        function set.fitwSpeeddome(obj,val)
            obj.p_fitFlagAtm(obj.p_idxwSpeeddome) = logical(val(:));
        end
        function out = get.LBwSpeeddome(obj)
            out = obj.p_LBAtm(obj.p_idxwSpeeddome);
        end
        function set.LBwSpeeddome(obj,val)
            obj.p_LBAtm(obj.p_idxwSpeeddome) = val(:);
        end
        function out = get.UBwSpeeddome(obj)
            out = obj.p_UBAtm(obj.p_idxwSpeeddome);
        end
        function set.UBwSpeeddome(obj,val)
            obj.p_UBAtm(obj.p_idxwSpeeddome) = val(:);
        end
        
        %% Get/Set wDirection
        function out = get.wDirectiondome(obj)
            out = obj.p_paramAtm(obj.p_idxwDirectiondome);
        end
        function set.wDirectiondome(obj,val)
            obj.p_paramAtm(obj.p_idxwDirectiondome) = val(:);
        end
        function out = get.fitwDirectiondome(obj)
            out = obj.p_fitFlagAtm(obj.p_idxwDirectiondome);
        end
        function set.fitwDirectiondome(obj,val)
            obj.p_fitFlagAtm(obj.p_idxwDirectiondome) = logical(val(:));
        end
        function out = get.LBwDirectiondome(obj)
            out = obj.p_LBAtm(obj.p_idxwDirectiondome);
        end
        function set.LBwDirectiondome(obj,val)
            obj.p_LBAtm(obj.p_idxwDirectiondome) = val(:);
        end
        function out = get.UBwDirectiondome(obj)
            out = obj.p_UBAtm(obj.p_idxwDirectiondome);
        end
        function set.UBwDirectiondome(obj,val)
            obj.p_UBAtm(obj.p_idxwDirectiondome) = val(:);
        end
        
        %% Get/Set fr0unsensed
        function out = get.fr0unsensed(obj)
            out = obj.p_paramAtm(obj.p_idxfr0unsensed);
        end
        function set.fr0unsensed(obj,val)
            obj.p_paramAtm(obj.p_idxfr0unsensed) = val(:);
        end
        function out = get.fitfr0unsensed(obj)
            out = obj.p_fitFlagAtm(obj.p_idxfr0unsensed);
        end
        function set.fitfr0unsensed(obj,val)
            obj.p_fitFlagAtm(obj.p_idxfr0unsensed) = logical(val(:));
        end
        function out = get.LBfr0unsensed(obj)
            out = obj.p_LBAtm(obj.p_idxfr0unsensed);
        end
        function set.LBfr0unsensed(obj,val)
            obj.p_LBAtm(obj.p_idxfr0unsensed) = val(:);
        end
        function out = get.UBfr0unsensed(obj)
            out = obj.p_UBAtm(obj.p_idxfr0unsensed);
        end
        function set.UBfr0unsensed(obj,val)
            obj.p_UBAtm(obj.p_idxfr0unsensed) = val(:);
        end
        
        %% Get/Set r0unsensed
        function out = get.r0unsensed(obj)
            out = obj.p_paramAtm(obj.p_idxfr0unsensed).^(-3/5);
        end
        function set.r0unsensed(obj,val)
            obj.p_paramAtm(obj.p_idxfr0unsensed) = (val(:)).^(-5/3);
        end
        function out = get.fitr0unsensed(obj)
            out = obj.p_fitFlagAtm(obj.p_idxfr0unsensed);
        end
        function set.fitr0unsensed(obj,val)
            obj.p_fitFlagAtm(obj.p_idxfr0unsensed) = logical(val(:));
        end
        
        %% Get/Set L0unsensed
        function out = get.L0unsensed(obj)
            out = obj.p_paramAtm(obj.p_idxL0unsensed);
        end
        function set.L0unsensed(obj,val)
            obj.p_paramAtm(obj.p_idxL0unsensed) = val(:);
        end
        function out = get.fitL0unsensed(obj)
            out = obj.p_fitFlagAtm(obj.p_idxL0unsensed);
        end
        function set.fitL0unsensed(obj,val)
            obj.p_fitFlagAtm(obj.p_idxL0unsensed) = logical(val(:));
        end
        function out = get.LBL0unsensed(obj)
            out = obj.p_LBAtm(obj.p_idxL0unsensed);
        end
        function set.LBL0unsensed(obj,val)
            obj.p_LBAtm(obj.p_idxL0unsensed) = val(:);
        end
        function out = get.UBL0unsensed(obj)
            out = obj.p_UBAtm(obj.p_idxL0unsensed);
        end
        function set.UBL0unsensed(obj,val)
            obj.p_UBAtm(obj.p_idxL0unsensed) = val(:);
        end

        %% Get/Set noise
        function out = get.noise(obj)
            out = obj.p_paramAtm(obj.p_idxnoise);
        end
        function set.noise(obj,val)
            obj.p_paramAtm(obj.p_idxnoise) = val(:);
        end
        function out = get.fitnoise(obj)
            out = obj.p_fitFlagAtm(obj.p_idxnoise);
        end
        function set.fitnoise(obj,val)
            obj.p_fitFlagAtm(obj.p_idxnoise) = logical(val(:));
        end
        function out = get.LBnoise(obj)
            out = obj.p_LBAtm(obj.p_idxnoise);
        end
        function set.LBnoise(obj,val)
            obj.p_LBAtm(obj.p_idxnoise) = val(:);
        end
        function out = get.UBnoise(obj)
            out = obj.p_UBAtm(obj.p_idxL0unsensed);
        end
        function set.UBnoise(obj,val)
            obj.p_UBAtm(obj.p_idxL0unsensed) = val(:);
        end
        
        %% %%%%%%%%%%%
        %Fitting procedure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fitX=testFitting(obj, atm, varargin)
            %% simple fitting
            
            %Checking inputs
            inputs = inputParser;
            inputs.addParameter('init', obj.fitParam, @isnumeric);
            inputs.addParameter('LB', obj.fitLB, @isnumeric);
            inputs.addParameter('UB', obj.fitUB, @isnumeric);
            inputs.parse(varargin{:});
            
            l_init = inputs.Results.init;
            l_LB = inputs.Results.LB;
            l_UB = inputs.Results.UB;
            
            obj.createMask;
            
            if ~obj.twoSteps % one Step
                mcov = obj.computeTestCov(atm);
                
                obj.p_paramAtm(obj.p_fitFlagAtm) = l_init;
                obj.p_LBAtm(obj.p_fitFlagAtm)    = l_LB;
                obj.p_UBAtm(obj.p_fitFlagAtm)    = l_UB;
                
                [fitX,~,~,~,obj.fitOutput] = lsqnonlin(@(x) obj.computeModCov(x)-mcov,...
                    obj.fitParam,...
                    obj.fitLB,...
                    obj.fitUB,...
                    obj.OPT);
                
            else % two Steps
                
                % first step
                % save init for GL
                glfr0 = obj.fr0(1);
                glL0  = obj.L0(1);
                
                % reset GL
                obj.fr0(1)    = 0;
                obj.fith(1)   = false;
                obj.fitfr0(1) = false;
                obj.fitL0(1)  = false;
                
                % S is slope in arcsec
                mcov = obj.computeTestCov(atm);
                
                lsqnonlin(@(x) obj.computeModCov(x)-mcov,...
                    obj.fitParam,...
                    obj.fitLB,...
                    obj.fitUB,...
                    obj.OPT);
                
                % second step
                obj.twoSteps = false;
                
                % reset flag and init
                obj.fith      = false;
                obj.h(1)      = 0;
                obj.fitfr0    = false;
                obj.fitfr0(1) = true;
                obj.fr0(1)    = glfr0;
                obj.fitL0     = false;
                obj.fitL0(1)  = true;
                obj.L0(1)     = glL0;
                
                mcov = obj.computeTestCov(atm);
                
                [fitX,~,~,~,obj.fitOutput] = lsqnonlin(@(x) obj.computeModCov(x)-mcov,...
                    obj.fitParam,...
                    obj.fitLB,...
                    obj.fitUB,...
                    obj.OPT);
                
                obj.fith = true;
                obj.fitr0 = true;
                obj.fitL0 = true;
                obj.twoSteps = true;
            end
        end
        
        function fitX=fitting(obj,S, varargin)
            %% simple fitting
            
            %Checking inputs
            inputs = inputParser;
            inputs.addParameter('init', obj.fitParam, @isnumeric);
            inputs.addParameter('LB', obj.fitLB, @isnumeric);
            inputs.addParameter('UB', obj.fitUB, @isnumeric);
            inputs.parse(varargin{:});
            
            l_init = inputs.Results.init;
            l_LB = inputs.Results.LB;
            l_UB = inputs.Results.UB;
            
            obj.p_paramAtm(obj.p_fitFlagAtm) = l_init;
            obj.p_LBAtm(obj.p_fitFlagAtm)    = l_LB;
            obj.p_UBAtm(obj.p_fitFlagAtm)    = l_UB;
            
            obj.createMask;
            
            if ~obj.twoSteps % one Step
                if ~isvector(S)
                    mcov = obj.computeObsCov(S);
                else
                    mcov = S(:)';
                end
                if ~obj.isMyLMA
                    [fitX,obj.fitResnorm,obj.fitResidual,~,obj.fitOutput] = lsqnonlin(@(x) obj.computeModCov(x)-mcov,...
                        obj.fitParam,...
                        obj.fitLB,...
                        obj.fitUB,...
                        obj.OPT);
                else
                    fitX = myLevenbergMarquardt(@(x) obj.computeModCov(x),...
                        obj.fitParam,...
                        mcov,...
                        0,...
                        'nMAX',obj.maxIter);
                end
                
            else % two Steps
                
                % first step
                % save init for GL
                glfr0 = obj.fr0(1);
                glL0  = obj.L0(1);
                
                % reset GL
                obj.fr0(1)    = 0;
                obj.fith(1)   = false;
                obj.fitfr0(1) = false;
                obj.fitL0(1)  = false;
                
                % S is slope in arcsec
                P = obj.covModel.groundLayerFilter;
                mcov = obj.computeObsCov(P*S);
                
                lsqnonlin(@(x) obj.computeModCov(x)-mcov,...
                    obj.fitParam,...
                    obj.fitLB,...
                    obj.fitUB,...
                    obj.OPT);

                % second step
                obj.twoSteps = false;
                
                % reset flag and init
                obj.fith      = false;
                obj.h(1)      = 0;
                obj.fitfr0    = false;
                obj.fitfr0(1) = true;
                obj.fr0(1)    = glfr0;
                obj.fitL0     = false;
                obj.fitL0(1)  = true;                
                obj.L0(1)     = glL0;
                
                mcov = obj.computeObsCov(S);
                
                fitX = lsqnonlin(@(x) obj.computeModCov(x)-mcov,...
                    obj.fitParam,...
                    obj.fitLB,...
                    obj.fitUB,...
                    obj.OPT);
                
                obj.fith = true;
                obj.fitr0 = true;
                obj.fitL0 = true;
                obj.twoSteps = true;
            end
        end
        
        function cov = computeObsCov(obj,S)
            %% compute observed covariance from slope
            
            nf = size(S,2);
            nSub = obj.wfs(1).lenslets.nLenslet;
            nm = size(S,1)/(2*obj.nGs);
            mask = obj.wfs(1).validLenslet;

            % base mask
            if obj.isMap
                n1 = 2;
                if obj.isXY
                    n1 = 3;
                end
                n2 = nchoosek(obj.nGs,2)+1;
            end
            
            % mode fitrering
            if obj.flagTT
                mode = [2 3];
                if obj.flagFocus
                    mode = [2 3 4];
                end
                zer = zernike(1:max(mode),...
                    'resolution',nSub,...
                    'pupil',double(mask));
                S = mat2cell(S,ones(2*obj.nGs,1)*nm,nf);
                
                z = sum(zer.xDerivative(mask,mode),2);
                for k=1:2:obj.nGs*2
                    S{k} = S{k} - z*inv(z'*z)*z'*S{k};
                end
                
                z = sum(zer.yDerivative(mask,mode),2);
                for k=2:2:obj.nGs*2
                    S{k} = S{k} - z*inv(z'*z)*z'*S{k};
                end
                
                S = cell2mat(S);
            end
            
            if ~obj.isMap
                cov = S*S'/nf;
            else
                S = mat2cell(S,ones(2*obj.nGs,1)*nm,nf);
                cov = cell(1,n2);

                % auto-correlation
                if obj.isAuto
                    cov{1} = zeros(n1*(nSub*2-1),(nSub*2-1));
                    for iGs = 1:obj.nGs
                        tmp = [tools.covMatrix2Map(S{2*iGs-1}*S{2*iGs-1}'/nf,nSub,nSub,'mask1',mask,'mask2',mask);...
                            tools.covMatrix2Map(S{2*iGs}*S{2*iGs}'/nf,nSub,nSub,'mask1',mask,'mask2',mask);];
                        if obj.isXY
                            tmp = [tmp;...
                                tools.covMatrix2Map(S{2*iGs-1}*S{2*jGs}'/nf,nSub,nSub,'mask1',mask,'mask2',mask);];
                        end
                        cov{1} = cov{1} + tmp;
                    end
                    cov{1} = cov{1}/obj.nGs;
                end
                
                % cross-correlation
                t = obj.isAuto;
                for iGs = 1:obj.nGs
                   for jGs = iGs+1:obj.nGs
                       t = t+1;
                       cov{t} = [tools.covMatrix2Map(S{2*iGs-1}*S{2*jGs-1}'/nf,nSub,nSub,'mask1',mask,'mask2',mask);...
                           tools.covMatrix2Map(S{2*iGs}*S{2*jGs}'/nf,nSub,nSub,'mask1',mask,'mask2',mask);];
                       if obj.isXY
                           cov{t} = [cov{t};...
                               tools.covMatrix2Map(S{2*iGs-1}*S{2*jGs}'/nf,nSub,nSub,'mask1',mask,'mask2',mask);];
                       end
                   end
                end
                cov = cell2mat(cov);
            end
            cov(~obj.mask) = [];
        end
        
        function cov = computeModCov(obj,x)
            %% compute model covariance from inputs
            if ~isempty(x)
                obj.p_paramAtm(obj.p_fitFlagAtm) = x;
            end
            if obj.isFlatL0
                   obj.L0 = obj.L0(1);
            end
            obj.covModel.h        = obj.h;
            obj.covModel.fr0      = obj.fr0;
            obj.covModel.L0       = obj.L0;
            obj.covModel.noiseVar = obj.noise;
            if ~obj.isMap
                cov = obj.covModel.getCovarianceMatrix;
            else
                cov = obj.covModel.getCovarianceMap;
            end
            cov(~obj.mask) = [];
            
            if false
                figure(101);
                subplot(1,3,1);
                imagesc(cov); axis equal tight;
                
                w = bsxfun(@(x,y) abs(x-y),obj.h,obj.h');
                w = min(w(w>0));
                
                subplot(1,3,2);
                barh(obj.h/1000,obj.fr0,500/w,'r');
                xlabel('fr0');
                ylabel('Altitude [km]');
                
                subplot(1,3,3);
                barh(obj.h/1000,obj.L0,500/w,'r');
                xlabel('L0');
                ylabel('Altitude [km]');
                
                drawnow;
            end
            
        end
        
        function cov = computeTestCov(obj,atm)
            %% compute model covariance from inputs
            inh  = zeros(atm.nLayer,1);
            infr0 = zeros(atm.nLayer,1);
            inL0 = zeros(atm.nLayer,1);
            for kLayer=1:atm.nLayer
                infr0(kLayer) = atm.r0^(-5/3)*atm.layer(kLayer).fractionnalR0;
                inL0(kLayer) = atm.layer(kLayer).layeredL0;
                inh(kLayer) = atm.layer(kLayer).altitude;
            end
            obj.covModel.h        = inh;
            obj.covModel.fr0      = infr0;
            obj.covModel.L0       = inL0;
            obj.covModel.noiseVar = 0;
            if ~obj.isMap
                cov = obj.covModel.getCovarianceMatrix;
            else
                cov = obj.covModel.getCovarianceMap;
            end
            cov(~obj.mask) = 0;
        end
        
        function createMask(obj)
            
            % base mask
            if ~obj.isMap
                ns = sum(obj.wfs(1).validLenslet(:));
                obj.mask = true(2*ns*obj.nGs);
            else
                nsub = obj.wfs(1).lenslets.nLenslet;
                n = 2*nsub-1;
                n1 = 2;
                if obj.isXY
                    n1 = 3;
                end
                n2 = nchoosek(obj.nGs,2)+obj.isAuto;
                obj.mask = true(n*n1,n*n2);
            end
            if obj.isMap
                m = obj.wfs(1).validLenslet;
                pmask = fftshift(ifft2(fft2(m,n,n).*conj(fft2(m,n,n))));
                pmask(pmask<.5)=0;
                pmask=logical(pmask);
                obj.mask = logical(obj.mask.*repmat(pmask,n1,n2));
            end
            
            % diagonal or central peak
            if obj.rmDiagonal
                if ~obj.isMap
                    obj.mask = logical(obj.mask - eye(size(obj.mask)));
                else
                    if obj.isAuto
                        tmp = mat2cell(obj.mask,ones(n1,1)*n,ones(n2,1)*n);
                        for t=1:n2
                            tmp{t,1}(nsub,nsub) = false;
                        end
                        obj.mask = cell2mat(tmp);
                    end
                end
            end
            
            % additional mask for map
            if obj.isMap && obj.isMask4Map
                tmp = mat2cell(obj.mask,ones(n1,1)*n,ones(n2,1)*n);
                
                % auto-correlation
                if obj.isAuto
                    for t=1:n1
                        w = 1+obj.boxWidth4mask;
                        tmp{t,1}(1:nsub-w,1:nsub-w) = false;
                        tmp{t,1}(1:nsub-w,nsub+w:end) = false;
                        tmp{t,1}(nsub+w:end,1:nsub-w) = false;
                        tmp{t,1}(nsub+w:end,nsub+w:end) = false;
                    end
                end
                
                k = obj.isAuto;
                for iGs = 1:obj.nGs
                    for jGs = iGs+1:obj.nGs
                        k=k+1;
                        vi = obj.gs(iGs).directionVector(1:2);
                        vj = obj.gs(jGs).directionVector(1:2);
                        ang = atan2(vi(2)-vj(2),vi(1)-vj(1));
                        tmap = false(n);
                        if -pi/4<=ang && ang<=pi/4
                            x=0:nsub-1;
                            y=tan(ang)*x;
                        elseif pi/4<ang && ang<=3/4*pi
                            y=0:nsub-1;
                            x=-tan(ang-pi/2)*y;
                        elseif -3/4*pi>=ang && ang>-pi/4
                            y=-(0:nsub-1);
                            x=-tan(ang+pi/2)*y;
                        else
                            x=-(0:nsub-1);
                            y=-tan(ang)*x;
                        end
                        x=round(x)+nsub;
                        y=round(y)+nsub;
                        
                        for p=1:nsub
                            xr = (x(p)-obj.boxWidth4mask):(x(p)+obj.boxWidth4mask);
                            yr = (y(p)-obj.boxWidth4mask):(y(p)+obj.boxWidth4mask);
                            
                            xr(xr<1)=[];
                            xr(xr>n)=[];
                            yr(yr<1)=[];
                            yr(yr>n)=[];
                            
                            if ~isempty(xr) && ~ isempty(yr)
                                tmap(yr,xr) = true;
                            end
                        end
                        for t=1:n1
                             tmp{t,k} = tmap;
                        end
                    end
                end
                obj.mask = cell2mat(tmp);
            end
        end
        
        function getHeight(obj,hresolution,hlimit)

            pairs = nchoosek(1:obj.nGs,2);
            npair = size(pairs,1);
            
            hmax = zeros(npair,1);
            dv = zeros(npair,1);
            nlayer = zeros(npair,1);
            nsub = obj.wfs(1).lenslets.nLenslet;
            dsub = obj.D/nsub;
            
            n = 2*nsub-1;
            m = obj.wfs(1).validLenslet;
            pmask = fftshift(ifft2(fft2(m,n,n).*conj(fft2(m,n,n))));
            pmask(pmask<1)=0;
            pmask=logical(pmask);
            
            for p = 1:npair
                                    
                vi = obj.gs(pairs(p,1)).directionVector(1:2);
                vj = obj.gs(pairs(p,2)).directionVector(1:2);
                dv(p) = sqrt((vi(1)-vj(1))^2+(vi(2)-vj(2))^2);
                ang = atan2(vi(2)-vj(2),vi(1)-vj(1));
                %if (-pi/4<=ang && ang<=pi/4) || (-3/4*pi>=ang && ang>-pi/4)
                if (-pi/4<=ang && ang<=pi/4) || (ang<-3*pi/4 || ang>3*pi/4)
                    rdsub = dsub*abs(cos(ang));
                else
                    rdsub = dsub*abs(sin(ang));
                end
                

                if hlimit==0
                    if -pi/4<=ang && ang<=pi/4
                        x=0:nsub-1;
                        y=tan(ang)*x;
                    elseif pi/4<ang && ang<=3/4*pi
                        y=0:nsub-1;
                        x=-tan(ang-pi/2)*y;
                    elseif -3/4*pi>=ang && ang>-pi/4
                        y=-(0:nsub-1);
                        x=-tan(ang+pi/2)*y;
                    else
                        x=-(0:nsub-1);
                        y=-tan(ang)*x;
                    end
                    x=round(x)+nsub;
                    y=round(y)+nsub;
                    for t=1:nsub
                        if ~isinf(obj.gs(1).height)
                            h = (t-1)*rdsub*obj.gs(1).height./(obj.gs(1).height*dv(p)+(t-1)*rdsub);
                        else
                            h = (t-1)*rdsub/dv(p);
                        end
                        if x(t)<1 && y(t)<1 || x(t)>n || y(t)>n
                            hmax(p) = h;
                            break;
                        elseif ~pmask(y(t),x(t))
                            hmax(p) = h;
                            break;
                        elseif hlimit<h
                            hmax(p) = h;
                        end
                    end
                end
            end
            obj.hmax = max(hmax);
            [~,ip] = max(dv);
            
            vi = obj.gs(pairs(ip,1)).directionVector(1:2);
            vj = obj.gs(pairs(ip,2)).directionVector(1:2);
            dv = sqrt((vi(1)-vj(1))^2+(vi(2)-vj(2))^2);
            ang = atan2(vi(2)-vj(2),vi(1)-vj(1));
            %if (-pi/4<=ang && ang<=pi/4) || (-3/4*pi>=ang && ang>-pi/4)
            if (-pi/4<=ang && ang<=pi/4) || (ang<-3*pi/4 || ang>3*pi/4)
                rdsub = dsub*abs(cos(ang));
            else
                rdsub = dsub*abs(sin(ang));
            end
            
            t = 0;
            H = [];
            while 1
                if ~isinf(obj.gs(1).height)
                    h = t*rdsub*obj.gs(1).height./(obj.gs(1).height*dv+t*rdsub);
                else
                    h = t*rdsub/dv;
                end

                if h < obj.hmax
                    t = t + 1;
                    H(t,1) = h;
                else
                    break;
                end
            end
            obj.nLayer = t;
            obj.h = H;
        end
        
        function plotAns(obj)
            w = bsxfun(@(x,y) abs(x-y),obj.h,obj.h');
            w = min(w(w>0));
            
            figure;
            subplot(1,2,1);
            barh(obj.h/1000,obj.fr0,500/w,'r');
            xlabel('fr0');
            ylabel('Altitude [km]');
            
            subplot(1,2,2);
            barh(obj.h/1000,obj.L0,500/w,'r');
            xlabel('L0');
            ylabel('Altitude [km]');
        end
    end
end

