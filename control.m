classdef control < handle
    
    % CONTROL Create a temporal controller object
    %
    % int = control('typeI',2,'controlMatrix',calibDm.M,'gainInt',0.5); creates a
    % single integrator ("typeI") with 2 steps delay and gain 0.5 implementing the recursion:
    % int.u = int.u - int.gainInt*controlMatrix*int.data; The object is called in a loop as follows:
    % int.data = wfs.slopes;
    % dm.coefs(:,1) = int.u;
    %
    %
    
    properties
        tag = 'CONTROL';
        type;                   % controller type: typeI, typeII, POLC, SALQG
        nDelay;                 % n-steps delay
        buffer;                 % buffer for accumulating the commands
        gainInt = 0.5;
        gainPolc = 0.7;
        gainDoubleInt;
        parmsDoubleInt;
        reconstructor;
        controlMatrix;
        interactionMatrix;
        fittingMatrix;
        u;
        samplingTime;
        partialDelayCoefs = [0 0 1];
        % commands listener
        commandsListener;
        %dataListener
        dataListener;
    end
    properties (Dependent)
        uInit;
    end
    
    properties(SetAccess=public, SetObservable=true)
        data;
    end
    
    properties (Access=private)
        p_u;
        p_uk;
        p_ukm1;
        p_uInit;
        
        C0;
        C1;
        D1;
        o1;
        o2;
        o3;
        log;
    end
    
    methods
        %% Constructor
        function obj = control(type,nDelay,varargin)
            p = inputParser;
            addRequired(p,'type',@ischar);
            addRequired(p,'nDelay',@isnumeric);
            addParameter(p,'gainInt',0.5, @isnumeric);
            addParameter(p,'gainPolc',0.7, @isnumeric);
            addParameter(p,'parmsDoubleInt',[0.02,0.0680 0.0279], @isnumeric);
            addParameter(p,'samplingTime',0.002, @isnumeric);
            addParameter(p,'reconstructor',[],@(x) isa(x,'slopesLinearMMSE') );
            addParameter(p,'controlMatrix',[], @isnumeric);
            addParameter(p,'interactionMatrix', [], @isnumeric);
            addParameter(p,'fittingMatrix',[], @isnumeric);
            addParameter(p,'commandInitVal',0, @isnumeric);
            parse(p,type,nDelay,varargin{:})
            
            obj.type                = p.Results.type;
            obj.nDelay              = p.Results.nDelay;
            obj.gainInt             = p.Results.gainInt;
            obj.gainPolc            = p.Results.gainPolc;
            obj.parmsDoubleInt      = p.Results.parmsDoubleInt;
            obj.samplingTime        = p.Results.samplingTime;
            obj.reconstructor       = p.Results.reconstructor;
            obj.controlMatrix       = p.Results.controlMatrix;
            obj.interactionMatrix   = p.Results.interactionMatrix;
            obj.fittingMatrix       = p.Results.fittingMatrix;
            obj.p_uInit             = p.Results.commandInitVal;
            
            
            obj.dataListener = addlistener(obj,'data','PostSet',...
                @(src,evnt) obj.updateController );
            obj.dataListener.Enabled = true;
            
            obj.log = logBook.checkIn(obj);
            
            
            if ~isempty(obj.controlMatrix)
                obj.p_u = zeros(size(obj.controlMatrix,1),1) + obj.p_uInit;
                obj.buffer = zeros(size(obj.controlMatrix,1), 3);
            elseif ~isempty(obj.interactionMatrix)
                obj.p_u = zeros(size(obj.interactionMatrix,2),1);
                obj.buffer = zeros(size(obj.interactionMatrix,2), 3);
            end
            if isempty(obj.fittingMatrix)
                obj.fittingMatrix = eye(size(obj.controlMatrix,1));
            end
            if strcmp(obj.type,'typeII')
                
                obj.gainDoubleInt   = obj.parmsDoubleInt(1);
                a                   = obj.parmsDoubleInt(2);
                Tlead               = obj.parmsDoubleInt(3);
                Ts                  = obj.samplingTime;
                
                obj.C0 = (1+2*Tlead/Ts)/(1+2*Tlead*a/Ts);
                obj.C1 = (1-2*Tlead/Ts)/(1+2*Tlead*a/Ts);
                obj.D1 = (1-2*Tlead*a/Ts)/(1+2*Tlead*a/Ts);
                
                obj.o1 = zeros(size(obj.controlMatrix,1),1);
                obj.o2 = obj.o1;
                obj.o3 = obj.o1;
                
                obj.p_uk = zeros(size(obj.controlMatrix,1),1);
                obj.p_ukm1 = obj.p_uk;
                obj.buffer = zeros(size(obj.fittingMatrix,1), obj.nDelay+1);
            end
        end
        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.log)
                checkOut(obj.log,obj);
            end
        end
        %% SETS AND GETS
        %% Get u
        function out = get.u(obj)
            out = dmCommandOverWfsExposure(obj);
            %out = obj.p_u;%obj.buffer(:,1);
        end
        %% set and get uInit
        function out = get.uInit(obj)
            out = obj.p_uInit;
        end
        function set.uInit(obj,val)
            obj.p_uInit = val;
            obj.p_u = val;
            obj.buffer = repmat(obj.p_u,1,obj.nDelay+1);
        end
        %% SET PARTIAL DELAY COEFFICIENTS
        function set.nDelay(obj, val)
            obj.nDelay = val;
            if mod(obj.nDelay,1) % i.e. is a non-integer
                if 1 < obj.nDelay < 2
                    b = obj.nDelay-1;
                    a = 1-b;
                    c = 0;
                elseif 2 < obj.nDelay < 3
                    a = 0;
                    c = obj.nDelay-2;
                    b = 1-c;
                end
                obj.partialDelayCoefs = [c, b, a];
            elseif obj.nDelay == 1
                obj.partialDelayCoefs = [0, 0, 1];
            elseif obj.nDelay == 2
                obj.partialDelayCoefs = [0, 1, 0];
            elseif obj.nDelay == 3
                obj.partialDelayCoefs = [1, 0, 0];
            end
        end
        %%
        function updateController(obj)
            obj = obj*obj.data;
            iterate(obj);
            fifo(obj);
        end
        %%
        function obj = mtimes(obj,wfs)
            % * Matrix multilplication
            if isa(wfs,'shackHartmann')
                obj.data = wfs.slopes(:);
            else
                obj.data = wfs(:);
            end
            switch obj.type
                case 'typeI'
                    obj.p_uk = obj.controlMatrix*obj.data;
                case 'typeII'
                    obj.p_uk = obj.controlMatrix*obj.data;
                case 'polc'
                    obj.p_uk = obj.reconstructor * ( bsxfun( @minus, obj.data, obj.interactionMatrix*obj.p_u ) );
            end
        end
        %%
        function iterate(obj)
            switch obj.type
                case 'typeI'
                    obj.p_u = obj.p_u - obj.gainInt*obj.fittingMatrix*obj.p_uk;
                case 'typeII'
                    % DOUBLE INTEGRATOR
                    obj.o1 = obj.gainDoubleInt*(obj.C0*obj.p_uk + obj.C1*obj.p_ukm1) - obj.D1*obj.o1;
                    obj.o2 = obj.o1 + obj.o2;
                    obj.o3 = obj.o2 + obj.o3;
                    obj.p_u = obj.fittingMatrix*obj.o3;
                    obj.p_ukm1 = obj.p_uk;
                case 'polc'
                    obj.p_u = (1-obj.gainPolc)*obj.p_u + obj.gainPolc*obj.fittingMatrix*obj.p_uk;
            end
        end
        
        
        %% FIFO with internal buffer of demands
        function fifo(obj)
            if size(obj.buffer,2) > 1
                obj.buffer = [obj.buffer(:,2:end) obj.p_u];
            else
                obj.buffer = obj.p_u;
            end
            %             for kDelay = 1:obj.nDelay
            %                 obj.buffer(:,kDelay) = obj.buffer(:,kDelay+1);
            %             end
            %             obj.buffer(:,obj.nDelay+1) = obj.p_u;
        end
        
        
        %% computation of the linear combination of DM commands seen over the WFS exposure
        function out = dmCommandOverWfsExposure(obj)
            % In this case the loop is simulated with a fractional delay. In the
            % previous example, the AO has a fixed 1-frame delay, which corresponds to
            % a null RTC latency, since if we consider the discrete samples
            % representing the mid-points in the exposure, the commands are computed
            % and applied directly at the end of the WFS exposure, such that at iteration k the measurement
            %
            % s_k =  phi_k - u_k-1
            %
            % The general case is with a delay greater than one, in which case
            %
            % s_k =  phi_k - a*u_k-1 - b*u_k-2 - c*u_k-3
            %
            % with
            % b = delay-1;      a = 1-b;        c = 0;      if 1 < delay < 2
            % a = 0;            c = delay-2;    b = 1-c;    if 2 < delay < 3
            
            
            out = sum(obj.buffer.*obj.partialDelayCoefs,2);
    
        
        end
    end
    
end