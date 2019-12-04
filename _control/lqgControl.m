classdef lqgControl < handle
    %{
     lqgControl Creates an LQG temporal controller object
     
     obj = lqgControl(samplingTime, nDelay, 'atmParms',[pole1, pole2 varAtm],'vibParms',[f0, deltaF, varVib])
     creates an lqgControl object with samplingTime temporal sampling,
     nDelay the number of integer multiple of samplingTime with
     'atmParms'=[pole1, pole2 varAtm] the poles of a continuous-time
     transfer function that approximates the temporal spectrum of the
     atmoospheric disturbance and 'vibParms' = [f0, deltaF, varVib] the
     central frequency, the FWHM and the total variance of the vibration
     peak.
    %}
    properties
        tag = 'lqgControl';
        samplingTime;           % discrete time-step
        vibParms;               % vibration parameters
        atmParms;
        
        continuousModel;
        discreteModel;
        
        nDelay;                 % n-steps delay
        noiseVar;               % measurement noise variance
        
        KalmanFilter;           % Kalman filter matrices after synthesis
        AR;                     % whether this is an AR model based
    end
    properties (Dependent)
        nVibPeaks;              % number of vibration peaks
    end
    
    properties(SetAccess=public, SetObservable=true)
    end
    
    properties (Access=private)
        log;
    end
    
    methods
        %% Constructor
        function obj = lqgControl(samplingTime,nDelay,varargin)
            
            charchk = {'char'};
            numchk = {'numeric'};
            nempty = {'nonempty'};
            p = inputParser;
            addRequired(p,'samplingTime',@isnumeric);
            addRequired(p,'nDelay',@isnumeric);
            addOptional(p,'atmParms',[1 1 1], @(x)validateattributes(x,numchk,nempty))
            addOptional(p,'vibParms',[1 1 1],@isnumeric);
            addOptional(p,'noiseVar',1,@isnumeric);
            addOptional(p,'AR',0,@isnumeric);
            parse(p,samplingTime,nDelay,varargin{:})
                        
            obj.atmParms = p.Results.atmParms;
            
            obj.samplingTime        = p.Results.samplingTime;
            obj.nDelay              = p.Results.nDelay;
            obj.vibParms            = p.Results.vibParms;
            obj.noiseVar            = p.Results.noiseVar;
            obj.log                 = logBook.checkIn(obj);
            obj.AR                  = p.Results.AR;
            
            % continuous-time model
            continuousTimeModel(obj);
            
            % discrete-time model
            discreteTimeModel(obj);
            
            %model stack
            modelStack(obj);
            
            %Kalman filter synthesis
            synthesis(obj);
            
        end
        %% Destructor
        function delete(obj)
            if ~isempty(obj.log)
                checkOut(obj.log,obj);
            end
        end
        
        %% SETS AND GETS
        %% Get and Set nVibPecks
        function out = get.nVibPeaks(obj)
            if ~isempty(obj.vibParms)
                out = size(obj.vibParms,1);
            else
                out = 0;
            end
            
        end
        
        %% CONTINUOUS-TIME MODELS
        %-------------------------------------------------
        %     CONTINUOUS TEMPORAL MODEL
        %-------------------------------------------------
        function continuousTimeModel(obj)
            T = obj.samplingTime;
            Ac = sparse(obj.nVibPeaks*3 + 3,obj.nVibPeaks*3 + 3);
            covMatContDriveNoise = zeros((obj.nVibPeaks+1)*3);
            
            % --- vibrations ---
            AcVib = {}; BcVib = {}; CcVib = {};
            for kVib = 1:obj.nVibPeaks
                w0v(kVib) = obj.vibParms(kVib,1)*2*pi; % central vibration frequency
                csi(kVib) = obj.vibParms(kVib,2)/sqrt(14); % FWHM converted to damping coefficient
                varVib(kVib)  = obj.vibParms(kVib,3);
                
                [AcVib{kVib}, BcVib{kVib}, CcVib{kVib}] = tf2ss(w0v(kVib)^2,[1 2*csi(kVib)*w0v(kVib) w0v(kVib)^2]);
                
                % ---> average phase within two sampling instants at sampling T
                %Avr     =  1/Ts*CcVib*(expm(AcVib*Ts)-eye(2))/(AcVib);
                
                % ---> State covariance matrix for the continuous temporal state-space model
                Sigmaxv  = [-AcVib{kVib}(1,2)*varVib(kVib)/(CcVib{kVib}(end))^2,0;0 varVib(kVib)/(CcVib{kVib}(end))^2];
                
                % ---> Drive Noise matrix for the continuous temporal state-space model
                Sigmav  = [2*AcVib{kVib}(1,1)*AcVib{kVib}(1,2)*varVib(kVib)/(CcVib{kVib}(end))^2,0;0,0];
                
                % ---> Covariance matrix of the continuous drive noise
                Sigmav  = -AcVib{kVib}*Sigmaxv - Sigmaxv*AcVib{kVib}'; % == Sigmav
                
                % ---> Covariance matrix of the discrete noise
                Sigmadv  = Sigmaxv - expm(T*AcVib{kVib})*Sigmaxv*expm(T*AcVib{kVib}'); % this Sigmad should also be equal to the noise variance of the discrete system computed from Sigmad = -Ad11*Sigmav*Ad11 + Sigmav (the discrete time Lyapunov equation for the system X_{k+1} = A X_k + \nu_k)
                
                
                
                Ac((1:2) + (kVib-1)*3, (1:2) + (kVib-1)*3) = AcVib{kVib};
                Ac(3+(kVib-1)*3, (1:2)) = 1/T*CcVib{kVib};
                
                covXvib   = [-AcVib{kVib}(1,2)*varVib(kVib)/(CcVib{kVib}(end))^2,0;0 varVib(kVib)/(CcVib{kVib}(end))^2];
                covMatContDriveNoise((1:2) + (kVib-1)*3, (1:2) + (kVib-1)*3) = -AcVib{kVib}*covXvib - covXvib*AcVib{kVib}';
            end
            
            obj.continuousModel.AcVib = AcVib;
            obj.continuousModel.BcVib = BcVib;
            obj.continuousModel.CcVib = CcVib;
            
            % --- atmospheric modes ---
            p1Atm = obj.atmParms(1);
            p2Atm = obj.atmParms(2);
            varAtm  = obj.atmParms(3);
            
            
            [AcAtm, BcAtm, CcAtm] = tf2ss(p1Atm*p2Atm,[1 p1Atm+p2Atm p1Atm*p2Atm]);
            
            %ssAtm  = ss(AcAtm, BcAtm, CcAtm, 0);
            
            % --- state noise covariance matrix, output has varAtm variance ---
            covXatm   = [-AcAtm(1,2)*varAtm/(CcAtm(end))^2,0;0 varAtm/(CcAtm(end))^2];
            covMatContDriveNoiseAtm    = -AcAtm*covXatm - covXatm*AcAtm';   % result that comes from the solution of a Lyapunov equation: Atur*Sigmav + Sigmav*Atur' + Sigmac = 0, for stationary processes
            
            obj.continuousModel.AcAtm = AcAtm;
            obj.continuousModel.BcAtm = BcAtm;
            obj.continuousModel.CcAtm = CcAtm;
            
            % concatenation of continuous models
            if isempty(kVib)
                kVib = 0;
            end
            Ac((1:2) + (kVib)*3, (1:2) + (kVib)*3) = AcAtm;
            Ac(3+(kVib)*3, (1:2) + (kVib)*3) = 1/T*CcAtm;
            covMatContDriveNoise((1:2) + (kVib)*3, (1:2) + (kVib)*3) = -AcAtm*covXatm - covXatm*AcAtm';
            
            
            obj.continuousModel.Ac = full(Ac);
            obj.continuousModel.covMatContDriveNoise = covMatContDriveNoise;
        end
        
        %% DISCRETE TIME MODEL (ATM+VIB+NOISE)
        function obj = discreteTimeModel(obj)
            % discrete-time vibration model
            discreteVibrationModel(obj);
            % discrete-time atmospheric model
            discreteTurbulenceModel(obj);
            % discrete-time state-noise model
            discreteDriveNoiseModel(obj);
        end
        
        %% DISCRETE TURBULENCE MODEL
        function discreteTurbulenceModel(obj)             % --- turbulence ---
            T = obj.samplingTime;
            AcAtm = obj.continuousModel.AcAtm;
            CcAtm = obj.continuousModel.CcAtm;
            %-------------------------------------------------
            %     DISCRETE TEMPORAL MODEL
            %-------------------------------------------------
            AdAtm = zeros(4);
            AdAtm(1:2,1:2) = expm(AcAtm*T);
            AdAtm(3,1:2)   = 1/T*CcAtm*(expm(AcAtm*T)-eye(2))/AcAtm;
            AdAtm(4,3) = 1;
            
            if obj.AR
                % AdAtm = [1.9960   -0.99600];
                p1 = 1;
                p2 = 1.0;
                
                z1 = exp(-p1*obj.samplingTime);
                z2 = exp(-p2*obj.samplingTime);
                
                Hz = zpk([], [z1 z2], 1,obj.samplingTime); % a normalization by the gain@0 may be needed here
                Hz = tf(Hz); % convert to tf
                pp = get(Hz,'den');
                pp = pp{1};
                A = -pp(2);
                B = -pp(3);
                AdAtm = [A B;1 0];
                AdAtm(2,1) = 1;
                AdAtm(3,1) = 1;
                AdAtm(4,4) = 0;
                AdAtm(4,2) = 1;
            end
            
            BdAtm   = zeros(4,1);
            
            if obj.nDelay == 1
                CdAtm   = [0 0 1 0];
            elseif obj.nDelay == 2
                CdAtm   = [0 0 0 1];
            end
            obj.discreteModel.AdAtm = AdAtm;
            obj.discreteModel.BdAtm = BdAtm;
            obj.discreteModel.CdAtm = CdAtm;
        end
        
        %% DISCRETE VIBRATION MODEL
        function discreteVibrationModel(obj) % --- vibration ---
            T = obj.samplingTime;
            AcVib = obj.continuousModel.AcVib;
            BcVib = obj.continuousModel.BcVib;
            CcVib = obj.continuousModel.CcVib;
            %-------------------------------------------------
            %     DISCRETE TEMPORAL MODEL
            %-------------------------------------------------
            AdVib = {}; BdVib = {}; CdVib = {};covXvib = {};
            for kVib = 1:obj.nVibPeaks
                w0v = obj.vibParms(kVib,1)*2*pi; % central vibration frequency
                csi = obj.vibParms(kVib,2)/sqrt(14); % FWHM converted to dampling coefficient
                varVib  = obj.vibParms(kVib,3);
                a1v = 2*exp(-csi*w0v*T)*cos(w0v*sqrt(1-csi^2)*T);
                a2v = -exp(-2*csi*w0v*T);
                covXvib{kVib} = varVib*(1-a1v^2 - a2v^2 -2*a1v*a2v*a1v/(1-a2v));%w0v^2*T^2*18.8/5;
                
                
                AdVib{kVib} = zeros(4);
                AdVib{kVib}(1:3,1:3) = [a1v a2v 0; 1 0 0; 1 0 0];
                AdVib{kVib}(4,3) = 1;
                BdVib{kVib} = [0 0 0 0]';
                if obj.nDelay == 1
                    CdVib{kVib}   = [0 0 1 0];
                elseif obj.nDelay == 2
                    CdVib{kVib}   = [0 0 0 1];
                end
            end
            obj.discreteModel.AdVib = AdVib;
            obj.discreteModel.BdVib = BdVib;
            obj.discreteModel.CdVib = CdVib;
            obj.discreteModel.covXvib = covXvib;
            
            
        end
        %% DISCRETE STATE DRIVE NOISE MODEL
        function discreteDriveNoiseModel(obj)
            T = obj.samplingTime;
            Ac = obj.continuousModel.Ac;
            covMatContDriveNoise = obj.continuousModel.covMatContDriveNoise;
            
            nt = length(Ac);
            [t,Result] = ode113(@lyapoupou,[0 0.5*T T],zeros(1,nt*nt),[],Ac,covMatContDriveNoise);
            covMatDiscreteDriveNoiseBuf = reshape(Result(3,:),nt,nt);
            [covMatDiscreteDriveNoise{1:obj.nVibPeaks+1, 1:obj.nVibPeaks+1}] = deal(zeros(4));
            if obj.AR
                covMatDiscreteDriveNoiseBuf(1:3,1:3) = 0;
                covMatDiscreteDriveNoiseBuf(1) = obj.atmParms(3);
            end
            for kVib = 1:obj.nVibPeaks+1
                covMatDiscreteDriveNoise{kVib,kVib}(1:3,1:3) = covMatDiscreteDriveNoiseBuf((kVib-1)*3+1:(kVib-1)*3+3,(kVib-1)*3+1:(kVib-1)*3+3);
            end
            
            % state drive noise compatible with user-input vibration
            % model
            for kVib = 1:obj.nVibPeaks
                covMatDiscreteDriveNoise{kVib,kVib} = zeros(4);
                covMatDiscreteDriveNoise{kVib,kVib}(1,1) = obj.discreteModel.covXvib{kVib};
            end
            covMatDiscreteDriveNoise = cell2mat(covMatDiscreteDriveNoise);
            
            obj.discreteModel.covMatDiscreteDriveNoise = covMatDiscreteDriveNoise;
            
            %% Discrete-time Lyapunov solver
            function Y=lyapoupou(t,X,A,Q)
                [m,n]=size(A);
                sigma = reshape(X,n,n);
                sigmap = A*sigma + sigma*A' + Q;
                Y = sigmap(:);
            end
        end
        
        %% MODEL STACK
        function modelStack(obj)
            [Ad{1:obj.nVibPeaks+1, 1:obj.nVibPeaks+1}] = deal(zeros(4));
            for kVib = 1:obj.nVibPeaks
                Ad{kVib,kVib}    = obj.discreteModel.AdVib{kVib};
                Bd{kVib,1}       = obj.discreteModel.BdVib{kVib};
                Cd{1,kVib}       = obj.discreteModel.CdVib{kVib};
            end
            Ad{obj.nVibPeaks+1, obj.nVibPeaks+1} = obj.discreteModel.AdAtm;
            Bd{obj.nVibPeaks+1,1} = obj.discreteModel.BdAtm;
            Cd{1,obj.nVibPeaks+1} = obj.discreteModel.CdAtm;
            
            Ad = cell2mat(Ad);
            Bd = cell2mat(Bd);
            Cd = cell2mat(Cd);
            
            % controls in ss A matrix
            Ad(obj.nVibPeaks*4 + 4+2, obj.nVibPeaks*4 + 4+1 ) = 1;
            Ad(obj.nVibPeaks*4 + 4+2, obj.nVibPeaks*4 + 4+2 ) = 0;
            Bd(obj.nVibPeaks*4 + 4+1) = 1;
            Bd(obj.nVibPeaks*4 + 4+2) = 0;
            Cd(obj.nVibPeaks*4 + 4+2) = 0;
            if obj.nDelay == 1
                Cd(obj.nVibPeaks*4 + 4+1) = -1;
            elseif obj.nDelay == 2
                Cd(obj.nVibPeaks*4 + 4+2) = -1;
            end
            
            obj.discreteModel.Ad = Ad;
            obj.discreteModel.Bd = Bd;
            obj.discreteModel.Cd = Cd;
        end
        
        %% CONTROLLER SYNTHESIS
        function synthesis(obj)
            T           = obj.samplingTime;
            Qn          = obj.discreteModel.covMatDiscreteDriveNoise;               % State-equation noise covariance matrix
            Qn          = padarray(Qn,[2 2],0, 'post');
            Rn          = obj.noiseVar;                         % Measurement noise covariance matrix
            Nn          = 0;                                % Covariance of state w.r.t. measurement noise
            
            % --- state noise ---
            nVibPeaks = obj.nVibPeaks;
            [m,n] = size(obj.discreteModel.Ad);
            Ns = sparse(m,n);
            for kVib = 1:nVibPeaks
                Ns((kVib-1)*4+1:(kVib-1)*4+3,(kVib-1)*4+1:(kVib-1)*4+3) = eye(3);
            end
            if isempty(kVib)
                kVib = 0;
            end
            Ns((kVib)*4+1:(kVib)*4+3, (kVib)*4+1:(kVib)*4+3) = eye(3);
            sysd       = ss(obj.discreteModel.Ad, [obj.discreteModel.Bd full(Ns)], obj.discreteModel.Cd, 0, T);  % Create a discrete space-state model
            
            [kest,L,P_kalman,M,Z] = kalman(sysd,Qn,Rn,Nn);% Kalman-filter solver
            obj.KalmanFilter.M = M;
            obj.KalmanFilter.ssd = sysd;
            obj.KalmanFilter.kest = kest;
            
            K = zeros(1,4*(nVibPeaks+1)+2);
            for kVib = 1:nVibPeaks
                K((kVib-1)*4+1) = -1;
            end
            if isempty(kVib), kVib = 0;end
            if obj.AR
                K(kVib*4+1) = -1;
            else
                K(kVib*4+3) = -1;
            end
            obj.KalmanFilter.K = K;
        end
        %% BODE DIAGRAM
        function bode(obj,nu, inputPsd)
            T = obj.samplingTime;
            z = tf('z',T);
            if ~exist('nu','var')
                nu = logspace(-3,log10(1/(2*T)),500-1);
            end
            % REGULATOR
            rlqg = lqgreg(obj.KalmanFilter.kest,obj.KalmanFilter.K);
            % FEEDBACK CONTROLLER
            sysCL = feedback(1, z^(-(obj.nDelay))*rlqg);
            %sysCL = rlqg;
            [MAG,  ~, W] = bode(sysCL,2*pi*nu);
            hold on
            semilogx(W/2/pi, 20*log10(squeeze(MAG)));
            %loglog(W/2/pi, squeeze(MAG).^2);
            hold on
            axis tight
            title('Bode diagrame for the closed loop system')
            xlabel('frequency, [Hz]')
            ylabel('magnitude, [dB]')
            
            % for verification, next 3 lines give exactly same response as
            % bode
            [tfLQG,floc, RTF, NTF] = LQG_controller_TF(obj.discreteModel.Ad,...
                obj.discreteModel.Bd,obj.discreteModel.Cd,-obj.nDelay,obj.KalmanFilter.M,obj.KalmanFilter.K,T,nu,'estimator');
            loglog(floc, 10*log10(abs(RTF).^2),'r--')
            if nargin == 3
                figure
                loglog(W/2/pi, inputPsd)
                hold on
                loglog(W/2/pi, squeeze(MAG).^2.*inputPsd)
            end
        end
        %% RESIDUAL PSD
        function out = residualPSD(obj,timeSeries,hOl)
            T = obj.samplingTime;
            z = tf('z',T);
            if nargin < 3
                hOl = 1;
            end
            
            % disturbance PSD
            [psd, nu] = myPSD(timeSeries, 1/obj.samplingTime,0);
            % REGULATOR
            rlqg = lqgreg(obj.KalmanFilter.kest,obj.KalmanFilter.K);
            % FEEDBACK CONTROLLER
            sysCL = feedback(1, z^(-(obj.nDelay))*hOl*rlqg);
            
            [MAG,  ~, W] = bode(sysCL, 2*pi*nu);
            semilogx(W/2/pi, 10*log10(psd));
            hold on
            MAG = squeeze(MAG);
            MAG(isnan(MAG)) = 0;
            semilogx(W/2/pi, 10*log10(MAG.^2));
            semilogx(W/2/pi, 10*log10(MAG.^2.*psd));
            legend('input disturbance','RTF','residual disturbance')
            
            pbaspect([ 1.618 1 1])
            title('Single-Sided Amplitude Spectrum of y(t)','fontsize',14)
            xlabel('Temporal frequency [Hz]','fontsize',14)
            ylabel('a.u. |.|^2/Hz, [dB]','fontsize',14)
            box on
            axis tight
            fprintf('-----------------------------------------------\n')
            fprintf(' ------------    SUMMARY    -------------------\n')
            fprintf('rms of input disturbance: %f units rms\n',sqrt(trapz(nu, psd)))
            fprintf('rms of residual disturbance: %f units rms\n',sqrt(trapz(nu, squeeze(MAG).^2.*psd)))
            fprintf('-----------------------------------------------\n')
            out = sqrt(trapz(nu, squeeze(MAG).^2.*psd));
            noiseCL = feedback(rlqg, z^(-(obj.nDelay)));
            [MAG,  ~, W] = bode(noiseCL, 2*pi*nu);
            out(2) = 2*trapz(nu, squeeze(MAG).^2)*T;
        end
        %% TEMPORAL SIMULATION
        function sim(obj,timeSeries)
            T = obj.samplingTime;
            z = tf('z',T);
            rlqg = lqgreg(obj.KalmanFilter.kest,obj.KalmanFilter.K);
            % FEEDBACK CONTROLLER
            sysCL = feedback(z^(-(obj.nDelay)), rlqg);
            
            
            [o, tsim] = lsim(sysCL, timeSeries);
            figure
            plot(tsim, timeSeries, tsim, o,'r')
            xlabel('Time, [s]')
            ylabel('Input units, [n/a]')
            title('Temporal simulation')
            legend('Input','Output')
            pbaspect([ 1.618 1 1])
            grid
        end
    end
end
