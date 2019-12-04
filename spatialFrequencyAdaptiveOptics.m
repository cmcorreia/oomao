classdef spatialFrequencyAdaptiveOptics < handle
    
    properties
        psf;
        psfDl;
        structureFunction;
        paramListener;
        psfRootListener;
    end
    
    properties (SetObservable=true)
        tel;
        atm;
        src;
        noiseVariance;
        nActuator;
        psfResolution;
        psfPixelScaleInMas;
        loopGain;     % loop gain
        samplingTime; % exposure time
        latency;      % delay
        psfRoot = 1;
        Rx; Ry;
        noiseGain;
        h1;h2;hn;
        rtf;ntf;atf;
        SxAv; SyAv;
        nTimes = 3;
        fx;fy;
        MV; % 1 if MV, 0 if LS
        nTh; % # of points for average over # windDirections [0...2/pi]
        modulation = 4;
        filter;
        binning;
        f_filter; %Spatial filter frequency
    end
    
    properties (Dependent)
        fc;
        fcInMas;
    end
    
    properties (Access=private)
        figHandle
    end
    
    methods
        %%
        function obj = spatialFrequencyAdaptiveOptics(tel, atm, ...
                nActuator, noiseVariance, loopGain, samplingTime, latency, resolution, filter, MV, varargin)
            expectedFilters = {'Rigaut','Fried','Hudgin','Anti-Alias','Southwell','hudgin', 'fried', 'southwell', 'rigaut', 'antialias', 'freischlad86','pyr','pyramid','Fauvarque'};
            defaultFilter = 'Rigaut';
            inputs = inputParser;
            addRequired(inputs,'tel',@(x) isa(x,'telescope') );
            addRequired(inputs,'atm',@(x) isa(x,'atmosphere') );
            addRequired(inputs,'nActuator', @isnumeric);
            addRequired(inputs,'noiseVariance', @isnumeric);
            addOptional(inputs,'loopGain', 0.5, @isnumeric);
            addOptional(inputs,'samplingTime', 1e-3, @isnumeric);
            addOptional(inputs,'latency', 1e-3, @isnumeric);
            addOptional(inputs,'resolution', 81, @isnumeric);
            addOptional(inputs,'filter',defaultFilter, @(x) any(validatestring(x,expectedFilters)));
            addOptional(inputs,'MV',1,@isnumeric);
            addOptional(inputs,'modulation',1,@isnumeric);
            addOptional(inputs,'nTimes',2,@isnumeric);
            addOptional(inputs,'nTh',1,@isnumeric);
            addOptional(inputs,'binning',1,@isnumeric);
            addOptional(inputs,'f_filter',[],@isnumeric);
            parse(inputs, tel, atm, ...
                nActuator, noiseVariance, loopGain, samplingTime, latency, resolution,filter, MV, varargin{:});
            
            obj.tel           = inputs.Results.tel;
            obj.atm           = inputs.Results.atm;
            obj.nActuator     = inputs.Results.nActuator;
            obj.noiseVariance = inputs.Results.noiseVariance;
            obj.loopGain      = inputs.Results.loopGain;
            obj.latency       = inputs.Results.latency;
            obj.samplingTime  = inputs.Results.samplingTime;
            obj.samplingTime  = inputs.Results.samplingTime;
            resolution        = inputs.Results.resolution;             % RESOLUTION AND SPATIAL-FREQUENCY VECTORS
            obj.filter        = inputs.Results.filter;
            obj.MV            = inputs.Results.MV;
            obj.nTh           = inputs.Results.nTh;
            obj.modulation    = inputs.Results.modulation;
            obj.nTimes        = inputs.Results.nTimes;
            obj.binning       = inputs.Results.binning;
            obj.psfResolution = resolution*obj.nTimes;
            obj.f_filter      = inputs.Results.f_filter*obj.fc;
            
            [obj.fx,obj.fy] = freqspace(resolution,'meshgrid');
            obj.fx = obj.fx*obj.fc + 1e-7;
            obj.fy = obj.fy*obj.fc + 1e-7;
            
            % RECONSTRUCTION FILTER
            if strcmp(filter,'pyr') || strcmp(filter,'pyramid')|| strcmp(filter,'Conan')|| strcmp(filter,'Fauvarque')
                pyramidReconstructionFilter(obj)
            else
                reconstructionFilter(obj);
            end
            % CONTROLLER
            if obj.loopGain < 0
                controller(obj,'lqg')
            else
                controller(obj,'int')
            end
            % DEFAULT SOURCE
            % Source
            obj.src = source('zenith',0*constants.arcsec2radian,'azimuth',0,'wavelength',photometry.H,'magnitude',12);
            
            obj.paramListener = ...
                addlistener(obj,...
                {'noiseVariance','nActuator','tel',...
                'atm','psfResolution','psfPixelScaleInMas',...
                'loopGain','samplingTime','latency','src'},...
                'PostSet',@obj.resetPsf);
            obj.psfRootListener = addlistener(obj,'psfRoot','PostSet',@obj.resetPsfScale);
            
        end
        %%
        function out = get.fc(obj)
            out = 1/obj.binning*0.5*(obj.nActuator-1)/obj.tel.D;
        end
        %%
        function out = get.fcInMas(obj)
            out = obj.fc*obj.atm.wavelength*constants.radian2mas;
        end
        %%
        %         function out = get.noise(obj)
        %
        %         end
        
        %%
        function pyramidReconstructionFilter(obj,fx_,fy_)
            nL      = obj.nActuator - 1;
            d       = obj.tel.D/nL;
            if nargin > 1
                f       = hypot(fx_,fy_);
            else
                f       = hypot(obj.fx,obj.fy);
            end
            Wn      = obj.noiseVariance/(2*obj.fc)^2*obj.MV;
            Wphi    = phaseStats.spectrum(f,obj.atm);
            u       = obj.fx;
            if strcmp(obj.filter,'pyr') || strcmp(obj.filter,'pyramid')|| strcmp(obj.filter,'Conan')
                % forward pyramid filter (continuous)
                umod = 1/(2*d)/(nL/2)*obj.modulation;
                Sx = zeros(size(u));
                idx = abs(u) > umod;
                Sx(idx) = 1i*sign(u(idx));
                idx = abs(u) <= umod;
                Sx(idx) = 2*1i/pi*asin(u(idx)/umod);
                Av = sinc(obj.binning*d*u).*sinc(obj.binning*d*u)';
                Sy = Sx.';
                obj.SxAv = Sx.*Av;
                obj.SyAv = Sy.*Av;
                
                %reconstruction filter
                AvRec = Av;
                SxAvRec = Sx.*AvRec;
                SyAvRec = Sy.*AvRec;
                
            elseif strcmp(obj.filter,'Fauvarque')
                %% transforms
                fft2c       = @(x) tools.fftcentre(x,1);
                ifft2c      = @(x) (tools.fftcentre(x,-1));
                %fftConv2    = @(x,y) ifft2c(fft2c(x).*(fft2c(y)));
                fftConv2    = @(x,y) conv2(x,y,'same');
                %% create a transparency mask
                
                nPx = size(obj.fx,2);
                % extend the mask to avoid Fourier periodicity
                if iseven(nPx)
                    nPx = nPx*2;
                else
                    nPx = 2*(nPx-1) + 1;
                end
                [fxn,fyn] = freqspace(nPx,'meshgrid');
                %vec = linspace(-1,1,nPx)*2;
                %[fxn, fyn] = meshgrid(vec,vec);
                
                l_alpha = ones(4,1);%*pi/2;
                
                % pyramid face transmitance
                mask  = heaviside(fxn).*heaviside(fyn);
                phase = -l_alpha(1).*(fxn+fyn);
                pym   = mask.*exp(1i.*phase);%
                m1    = pym;
                
                mask  = heaviside(fxn).*heaviside(-fyn);
                phase = -l_alpha(1).*(fxn-fyn);
                m2    = mask.*exp(1i.*phase);
                pym   =  m2 + pym;
                
                mask  = heaviside(-fxn).*heaviside(fyn);
                phase = -l_alpha(1).*(-fxn+fyn);
                m3    = mask.*exp(1i.*phase);
                pym   = m3 + pym;
                
                mask  = heaviside(-fxn).*heaviside(-fyn);
                phase = -l_alpha(1).*(-fxn-fyn);
                m4    = mask.*exp(1i.*phase);
                pym   = m4 + pym;
                %% O. Fauvarque's slopesMaps model
                %oversampling by a factor 2 to avoid circularity + aliasing
                os = (nPx-1)/(obj.nActuator-1);
                ff = hypot(fxn,fyn);
                thePupil = utilities.piston(obj.nActuator-1,nPx);
                weightSignal = ifft2c(thePupil.*besselj(0,obj.modulation*os*2*ff));
                weightSignal = weightSignal/sum(real(weightSignal(:)));
                
                filterFunc = fftConv2(m3,m2.*weightSignal)-fftConv2(m2,m3.*weightSignal)+...
                    fftConv2(m1,m4.*weightSignal)-fftConv2(m4,m1.*weightSignal);
                Sx = 2*1i*(real(filterFunc));%.*sinc(fxn).*sinc(fyn);
                
                % crop to the original size
                nPts = size(obj.fx,2);
                offSet = (nPts-1)/2;%-2;
                Sx = Sx(offSet + (1:nPts),offSet + (1:nPts))*2;

                Sy = Sx.';
                obj.SxAv = Sx;
                obj.SyAv = Sy;
                
                %reconstruction filter
                AvRec = 1;
                SxAvRec = Sx;
                SyAvRec = Sy;
            end
            
            
            % --------------------------------------
            %   MMSE filter
            % --------------------------------------
            gPSD = abs(Sx.*AvRec).^2 + abs(Sy.*AvRec).^2 + Wn./Wphi +1e-7;
            obj.Rx = conj(SxAvRec)./gPSD;
            obj.Ry = conj(SyAvRec)./gPSD;
        end
        %% Shack-Hartmann filter and reconstructor
        function reconstructionFilter(obj,fx_,fy_)
            d = obj.tel.D/(obj.nActuator-1);
            N = size(obj.fx,1);
            Ts = obj.samplingTime;
            td = obj.latency;
            [vx,vy] = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            Wn = obj.noiseVariance/(2*obj.fc)^2*obj.MV;
            if nargin > 2
                f  = hypot(fx_,fy_);
            else
                f  = hypot(obj.fx,obj.fy);
                fx_ = obj.fx;
                fy_ = obj.fy;
            end
            Wphi = phaseStats.spectrum(f,obj.atm);
            
            Sx = 1i*2*pi*fx_*d ;
            Sy = 1i*2*pi*fy_*d ;
            
            timeAverage = 0;
            for kLayer = 1:obj.atm.nLayer
                timeAverage     = timeAverage + ...
                    obj.atm.layer(kLayer).fractionnalR0*...
                    sinc(vx(kLayer)*Ts*fx_).*sinc(vy(kLayer)*Ts*fy_).*exp(2*1i*pi*(fx_*vx(kLayer)*td+fy_*vy(kLayer)*td));
            end
            timeAverage = 1;
            Av              = sinc(d*fx_).*sinc(d*fy_).*exp(1i*pi*d*(fx_+fy_)).*timeAverage;
            obj.SxAv = Sx.*Av;
            obj.SyAv = Sy.*Av;
            
            % average the temporal shift operators over all wind directions
            v = sqrt(vx.^2 + vy.^2);
            nTh_ = obj.nTh;
            thetaWind = linspace(0, 2*pi-2*pi/nTh_,nTh_) + atan(vy/vx);
            cubeInt = [];
            for kLayer = 1:obj.atm.nLayer
                for kTheta = 1:nTh_
                    kvx = v(kLayer)*cos(thetaWind(kTheta));
                    kvy = v(kLayer)*sin(thetaWind(kTheta));
                    cubeInt(:,:,kTheta)  = sinc(kvx*Ts*fx_).*sinc(kvy*Ts*fy_).*exp(2*1i*pi*(fx_*kvx*td+fy_*kvy*td));%;
                end
                timeAverageThetaWind(:,:,kLayer) = obj.atm.layer(kLayer).fractionnalR0*sum(cubeInt,3)/nTh_;
            end
            
            timeAverageThetaWind = sum(timeAverageThetaWind,3);
            timeAverageThetaWind = 1;
            clear cubeInt
            AvRec           = sinc(d*fx_).*sinc(d*fy_).*exp(1i*pi*d*(fx_+fy_)).*timeAverageThetaWind;
            SxAvRec = Sx.*AvRec;
            SyAvRec = Sy.*AvRec;
            %Av              = sinc(d*fx).*sinc(d*fy).*exp(1i*pi*d*(fx+fy)).*timeAverage;
            filter_ = obj.filter;
            if strcmp(filter_, 'hudgin') || strcmp(filter_, 'Hudgin')
                % --------------------------------------
                %   HUDGIN filter
                % --------------------------------------
                gPSD = 4*(sin(pi*d*fx_).^2 + sin(pi*d*fy_).^2) + Wn./Wphi;
                %                 % USED IN THE JOSA-A 2014 PAPER
                %                 Rx = (exp(-2*pi*1i*d*fx)-1)./gPSD.*exp(-2*1i*pi*d/8*(fx+fy));
                %                 Ry = (exp(-2*pi*1i*d*fy)-1)./gPSD.*exp(-2*1i*pi*d/8*(fx+fy));
                
                % modified Hudgin a-la GPI (after Poyneer's remarks)
                obj.Rx = (exp(-2*pi*1i*d*fx_)-1)./gPSD.*exp(-2*1i*pi*d/2*(fy_));
                obj.Ry = (exp(-2*pi*1i*d*fy_)-1)./gPSD.*exp(-2*1i*pi*d/2*(fx_));
                
            elseif strcmp(filter_, 'fried') || strcmp(filter_, 'Fried')
                % --------------------------------------
                %   FRIED filter
                % --------------------------------------
                gPSD = 8*((sin(pi*d*fx_).^2).*(cos(pi*d*fy_).^2) + (sin(pi*d*fy_).^2).*(cos(pi*d*fx_).^2)) + 2*Wn./Wphi;
                Rx = (exp(-2*pi*1i*d*fx_)-1).*(exp(-2*pi*1i*d*fy_)+1)./gPSD;
                Ry = (exp(-2*pi*1i*d*fy_)-1).*(exp(-2*pi*1i*d*fx_)+1)./gPSD;
                
                % local waffle filtering
                LWR = 1/4*(3 + exp(-2*pi*1i*d*fy_) + exp(-2*pi*1i*d*fx_) -   exp(-2*pi*1i*d*(fx_ + fy_)));
                
                obj.Rx = Rx.*LWR;
                obj.Ry = Ry.*LWR;
            elseif strcmp(filter_, 'southwell') || strcmp(filter_, 'Southwell')
                % --------------------------------------
                %   SOUTHWELL filter
                % --------------------------------------
                gPSD = 4*(sin(pi*d*fx_).^2 + sin(pi*d*fy_).^2) + Wn./Wphi;
                obj.Rx = (-1i*sin(2*pi*d*fx_))./gPSD.*exp(-2*1i*pi*d/2*(fx_+fy_));
                obj.Ry = (-1i*sin(2*pi*d*fy_))./gPSD.*exp(-2*1i*pi*d/2*(fx_+fy_));
                
            elseif strcmp(filter_, 'rigaut') || strcmp(filter_, 'Rigaut')
                % --------------------------------------
                %   EXACT filter
                % --------------------------------------
                gPSD = abs(Sx.*AvRec).^2 + abs(Sy.*AvRec).^2 + Wn./Wphi;
                obj.Rx = conj(SxAvRec)./gPSD;
                obj.Ry = conj(SyAvRec)./gPSD;
                
            elseif strcmp(filter_, 'Anti-Alias') || strcmp(filter_, 'antialias')
                if obj.atm.nLayer > 1
                    r0 = [obj.atm.r0 obj.atm.layer.fractionnalR0];
                else
                    r0 = obj.atm.r0;
                end
                [~, WSx, WSy] = mk_aliased_PSD(fx_,2,d, r0,obj.atm.L0,0, Ts, 0, vx, vy);
                gPSD = WSx + WSy + Wn;
                obj.Rx = conj(SxAvRec).*Wphi./gPSD;
                obj.Ry = conj(SyAvRec).*Wphi./gPSD;
                
            elseif strcmp(filter_, 'freischlad86')
                gPSD = (4*(1-cos(2*pi*d*fx_).*cos(2*pi*d*fy_)));
                obj.Rx = (exp(-2*pi*1i*d*(fx_+fy_))-1)./gPSD;
                obj.Ry = (exp(-2*pi*1i*d*fx_)-exp(-2*pi*1i*d*fy_))./gPSD;
            else
                fprintf('No filter chosen. Please choose "hudgin", "fried", "southwell", "Rigaut", "Anti-Alias" or "freischlad86".')
            end
            obj.Rx(isnan(obj.Rx)) = 0;
            obj.Ry(isnan(obj.Ry)) = 0;
            
            % Set central point (i.e. kx=0,ky=0) to zero
            obj.Rx(ceil((N+1)/2),ceil((N+1)/2)) = 0;
            obj.Ry(ceil((N+1)/2),ceil((N+1)/2)) = 0;
        end
        %%
        function  controller(obj,type,fx_,fy_,varargin)
            if nargin < 3
                fx_ = obj.fx;
                fy_ = obj.fy;
            end
            Ts = obj.samplingTime;
            td = obj.latency;
            delay = floor(td/Ts);
            [vx,vy] = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            z = tf('z',Ts);
            %floc = linspace(1e-3,0.5/Ts,500);
            floc = logspace(-3,log10(0.5/Ts),500);
            tPSD = phaseStats.temporalSpectrum(floc,obj.atm);
            if strcmp(type,'int')
                %--------------------------------------------------------------------------
                % INTEGRATOR
                hInt = obj.loopGain/(1-z^(-1));
                rtfInt       = 1/(1+hInt*z^(-delay));
                atfInt = hInt*z^(-delay)*rtfInt;%/(1+hInt*z^(-delay));
                if obj.loopGain == 0 || isempty(obj.loopGain)
                    ntfInt = tf(1);
                else
                    ntfInt = atfInt/z;%=hInt*z^(-delay)/z/(1+hInt*z^(-delay)); % = ntfInt = z^(-1)*gloop*hInt*rtfInt;
                end
                magnitudeNtfInt = bode(ntfInt, 2*pi*floc);
                integratorNoisePropFactor = trapz(floc, squeeze(magnitudeNtfInt).^2)*2*Ts;
                
                obj.rtf = rtfInt;
                obj.atf = atfInt;
                obj.ntf = ntfInt;
                obj.noiseGain = integratorNoisePropFactor;
                
                
                %nTh = obj.nTh;
                [obj.h1, obj.h2] = integralOverWindDirections(obj,atfInt,fx_,fy_);
                [obj.hn] = integralOverWindDirections(obj,ntfInt,fx_,fy_);
                
                if exist('hTFs','var')
                    figure(hTFs);
                else
                    hTFs = figure;
                end
                magRtf = bode(obj.rtf, 2*pi*floc);
                magNtf = bode(obj.ntf, 2*pi*floc);
                magAtf = bode(obj.atf, 2*pi*floc);
                
                semilogx(floc, 10*log10(tPSD),'k',...
                    floc, 10*log10(squeeze(magRtf).^2),'g--',...
                    floc, 10*log10(squeeze(magNtf).^2),'r--',...
                    floc, 10*log10(squeeze(magAtf).^2),'b-.')
                title('LOOP TRANSFER FUNCTIONS','fontweight','bold','fontsize',18)
                xlabel('temporal frequency, [Hz]', 'fontsize',18)
                ylabel('[dB]', 'fontsize',18)
                axis tight, grid on
                legend('uncorrected phase','signal rejection','noise rejection','aliasing rejection')
                
                fprintf('Integrator residual: %f\n',trapz(floc, tPSD'.*squeeze(magRtf).^2)*2 + integratorNoisePropFactor*obj.noiseVariance)
                
                
                % CC Comment: if the integral controller is employed in a loop with no
                % delay then the INTnoisePropFactor = g/(2-g), Flicker's Eq. (39) in "Analytical evaluations of closed-loop adaptive optics spatial power spectral densities"
                % closedLoopRTF_Intz       = 1/(1+z^(-0)*gloop*Hint);
                % NTF_Intz = z^(-0)*gloop*Hint*closedLoopRTF_Intz;
                % MAG = bode(NTF_Intz, 2*pi*floc);
                % INTnoisePropFactor = trapz(floc, squeeze(MAG).^2)*2*Ts
                
                
            elseif strcmp(type,'lqg')
                %--------------------------------------------------------------------------
                % LQG
                
                nPts    = size(fx_,1);
                d       = obj.tel.D/(obj.nActuator-1);
                [vx,vy] = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
                nLayer  = obj.atm.nLayer;
                for kLayer = 1:nLayer
                    nu = -vx(kLayer)*fx_ + -vy(kLayer)*fy_;
                    Atur{kLayer}    = exp(-1i*pi*2*Ts*nu); % the minus sign is because temporal prediction forward in time is shift backwards in space
                    %                     Alias{kLayer} = Atur{kLayer};
                    %                     m = 1;%obj.nTimes*0;
                    %                     for mi = -m:m
                    %                         for ni = -m:m
                    %                             if mi~=0 || ni~=0
                    %                                 fm = fx - mi/d;
                    %                                 fn = fy - ni/d;
                    %                                 nu = -vx(kLayer)*fm + -vy(kLayer)*fn;
                    %                                 Alias{kLayer}    = Alias{kLayer}.* exp(-1i*pi*2*Ts*nu); % the minus sign is because temporal prediction forward in time is shift backwards in space
                    %                             end
                    %                         end
                    %                     end
                end
                
                l = ones(1,nLayer);
                a = 0.9999;
                phaseSpatialSpectrum = phaseStats.spectrum(hypot(fx_,fy_),obj.atm);
                e = phaseSpatialSpectrum*(1- a^2);
                
                
                Wa = mk_aliased_PSD(fx_, obj.nTimes, d, obj.atm.r0, obj.atm.L0, 1, obj.tel.samplingTime);
                W0 = mk_aliased_PSD(fx_, 0, d, obj.atm.r0, obj.atm.L0, 1, obj.tel.samplingTime);
                aliasingPsd = Wa - W0; % aliasing PSD
                b = 0.999;
                eAliasing = aliasingPsd*(1-b^2);
                fractionalR0 = [obj.atm.layer.fractionnalR0];
                nLayer = obj.atm.nLayer;
                noiseVariance_ = obj.noiseVariance;
                samplingTime_ = obj.samplingTime;
                fc_ = obj.fc;
                for kx = 1:nPts
                    for ky = 1:nPts
                        A{kx,ky} = zeros(nLayer + 5);
                        B{kx,ky} = zeros(nLayer + 5,1);
                        C{kx,ky} = zeros(1,nLayer + 5);
                        for kLayer = 1:nLayer
                            A{kx,ky}(kLayer,kLayer) = a*Atur{kLayer}(kx,ky);
                            A{kx,ky}(nLayer+1,kLayer) = a*Atur{kLayer}(kx,ky);
                        end
                        %A{kx,ky}(nLayer+1,1:nLayer) = l;
                        A{kx,ky}(nLayer+2,nLayer+1) = 1;
                        A{kx,ky}(nLayer+3,nLayer+2) = 1;
                        A{kx,ky}(nLayer+5,nLayer+4) = 1;
                        
                        B{kx,ky}(nLayer+4) = 1;
                        
                        r = obj.Rx(kx,ky)*obj.SxAv(kx,ky) + obj.Ry(kx,ky)*obj.SyAv(kx,ky);
                        C{kx,ky}(nLayer+3) = r;
                        C{kx,ky}(nLayer+5) = -r;
                        %C{kx,ky}(nLayer+3) = 1;
                        %C{kx,ky}(nLayer+5) = -1;
                        
                        %                         C{kx,ky}(1,nLayer+3) = obj.SxAv(kx,ky);
                        %                         C{kx,ky}(2,nLayer+3) = obj.SyAv(kx,ky);
                        %                         C{kx,ky}(1,nLayer+5) = -obj.SxAv(kx,ky);
                        %                         C{kx,ky}(2,nLayer+5) = -obj.SyAv(kx,ky);
                        % KF synthesis
                        V = [eye(nLayer+1);zeros(4,nLayer+1)];
                        sys = ss(A{kx,ky}, [B{kx,ky} V], C{kx,ky}, 0, Ts);
                        
                        
                        QN = diag([fractionalR0*e(kx,ky) 0*1e-3]); % this does not seem to provide the MV results (24 Feb 2016)
                        RN = noiseVariance_/(fc_*2)^2*eye(1);
                        
                        QN = diag([fractionalR0*mean(e(:)) 1e-3]);
                        %RN = 0.01;
                        [KEST,Linf,P,Minfi{kx,ky}] = kalman(sys,QN,RN,0);
                        
                        Kinf = zeros(1,nLayer+5);
                        Kinf(nLayer+1) = -1;
                        % ------------------------------------------------
                        %[~,floc, ~, NTFi] = LQG_controller_TF(A{kx,ky},B{kx,ky},C{kx,ky},-1,Minfi{kx,ky},Kinf,samplingTime_,[],'estimator');
                        %lqgNoisePropFactor(kx,ky) = trapz(floc, abs(NTFi))*Ts;
                        % updated previous 2 lines on 17th Jan 17
                        [~,floc, RTFi, NTFi] = LQG_controller_TF(A{kx,ky},B{kx,ky},C{kx,ky},-2-1i,Minfi{kx,ky},Kinf,samplingTime_,[],'estimator');
                        loglog(floc(end/2+1:end), abs(RTFi(end/2+1:end)),'b')
                        hold on
                        loglog(floc(end/2+1:end),abs(RTFi(1:end/2)),'r')
                        lqgNoisePropFactor(kx,ky) = ...
                            trapz(floc(end/2+1:end), abs(NTFi(end/2+1:end)))*Ts-...
                            trapz(floc(1:end/2), abs(NTFi(1:end/2)))*Ts;
                        % ------------------------------------------------
                    end
                end
                obj.noiseGain = lqgNoisePropFactor;
                %nTh = obj.nTh;
                [obj.h1, obj.h2, obj.hn] = integralOverWindDirectionsDKF(obj,A, B, C, Minfi,fx_,fy_);
                
            end
            
            
            
        end
        
        function [h1, h2] = integralOverWindDirections(obj,RTF,fx_,fy_)
            if nargin < 3
                fx_ = obj.fx;
                fy_ = obj.fy;
            end
            %--------------------------------------------------------------------------
            [vx,vy]     = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            nPts        = size(fx_,1);
            nTh_        = obj.nTh;
            thetaWind   = linspace(0, 2*pi-2*pi/nTh_,nTh_);
            costh       = cos(thetaWind);
            fr0         = [obj.atm.layer.fractionnalR0];
            h1          = zeros(nPts); h2 = zeros(nPts);
            for kLayer = 1:obj.atm.nLayer
                for iTheta = 1:nTh_
                    fi = -vx(kLayer)*fx_*costh(iTheta) + -vy(kLayer)*fy_*costh(iTheta);
                    idx = abs(fi) <1e-7;
                    fi(idx) = 1e-8.*sign(fi(idx));
                    [MAG, PH] = bode(RTF, 2*pi*fi(:));
                    MAG = reshape(MAG,[nPts,nPts]);
                    MAG(fi == 0) = 1;
                    PH = reshape(PH,[nPts,nPts]);
                    h2buf(:,:,iTheta) = abs(MAG.*exp(1i*PH/180*pi)).^2;%*2*pi;
                    h1buf(:,:,iTheta) = MAG.*exp(1i*PH/180*pi);%*2*pi;
                end
                h1 = h1 + fr0(kLayer)*sum(h1buf,3)/nTh_;
                h2 = h2 + fr0(kLayer)*sum(h2buf,3)/nTh_;
            end
        end
        
        function [h1, h2, hn] = integralOverWindDirectionsDKF(obj,A, B, C, Minfi,fx,fy)
            if nargin < 3
                fx = obj.fx;
                fy = obj.fy;
            end
            %--------------------------------------------------------------------------
            [vx,vy]     = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            nPts        = size(fx,1);
            nTh_        = obj.nTh;
            thetaWind   = linspace(0, 2*pi-2*pi/nTh_,nTh_);
            costh       = cos(thetaWind);
            fr0         = [obj.atm.layer.fractionnalR0];
            h1          = zeros(nPts); h2 = zeros(nPts);hn = zeros(nPts);
            Kinf        = zeros(1,obj.atm.nLayer+5);
            Kinf(obj.atm.nLayer+1) = -1;
            %fx = (fx + obj.fc)/(2*obj.fc); % convert to discrete meaningful frequencies [0:N-1]/N
            %fy = (fy + obj.fc)/(2*obj.fc);
            fx = obj.fx; fy = obj.fy;
            for kLayer = 1:obj.atm.nLayer
                for iTheta = 1:nTh_
                    nu = (-vx(kLayer)*fx*costh(iTheta) + -vy(kLayer)*fy*costh(iTheta));
                    for kx = 1:nPts
                        parfor ky = 1:nPts
                            [~,~, RTF(kx,ky), NTF(kx,ky), olRTF(kx,ky)] = ...
                                LQG_controller_TF(A{kx,ky},B{kx,ky},C{kx,ky},-2-i,Minfi{kx,ky},Kinf,obj.samplingTime,nu(kx,ky),'estimator');
                        end
                    end
                    MAG = abs(olRTF);
                    PH = angle(olRTF);
                    
                    h2buf(:,:,iTheta) = abs(MAG.*exp(1i*PH/180*pi)).^2;%*2*pi;
                    h1buf(:,:,iTheta) = MAG.*exp(1i*PH/180*pi);%*2*pi;
                    
                    MAG = abs(NTF);
                    PH = angle(NTF);
                    hnbuf(:,:,iTheta) = MAG.*exp(1i*PH/180*pi);%*2*pi;
                end
                h1 = h1 + fr0(kLayer)*sum(h1buf,3)/nTh_;
                h2 = h2 + fr0(kLayer)*sum(h2buf,3)/nTh_;
                hn = hn + fr0(kLayer)*sum(hnbuf,3)/nTh_;
            end
        end
        
        function out = fittingPSD(obj,fx,fy,flagFitting)
            %% FITTINGPSD Fitting error power spectrum density
            if nargin < 2
                fx = obj.fx;
                fy = obj.fy;
            end
            if nargin<4
                flagFitting = 'square';
            end
            
            [fx,fy] = freqspace(size(fx,1)*obj.nTimes,'meshgrid');
            fx = fx*obj.fc*obj.nTimes;
            fy = fy*obj.fc*obj.nTimes;
            out   = zeros(size(fx));
            if strcmp(flagFitting,'square')
                index  = abs(fx)>obj.fc | abs(fy)>obj.fc;
                
            else
                index  = hypot(fx,fy) > obj.fc;
            end
            f     = hypot(fx(index),fy(index));
            out(index) ...
                = phaseStats.spectrum(f,obj.atm);
            out = out.*pistonFilter(obj,hypot(fx,fy));
        end
        
        function out = noisePSD(obj,fx_,fy_)
            %% NOISEPSD Noise error power spectrum density
            if nargin < 2
                fx_ = obj.fx;
                fy_ = obj.fy;
            end
            fc = obj.fc;
            out   = zeros(size(fx_));
            if obj.noiseVariance>0
                index = ~(abs(fx_)>fc | abs(fy_)>fc) & hypot(fx_,fy_)>0;
                %                f     = hypot(fx(index),fy(index));
                %                 out(index) = obj.noiseVariance./...
                %                     ( 2*pi*f.*tools.sinc(0.5*fx(index)/fc).*tools.sinc(0.5*fy(index)/fc)).^2;
                out(index) = obj.noiseVariance/(2*obj.fc)^2.*(abs(obj.Rx(index)).^2 + abs(obj.Ry(index)).^2);
                %out(index) = obj.noiseVariance*(abs(obj.Rx(index)).^2 + abs(obj.Ry(index).^2)); % See JM Conan thesis, Annex A
                if obj.loopGain
                    
                    % Next 2 lines account for noise propagation twice
                    %E = abs(obj.hn).^2;
                    %out = out.*E.*pistonFilter(obj,hypot(fx_,fy_)).*obj.noiseGain;
                    
                    %out = out.*pistonFilter(obj,hypot(fx_,fy_)).*obj.noiseGain*5*5;
                    % this 5x5 factor was initially used to check the
                    % pyramid simulations but is deprecated
                    
                    out = out.*pistonFilter(obj,hypot(fx_,fy_)).*obj.noiseGain;
                else
                    out = out.*pistonFilter(obj,hypot(fx_,fy_));
                end
                
            end
        end
        function [out, out1] = aliasingPSD(obj,fx,fy,flagFitting)
            if nargin < 2
                fx = obj.fx;
                fy = obj.fy;
            end
            if nargin < 4
                flagFitting = 'square';
            end
            out   = zeros(size(fx));            
            index = ~(abs(fx)>obj.fc | abs(fy)>obj.fc);% & hypot(fx,fy)>1/obj.tel.D;           
            fx      = fx(index);
            d       = obj.tel.D/(obj.nActuator-1);
            Ts      = obj.samplingTime;
            td      = obj.latency;
            [vx,vy] = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            side    = sqrt(length(fx));
            
            fx = reshape(fx, side,side);
            localRx = reshape(obj.Rx(index),side,side);
            localRy = reshape(obj.Ry(index),side,side);
            if obj.loopGain == 0 || isempty(obj.loopGain)
                tf = [];
            else
                tf = reshape(obj.h1(index),side,side);
            end
            if obj.atm.nLayer > 1
                r0 = [obj.atm.r0 obj.atm.layer.fractionnalR0];
            else
                r0 = obj.atm.r0;
            end
            [PSD_ALIASING] = mk_propagated_aliased_PSD(fx, obj.nTimes, d, r0, obj.atm.L0,...
                localRx, localRy, vx, vy, Ts, td, obj.tel.D, tf,obj.filter,obj.f_filter);
            
            if strcmp(flagFitting,'circular')
                PSD_ALIASING = PSD_ALIASING.*(hypot(fx,fy)<=obj.fc);
            end
            
            out(index) = PSD_ALIASING;
            out1(index) = PSD_ALIASING;
        end
        function out = aliasingPSDold(obj,fx,fy)
            %% ALIASINGPSD Aliasing error power spectrum density
            
            fc    = obj.fc;
            out   = zeros(size(fx));
            index = ~(abs(fx)>fc | abs(fy)>fc);% & hypot(fx,fy)>1/obj.tel.D;
            pf = pistonFilter(obj,hypot(fx,fy));
            fx     = fx(index);
            fy     = fy(index);
            %             figure,imagesc(ind),pause
            [fo,f] = cart2pol(fx,fy);
            n = 5;
            al = 0;
            nv = -n:n;
            for l=nv
                flx = fx - 2*l*fc;
                for m=nv
                    if (l~=0) && (m~=0)
                        fmy = fy - 2*m*fc;
                        flm = hypot(flx,fmy);
                        al = al + 0.25*sin(2*fo).^2.*(fx./fmy+fy./flx).^2.*...
                            phaseStats.spectrum(flm,obj.atm);
                    end
                end
            end
            l=0;
            flx = fx - 2*l*fc;
            ind = flx==0;
            nind = ~ind;
            nv(nv==0) = [];
            for m=nv
                fmy = fy - 2*m*fc;
                flm = hypot(flx,fmy);
                al(ind) = al(ind) + phaseStats.spectrum(flm(ind),obj.atm);
                al(nind) = al(nind) + 0.25*sin(2*fo(nind)).^2.*(fx(nind)./fmy(nind)+fy(nind)./flx(nind)).^2.*...
                    phaseStats.spectrum(flm(nind),obj.atm);
            end
            m = 0;
            fmy = fy - 2*m*fc;
            ind = fmy==0;
            nind = ~ind;
            for l=nv
                flx = fx - 2*l*fc;
                flm = hypot(flx,fmy);
                al(ind) = al(ind) + phaseStats.spectrum(flm(ind),obj.atm);
                al(nind) = al(nind) + 0.25*sin(2*fo(nind)).^2.*(fx(nind)./fmy(nind)+fy(nind)./flx(nind)).^2.*...
                    phaseStats.spectrum(flm(nind),obj.atm);
            end
            if obj.loopGain
                out(index) =  al.*averageClosedLoopAliasing(obj,fx,fy);
            else
                out(index) =  al;
            end
            out = out.*pf;
        end
        
        function out = servoLagPSD(obj,fx,fy)
            %% SERVOLAGPSD Servo-lag power spectrum density
            if nargin < 2
                fx = obj.fx;
                fy = obj.fy;
            end
            fc    = obj.fc;
            out   = zeros(size(fx));
            index = ~(abs(fx)>=fc | abs(fy)>=fc);
            pf = pistonFilter(obj,hypot(fx,fy));
            fx     = fx(index);
            fy     = fy(index);
            
            %out(index) = phaseStats.spectrum(hypot(fx,fy),obj.atm).*averageClosedLoopRejection(obj,fx,fy);
            F = (obj.Rx(index).*obj.SxAv(index) + obj.Ry(index).*obj.SyAv(index));
            if obj.loopGain == 0 || isempty(obj.loopGain)
                out(index) = abs(1-F).^2.*phaseStats.spectrum(hypot(fx,fy),obj.atm);
            else
                out(index) = (1 + abs(F).^2.*obj.h2(index) - F.*obj.h1(index) - conj(F.*obj.h1(index))) ...
                    .*phaseStats.spectrum(hypot(fx,fy),obj.atm);
            end
            out = pf.*real(out);
        end
        
        function out = anisoServoLagPSD(obj,fx,fy,iSrc)
            %% SERVOLAGPSD Servo-lag power spectrum density
            if nargin < 2
                fx = obj.fx;
                fy = obj.fy;
            end
            if(nargin<4)
                iSrc=1;
            end
            fc    = obj.fc;
            out   = zeros(size(fx));
            index = ~(abs(fx)>=fc | abs(fy)>=fc);
            pf = pistonFilter(obj,hypot(fx,fy));
            fx     = fx(index);
            fy     = fy(index);
            if ~isempty(obj.src)
                zLayer = [obj.atm.layer.altitude];
                fr0    = [obj.atm.layer.fractionnalR0];
                A = zeros(size(fx));
                for kLayer=1:obj.atm.nLayer
                    red = 2*pi*zLayer(kLayer)*...
                        ( fx*obj.src(iSrc).directionVector(1) + fy*obj.src(iSrc).directionVector(2) );
                    A  = A + fr0(kLayer)*exp(1i*red);
                end
            else
                A = ones(size(fx));
            end
            %out(index) = phaseStats.spectrum(hypot(fx,fy),obj.atm).*averageClosedLoopRejection(obj,fx,fy);
            F = (obj.Rx(index).*obj.SxAv(index) + obj.Ry(index).*obj.SyAv(index));
            if obj.loopGain == 0 || isempty(obj.loopGain)
                out(index) = abs(1-F).^2.*phaseStats.spectrum(hypot(fx,fy),obj.atm);
            else
                out(index) = (1 + abs(F).^2.*obj.h2(index) - F.*obj.h1(index).*A(index) - conj(F.*obj.h1(index).*A(index))) ...
                    .*phaseStats.spectrum(hypot(fx,fy),obj.atm);
            end
            out = pf.*real(out);
        end
        
        function out = anisoplanatismPSD(obj,fx,fy,iSrc,flagFilterPiston)
            if nargin < 2
                fx = obj.fx;
                fy = obj.fy;
            end
            if nargin <5
                flagFilterPiston = true;
            end
            %% ANISOPLANATISM Anisoplanatism power spectrum density
            if(nargin<4)
                iSrc=1;
            end
            zLayer = [obj.atm.layer.altitude];
            fr0    = [obj.atm.layer.fractionnalR0];
            A = zeros(size(fx));
            for kLayer=1:obj.atm.nLayer
                red = 2*pi*zLayer(kLayer)*...
                    ( fx*obj.src(iSrc).directionVector(1) + fy*obj.src(iSrc).directionVector(2) );
                A  = A + 2*fr0(kLayer)*( 1 - cos(red) );
            end     
            pf = 1;
            if flagFilterPiston
                pf = pistonFilter(obj,hypot(fx,fy));
            end
            out = pf.*A.*phaseStats.spectrum(hypot(fx,fy),obj.atm);
        end
        
        function out = powerSpectrumDensity(obj,fx,fy)
            %% POWERSPECTRUMDENSITY AO system power spectrum density
            if nargin < 2
                fx = obj.fx;
                fy = obj.fy;
            end
            fc    = obj.fc;
            
            [fxExt,fyExt] = freqspace(size(fx,1)*obj.nTimes,'meshgrid');
            fxExt = fxExt*obj.fc*obj.nTimes;
            fyExt = fyExt*obj.fc*obj.nTimes;
            out   = zeros(size(fxExt));
            index = abs(fxExt)<fc & abs(fyExt)<fc;
            
            out(index) = noisePSD(obj,fx,fy) + ...
                aliasingPSD(obj,fx,fy) + ...
                anisoServoLagPSD(obj,fx,fy);
            
            out = out + fittingPSD(obj,fx,fy);
            
        end
        
        function out = varFitting(obj)
            fc  = obj.fc;
            a = phaseStats.variance(obj.atm);
            b = integral2( @(fx,fy) phaseStats.spectrum( hypot(fx,fy) , obj.atm ) , ...
                -fc,fc,-fc,fc);
            out = a - b;
        end
        
        
        
        function out = varAliasing(obj)
            %fc  = obj.fc;
            %out = quad2d( @(fx,fy) aliasingPSD(obj,fx,fy),-fc,fc,-fc,fc);
            fx_ = obj.fx;
            fy_ = obj.fy;
            out = trapz(fy_(:,1),trapz(fx_(1,:),aliasingPSD(obj,fx_,fy_),2));
        end
        
        
        function out = varServoLag(obj)
            %fc  = obj.fc;
            %out = quad2d( @(fx,fy) servoLagPSD(obj,fx,fy),-fc,fc,-fc,fc);
            fx_ = obj.fx;
            fy_ = obj.fy;
            out = trapz(fy_(:,1),trapz(fx_(1,:),anisoServoLagPSD(obj,fx_,fy_),2));
        end
        
        function out = varNoise(obj)
            %fc  = obj.fc;
            %out = quad2d( @(fx,fy) noisePSD(obj,fx,fy),-fc,fc,-fc,fc);
            fx_ = obj.fx;
            fy_ = obj.fy;
            out = trapz(fy_(:,1),trapz(fx_(1,:),noisePSD(obj,fx_,fy_),2));
        end
        
        function out = varAo(obj,fLim)
            %out = quad2d( @(fx,fy) powerSpectrumDensity(obj,fx,fy),-fLim,fLim,-fLim,fLim);
            out = obj.varAliasing + obj.varFitting + obj.varNoise + obj.varServoLag;
        end
        
        function varargout = image(obj,resolution,pixelScaleInMas, telOtf)
            %% IMAGE Point spread function
            wvlRatio = obj.atm.wavelength/obj.src.wavelength;
            
            fprintf('Computing image plane ...\n')
            %obj.psfResolution = resolution*obj.nTimes;
            %obj.psfPixelScaleInMas = pixelScaleInMas;
            
            pixelScale = pixelScaleInMas*1e-3*constants.arcsec2radian/obj.src.wavelength;
            
            %[fx_,fy_] = freqspace(resolution,'meshgrid');
            %fx_ = pixelScale*fx_*resolution/2;
            %fy_ = pixelScale*fy_*resolution/2;
            
            psd = powerSpectrumDensity(obj)*(wvlRatio)^2;
            %obj.atm.wavelength = obj.src.wavelength;
            
            
            sf  = fft2(fftshift(psd))*pixelScale^2;
            sf  = 2*fftshift( sf(1) - sf );
            obj.structureFunction = sf;
            
            [rhoX,rhoY] = freqspace(resolution*obj.nTimes,'meshgrid');
            rhoX = 0.5*rhoX/pixelScale;
            rhoY = 0.5*rhoY/pixelScale;
            if all(abs(rhoX)<obj.tel.D)
                warning('oomao.fourierAdaptiveOptics.image','OTF is truncated: increase the resolution or decrease the pixel scale')
            end
            
            [u,v] = freqspace(resolution*obj.nTimes,'meshgrid');
            fftPhasor = exp(1i.*pi.*(u+v)*0.5);
            %             rho  = hypot(rhoX,rhoY);
            if ~exist('telOtf','var')
                telOtf = otf(obj.tel, rhoX+1i*rhoY);
            end
            obj.psfDl = real(ifftshift(ifft2(ifftshift(telOtf.*fftPhasor))))/pixelScale^2;
            obj.psfDl = obj.psfDl/obj.tel.area;
            
            thisOtf = telOtf.*exp(-0.5*sf);
            obj.psf = real(ifftshift(ifft2(ifftshift(thisOtf.*fftPhasor))))/pixelScale^2;
            obj.psf = obj.psf/obj.tel.area;
            
            alpha = pixelScaleInMas*((0:resolution*obj.nTimes)-resolution*obj.nTimes/2);%(-resolution+1:2:resolution-1)/2;
            if isempty(obj.figHandle)
                obj.figHandle = figure;
            end
            if isvalid(obj.figHandle)
                figure(obj.figHandle)
            else
                figure
            end
            imagesc(alpha,alpha,log10(abs(obj.psf).^(1/obj.psfRoot)))
            if any(abs(alpha)>obj.fcInMas)
                u = [-1 1 1 -1 -1]*obj.fcInMas/wvlRatio;
                v = [-1 -1 1 1 -1]*obj.fcInMas/wvlRatio;
                line(u,v,'color','r','linestyle',':')
            end
            axis xy equal tight
            xlabel('X axis [mas]')
            ylabel('Y axis [mas]')
            title(sprintf('Strehl: %4.2f%',sum(thisOtf(:))./sum(telOtf(:))*100))
            
            if nargout>0
                varargout{1} = obj.psf;
                varargout{2} = abs(sum(thisOtf(:))./sum(telOtf(:)));
            end
            %             pxScaleAtNyquist = 0.25/tel.D;
            %             imgLens = lens;
            %             imgLens.nyquistSampling = ...
            %                 pxScaleAtNyquist/pixelScale;
            
            obj.atm.wavelength = 0.5e-6;
        end
        
        function resetPsf(obj,varargin)
            if ~isempty(obj.psf)
                image(obj,obj.psfResolution/obj.nTimes,obj.psfPixelScaleInMas)
            end
        end
        
        function resetPsfScale(obj,varargin)
            if ~isempty(obj.psf)
                figure(obj.figHandle)
                h = findobj(obj.figHandle,'type','image');
                set(h,'Cdata',obj.psf.^(1/obj.psfRoot))
            end
        end
        
        function out = pistonFilter(obj,f)
            red = pi*obj.tel.D*f;
            out = 1 - 4*tools.sombrero(1,red).^2;
            
        end
        
        function out = closedLoopRejection(obj,nu)
            %% CLOSEDLOOPREJECTION Closed loop rejection transfer function
            
            out   = zeros(size(nu));
            index = nu~=0;
            nu   = nu(index);
            red = obj.loopGain.*tools.sinc(nu.*obj.samplingTime)./(2*pi*nu*obj.samplingTime);
            out(index) = ....
                1./( 1 + red.^2 - 2.*red.*sin( 2*pi*nu*(obj.samplingTime+obj.latency) ) );
            
        end
        
        function E = averageClosedLoopRejection(obj,fx,fy)
            %% AVERAGECLOSEDLOOPREJECTION Atmosphere average closed loop rejection transfer function
            
            E = averageRejection(obj,fx,fy,@(nu)closedLoopRejection(obj,nu));
        end
        
        function out = closedLoopAliasing(obj,nu)
            %% CLOSEDLOOPALIASING Closed loop aliasing transfer function
            
            out   = ones(size(nu));
            index = nu~=0;
            nu   = nu(index);
            red = obj.loopGain.*tools.sinc(nu.*obj.samplingTime)./(2*pi*nu*obj.samplingTime);
            out(index) = ....
                red.^2./( 1 + red.^2 - 2.*red.*sin( 2*pi*nu*(obj.samplingTime+obj.latency) ) );
            
        end
        
        function E = averageClosedLoopAliasing(obj,fx,fy)
            %% AVERAGECLOSEDLOOPALIASING Atmosphere average closed loop aliasing transfer function
            
            E = averageRejection(obj,fx,fy,@(nu)closedLoopAliasing(obj,nu));
        end
        
        function out = closedLoopNoise(obj,nu)
            %% CLOSEDLOOPALIASING Closed loop noise transfer function
            
            out   = ones(size(nu));
            index = nu~=0;
            nu   = nu(index);
            red = obj.loopGain.*tools.sinc(nu.*obj.samplingTime)./(2*pi*nu*obj.samplingTime);
            out(index) = (red./tools.sinc(nu.*obj.samplingTime)).^2./...
                ( 1 + red.^2 - 2.*red.*sin( 2*pi*nu*(obj.samplingTime+obj.latency) ) );
            
        end
        
        function E = averageClosedLoopNoise(obj,fx,fy,rtf)
            %% AVERAGECLOSEDLOOPNOISE Atmosphere average closed loop noise transfer function
            
            %E = averageRejection(obj,fx,fy,@(nu)closedLoopNoise(obj,nu));
            E = averageRejection(obj,fx,fy,@(nu)closedLoopFilter(obj,nu,rtf));
        end
        
        function E = averageRejection(obj,fx,fy,fun,tf)
            [vx,vy] = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            fr0     = [obj.atm.layer.fractionnalR0];
            E = zeros(size(fx));
            for kLayer=1:obj.atm.nLayer
                nu = fx*vx(kLayer) + fy*vy(kLayer);
                E  = E + fr0(kLayer)*fun(nu,tf);
            end
            
        end
    end
    
    %%
    %**************************************************************
    %**************************************************************
    % EXAMPLES
    %**************************************************************
    %**************************************************************
    methods (Static)
        
        function fao = psdDemo
            fao = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',0,'Rigaut',1);
        end
        function pyramidDemo
            spherePyr = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',0,'pyr',1);
        end
        %--------------------------------------------------------------
        % JOSA-A PAPER: Modelling astronomical adaptive optics performance with temporally filtered Wiener reconstruction of slope data
        %--------------------------------------------------------------
        %% section "Limiting post-coronagraphic contrast under predictive closed-loop control"
        function limitingContrast
            mags = [0 5 10 12 14 16];
            for m = 1:length(mags)
                close all
                fao = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',0,'Rigaut',1,mags(m));
                varFit = fao.varFitting;
                varNoise = fao.varNoise;
                varAliasing = fao.varAliasing;
                varAnisoServoLag = fao.varServoLag;
                varTotal = varFit + varNoise + varAliasing + varAnisoServoLag;
                SRint(m) = 100*exp(-varTotal*(0.5/1.65)^2);
                cInt(:,:,m) = fao.powerSpectrumDensity(fao.fx,fao.fy);
                psfInt(:,:,m) = fao.psf;
                psfDlInt(:,:,m) = fao.psfDl;
                psdResolution = length(fao.fx);
                pixelScale  = fao.atm.wavelength*(fao.nActuator-1)/fao.tel.D/psdResolution*constants.radian2arcsec;
                o = linspace(-1,1,psdResolution)*pixelScale*(1.65/0.5)*psdResolution/2;
                psdTotalIntLS = anisoServoLagPSD(fao) + aliasingPSD(fao) + noisePSD(fao);
                deltaF = fao.fy(2) - fao.fy(1);
                psdTotalIntLS = psdTotalIntLS*deltaF^2;
                psdTotalIntLS(end/2+0.5, end/2+0.5) = nan;
                
                %                 close all
                %                 fao = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',1,'Rigaut',1,mags(m));
                %                 varFit = fao.varFitting;
                %                 varNoise = fao.varNoise;
                %                 varAliasing = fao.varAliasing;
                %                 varAnisoServoLag = fao.varServoLag;
                %                 varTotal = varFit + varNoise + varAliasing + varAnisoServoLag;
                %                 SRmv(m) = 100*exp(-varTotal*(0.5/1.65)^2);
                %                 cMv(:,:,m) = fao.powerSpectrumDensity(fao.fx,fao.fy);
                %                 psfMv(:,:,m) = fao.psf;
                %                 psfDlMv(:,:,m) = fao.psfDl;
                
                
                close all
                faoDkf = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',1,'Rigaut',-1,mags(m));
                varFit = faoDkf.varFitting;
                varNoise = faoDkf.varNoise;
                varAliasing = faoDkf.varAliasing;
                varAnisoServoLag = faoDkf.varServoLag;
                varTotal = varFit + varNoise + varAliasing + varAnisoServoLag;
                SRdkf(m) = 100*exp(-varTotal*(0.5/1.65)^2);
                cDkf(:,:,m) = faoDkf.powerSpectrumDensity(faoDkf.fx,faoDkf.fy);
                psfDkf(:,:,m) = faoDkf.psf;
                psfDlDkf(:,:,m) = faoDkf.psfDl;
                psdTotalDkf = anisoServoLagPSD(faoDkf) + aliasingPSD(faoDkf) + noisePSD(faoDkf);
                deltaF = faoDkf.fy(2) - faoDkf.fy(1);
                psdTotalDkf = psdTotalDkf*deltaF^2;
                psdTotalDkf(end/2+0.5, end/2+0.5) = nan;
                
                
                % contrast ratio image
                figure
                imagesc(o, o, log10(psdTotalDkf./psdTotalIntLS*SRint(m)/SRdkf(m)))
                hold on
                contour(o, o, log10(psdTotalDkf./psdTotalIntLS*SRint(m)/SRdkf(m))+0.5,10)
                title(['post-coronagraphic contrast ratio, m_v=' num2str(mags(m))] ,'fontweight','bold','fontsize',24)
                xlabel('[arcsec]','fontsize',24,'fontweight','bold')
                ylabel('[arcsec]','fontsize',24,'fontweight','bold')
                %colormap('bone')
                load mycmap.mat % contrast colormap with
                colormap(mycmap)
                set(gca,'FontSize',18)
                colorbar
                axis square
                %print(['sphereConstrastEnhancementMag' num2str(mags(m))],'-depsc')
                savefig(['sphereContrastEnhancementMag' num2str(mags(m))]);
            end
        end
        %% section: "Keck?s NIRC2 performance and post- coronagraphic contrast enhancement"
        function keckPaperResults
            fao1 = spatialFrequencyAdaptiveOptics.psdPaper('keck','keck',0,'Rigaut',1);
            fao2 = spatialFrequencyAdaptiveOptics.psdPaper('keck','keck',1,'Rigaut',1);
            fao3 = spatialFrequencyAdaptiveOptics.psdPaper('keck','keck',1,'Anti-Alias',1);
            fao4 = spatialFrequencyAdaptiveOptics.psdPaper('keck','keck',0,'Rigaut',-1);
            fao5 = spatialFrequencyAdaptiveOptics.psdPaper('keck','keck',1,'Rigaut',-1);
            fao6 = spatialFrequencyAdaptiveOptics.psdPaper('keck','keck',1,'Anti-Alias',-1);
            c1 = fao1.powerSpectrumDensity(fao1.fx,fao1.fy);
            c4 = fao4.powerSpectrumDensity(fao1.fx,fao1.fy);
            c6 = fao6.powerSpectrumDensity(fao1.fx,fao1.fy);
            
            resolution = 41;
            pixelScale  = fao1.src.wavelength*(fao1.nActuator-1)/fao1.tel.D/resolution;
            pixelScaleInMas = constants.radian2arcsec*1e3*pixelScale;
            alpha       = constants.radian2arcsec*pixelScale*((0:resolution*fao1.nTimes)-resolution*fao1.nTimes/2); % arcsec on sky :: pixelScale = n lambda/D /nPts, n=nLenslets=fc*2*D
            
            imagesc(alpha, alpha,log10(c6./c1))
            hold on
            contour(alpha(1:end-1), alpha(1:end-1),log10(c6./c1)+0.1)
            title('Contrast ratio map','fontsize',22,'fontweight','bold')
            xlabel('[arcsec]','fontsize',22), ylabel('[arcsec]','fontsize',22)
            colorbar
            xlim([-0.3 0.3]),ylim([-0.3 0.3])
            set(gca,'FontSize',18)
            figure
            c = c6./c4;
            imagesc(alpha, alpha,log10(c))
            hold on
            contour(alpha(1:end-1), alpha(1:end-1),log10(c)+0.1)
            title('Contrast ratio map','fontsize',22,'fontweight','bold')
            xlabel('[arcsec]','fontsize',22), ylabel('[arcsec]','fontsize',22)
            colorbar
            xlim([-0.3 0.3]),ylim([-0.3 0.3]), zlim([-3 1.1])
            set(gca,'FontSize',18)
        end
        
        %%
        function sphere % comparing to Jolissaint EOS 2010 - removed from revised paper
            sphereInt = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',0,'Rigaut',1);
            figure(4)
            set(gca','YLim',[1e-0 10e4])
            set(gca,'XLim',[-1 1]*5)
            
            sphereMV = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',1,'Rigaut',1);
            figure(4)
            set(gca','YLim',[1e-0 10e4])
            set(gca,'XLim',[-1 1]*5)
            
            sphereDKF = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',1,'Rigaut',-1);
            figure(4)
            set(gca','YLim',[1e-0 10e4])
            set(gca,'XLim',[-1 1]*5)
            
            sphereAADKF = spatialFrequencyAdaptiveOptics.psdPaper('vlt9layers','sphere',1,'Anti-Alias',-1);
            figure(4)
            set(gca','YLim',[1e-0 10e4])
            set(gca,'XLim',[-1 1]*5)
        end
        %%
        function eltScaoResults % - removed from revised paper
            eltInt = spatialFrequencyAdaptiveOptics.psdPaper('elt','hScaoElt',0,'Rigaut',1);
            figure
            subplot(1,3,1)
            plot([eltInt.atm.layer.fractionnalR0], [eltInt.atm.layer.altitude]/1e3,'k-o')
            ylabel('altitude [km]')
            xlabel('fractional r_0')
            title('C_n^2 profile')
            axis tight
            subplot(1,3,2)
            plot([eltInt.atm.layer.windSpeed], [eltInt.atm.layer.altitude]/1e3,'ko')
            xlabel('wind speed [m/s]')
            axis tight
            title('wind profile')
            subplot(1,3,3)
            plot( [eltInt.atm.layer.windDirection], [eltInt.atm.layer.altitude]/1e3,'ko')
            xlabel('wind direction [rad]')
            axis tight
            title('wind direction')
            
            eltDKF = spatialFrequencyAdaptiveOptics.psdPaper('elt','hScaoElt',1,'Anti-Alias',-1);
            
            % contrast improvement
            cInt = eltInt.powerSpectrumDensity(eltInt.fx,eltInt.fy);
            cDKF = eltDKF.powerSpectrumDensity(eltInt.fx,eltInt.fy);
            
            resolution = 149;
            pixelScale  = eltInt.src.wavelength*(eltInt.nActuator-1)/eltInt.tel.D/resolution;
            pixelScaleInMas = constants.radian2arcsec*1e3*pixelScale;
            alpha       = constants.radian2arcsec*pixelScale*((0:resolution*eltInt.nTimes)-resolution*eltInt.nTimes/2); % arcsec on sky :: pixelScale = n lambda/D /nPts, n=nLenslets=fc*2*D
            
            figure
            imagesc(alpha, alpha,log10(cDKF./cInt), [-1 0.5])
            title('Contrast ratio map','fontsize',16,'fontweight','bold')
            xlabel('[arcsec]','fontsize',16), ylabel('[arcsec]','fontsize',14)
            colorbar
            xlim([-0.35 0.35]),ylim([-0.35 0.35]), colormap bone
        end
        %% psdPaper instantiates the ATM, TEL, SRC and plots the PSDs, PSFs and post-coronagraphic PSFs and WFE for the individual error terms
        function fao = psdPaper(atmCase,aoSystem,mvFlag, filterType,CTYPE,mag,T, wfsWvl,noiseVar)
            % THIS EXAMPLE COMPUTES THE RESULTS TO BE GATHERED IN THE JOSA-A PAPER
            %--------------------------------------------------------------
            % ATMOSPHERE TYPES
            %--------------------------------------------------------------
            if ~exist('atmCase','var')
                atmCase = 'keck';
            end
            switch atmCase
                case 'mono'
                    % 1/ mono-layer atmosphere
                    d_atm = atmosphere(photometry.V0,15e-2,25,'windSpeed',15,'windDirection',pi,'altitude',0e3);
                case 'Jolissaint'
                    % 2/ Jolissaint EOS 2010 paper atmosphere
                    r0 = 15.5e-2;        % Fried parameter
                    L0 = 25;             % Outer scale
                    
                    % Altitude of the 9 different layers
                    alt = [42 140 281 562 1125 2250 4500 9000 18000];
                    % Wind speed of the different layers
                    wS = [15 13 13 9 9 15 25 40 21];
                    % Wind direction for the different layers
                    wD = [38 34 54 42 57 48 -102 -83 -77]*pi/180;
                    % Turbulent energy in each layer
                    fR0 = [53.28 1.45 3.5 9.57 10.83 4.37 6.58 3.71 6.71]/100;
                    % Atmosphere
                    d_atm = atmosphere(photometry.V0,r0,L0,'fractionnalR0',fR0,...
                        'altitude',alt,'windSpeed',wS,'windDirection',wD);
                    % Source
                    src = source('zenith',3*constants.arcsec2radian,'azimuth',pi/4*0,'wavelength',photometry.H,'magnitude',12);
                case 'vlt9layers'
                    % 2/ Jolissaint EOS 2010 paper atmosphere
                    r0 = 15.5e-2;        % Fried parameter
                    L0 = 25;             % Outer scale
                    
                    % Altitude of the 9 different layers
                    alt = [42 140 281 562 1125 2250 4500 9000 18000];
                    % Wind speed of the different layers
                    wS = [15 13 13 9 9 15 25 40 21];
                    % Wind direction for the different layers
                    wD = [38 34 54 42 57 48 -102 -83 -77]*pi/180;
                    % Turbulent energy in each layer
                    fR0 = [53.28 1.45 3.5 9.57 10.83 4.37 6.58 3.71 6.71]/100;
                    % Atmosphere
                    d_atm = atmosphere(photometry.V0,r0,L0,'fractionnalR0',fR0,...
                        'altitude',alt,'windSpeed',wS,'windDirection',wD);
                case 'multi'
                    % 4/ General Multi-layer atmosphere
                    fractionalR0    = [0.596, 0.224, 0.180]; %0.7,0.25,0.05];
                    altitude        = [0e3,5.5e3,11e3];
                    windSpeed       = [15,10,20];
                    windDirection   = [0,pi/4,pi];
                    
                    d_atm = atmosphere(photometry.V0,0.155,25,...
                        'fractionnalR0',fractionalR0,'altitude',altitude,...
                        'windSpeed',windSpeed,'windDirection',windDirection);
                case 'keck'
                    % 5/ Keck parameters from KAON 716
                    fractionalR0    = [0.517, 0.119, 0.063, 0.061, 0.105, 0.081, 0.054];
                    altitude        = [0,0.5, 1, 2, 4, 8, 16]*1e3;
                    windSpeed       = [6.7, 13.9, 20.8, 29.0, 29.0, 29.0, 29.0];
                    windDirection   = [0,pi/3,-pi/3, -pi, -4/3*pi, -pi/6, pi/8];
                    
                    d_atm = atmosphere(photometry.V0,0.16,75,...
                        'fractionnalR0',fractionalR0,'altitude',altitude,...
                        'windSpeed',windSpeed,'windDirection',windDirection);
                case 'elt'
                    profile = fitsread('profil_turbulent_eso.fits');
                    fractionalR0 = profile(:,4)/100;
                    r0 = 0.151;
                    fractionalR0 = fractionalR0/sum(fractionalR0);
                    altitude = profile(1:35,2);
                    windSpeed            = profile(1:35,3);
                    windDirection        = profile(1:35,9);
                    d_atm = atmosphere(photometry.V0,r0,50,...
                        'fractionnalR0',fractionalR0,'altitude',altitude,...
                        'windSpeed',windSpeed,'windDirection',windDirection);
                    
            end
            %d_atm.wavelength = photometry.H;
            
            %--------------------------------------------------------------
            % AO SYSTEM TYPES
            %--------------------------------------------------------------
            if ~exist('aoSystem','var')
                aoSystem = 'keck';
            end
            if ~exist('mvFlag','var')
                mvFlag = 1; % LS = 0, MV = 1
            end
            if ~exist('filterType','var')
                filterType = 'Rigaut';%'Anti-Alias';%
            end
            if ~exist('CTYPE','var')
                CTYPE = 1; % -1: LQG (or DKF), 0: open-loop, g \in \Re+: integrator with gain g
            end
            if ~exist('mag','var')
                mag = 12; %
            end
            if ~exist('wfsWvl','var') || isempty(wfsWvl)
                wfsWvl = 'R';
            end
            if ~exist('noiseVar','var')
            % SH WFS NOISE :: (values computed in caseStudy.m for SPHERE case, RON=1, CoG centroiding)
            mags = [0     5     8    10    12    14    16    18];
            switch wfsWvl
                case 'V0'
                    noiseVarVals = [3.20143552585353e-05,0.00321477083563130,0.0541207123617045,0.454830181340785,7.38239178717385,226.229877695295,8579.41470289181,338858.722400154];
                case 'R'
                    noiseVarVals = [2.64116486150240e-05,0.00265024120671990,0.0441610929362495,0.355787572846906,5.31626120132333,155.817713306330,5850.97376496079,230708.969764456];
                case 'H'
                    noiseVarVals = [9.60511477815319e-05,0.00972513256615308,0.182663218109808,2.17269113661420,54.3222204771618,1959.60125042775,76732.2540795827,3046684.26966993];
            end
            noiseVar = interp1(mags, noiseVarVals, mag);
            end
            switch aoSystem
                case 'general'
                    %fao = fourierAdaptiveOptics(d_tel,d_atm,10,0*.1,0.5,1e-3,1e-3);
                    d_tel = telescope(8);
                    resolution = 81; % needs be at least 2*nLenlets+1
                    T = 2e-3;
                    tau = T + 0.8e-3;
                    %noiseVar = 7.6;% (see caseStudy.m) before:.0045*25*0;
                    nTh_ = 1;
                    fao = spatialFrequencyAdaptiveOptics(d_tel,d_atm,41,noiseVar,CTYPE*0.5,T,tau,resolution,filterType,mvFlag);
                case 'sphere'
                    %--------------------------------------------------------------
                    % VLT-SPHERE-like
                    d_tel = telescope(8);
                    resolution = 81; % needs be at least 2*nLenlets+1
                    %noiseVar = 7.6;% (see caseStudy.m) for a mag=12
                    %noiseVar = 1;
                    if ~exist('T','var')
                        T = 0.00066;
                        tau = 2.8e-3;
                    else
                        tau = 1.5*T;
                    end
                    nAct = 41;
                    wfs = shackHartmann(nAct-1,(nAct-1)*6,0.5);
                    wfs.lenslets.throughput = 0.31;
                    wfs.camera.exposureTime = T;
                    switch wfsWvl
                        case 'V0'
                            wfs.camera.readOutNoise = 0.3;
                            ngs = source('wavelength',photometry.V0);
                        case 'R'
                            wfs.camera.readOutNoise = 0.3;
                            ngs = source('wavelength',photometry.R);
                        case 'H'
                            wfs.camera.readOutNoise = 1;
                            ngs = source('wavelength',photometry.H);
                    end
                    
                    ngs.magnitude = mag;
                    thNoiseVar = wfs.theoreticalNoise(d_tel, d_atm, ngs, ngs,'centroidingAlgorithm','cog','emccd',1);
                    if isempty(noiseVar)
                        noiseVar = thNoiseVar(1);
                    end
                    %noiseVar = 0; %debug
                    if strcmp(filterType,'pyr') || strcmp(filterType,'Fauvarque')
                        %noiseVar = noiseVar/pi^2; % commented since the
                        %user now defines this input :: ccorreia 11 Oct 2019
                        p_modulation = 0;
                        fao = spatialFrequencyAdaptiveOptics(d_tel,d_atm,nAct,noiseVar,CTYPE*0.3,T,tau,resolution,filterType,mvFlag,'modulation',p_modulation);
                    else
                        fao = spatialFrequencyAdaptiveOptics(d_tel,d_atm,nAct,noiseVar,CTYPE*0.3,T,tau,resolution,filterType,mvFlag);
                    end
                case 'keck'
                    %--------------------------------------------------------------
                    % Keck-like
                    resolution = 41; % needs be at least 2*nLenlets+1
                    
                    pupilKeck   = fitsread('keckPupil_460px.fits');
                    % interpolate to nRes resolution
                    xi = linspace(0,1,length(pupilKeck));
                    [Xi, Yi] = meshgrid(xi, xi);
                    xo= linspace(0,1,resolution);
                    [Xo, Yo] = meshgrid(xo, xo);
                    pupilKeck = interp2(Xi, Yi, pupilKeck, Xo, Yo);
                    pupilKeck = pupilKeck./max(pupilKeck(:));
                    pupilKeck(pupilKeck < 0.05) = 0;
                    d_tel = telescope(11.25);
                    d_tel.pupil = pupilKeck;
                    Ts = 2e-3;
                    tau = 2*Ts;
                    noiseVar = 0.0847; % mag 8, ron=4.5
                    %noiseVar = 3.7;% mag 12, ron=4.5
                    nTh_ = 1;
                    fao = spatialFrequencyAdaptiveOptics(d_tel,d_atm,20,noiseVar,CTYPE*0.5,Ts,tau,resolution,filterType,mvFlag);
                    
                    % create Keck pupil OTF
                    pupilKeck   = fitsread('keckPupil_460px.fits');
                    pupilKeck   = fitsread('keckPupil.fits'); %794x794 squared, aligned pupil (no rotation as former)
                    % interpolate to nRes resolution
                    xi = linspace(0,1,length(pupilKeck));
                    [Xi, Yi] = meshgrid(xi, xi);
                    xo= linspace(0,1,resolution*fao.nTimes);
                    [Xo, Yo] = meshgrid(xo, xo);
                    pupilKeck = interp2(Xi, Yi, pupilKeck, Xo, Yo);
                    pupilKeck = pupilKeck./max(pupilKeck(:));
                    pupilKeck(pupilKeck < 0.05) = 0;
                    otfTel = xcorr2(pupilKeck,pupilKeck)/(sum(pupilKeck(:)));
                    xii = linspace(0,1,length(otfTel));
                    [Xii, Yii] = meshgrid(xii, xii);
                    otfTel = interp2(Xii,Yii,otfTel, Xo, Yo);
                    %--------------------------------------------------------------
                case 'hScaoElt'
                    %--------------------------------------------------------------
                    % Harmoni-SCAO-like
                    d_tel = telescope(37);
                    Ts = 2e-3;
                    tau = 2*Ts;
                    resolution = 2*74+1; % needs be at least 2*nLenlets+1
                    noiseVar = 0.0524; % See Keck
                    nTh_ = 1;
                    fao = spatialFrequencyAdaptiveOptics(d_tel,d_atm,74,noiseVar,CTYPE*0.5,Ts,tau,resolution,filterType,mvFlag);
            end
            
            
            if strcmp(atmCase,'Jolissaint')
                fao.src = src;
            end
            %--------------------------------------------------------------
            % VAR DEFINITIONS
            %--------------------------------------------------------------
            % compute pixelScale from resolution and correctable spatial BW
            pixelScale  = fao.atm.wavelength*(fao.nActuator-1)/fao.tel.D/resolution;
            pixelScaleInMas = constants.radian2arcsec*1e3*pixelScale;
            alpha       = constants.radian2arcsec*pixelScale*((0:resolution)-resolution/2); % arcsec on sky :: pixelScale = n lambda/D /nPts, n=nLenslets=fc*2*D
            fx          = fao.fx;
            fy          = fao.fy;
            alpha       = fx(1,:);
            
            %--------------------------------------------------------------
            % FACE-ON PSD PLOTS
            %--------------------------------------------------------------
            figure
            subplot(2,2,1)
            CLIM = ([0 1.0]);
            imagesc(alpha*fao.nTimes,alpha*fao.nTimes,fittingPSD(fao,fx,fy),CLIM), title('fitting PSD','fontweight','bold','fontsize',20)
            colorbar
            ylabel('\kappa_y [m^{-1}]','fontsize',20,'fontweight','bold')
            %xlabel('\kappa_x [m^{-1}]','fontsize',20,'fontweight','bold')
            set(gca,'FontSize',18)
            subplot(2,2,2)
            CLIM = ([0 0.25]);
            imagesc(alpha,alpha,noisePSD(fao,fx,fy),CLIM),  title('noise PSD','fontweight','bold','fontsize',20)
            colorbar
            %ylabel('\kappa_y [m^{-1}]','fontsize',14,'fontweight','bold')
            %xlabel('\kappa_x [m^{-1}]','fontsize',14,'fontweight','bold')
            set(gca,'FontSize',18)
            subplot(2,2,3)
            CLIM = ([0 1.0]);
            imagesc(alpha,alpha,aliasingPSD(fao,fx,fy),CLIM),  title('aliasing PSD','fontweight','bold','fontsize',20)
            colorbar
            ylabel('\kappa_y [m^{-1}]','fontsize',20,'fontweight','bold')
            xlabel('\kappa_x [m^{-1}]','fontsize',20,'fontweight','bold')
            set(gca,'FontSize',18)
            subplot(2,2,4)
            %CLIM = ([0 2.5]);
            imagesc(alpha,alpha,anisoServoLagPSD(fao,fx,fy)),  title('aniso-servo-lag PSD','fontweight','bold','fontsize',20)
            colorbar
            %ylabel('\kappa_y [m^{-1}]','fontsize',20,'fontweight','bold')
            xlabel('\kappa_x [m^{-1}]','fontsize',20,'fontweight','bold')
            set(gca,'FontSize',18)
            colormap('hot')
            %varFit = trapz(fy(:,1),trapz(fx(1,:),fittingPSD(fao,fx,fy),2))
            varFit = fao.varFitting;
            
            %varNoise = trapz(fy(:,1),trapz(fx(1,:),noisePSD(fao,fx,fy),2));
            varNoise = fao.varNoise;
            
            %varAliasing = trapz(fy(:,1),trapz(fx(1,:),aliasingPSD(fao,fx,fy),2));
            varAliasing = fao.varAliasing;
            
            %varAnisoServoLag = trapz(fy(:,1),trapz(fx(1,:),anisoServoLagPSD(fao,fx,fy),2));
            varAnisoServoLag = fao.varServoLag;
            
            varTotal = varFit + varNoise + varAliasing + varAnisoServoLag;
            
            psdTotal = anisoServoLagPSD(fao,fx,fy) + aliasingPSD(fao,fx,fy) + noisePSD(fao,fx,fy);
            figure
            imagesc(alpha,alpha,log10(psdTotal)),  title('in-band redidual PSD','fontsize',18)
            colormap('hot')
            ylabel('\kappa_y [m^{-1}]','fontsize',14,'fontweight','bold')
            xlabel('\kappa_x [m^{-1}]','fontsize',14,'fontweight','bold')
            colorbar
            %--------------------------------------------------------------
            % ERROR BUDGET PLOTS
            %--------------------------------------------------------------
            rad2nmRms     = @(x) sqrt(x)/2/pi*500;
            fprintf('RMS [nm]: \n')
            fprintf('High frequency:     %g \n',rad2nmRms(varFit))
            fprintf('Aliasing:           %g \n',rad2nmRms(varAliasing))
            fprintf('Servo lag:          %g \n',rad2nmRms(varAnisoServoLag))
            fprintf('Noise:              %g \n',rad2nmRms(varNoise))
            fprintf('Total w/o fitting:  %g \n',rad2nmRms(varNoise + varAnisoServoLag + varAliasing))
            fprintf('Total:              %g \n',rad2nmRms(varNoise + varAnisoServoLag + varAliasing + varFit))
            fprintf('--------------------------------\n')
            fprintf('Fitting theoretical:%g \n',sqrt(0.23*(d_tel.D/(fao.nActuator-1)/d_atm.r0)^(5/3)*(0.5e-6/2/pi)^2*1e18))
            fprintf('Aliasing theoretical:%g \n',sqrt(0.073*(d_tel.D/(fao.nActuator-1)/d_atm.r0)^(5/3)*(0.5e-6/2/pi)^2*1e18))
            fprintf('--------------------------------\n')
            fprintf('Strehl-ratio (sci-band) through Marechal approximation: %1.2f %% \n', 100*exp(-varTotal*(0.5/1.65)^2))
            
            %--------------------------------------------------------------
            % SLAB PSD PLOTS
            %--------------------------------------------------------------
            %
            % * High frequency (or fitting) error
            % * Noise error
            % * Aliasing error
            % * Servo lag error
            %
            
            % Scale factor to convert from PSD in rad^2/(m^(-2)) to nm^2/(m^(-2))
            sc = fao.atm.wavelength^2/(2*pi)^2*(1e9)^2;
            
            % Calculating the 4 PSDs
            wHF = sc*fittingPSD(fao,fx,fy);
            wNS = sc*noisePSD(fao,fx,fy);
            wAL = sc*aliasingPSD(fao,fx,fy);
            %wASo = sc*(servoLagPSD(fao,fx,fy) + anisoplanatismPSD(fao,fx,fy));
            wAS = sc*(anisoServoLagPSD(fao,fx,fy));
            
            %%PSDs: x profile
            % Here we plot a profile of the PSD, for fy = 0.  This plot should be
            % compared with Fig 4 in vlt9layers 2010.
            
            % Frequency vector for HF PSD
            [fxHF,fyHF] = freqspace(resolution*fao.nTimes,'meshgrid');
            fxHF = fxHF*fao.fc*fao.nTimes;
            fyHF = fyHF*fao.fc*fao.nTimes;
            
            % Plotting different PSDs
            figure
            wAS(ceil(resolution/2),ceil(resolution/2)) = 0;
            semilogy(fx(1,:),wNS(ceil(resolution/2),:),'k',fx(1,:),wAL(ceil(resolution/2),:),'k--',...
                fx(1,:),wAS(ceil(resolution/2),:),'k-+','LineWidth',2)
            hold on;
            semilogy(fxHF(1,:),wHF(ceil(resolution*fao.nTimes/2),:),'k:','LineWidth',2)
            set(gca,'YLim',[10^(-1) 10^6])
            legend('noise','aliasing','aniso-servo-lag','high-frequency')
            set(gca,'XLim',[-1 1]*fao.fc*2)
            set(gca,'YLim',[10e-1 10e4])
            title('wave-front error PSD','fontweight','bold','fontsize',20)
            ylabel('[nm^2/m]','fontsize',18,'fontweight','bold')
            xlabel('\kappa_x [m^{-1}]','fontsize',18,'fontweight','bold')
            grid on
            pbaspect([ 1.618 1 1])
            xtick = get(gca,'XTick');
            set(gca,'FontSize',18)
            set(gca,'XTick', xtick)
            %%
            %--------------------------------------------------------------
            % PSF AND PSD (=IDEAL POST-CORONO IMAGES) PLOTS
            %--------------------------------------------------------------
            %fao.atm.wavelength = photometry.H;
            if exist('otfTel','var')
                fao.image(resolution, pixelScaleInMas*(1.65/0.5),otfTel)
            else
                fao.image(resolution, pixelScaleInMas*(1.65/0.5));
            end
            
            if iseven(size(fao.psf,1))
                o = linspace(-1,1,fao.psfResolution)*pixelScaleInMas*(1.65/0.5)/1000*fao.psfResolution/2;
            end            
            
            
            % radially-averaged PSF
            figure(100)
            r = radial(fao.psf);

            semilogy(o,[r(end:-1:1) r],'LineWidth',2)
            hold on
            title('Radially-averaged PSF','fontweight','bold','fontsize',16)
            xlabel('[arcsec]','fontsize',14,'fontweight','bold')
            ylabel('[a.u.]','fontsize',14,'fontweight','bold')
            axis tight
            grid on
            pbaspect([ 1.618 1 1])
            
            % face-on PSD/post-corono image
            psd = fao.powerSpectrumDensity(fx,fy);
            df = fy(2) - fy(1);
            psd = psd*df^2; % this normalisation assumes SR=1
            figure
            imagesc(o,o,log10(psd)) 
            title('post-coronagraphic image','fontweight','bold','fontsize',22)
            xlabel('[arcsec]','fontsize',20,'fontweight','bold')
            ylabel('[arcsec]','fontsize',20,'fontweight','bold')
            colormap('bone'), axis square, set(gca,'FontSize',18)
            ytick = get(gca,'YTick');
            set(gca,'XTick',ytick);
            % radially-averaged PSD/post-corono image
            figure(200)
            psd = radial(psd);% this normalisation assumes SR=1
            semilogy(o,[psd(end:-1:1) psd])
            title('Post-coronographic radial profile','fontweight','bold','fontsize',16)
            xlabel('[arcsec]','fontsize',20,'fontweight','bold')
            hold on
            grid on
            pbaspect([ 1.618 1 1])
            axis tight
            
        end
        
    end
    
end
