classdef slopesGenerator < handle
    
    properties (SetObservable=true)
        tel;
        atm;
        wfs;
        ast;
        nFrames;
        nSrc;
        slopes;
        PSDsx;
        PSDsy;
        checkFlag;
    end
    
    methods
        
        %CONSTRUCTOR
        function obj = slopesGenerator(tel,atm,wfs,ast,varargin)
            
            
             %Checking inputs
            inputs = inputParser;
            inputs.addRequired('tel',@(x) isa(x,'telescope') );
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann'));
            inputs.addRequired('ast',@(x) isa(x,'source'));
            inputs.addParameter('nFrames',1024, @isnumeric );
            inputs.addParameter('checkFlag',0, @isnumeric );
            inputs.parse(tel,atm,wfs,ast,varargin{:});
             %Instantiation
            obj.tel       = tel;
            obj.atm       = atm;
            obj.wfs       = wfs;
            obj.ast       = ast;
            obj.nFrames   = inputs.Results.nFrames;
            obj.checkFlag = inputs.Results.checkFlag;
            obj.nSrc      = length(obj.ast);
        end
        
        %% %%%%%%%%%%%%%%%%%
        % GENERATOR FUNCTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = generateOpenLoopSlopes(obj)           
            %Linking atmosphere to telescope
            obj.tel   = obj.tel + obj.atm;
            s         = zeros(obj.wfs.nSlope,obj.nSrc,obj.nFrames);
            nL        = obj.wfs.lenslets.nLenslet;
            if obj.checkFlag               
                obj.PSDsx = zeros(nL,nL,obj.nSrc);
                obj.PSDsy = zeros(nL,nL,obj.nSrc);
                wS        = mywindow('hann',nL);
                s2        = sum(sum(wS.^2));
                xmap      = zeros(nL,nL,obj.nSrc);
                ymap      = zeros(nL,nL,obj.nSrc);
            end

            hwait = waitbar(0,'Generating slopes...');
            for k = 1:obj.nFrames
                +obj.tel;
                %Propagate
                obj.ast   = obj.ast.*obj.tel*obj.wfs;
                %Getting slopes
                s(:,:,k)     = obj.wfs.slopes;
                % Getting PSDs
                if obj.checkFlag                    
                    xmap(obj.wfs.validLenslet,:) = s(1:end/2,:,k);
                    ymap(obj.wfs.validLenslet,:) = s(1+end/2:end,:,k);
                    obj.PSDsx = obj.PSDsx + fftshift(abs(fft2(wS.*xmap)).^2);
                    obj.PSDsy = obj.PSDsy + fftshift(abs(fft2(wS.*ymap)).^2);
                end
                waitbar(k/obj.nFrames);
            end
            close(hwait);
            % Normalization
            if obj.checkFlag
                obj.PSDsx = obj.PSDsx./s2/obj.nFrames;
                obj.PSDsy = obj.PSDsy./s2/obj.nFrames;
            end
            
            % Putting slopes into arcsec
            D          = obj.tel.D;
            d          = D/nL;
            ps         = obj.ast(1).wavelength/(2*d*obj.wfs.lenslets.nyquistSampling);
            obj.slopes = s.* ps*constants.radian2arcsec.*obj.wfs.slopesUnits;
            obj.slopes = reshape(obj.slopes,obj.wfs.nSlope*obj.nSrc,obj.nFrames);
        end
        
        function checkingSlopesStats(obj)
                                               
            % Getting the PSD from FOURIER
            nL   = obj.wfs.lenslets.nLenslet;
            d    = obj.tel.D/nL;
            nPx  = obj.wfs.camera.resolution(1);
            nPxL = nPx/nL;

            Dext = obj.tel.D/nL*nL;
            kx = 1/(Dext)*(-floor(nL/2):ceil(nL/2-1));
            [~,WSx,WSy,WSxA,WSyA] = fourierReconstructor.createPhasePSD(kx,nPxL-1,d,nPxL...
                ,obj.atm.r0,obj.atm.L0);
            
            n0  = ceil((nL+1)/2);
            idx = find(abs(kx)>abs(kx(n0+1)));
            WSx = WSx/sum(sum(WSx(idx,idx)));
            i=1;
            PSD = obj.PSDsx(:,:,i)/sum(sum(obj.PSDsx(idx,idx,i)));
            
            ps  = obj.ast(1).wavelength/(2*d*obj.wfs.lenslets.nyquistSampling);
            PSD = PSD/( ps*constants.radian2arcsec.*obj.wfs.slopesUnits);
            % DISPLAY
            figure
            loglog(kx,PSD(n0,:))
            hold on;
            loglog(kx,WSx(n0,:))
            legend('Simulation','Fourier');
        end
        
        function obj = generateFakeSlopes(obj)
            
            tel  = obj.tel;
            atm  = obj.atm;
            ngs  = obj.gs;
            wfs  = obj.wfs;
            %Parameters
            nL   = obj.wfs.lenslets.nLenslet;
            D    = obj.tel.D;
            d    = D/nL;
            nPx  = obj.tel.resolution;          
            fS   = 1/D;                     
            %Spatial frequencies 
            kx      = 1/(D)*(-floor(nPx/2):ceil(nPx/2-1));
            [KX,KY] = meshgrid(kx,kx);
            PSD     = phaseStats.spectrum(hypot(KX,KY),obj.atm);%rd^2/m^2
            
            hwait = waitbar(0,'Generating slopes...');
            for k = 1:obj.nFrames
                % Add random noise to phase and amplitude of each frequency component
                n = fft2(randn(nPx,nPx) + 1i*randn(nPx,nPx))./nPx;
                % Use inverse FFT to calculate the phase screen for the total phase and
                phi(:,:,k) = real(fS*ifft2(fftshift(sqrt(PSD)).*n)).*nPx^2;
                phi(:,:,k) = phi(:,:,k) - mean(phi(:,k)).*ones(nPx);
                %Propagate
                ngs        = ngs.*tel;
                ngs.phase  = phi(:,:,k);
                ngs        = ngs*wfs;
                %Getting slopes
                s(:,k)     = wfs.slopes;
                waitbar(k/obj.nFrames);
            end
            close(hwait);
            %Putting slopes into arcsec
            ps         = ngs.wavelength/(2*d*wfs.lenslets.nyquistSampling);
            obj.slopes = s.* ps*constants.radian2arcsec.*wfs.slopesUnits;
        end
                   
    end
    
    
    
end
   