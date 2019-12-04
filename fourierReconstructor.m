classdef fourierReconstructor < handle
    %% fourierReconstructor 
    % Wavefront reconstructor computed in the fourier domain, using FFTs to
    % reconstruct the wave-front and get the DM commands
    %
    % Required inputs: wfs, tel, atm
    %
    % Optional inputs       | Defaults
    %----------------------------------------------------------------------
    % dmFunction            | [] (i.e. no DM)
    % stroke                | 0 (must specify a stroke if dmFunction is
    %                       | defined)
    % extension             | 3
    % extensionMethod       | 'Gerchberg'
    % GerchbergFilter       | {}        (if method =/= 'Gerchberg')
    %                       | 'filter'  (i.e. same as recon filter)
    % GerchbergIterations   | 6
    % filter                | 'exact'
    % modeRemoval           | [0 0 0]
    %                       | A: global waffle
    %                       | B: local waffle
    %                       | C: piston
    % interpolationFactor   | nPx/nLenslet (i.e. goes from measurement to
    %                       | full resolution.
    % 
    
    
    
    properties
        tag = 'fourierReconstructor';
        % Wavefront sensor (Shack-Hartmann/ pyramid object)
        wfs;
        % Telescope object
        tel;
        % Atmosphere object
        atm;
        % Parameters (resolution, no. of pixels etc.)
        params;
        % Method for extending data
        extensionMethod;
        % Extension of data (no. of pixels)
        extension;
        % Enforce periodicity when slopes are extended (default = on)
        periodicCond;
        % Filter used for reconstrction
        filter;
        % Coefficient for noise within filter (i.e. gamma in literature)
        kNoise;
        % Variance of the noise.
        noiseVar;
        % Coefficients for aliasing (i.e. factors of aliased PSD included
        % in the filter)
        kAliasing;
        % Shifts to increase corespondence between model and measurement.
        % Should be a vector [xshift yshift]
        shifts;
        % Mode removal
        modeRemoval;
        doModeRemoval;
        % Shannon Interpolation factor (for the interpolation of low
        % resolution phase to high resolution phase)
        interpolationFactor
        % Valid actuators (for padded data)
        validActuator;
        % Valid lenslets (for padded data)
        validLenslet;
        % Number of data points for reconstructed phase
        N;
        % Filters (Fourier domain) for x and y slopes
        RX;
        RY;
        % wfs response to spatial frequencies
        GX;
        GY;
        % FFT filter coresponding to DM influence function
        dm;
        % Gerchberg filter used for extension (if Gerchberg is used as
        % extension method (should this just be same as filter?)
        Gerchberg;
        % Mask for Gerchberg filter
        GerchbergMask;
        
        % N.B. add dependent variables at some point (ones which will
        % update)
        
        % Orginal slopes
        xSlopes;
        ySlopes;
        % Extended slopes
        Sx;
        Sy;
        % Convert Pyramid signals to SH-like signals
        pym2sh;
        % modulation in units of lambda/D
        modulation;
        
    end
    
    methods
        
        %% Constructor
        function obj = fourierReconstructor(wfs,tel,atm,varargin)
           
            inputs = inputParser;
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann') || isa(x,'geomGrad')...
                || isa(x,'pyramid'));
            inputs.addRequired('tel',@(x) isa(x,'telescope'));
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addOptional('dm',[],@(x) isa(x,'deformableMirror'));
            inputs.addOptional('extension',3,@(x) isnumeric(x) || ischar(x));
            inputs.addOptional('extensionMethod','simple',@ischar);
            inputs.addOptional('periodicCond',1,@isnumeric);
            inputs.addOptional('GerchbergFilter',[],@ischar);
            inputs.addOptional('GerchbergIterations',6,@isnumeric);
            inputs.addOptional('GerchbergMask',[],@isnumeric);
            inputs.addOptional('filter','exact',@ischar);
            inputs.addOptional('kNoise',[],@isnumeric);
            inputs.addOptional('noiseVariance',0,@isnumeric);
            inputs.addOptional('kAliasing',[],@isnumeric);
            inputs.addOptional('shifts',[],@isnumeric);
            inputs.addOptional('modeRemoval',[],@isnumeric);
            inputs.addOptional('interpolationFactor',[],@isnumeric);
            inputs.addOptional('convertPym2SH',0,@isnumeric);
            inputs.addOptional('modulation',0,@isnumeric);
            
            inputs.parse(wfs,tel,atm,varargin{:});
            
            % Required parameters
            obj.wfs                 = inputs.Results.wfs;
            obj.tel                 = inputs.Results.tel;
            obj.atm                 = inputs.Results.atm;
            
            % Parameters
            if isa(obj.wfs,'pyramid')
                obj.params.nLenslet = inputs.Results.wfs.nLenslet;
            else
                obj.params.nLenslet = inputs.Results.wfs.lenslets.nLenslet;
            end
            if isempty(inputs.Results.dm)
                obj.params.vActuator    = inputs.Results.wfs.validActuator;
            else
                obj.params.vActuator    = inputs.Results.dm.validActuator;
            end
            obj.params.vLenslet     = inputs.Results.wfs.validLenslet;
            obj.params.nPx          = inputs.Results.tel.resolution;
            obj.params.r0           = inputs.Results.atm.r0;
            obj.params.L0           = inputs.Results.atm.L0;
            obj.params.d            = inputs.Results.tel.D/obj.params.nLenslet;
            obj.extensionMethod     = inputs.Results.extensionMethod;
            obj.extension           = inputs.Results.extension;
            obj.periodicCond        = inputs.Results.periodicCond;
            obj.filter              = inputs.Results.filter;
            obj.noiseVar            = inputs.Results.noiseVariance;
            obj.shifts              = inputs.Results.shifts;
            obj.modeRemoval         = inputs.Results.modeRemoval;
            obj.pym2sh              = inputs.Results.convertPym2SH;
            obj.modulation          = inputs.Results.modulation;
            if isempty(inputs.Results.interpolationFactor)
                obj.interpolationFactor = obj.params.nPx/obj.params.nLenslet;
            else
                obj.interpolationFactor = inputs.Results.interpolationFactor;
            end
            
            % If we have an even number of lenslets this coresponds to an
            % odd number of actuator points.  We want to always use an even
            % number of points as it makes filter out frequencies such as
            % the waffle easier.
            
            % If there is DM reconstruct at actuator points, otherwise
            % reconstruct at measurement (lenslet) points.
            % N.B. for now ony consider DM with Fried geometry (i.e.
            % nLenslet+1 actuators).
            if isempty(inputs.Results.dm)
                n0                  = obj.params.nLenslet;
            else
                n0                  = obj.params.nLenslet+1;
            end
            
            
            % Always use even number of points (so if original is odd
            % we must pad by extra point).
            if n0/2~=round(n0/2)
                % Odd original griding
                obj.N               = n0 + 1 + 2*obj.extension;
                % First pad with 1 additional point (added to the
                % reconstruction grid).
                if n0==obj.params.nLenslet % I.e. no DM
                    obj.validLenslet    = logical([obj.params.vLenslet zeros(n0,1) ;zeros(1,n0+1);]);
                    obj.validActuator   = obj.validLenslet;
                else % With DM
                    obj.validLenslet    = logical([zeros(1,n0+1); zeros(n0-1,1) obj.params.vLenslet zeros(n0-1,1); zeros(1,n0+1)]);
                    obj.validActuator   = logical([zeros(1,n0+1); zeros(n0,1) obj.params.vActuator]);
                end
                
                % Then pad to extension
                obj.validActuator   = logical(fourierReconstructor.padarray(obj.validActuator,obj.extension));
                obj.validLenslet    = logical(fourierReconstructor.padarray(obj.validLenslet,obj.extension));
                
            else
                obj.N               = n0 + 2*obj.extension;
                % No additional padding for actuator array (we already have
                % an even number of points).
                
                if n0==obj.params.nLenslet % I.e. no DM
                    obj.validLenslet    = obj.params.vLenslet;
                    obj.validActuator   = obj.validLenslet;
                else
                    % Extend lenslets to no. of actuator points (if we have
                    % a DM)
                    obj.validLenslet    = logical([obj.params.vLenslet zeros(n0-1,1) ;zeros(1,n0);]);
                    obj.validActuator   = obj.params.vActuator;
                end
                
                % Pad actuator and lenslet arrays
                obj.validActuator   = logical(fourierReconstructor.padarray(obj.validActuator,obj.extension));
                obj.validLenslet    = logical(fourierReconstructor.padarray(obj.validLenslet,obj.extension));
                
            end
            
            if isempty(inputs.Results.GerchbergMask)
                obj.GerchbergMask = obj.validLenslet;
            else
                if n0/2~=round(n0/2)
                    if n0==obj.params.nLenslet % I.e. no DM
                        obj.GerchbergMask = logical([inputs.Results.GerchbergMask zeros(n0,1); zeros(1,n0+1)]);
                    else
                        obj.GerchbergMask = logical([zeros(1,n0+1); zeros(n0-1,1) inputs.Results.GerchbergMask zeros(n0-1,1); zeros(1,n0+1)]);
                    end
                else
                    if n0==obj.params.nLenslet
                        obj.GerchbergMask = logical(inputs.Results.GerchbergMask);
                    else
                        obj.GerchbergMask = logical([inputs.Results.GerchbergMask zeros(n0-1,1) ;zeros(1,n0);]);
                    end
                end
                obj.GerchbergMask = logical(fourierReconstructor.padarray(obj.GerchbergMask,obj.extension));
            end
                
            
            % If a noise variance is specified and no value is given for
            % the noise coefficient, it is set to 1.
            if isempty(inputs.Results.kNoise) && obj.noiseVar~=0
                obj.kNoise          = 1;
            else
                obj.kNoise          = inputs.Results.kNoise;
            end
                
            % If the aliasing filter is specified but no aliasing
            % coefficient is given it is set to 1.
            if isempty(inputs.Results.kAliasing) && strcmp(obj.filter,'antialias')
                obj.kAliasing       = 1;
            else
                obj.kAliasing       = inputs.Results.kAliasing;
            end
            
            [RX,RY,GX,GY] = fftFilters(obj);
            obj.RX                  = RX;
            obj.RY                  = RY;
            obj.GX                  = GX;
            obj.GY                  = GY;

            if strcmp(obj.extensionMethod,'Gerchberg')
                if isempty(inputs.Results.GerchbergFilter)
                    obj.Gerchberg.filter = obj.filter;
                else
                    obj.Gerchberg.filter = inputs.Results.GerchbergFilter;
                end
                obj.Gerchberg.nIter      = inputs.Results.GerchbergIterations;
            end
            
            if ~isempty(obj.modeRemoval) && (obj.modeRemoval(1)==1 || ...
                obj.modeRemoval(2)==1 || obj.modeRemoval(3)==1)
                obj.doModeRemoval = 1;
            else
                obj.doModeRemoval = 0;
            end
            
            if ~isempty(inputs.Results.dm)
                obj.dm.function         = inputs.Results.dm.modes; 
                obj.dm.modes            = inputs.Results.dm.modes.modes;
                obj.dm.filter           = fourierReconstructor.dmFilter(obj.dm.function,...
                    obj.extension,obj.params.nLenslet,obj.params.nPx);                
                obj.dm.compensation     = 1;
            else
                obj.dm.compensation     = 0;
            end
            
            % Set up empty matrices to store x/y slopes
            obj.xSlopes = zeros(obj.params.nLenslet,obj.params.nLenslet);
            obj.ySlopes = zeros(obj.params.nLenslet,obj.params.nLenslet);
            
        end    
        

        %% Wavefront reconstruction
        function [Phi,DMcoeffs,PhiDM,PhiHR] = mtimes(obj,wfs)
            
            % Returns:
            % Phi:      low resolution phase (nLenslet for no DM,
            %           nLenslet+1 with DM).
            % DMcoeffs: dm coefficients for compensation
            % PhiDM:    nPx resolution phase coresponding to phase of the
            %           dm.
            % PhiHR:    high resolution phase equivalent to Phi, where the
            %           resolution is determined by the interpolation
            %           factor
            
            % To Do: add mode removal
            
            % Can use either SH object or the wfs slopes themselves
            if isa(wfs,'shackHartmann') || isa(wfs,'pyramid')
                obj.xSlopes(obj.params.vLenslet) = wfs.slopes(1:end/2);
                obj.ySlopes(obj.params.vLenslet) = wfs.slopes(end/2+1:end);
            else
                obj.xSlopes(obj.params.vLenslet) = wfs(1:end/2);
                obj.ySlopes(obj.params.vLenslet) = wfs(end/2+1:end);
            end
           
            N = obj.N;
            nL = obj.params.nLenslet;
            mask = obj.validActuator;
            
            Sx = obj.xSlopes;
            Sy = obj.ySlopes;
            % Pad slopes to number of actuators points, with extension 
            % (always pad to even number)
            if nL/2==round(nL/2)
                % Extend slopes to the same size as the no. of actuator
                % points (i.e. recreated phase points)
                if obj.dm.compensation==1
                    Lext = obj.extension+1;
                else
                    Lext = obj.extension;
                end
                Sx = fourierReconstructor.padarray(Sx,Lext);
                Sy = fourierReconstructor.padarray(Sy,Lext);
            else
                % Extend slopes to the same as the number of actuator
                % points
                Sx = [Sx zeros(nL,1); zeros(1,nL+1)];
                Sy = [Sy zeros(nL,1); zeros(1,nL+1)];
                Sx = fourierReconstructor.padarray(Sx,obj.extension);
                Sy = fourierReconstructor.padarray(Sy,obj.extension);
            end
            
            
            % Extend slopes to edges of grid (with given method)
            [Sx,Sy] = extendSlopes(Sx,Sy,obj);

            % Store extended slopes
            obj.Sx = Sx;
            obj.Sy = Sy;
            
            % FFT of slopes
            fft_Sx = 1/N^2 * fftshift(fft2(Sx,N,N));
            fft_Sy = 1/N^2 * fftshift(fft2(Sy,N,N));
            
            % Use filters to create FFT of phase
            fft_Phi = obj.RX.*fft_Sx + obj.RY.*fft_Sy;
            % Remove piston
            fft_Phi(floor(N/2+1),floor(N/2+1)) = 0;
            
            % Modal removal (in Fourier space)
            if obj.doModeRemoval==1
                % Remove selected mdoes
                fft_Phi = fourierReconstructor.removeModes(fft_Phi,N,ones(N,N),obj.modeRemoval,'fourier');
%                 Phi = (N^2 * ifft2(ifftshift(fft_Phi)));
%                 Phi = fourierReconstructor.removeModes(Phi,N,ones(N,N),obj.modeRemoval,'direct');
%                 fft_Phi = 1/N^2 * fftshift(fft2(Phi));
            end   
            
            % Inverse Fourier transform to get LR phase
            Phi = real(N^2 * ifft2(ifftshift(fft_Phi)));
            
            % Shannon interpolation to return high resolution phase
            if obj.dm.compensation == 0
                nSh = N * obj.interpolationFactor;
            else
                nSh = (N-1) * obj.interpolationFactor;
            end
            fft_PhiHR = zeros(nSh,nSh);
            % Coordinates for central frequency (high res)
            nSh0 = floor(nSh/2+1);
            % Central FFT is equivalent to FFT of phi
            fft_PhiHR(nSh0-floor(N/2):nSh0+ceil(N/2-1),nSh0-floor(N/2):nSh0+ceil(N/2-1)) = fft_Phi;
                        
            % Use inverse FFT to get high resolution phase
            PhiHR = nSh^2 * ifft2(fftshift(fft_PhiHR),nSh,nSh);
            PhiHR = real(PhiHR);
            
            % Impose pupil (i.e. only return valid measurements)
            Phi = mask.*Phi;
                        
            % Remove extended results
            if nL/2==round(nL/2) && obj.dm.compensation==1
                Phi = Phi(2+obj.extension:N-obj.extension,2+obj.extension:N-obj.extension);
            else
                Phi = Phi(1+obj.extension:N-obj.extension,1+obj.extension:N-obj.extension);
            end
            
            % And remove extended results for high resolution phase
            M = obj.extension*obj.interpolationFactor;
            if nL/2==round(nL/2) && obj.dm.compensation==1
                PhiHR = PhiHR(1+obj.interpolationFactor+M:nSh-M,1+obj.interpolationFactor+M:nSh-M);
            else
                PhiHR = PhiHR(1+M:nSh-M,1+M:nSh-M);
            end
            
            % Calculate dm coefficients
            DMcoeffs = [];
            PhiDM = [];
            if obj.dm.compensation==1
               
                % To do: check DM filter -> understand and check 'aliasing
                % effect' is properly taken into account.
                coeffs = N^2 * ifft2(ifftshift(fft_Phi).*obj.dm.filter,N,N);
                % Take real part and calibrate with stroke (and factor of 2
                % to take reflection from dm into account)
                coeffs = obj.atm.wavelength/2 * real(coeffs);
                DMcoeffs = coeffs(obj.validActuator);
             
                % Equivalent phase of the deformable mirror with the given
                % coefficients
                PhiDM = obj.dm.modes * DMcoeffs;
                PhiDM = reshape(PhiDM,obj.params.nPx,obj.params.nPx);
            end
            
            % Scale phases
            Phi = pi*Phi;
            PhiHR = pi*PhiHR;
            
        end
        
        %% Calculation of filters
        function [RX, RY, GX, GY] = fftFilters(obj)

            % ------------HEADER-----------------
            % N.B. Need to edit comments etc.
            %
            % Filters:  
            %
            % Objective         ::  Implement several filters to reconstruct (LSR & MVR) the
            %                       phase in Fourier domain
            % INPUT VARS
            % filter            ::  the type of filter
            %                       'exact'
            %                       'Hudgin'
            %                       'modifiedHudgin'
            %                       'Freischlad86'
            %                       'Poyneer02'
            % N                 ::  The number of computational points
            % hshift, vshift    ::  The horizontal and vertical shifts
            %
            % OUTPUT VARS
            % CX                ::  Filter for the X slopes
            % CY                ::  Filter for the Y slopes
            % div               ::  The denominator, div = abs(CX)^2+abs(CY)^2
            % CXo               ::  Filter for the X slopes without being divided by
            %                       div
            % CYo               ::  Filter for the Y slopes without being divided by
            %                       div
            %Created by         ::  Carlos Correia
            %Creation date      ::  16/02/2007
            %Change Record:     ::  CCorreia, July09
            %                       - Added two output variables, to have the x and y
            %                       filters without being divided by the denominator
            %                       CCorreia, Jun09
            %                       - horizontal and vertical shifts were incorrectly computed. Added the
            %                       frequency vectors fx and fy with the correct numbering to allow for
            %                       fractional shifts

            %                       JTeixeira, July14
            %                       - Added option to Wiener filter;
            %                       - input parameter changed:
            %                           params :: structure containing parameters:
            %                              . filter
            %                              . N
            %                              . r0 & L0 to estimate PSD's
            %                              . kn & ka to the regularization parameter
            %                                coeficients, kn for noise, ka for aliasing;
            %
            %
            %       CORRECTIONS TO DO: Extra-alignment in Hudgin & Southwell is
            %                              incorrect, each filters should be shifted in
            %                              both directions x and y.

            % ------------HEADER END----------------

            N = obj.N;
            d = obj.params.d;
            % No. of pixels per lenslet
            nPxL = (obj.params.nPx/obj.params.nLenslet);
            % Averaging length (d-size of 1 pixel)
            d_av = d - d/nPxL;
            
            if isempty(obj.shifts)
                xshift = 0;
                yshift = 0;
            else
                xshift = obj.shifts(1);
                yshift = obj.shifts(2);
            end

            % Diameter over defined phase points (plus extension)
            D = N*d;
            
            % Spatial frequencies
            kx = 1/D * (-floor(N/2):ceil(N/2-1));

            ky = kx;
            [KX,KY] = meshgrid(kx,ky);


            % PSDs for noise (i.e. for Wiener filters)
            if ~isempty(obj.kNoise)
                
                % Calculate phase PSD (in band -> no aliasing)
                m = 0;
                Wphi = fourierReconstructor.createPhasePSD(kx,m,d,nPxL,obj.params.r0,obj.params.L0);
                
                % Calculate white noise PSD
                Wnoise = fourierReconstructor.createNoisePSD(kx,obj.noiseVar);
                
                % Noise PSD normalised against phase PSD
                Rnoise = obj.kNoise * Wnoise./Wphi;
            else
                Rnoise = 0;
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %--------------  filters   ----------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(obj.filter,'Hudgin')|| strcmp(obj.filter,'hudgin')
                % -------------------------------------------
                % ---------------------------------
                %----------  modified Hudgin  -----
                % ---------------------------------

                % Gx/Gy
                GX = exp(2*1i*pi*d*KX)-1;
                GY = exp(2*1i*pi*d*KY)-1;

                % Additional shift to increase correlation between model 
                % and SH.  User defined shifts more variable (i.e. can 
                % change depending on where phase is reconstructed and
                % number of pixels per lenslet).
                Ex = exp(-1i*pi*d*KY).*exp(1i*pi*d*(xshift*KX+yshift*KY));
                Ey = exp(-1i*pi*d*KX).*exp(1i*pi*d*(xshift*KX+yshift*KY));
                %E = 1;
                
                % Calculate reconstruction filters
                RX = conj(GX).*Ex./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*Ey./(abs(GX).^2+abs(GY).^2+Rnoise);
                
            elseif strcmp(obj.filter,'exact') || strcmp(obj.filter,'Exact')
                % --------------------------------
                %----------   Exact  --------------

                GX = 2*1i*pi*d*KX .* sinc(d_av*KX) .* sinc(d*KY) ...
                    .* exp(1i*pi*d*(KX+KY));
                GY = 2*1i*pi*d*KY .* sinc(d*KX) .* sinc(d_av*KY) ...
                    .* exp(1i*pi*d*(KX+KY));
                
                % Shifts for better corespondence between with SH.  For
                % exact filter will generally only be used if we are
                % changing where reconstruction happens (i.e. not on the
                % actuators).
                E = exp(1i*pi*d*(xshift*KX+yshift*KY));
                
                % Calculate reconstruction filters
                RX = conj(GX).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*E./(abs(GX).^2+abs(GY).^2+Rnoise);


            elseif strcmp(obj.filter,'Freischlad86')
                % ---------------------------------------
                %----------  Freischlad86  --------------
                % ---------------------------------------
                
                % N.B. doesn't work for dm coeffs????
                
                GX = exp(1i*2*pi*d*KX) - 1;
                GY = exp(1i*2*pi*d*KY) - 1;
                
                E = exp(1i*pi*d*(xshift*KX+yshift*KY));
                
                % Calculate reconstruction filters
                RX = conj(GX).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                
            elseif    strcmp(obj.filter,'Poyneer02') || strcmp(obj.filter,'Fried')|| strcmp(obj.filter,'fried')
                % -------------------------------------------
                % -------------------------------------------
                %----------  Poyneer02, exp. 13  ------------
                % -------------------------------------------
                % theoretical filter for the Fried geometry (See Poyneer02, 
                % exp. 13) THE APPLICATION OF LISA'S FILTER TO THE FRIED 
                % GEOMETRY AFTER LINEAR RECOMBINATION OF GRIDS (EXP. 18 AND
                % 19) GIVE VERY ACCURATE RESULTS, CLOSE OR EVEN BETTER THAN
                % THOSE OBTAINED WITH FREISCHLAD'S FILTER...

                % TO DO: add back in LISA's addition to this filter (waffle
                % removal) done
                % TO DO: check shifts

                GX = (exp(2*1i*pi*d*KX)-1) * (1/2) .* (1+exp(2*1i*pi*d*KY));
                GY = (exp(2*1i*pi*d*KY)-1) * (1/2) .* (1+exp(2*1i*pi*d*KX));
                
                % Shift
                E = exp(1i*pi*d*(xshift*KX+yshift*KY));
                
                % Calculate reconstruction filters
                RX = conj(GX).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                
                if N/2==round(N/2)
                    % If even number of points remove waffle mode
                    RX(1,1) = 0;
                    RY(1,1) = 0;
                end

            elseif strcmp(obj.filter,'DCT')
        %         % TO DO: what is this filter?
        %         % re-write in my format
        %         for k=1:N
        %             for l=1:N
        %                 CX(k,l) =  (exp(-imath*pi*(k-1)/N)-1); %2* (sin(pi*(k-1)/(2*N)));
        %                 CY(k,l) = (exp(-imath*pi*(l-1)/N)-1); %2* (sin(pi*(l-1)/(2*N)));%
        %                 div(k,l) = (4*((sin(pi*(k-1)/(2*N)))^2 + (sin(pi*(l-1)/(2*N)))^2))^(-1);
        %             end
        %         end


            elseif strcmp(obj.filter,'Southwell')|| strcmp(obj.filter,'southwell')
                %----------  Southwell_80  --------------

                GX = exp(2*1i*pi*d*KX)-1;
                GY = exp(2*1i*pi*d*KY)-1;
                
                % Shift and slope averaging
                E = exp(-1i*pi*d*(KX+KY)).*exp(1i*pi*d*(xshift*KX+yshift*KY));
                SX = 1/2*(1+exp(-2*1i*pi*d*KX));
                SY = 1/2*(1+exp(-2*1i*pi*d*KY));
                
                % Calculate reconstruction filters
                RX = conj(GX).*conj(SX).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*conj(SY).*E./(abs(GX).^2+abs(GY).^2+Rnoise);


            elseif strcmp(obj.filter,'antialias')
                % -------------------------------------------
                % -------------------------------------------
                % --------------------------------------
                %   ANTI-ALIASING filter
                % --------------------------------------
                
                GX = 2*1i*pi*d*KX .* sinc(d_av*KX) .* sinc(d*KY) ...
                    .* exp(1i*pi*d*(KX+KY));
                GY = 2*1i*pi*d*KY .* sinc(d*KX) .* sinc(d_av*KY) ...
                    .* exp(1i*pi*d*(KX+KY));


                % PSD of the phase (and X and Y slopes) in band (no
                % aliasing)
                [Wphi,WSx0,WSy0] = fourierReconstructor.createPhasePSD(kx,0,d,nPxL,obj.params.r0,obj.params.L0);

                % Fold in sould be N-1 where N is the multiple of nLenslet
                % for the full resolution (in simulation)
                m = nPxL-1;
                [~,WSx,WSy] = fourierReconstructor.createPhasePSD(kx,m,d,nPxL,obj.params.r0,obj.params.L0);
                % The slope PSDs with only the aliasing component
                RAx = obj.kAliasing.*(WSx-WSx0)./Wphi;
                RAy = obj.kAliasing.*(WSy-WSy0)./Wphi;
                
                % Shift
                E = exp(1i*pi*d*(xshift*KX+yshift*KY));
                
                % Calculate reconstruction filters
                RX = conj(GX).*E./(abs(GX).^2+abs(GY).^2+Rnoise+RAx);
                RY = conj(GY).*E./(abs(GX).^2+abs(GY).^2+Rnoise+RAy);
                
            elseif strcmp(obj.filter,'pyramid')
                
                sgnKX = KX./abs(KX);
                sgnKX(isnan(sgnKX)) = 1;
                sgnKY = KY./abs(KY);
                sgnKY(isnan(sgnKY)) = 1;
                
                GX = 1i*sgnKX./sinc(1*d/nPxL*KX);%./sinc(d/nPxL*KY);
                GY = 1i*sgnKY./sinc(1*d/nPxL*KY);%./sinc(d/nPxL*KY);
                
                % Shift
                E = exp(1i*pi*d*(xshift*KX+yshift*KY));
                
                RX = conj(GX).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                
                %RX = 1./GX;
                %RY = 1./GY;
                
                % NEW CODE TO INCLUDE MODULATION 
                nL = obj.wfs.nLenslet;
                m = obj.modulation;
                umod = 1/(2*d)/(nL/2)*m;
                idx = abs(kx) > umod;
                F = zeros(size(kx));
                F(idx) = 1i*sign(kx(idx));
                idx = abs(kx) <= umod;
                F(idx) = 2*1i/pi*asin(kx(idx)/umod);
                F(isnan(F)) = 0;
                GX = ones(N,1)*F;
                GY = GX.';
                
                % Pixel averaging
                GX = GX.*sinc(KX*d).*sinc(KY*d);
                GY = GY.*sinc(KX*d).*sinc(KY*d);
                
                % Shift by half a lenslet (for fried geometry)
                GX = GX.*exp(1i*pi*d*(KX+KY));
                GY = GY.*exp(1i*pi*d*(KX+KY));

                RX = conj(GX).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RY = conj(GY).*E./(abs(GX).^2+abs(GY).^2+Rnoise);
                RX(isnan(RX)) = 0;
                RY(isnan(RY)) = 0;
                
            elseif strcmp(obj.filter,'custom')
                
                fprintf('Custom filter chosen, measuring filter ...\n');
                [RX,RY,GX,GY] = measureFilter(obj);
                
            else
                
                fprintf('No valid filter chosen ...\n');
                
            end

            % Convert slopes from Pyramid to shack-hartmann like slopes
            if obj.pym2sh==1
                
                sgnKX = KX./abs(KX);
                sgnKX(isnan(sgnKX)) = 1;
                sgnKY = KY./abs(KY);
                sgnKY(isnan(sgnKY)) = 1;
                
                gspX = 2*1i*pi*d*KX./(1i*sgnKX);
                gspY = 2*1i*pi*d*KY./(1i*sgnKY);
                RX = RX.*gspX;
                RY = RY.*gspY;
            end
            
            
            % Remove piston
            RX(floor(N/2+1),floor(N/2+1)) = 0;
            RY(floor(N/2+1),floor(N/2+1)) = 0;

        end
        
        %% Function to measure wfs filter
        function [Rx,Ry,Gx,Gy] = measureFilter(obj,ngs,amp,display)
            
            % To do: check capability for extended slopes
            % To do: change filter to go from slope measurements to
            % actuator phase (nLenslet -> nLenslet +1).
            
            %obj.filter = 'measured';
            
            % Number of points
            %N = obj.N;
            N = obj.params.nLenslet;
            nPx = obj.params.nPx;
            wfs = obj.wfs;
            tel = obj.tel;
            
            % x (n) and y (m) vectors in indices (high resolution phase)
            [n,m] = meshgrid(0:nPx-1,0:nPx-1);
            
            %[n,m] = meshgrid(0.5:nPx-0.5,0.5:nPx-0.5); % this sampling
            %vector samples the modes at half-pixel giving correct WFS
            %response to waffle modes (only cos, no sin). However, use
            %other definition for compatibility with FFT
            
            
            % Define x (l) and y (k) frequency vectors.
            l = -floor(N/2):floor((N-1)/2);
            k = -floor(N/2):floor((N-1)/2);
            
            % Central frequency coordinate
            iF0 = ceil((N+1)/2);
            
            % Number of points corresponding to high resolution phase for
            % actuator points (i.e. nPx is high resolution points for phase
            % over number of lenslets)
            %nM = nPx*N/(N-1);
            nM = nPx*N/N;
            
            Gx = zeros(N,N);
            Gy = zeros(N,N);
            
            % Initiate slope vectors
            Sx = zeros(N,N);
            Sy = zeros(N,N);
            
            % Cycle through each spatial frequency, propagating this
            % through the wfs
            if display
                figure
            end
%             for i=1:N
%                 for j=1:N
%                     C_kl = Dkl/N * cos(2*pi/nM*(-l(i)*n-k(j)*m));
%                     %imagesc(C_kl)
%                     imagesc(fftshift(real(fft2(C_kl))))
%                     pause(0.2)
%                     drawnow
%                     
%                 end
%             end
            for i=1:N
                for j=1:N
                    
                    
                    % Normalisation factor
                    Dkl = sqrt(2);
                    
                    % Eigen function (cos term)
                    C_kl = Dkl/N * cos(2*pi/nM*(-l(i)*n-k(j)*m));
                    %sum(C_kl(:).*C_kl(:))
                    % Propagate ngs through telescope and add custom phase.
                    % Then propagate to wfs
                    ngs = ngs.*+tel;
                    ngs.phase = amp*C_kl;
                    ngs = ngs*wfs;
                    
                    % Select the valid measurement of the slopes
                    Sx(wfs.validLenslet) = wfs.slopes(1:end/2);
                    Sy(wfs.validLenslet) = wfs.slopes(end/2+1:end);  
                    
                    % Take fourier transform of slopes
                    fftCx = (1/N)*fftshift(fft2((Sx),N,N));
                    fftCy = (1/N)*fftshift(fft2((Sy),N,N));
                    
                    % Eigen function (sine term)
                    S_kl = Dkl/N * sin(2*pi/nM*(-l(i)*n-k(j)*m));
                    
                    % Propagate ngs through telescope and add custom phase.
                    % Then propagate to wfs
                    ngs = ngs.*+tel;
                    ngs.phase = amp*S_kl;
                    ngs = ngs*wfs;
                    
                    % Select the valid measurement of the slopes
                    Sx(wfs.validLenslet) = wfs.slopes(1:end/2);
                    Sy(wfs.validLenslet) = wfs.slopes(end/2+1:end);  
                    
                    % Take fourier transform of slopes
                    fftSx = (1/N)*fftshift(fft2((Sx),N,N));
                    fftSy = (1/N)*fftshift(fft2((Sy),N,N));
                    
                    
                    fftEx = fftCx-1i*fftSx;
                    fftEy = fftCy-1i*fftSy;

                    
                    if display
                        imagesc([abs(fftEx) abs(fftEy)])
                        %imagesc(wfs.slopesMap)
                        axis tight equal
                        colorbar()
                        drawnow;
                    end
                    
                    % Record response for given frequency
                    Gx(j,i) = Dkl/amp/2 * fftEx(j,i);
                    Gy(j,i) = Dkl/amp/2 * fftEy(j,i);
                    
                end
            end
            
            Gx = pi*Gx;
            Gy = pi*Gy;
            
            % Calculate filter from response functions
            Rx = conj(Gx)./(abs(Gx).^2+abs(Gy).^2);
            Ry = conj(Gy)./(abs(Gx).^2+abs(Gy).^2);
            
            % Zero the 0th frequency component
            Rx(iF0,iF0)=0;
            Ry(iF0,iF0)=0;
            
            obj.GX = Gx;
            obj.GY = Gy;
            obj.RX = Rx;
            obj.RY = Ry;
            
            
        end
        
        %%
        function [Sx,Sy] = extendSlopes(Sx,Sy,obj)

            % Extension methods             | Description
            %-----------------------------------------------------------
            % 'simple'/'hudgin'    | 
            %
            % 'FriedClosedExtension'/       |
            % 'mesh2x2'                     |
            %
            % 
            mask = obj.validLenslet;

            if strcmp(obj.extensionMethod,'simple')
                %----------------------------------------------------------------------
                % Extension method: simple
                % Extends the slopes along x (for x slopes) and along y (for y slopes)
                % using the last value in the valid apperture
                %----------------------------------------------------------------------
                                
                for i = 1:length(mask)
                    % Extend x slopes up and down out of aperture
                    idx = find(mask(:,i));
                    Sx(1:min(idx)-1,i) = Sx(min(idx),i);
                    Sx(max(idx)+1:end,i) = Sx(max(idx),i);
                    
                    % Extend y slopes left and right out of aperture
                    idy = find(mask(i,:));
                    Sy(i,1:min(idy)-1) = Sy(i,min(idy));
                    Sy(i,max(idy)+1:end) = Sy(i,max(idy));
                end

            elseif strcmp(obj.extensionMethod,'FriedClosedExtension') ||strcmp(obj.extensionMethod,'mesh2x2')
                %-------------------------------------
                % Extension method for the FRIED geometry
                %-------------------------------------
                N = length(mask);
                masktmp = mask;
                surrounded = 1;
                %for debug purposes
                xs = Sx;
                ys = Sy;
                while surrounded > 0
                    surrounded = 0;
                    for i = 2:length(mask)-1
                        % --- X slopes ---
                        idx = find(masktmp(i,:));
                        %solve conflicts - a point that is surronded by three other inside
                        %the aperture
                        if ~isempty(idx) %if subapertures have been found!
                            point = min(idx) - 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                            point = max(idx) + 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                        end
                    end
                    %enlarge the new masktmp to account for the conflict-free
                    %subapertures
                    maskidx = find(Sx+Sy);
                    masktmp(maskidx) = 1;
                end
                for i = 2:length(masktmp)-1
                    % --- Y slopes ---
                    idx = find(masktmp(:,i));
                    Sx(1:min(idx)-1,i) = Sx(min(idx),i);
                    Sx(max(idx)+1:end,i) = Sx(max(idx),i);
                    %create the waffle-like pattern in the y-slope direction
                    if ~isempty(idx)
                        pos = 1:min(idx)-1;
                        Sy(1:min(idx)-1,i) = Sy(min(idx),i)*pattern(length(pos),'inv')';
                        pos = max(idx)+1:N;
                        Sy(max(idx)+1:end,i) = Sy(max(idx),i)*pattern(length(pos))';
                    end

                    % --- Y slopes ---
                    idy = find(masktmp(i,:));
                    Sy(i,1:min(idy)-1) = Sy(i,min(idy));
                    Sy(i,max(idy)+1:end) = Sy(i,max(idy));
                    %create the waffle-like pattern in the x-slope direction
                    if ~isempty(idy)
                        pos = 1:min(idy)-1;
                        Sx(i,1:min(idy)-1) = Sx(i,min(idy))*pattern(length(pos),'inv');
                        pos = max(idy)+1:N;
                        Sx(i,max(idy)+1:end) = Sx(i,max(idy))*pattern(length(pos));
                    end
                end

                maskidx = find(Sx+Sy);
                masktmp(maskidx) = 1;
                surrounded = 1;
                while surrounded > 0
                    surrounded = 0;
                    for i = 1:length(mask)
                        % --- X slopes ---
                        idx = find(masktmp(i,:));
                        %solve conflicts - a point that is surronded by three other inside
                        %the aperture
                        if ~isempty(idx) & length(idx) < N%if subapertures have been found!
                            point = min(idx) - 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                            point = max(idx) + 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                        end
                    end
                    %enlarge the new masktmp to account for the conflict-free
                    %subapertures
                    maskidx = find(Sx+Sy);
                    masktmp(maskidx) = 1;
                end

            elseif strcmp(obj.extensionMethod,'HudginSimpleExtension') || strcmp(obj.extensionMethod,'HudginExtension')
                %-------------------------------------
                % Extension method for the HUDGIN geometry
                %   (takes into account three-slopes subapertures)
                %-------------------------------------
                for i = 2:length(mask)-1
                    idx = find(mask(:,i));
                    Sx(1:min(idx)-1,i) = Sx(min(idx),i);

                    idy = find(mask(i,:));
                    Sy(i,1:min(idy)-1) = Sy(i,min(idy));
                    if (mask(max(idx),i+1))%if there is a right-neighbour
                        if isempty(idx)
                            idx = 1;
                        end
                        Sx(max(idx)+1:end,i) = Sx(max(idx),i) - Sy(max(idx),i) + Sy(max(idx),i+1);
                    else
                        Sx(max(idx)+1:end,i) = Sx(max(idx),i);
                    end
                    if (mask(i+1,max(idy)))%if there is a bottom-neighbour
                        if isempty(idy)
                            idy = 1;
                        end
                        Sy(i,max(idy)+1:end) = Sy(i,max(idy)) - Sx(i,max(idy)) + Sx(i+1,max(idy));
                    else
                        Sy(i,max(idy)+1:end) = Sy(i,max(idy));
                    end
                end

            elseif strcmp(obj.extensionMethod,'FriedSkewedExtension')
                %-------------------------------------
                % Extension method for the Fried geometry...
                %-------------------------------------
                [B,d] = spdiags(Sy);
                for i=1:size(B,2)
                    idx = find(B(:,i));
                    B(1:min(idx)-1,i) = B(min(idx),i);
                    B(max(idx)+1:end,i) = B(max(idx),i);
                end
                %   --- re-shape the B matrix so it resembles a circle again and not
                %   an elipse ---
                %     shift = -27:27;
                %     for i=1:size(B,2)
                %         B(:,i) = circshift(B(:,i), floor(-shift(i)/2));
                %     end
                Sy = spdiags(B,d,Sy);
                Sy = full(Sy);
                circ = 1:2:length(mask);
                circ = [circ circ];
                for i = 1:length(mask)
                    Sx(i,:) = Sx(i,length(Sx):-1:1);
                end
                [B,d] = spdiags(Sx);
                for i=1:size(B,2)
                    idx = find(B(:,i));
                    B(1:min(idx)-1,i) = B(min(idx),i);
                    B(max(idx)+1:end,i) = B(max(idx),i);
                end
                Sx = spdiags(B,d,Sx);
                Sx = full(Sx);
                for i = 1:length(mask)
                    Sx(i,:) = Sx(i,length(Sx):-1:1);
                end

            elseif strcmp(obj.extensionMethod,'Gerchberg')
                %-------------------------------------
                % Gerchberg
                %-------------------------------------
                
                [Sx, Sy] = GerchbergMethod(Sx,Sy,obj);
                
                
            elseif strcmp(obj.extensionMethod,'edge')
                %-------------------------------------
                % Edge correction (Poyneer05)
                % Requires at least 1 point of padding to work
                % N.B. produce error message if no padding
                %-------------------------------------
                for i = 1:length(mask)
                    idx = find(mask(i,:));
                    tmp = -1/2*sum(Sx(i,idx));
                    Sx(i,min(idx)-1) = tmp;
                    Sx(i,max(idx)+1) = tmp;

                    idy = find(mask(:,i));
                    tmp = -1/2*sum(Sy(idx,i));
                    Sy(min(idy)-1,i) = tmp;
                    Sy(max(idy)+1,i) = tmp;
                end
                
            elseif strcmp(obj.extensionMethod,'new')
                % Also requires at least ext = 1
                % Same idea as simple extension, but in this case the x
                % slopes are extended in x and the y slopes are extended in
                % y
                for i = 1:length(mask)
                    % Extending x slopes left to right
                    idx = find(Sx(i,:));
                    Sx(i,1:min(idx)-1) = Sx(i,min(idx));
                    Sx(i,max(idx)+1:end) = Sx(i,max(idx));
                    % Extending y slopes top to bottom
                    idy = find(Sy(:,i));
                    Sy(1:min(idy)-1,i) = Sy(min(idy),i);
                    Sy(max(idy)+1:end,i) = Sy(max(idy),i);
                end
            elseif strcmp(obj.extensionMethod,'none')
                %fprintf('Slopes not extended ...\n')
            else
                fprintf('No extension method chosen ...\n')
            end

            % EDIT 15/10/15
            % Testing if we should use a window for cyclic conditions,
            % rather than the sum method.  Also enforce in x and y for both
            % slopes.  So use square window with 1 in pupil area.  May need
            % large extension for this???
%             dEdge = 2*obj.extension;
%             fracCentre = 1 - dEdge/obj.N;
%             % Square Tukey window
%             w = mywindow('squareTukey',obj.N,fracCentre);
%             Sx = w(:,:,2).*Sx;
%             Sy = w(:,:,1).*Sy;
            
            if obj.periodicCond==1
                % Enforcing cyclic condition on slopes
                Sx(:,end) = -sum(Sx(:,1:end-1)');
                Sy(end,:) = -sum(Sy(1:end-1,:));
            end
        end
        
        %% Gerchberg extension method
        function [Sx,Sy] = GerchbergMethod(Sx,Sy,obj)
            
            % Number of iterations
            nIter = obj.Gerchberg.nIter;
            N = obj.N;
            ext = obj.extension;

            % Defining size of slopes measurement data
            mask = obj.GerchbergMask;
            [mx,my] = size(mask);

            % Inverse mask (for use later)
            imask = 1-mask;

            % Valid slopes measurements
            idx = find(mask);

            % Prepare slopes: subtract average
%             mSx = mean(Sx(idx));
%             Sx = mask.*(Sx-mSx);
%             mSy = mean(Sy(idx));
%             Sy = mask.*(Sy-mSy);

            % Estimated slopes initially set to original slopes
            Sx_est = Sx;
            Sy_est = Sy;

            % Calculate filters
            RX = obj.RX;
            RY = obj.RY;
            GX = obj.GX;
            GY = obj.GY;
            absG2 = abs(GX).^2+abs(GY).^2;
            absG2(absG2<1e-6)=1e-6;
            
            for n=1:nIter

                % Add estimated slopes to measured ones
                % I.e. area within valid pupil is measured slopes, outside this
                % area we have the estimated slopes
                Sx_est = Sx.*mask + Sx_est.*imask;
                Sy_est = Sy.*mask + Sy_est.*imask;
                
                % FFT of new slopes
                fft_Sx = fftshift(fft2(Sx_est));
                fft_Sy = fftshift(fft2(Sy_est));
        
                % Estimate FFT slopes by filtering to produce FFT of phase (using
                % our reconstructors) and then use the reverse process (s = G*phi)
                % to get an estimate of the slopes (in Fourier domain)
%                 fft_Sx_est = GX.*(RX.*fft_Sx + RY.*fft_Sy);
%                 fft_Sy_est = GY.*(RX.*fft_Sx + RY.*fft_Sy);
                fft_Sx_est = GX.*(conj(GX).*fft_Sx + conj(GY).*fft_Sy)./absG2;
                fft_Sy_est = GY.*(conj(GX).*fft_Sx + conj(GY).*fft_Sy)./absG2;

                % Now get estimate of slopes using inverse FFT
                Sx_est = ifft2(ifftshift(fft_Sx_est));
                Sy_est = ifft2(ifftshift(fft_Sy_est));

            end

            % Final slopes are real part of estimated slopes
            Sx = real(Sx.*mask + Sx_est.*imask);
            Sy = real(Sy.*mask + Sy_est.*imask);

        end
        
        
        
        
                   
            
    end

    
    methods (Static)
        
        %% Noise PSD
        function [Wnoise] = createNoisePSD(k,noiseVar)
            
            % Creates noise PSD from noise variance, assuming white noise
            % (constant noise over all frequencies).
            
            N = length(k);
            Wnoise = noiseVar/N^2 * ones(N,N);
            
        end
        
        %% Phase PSD
        function [Wphi,WSx,WSy,WSxA,WSyA] = createPhasePSD(k,m,d,nPxL,r0,L0)
            
            % creates PSD of phase over given frequencies.  Option to
            % include aliasing
            % k:    frequency vector
            % m:    factor dictating maximum frequency to be aliased (i.e.
            %       how many times spectrum is folded).
            % d:    sampling period, max(k) = 1/(2*d).
            % nPxL: pixels per lenslet
            % r0,L0:atmospheric parameters
            %
            % N.B. no piston removal, can be added if needed.
            
            [kx,ky] = meshgrid(k,k);
            
            Wphi = 0;
            WSx = 0;
            WSy = 0;
            WSxA = 0;
            WSyA = 0;
            % Averaging d for sinc function
            d_av = d - d/nPxL;
            
            for mi = -floor(m/2):ceil(m/2)
                for ni = -floor(m/2):ceil(m/2)
                   
                    % Frequency range
                    km = [kx(:,k<0)+mi/d kx(:,k>=0)-mi/d];
                    kn = [ky(k<0,:)+ni/d; ky(k>=0,:)-ni/d];
                    
                    % Phase PSD 
                    % For current range of frequencies
                    Wmn = (abs(km.^2 + kn.^2) + (1/L0)^2).^(-11/6);
                    % Adding to total PSD
                    Wphi = Wphi + Wmn;
                    
                    % Slopes PSD
                    % Averaging term
                    Avx = abs(sinc(d_av*km).*sinc(d*kn).*exp(1i*pi*d_av*(km+kn))).^2;
                    Avy = abs(sinc(d*km).*sinc(d_av*kn).*exp(1i*pi*d_av*(km+kn))).^2;
                    % Total slopes PSD
                    WSx = WSx + Wmn.*abs(1i*2*pi*d*km).^2.*Avx;
                    WSy = WSy + Wmn.*abs(1i*2*pi*d*kn).^2.*Avy;
                    
                    % Aliased part
                    if ~(mi==0 && ni==0)
                        WSxA = WSxA + Wmn.*abs(1i*2*pi*d*km).^2.*Avx;
                        WSyA = WSyA + Wmn.*abs(1i*2*pi*d*kn).^2.*Avy;
                    end
                    
                end
            end
            
            % Normalisation factors (~0.0229*r0^(-5/3))
            Wphi = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*Wphi;
            WSx = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*WSx;
            WSy = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*WSy;
            WSxA = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*WSxA;
            WSyA = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*WSyA;
        end
        
        %% Windows
        function [w] = window(wname,N,vargin)
            % Windows for measuring PSD of variable defined over 2D grid 
            % (i.e. to smooth edges of the data)

            [X,Y] = meshgrid(1:N,1:N);

            n0 = (N+1)/2;

            R = sqrt((X-n0).^2 + (Y-n0).^2);

            w = zeros(N,N);
            [idx] = find(R<=N/2);

            if strcmp(wname,'barthann')
                w(idx) = 0.62 - 0.48*abs(R(idx)/N) + 0.38*cos(2*pi*R(idx)/N);
            elseif strcmp(wname,'bartlett')
                w(idx) = 1-2*R(idx)/N;
            elseif strcmp(wname,'bohman')
                rho = R/(N/2);
                w(idx) = (1-rho(idx)).*cos(pi*rho(idx)) + (1/pi)*sin(pi*rho(idx));
            elseif strcmp(wname,'flattop')
                a0 = 1-0.21557895;
                a1 = 0.41663158;
                a2 = 0.277263158;
                a3 = 0.083578947;
                a4 = 0.006947368;
                w(idx) = a0 + a1*cos(2*pi*R(idx)/N) - a2*cos(4*pi*R(idx)/N) ...
                    + a3*cos(6*pi*R(idx)/N) - a4*cos(8*pi*R(idx)/N);
            elseif strcmp(wname,'hamming')
                w(idx) = 0.46 + 0.46*cos(2*pi*R(idx)/N);
            elseif strcmp(wname,'hann')
                w(idx) = 0.5 + 0.5*cos(2*pi*R(idx)/N);
            elseif strcmp(wname,'nuttall')
                a0 = 1-0.3635819;
                a1 = 0.4891775;
                a2 = 0.1365995;
                a3 = 0.0106411;
                w(idx) = a0 + a1*cos(2*pi*R(idx)/N) - a2*cos(4*pi*R(idx)/N) ...
                    + a3*cos(6*pi*R(idx)/N);
            elseif strcmp(wname,'parzen')
                % ???? Is this the right equation?
                w(idx) = 1-(R(idx)/(N/2)).^3;
            elseif strcmp(wname,'tukey')
                r_tukey = vargin*N/2;
                w(idx) = 1/2 * (cos(2*pi/(2*(N/2-r_tukey))*(R(idx)-r_tukey))+1);
                i = find(R<=r_tukey);
                w(i) = 1;
            end
    
    
        end
        
        
        %% Remove modes
        function [Phi] = removeModes(Phi,N,pupil,modes,method)
            
            % methods = 'direct'  (in spatial domain)
            %           'fourier' (spatial frequency domain)
            % N.B. if the 'fourier' method is used then the phase (Phi)
            % must be passed to the function in fourier space (and will
            % also be returned in fourier space).  This phase should be
            % fftshifted to have the 0 frequency component in the centre of
            % the matrix/vector
            %
            % N: number of measurements
            % pupil: logical matrix describing valid elements (i.e. inside
            % the pupil)
            %
            % N.B: only compatible for waffle when Phi has an even number
            % of points.  I.e. we remove the waffle frequency at (1,1),
            % which must be an eigenmode of the measurement space. a vector
            % of 1,-1 across a grid is only an eigenmode over an even
            % number of points
            
            if strcmp(method,'direct')
                % Global waffle
                if modes(1)==1
                    % Waffle mode
                    waffle = ones(N,N);
                    for i=1:N
                        if mod(i,2) == 0
                            waffle(1:2:N,i) = -1;
                        else
                            waffle(2:2:N,i) = -1;
                        end
                    end
                    
                    % Waffle coefficient in Phi (calculate using inner
                    % product)
                    cw = sum(Phi(pupil).*waffle(pupil)) / sum(waffle(pupil).*waffle(pupil));
                    Phi = Phi - cw*waffle;
                end
                                
                % Local waffle removal
                if modes(2)==1
                    Phi = conv2(Phi,[-0.25 0.25; 0.25 0.75],'same');
                    %Phi = Phi + (1/4)*conv2(Phi,[-1 1; 1 -1],'same');
                end
                                
                % Piston removal
                if modes(3)==1
                    Phi(pupil) = Phi(pupil) - mean(Phi(pupil(:)));
                end
                
            elseif strcmp(method,'fourier')
                % Both the global and local waffle modes can be filtered in
                % the Fourier space (no need to go to direct space); 
                % N.B.: As has been outlined in Pyoneer02 (OSA paper) the 
                % removal in Fourier space is not as much efective as in 
                % direct space due to artifacts originated by the extension 
                % method.
                
                % Global waffle
                if modes(1)==1
                    
                    % Waffle frequency (assuming an even number of points
                    % in the grid)
                    Phi(1,1) = 0;
%                     Phi(1,N/2+1) = 0;
%                     Phi(N/2+1,1) = 0;
                end
                
                % Local waffle
                if modes(2)==1
                    % the Fourier transform of the direct-space filter is used (See Poyneer
                    % 2002 (SPIE conf. paper) for the original filter)
                    % Local waffle mode in direct space
                    Wmode = zeros(N,N);
                    Wmode(floor(N/2+1)-1:floor(N/2+1),floor(N/2+1)-1:floor(N/2+1)) = [-0.25 0.25; 0.25 0.75];
                    Phi = Phi.*ifftshift(fft2(fftshift(Wmode)));
%                     Wmode(floor(N/2+1)-1:floor(N/2+1),floor(N/2+1)-1:floor(N/2+1)) = [-1 1; 1 -1];
%                     Phi = Phi./ifftshift((fft2(fftshift(Wmode))));
                    

                end
                
                % Global piston
                if modes(3)==1
                    Phi(floor(N/2+1),floor(N/2+1)) = 0;
                end
                
            end
        end
        
        %% Pad data
        function mat = padarray(mat,vv)
            % Function to pad array
            if length(vv) == 1
                x = vv;
                y = vv;
            elseif length(vv) == 2
                x = vv(1);
                y = vv(2);
            else
                fprintf('Bad input')
            end

            % Padding with zeros
            mat = [zeros(size(mat,1),x) mat zeros(size(mat,1),x)];
            mat = [zeros(y,size(mat,2)); mat; zeros(y,size(mat,2))];
        end

        %% DM filters in fourier domain
        function out = dmFilter(dmFunction,ext,nLenslet,nPx)

            % High resolution filter: then fold in 
            nLPx = nPx/nLenslet;
            x = 0:1/nLPx:nLenslet+1-1/nLPx;
            % Centre at correct index
            x = x-x(ceil((nPx+nLPx+1)/2));
            % Indices for non-zero values
            if isa(dmFunction,'influenceFunction')
                idx = x>= -dmFunction.bezier(end,1) & x<= dmFunction.bezier(end,1);
            else
                idx = x>= -dmFunction.gaussian(end,1) & x<= dmFunction.gaussian(end,1);
            end
            
            % 1D individual influence function (in direct space)
            wv = zeros(nPx+nLPx,1);
            wv(idx) = ppval(dmFunction.splineP,x(idx));
            
            % 2D individual IF
            h = wv*wv';
            
                        
            % DM filter ---------------------------------------------------
            % Pad with zeros
            if nLenslet/2~=round(nLenslet/2)
                if nLPx/2==round(nLPx/2)
                    h = fourierReconstructor.padarray(h,nLPx/2);
                else
                    h = [h zeros(nPx+nLPx,1); zeros(1,nPx+nLPx)];
                    h = fourierReconstructor.padarray(h,(nLPx-1)/2);
                end
            else
                h = fourierReconstructor.padarray(h,nLPx);
            end
            hh = fourierReconstructor.padarray(h,ext*nLPx);
            
%             size(hh)
%             
%             figure
%                 imagesc(hh)
%                 drawnow;
            
            hhFFT = fftshift(fft2(ifftshift(hh)));           
            
            N = 2*ceil((nLenslet+1+2*ext)/2); 
            n0 = ceil((N*nLPx+1)/2);
            idx = (n0-ceil((N-1)/2)):(n0+floor((N-1)/2));
            
            % Create folded DM filter
            hhFFT_fold = zeros(N,N);
            % For no fold mm=0;
            mm = nLPx-1;
            
            for n = -floor(mm/2):ceil(mm/2)
                for m = - floor(mm/2):ceil(mm/2)
                    
                    % Fold in
                    in = [idx(idx<n0)+n*N idx(idx>=n0)-n*N];
                    im = [idx(idx<n0)+m*N idx(idx>=n0)-m*N];
                    hhFFT_fold = hhFFT_fold + hhFFT(im,in);
                    
                end
            end
            
            hhFFT_fold = N^2/(N*nLPx)^2*ifftshift(hhFFT_fold);
            
            % DM filter
            out = conj(hhFFT_fold)./abs(hhFFT_fold).^2;


        end

    end
    
end
        
function p = pattern(length,arg)
    p = ones(1,length);
    for pi = 1:2:length
        p(pi) = -1;
    end
    if nargin > 1
        p = p(end:-1:1);
    end
end

function [xslopes,yslopes,nn] = solve_conflicts(i,point,N,mask,xslopes,yslopes)
    localmask = mask(max(i-1,1):min(i+1,N),max(point-1,1):min(point+1,N));
    nneighbours = nnz(localmask);
    if nneighbours > 2 & nnz(sum(localmask,2) > 0) > 1 & nnz(sum(localmask,1) > 0) > 1% the conditions stands for: more than two neighbours,
        %with at least one at the left/right and one above/below
        nn = 1;
        if i < N/2+1 && point < N/2+1 %1st quadrant
            Q = 1;
        elseif i < N/2+1 && point > N/2+1 %2nd quadrant
            Q = 2;
        elseif i > N/2+1 && point < N/2+1%3rd quadrant
            Q = 3;
        elseif i > N/2+1 && point > N/2+1%4th quadrant
            Q = 4;
        end
        if Q == 1
            xslopes(i,point) = -xslopes(i,point+1) + xslopes(i+1,point) + xslopes(i+1,point+1);
            yslopes(i,point) = -yslopes(i+1,point) + yslopes(i,point+1) + yslopes(i+1,point+1);
        elseif Q == 2
            xslopes(i,point) = -xslopes(i,point-1) + xslopes(i+1,point) + xslopes(i+1,point-1);
            yslopes(i,point) = -yslopes(i+1,point) + yslopes(i+1,point-1) + yslopes(i,point-1);
        elseif Q == 3
            xslopes(i,point) = xslopes(i-1,point) + xslopes(i-1,point+1) - xslopes(i,point+1);
            yslopes(i,point) = -yslopes(i-1,point) + yslopes(i-1,point+1) + yslopes(i,point+1);
        elseif Q == 4
            xslopes(i,point) = -xslopes(i,point-1) + xslopes(i-1,point-1) + xslopes(i-1,point);
            yslopes(i,point) = -yslopes(i-1,point) + yslopes(i-1,point-1) + yslopes(i,point-1);
        end
    else
        nn = 0;
    end
end
