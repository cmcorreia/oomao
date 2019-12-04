% Modification Log
%
% 12/12/2016 Yoshito Ono
% Implement functions to compute a map directly.
%   getCovarianceMapSingleLayer
%   getCovarianceMapSingleLayerFFT
%
% 12/12/2016 Yoshito Ono
% I found and fixed a bug in function of "dphi".
%
% 17/11/2016 Yoshito Ono
% Fixed bug : getSubapertureCoordinates
%   Before : x     = linspace(-D/2.,D/2., nL);
%   After  : d     = D/nL;
%            x     = linspace(-D/2.+d/2.,D/2.-d/2., nL);
% The slope is defined at the center of sub-aperture.
%
% 17/11/2016 Yoshito Ono
% I added FFT slope covariance Model.
%    "getCovarianceMatrixSingleLayer" for full cavriance matrix
%    "getCovarianceMapSingleLayer" for covariance map
% You can swith models by changing a flag of "isFFT".
% By default, "isFFT" is set to false and the Hudin model is used.
%
% 17/11/2016 Yoshito Ono
% I removed some inputs from several functions in this class like
%   Before : function covMat = getCovarianceMatrix(obj,fr0,h,L0,xSrc,ySrc,G,th,dX,dY,s,track,flagTT)
%   After  : function covMat = getCovarianceMatrix(obj)
% because all inputs are already saved in the object (obj.h, obj,fr0, etc.)
%
% 27/10/2016 Olivier Martin
% I added some functions to compute the covariance matrix of the residual error
% a MMSE reconstructor (slopes to slopes) and the wavefront tomographic
% error

% 26/10/2016  Yoshito Ono
% I add "inputPaser" function in the constructor to
% check inputs and to set default values for additional
% inputs. Now, you can call this class as
% >> model = slopesCovarianceModel(tel,atm,wfs,gs);
% In case you need additional inputs (e.g. mis-registration)
% you should call this funtion like
% >> model = slopesCovarianceModel(tel,atm,wfs,gs,'misMat',misMat);
%
% I also fix "Mis-registration" and "Additonal tip/tilt" in the
% constructor to set default values if there is no inputs for them.
%
%
classdef slopesCovarianceModel < handle
    
    properties (SetObservable = true)
        %Sub-classes
        tel;
        atm;
        wfs;
        gs;
        %Parameters
        ns                  %# slopes
        nwfs;               %#WFS
        nmeas;              %#Measurement = ns*nwfs
        pitch;
        parameters;
        %Atmosphere
        fr0;
        h;
        L0;
        %Sources
        xSrc;
        ySrc;
        %Mis-registration
        G;
        th;
        dX;
        dY;
        s;
        track;
        %Outputs
        covMat;
        %Flags
        flagTT;
        flagFocus;
        rmDiagonal;
        twoSteps;
        isFFT = false;
        isXY = false;
        isAuto = false;
        isMap = false;
        %Fitting procedure
        initAtm;
        initSys;
        idxFitAtm;
        idxFitSys;
        %Tomographic error
        R;
        Css;
        Csg;
        Cgs;
        Cgg;
        Cee;
        Cnn;
        Projector;
        ProjectorUnits;
        wfeTomo;
        
        NF = 512;
        sf = 4;
    end
    
    properties (Dependent)
        noiseVar;
    end
    
    properties (Access=private)
        p_noiseVar;
    end
    
    
    methods
        
        function obj = slopesCovarianceModel(tel,atm,wfs,gs,varargin)
            
            % Model description
            % covMat = f(fr0,h,L0,xSrc,ySrc,G,th,dX,dY,tracking,s,flagTT), with:
            % fr0 = [r0_k^(-5/3.)] : fractionnal r0^(-5/3.)         [m^(-5/3)]
            % h   = [h_k]          : altitude                       [m]
            % L0  = [L0_k]         : outer scale profile            [m]
            % xSrc                 : x locations of sources         [rad]
            % ySrc                 : y location of sources          [rad]
            % G                    : relative WFS magnification     [unitless]
            % th                   : relative WFS rotation          [deg]
            % dX                   : relative pupil x-shift         [m]
            % dY                   : relative pupil y-shift         [m]
            % tracking             : xx,yy,xy covariance error      [arcsec^2]
            %                        due to tel. tracking error
            %s                     : relative centroid gain         [untiless]
            %flagTT                : tip-tilt removed if flagTT == 1
            
            
            % Check Inputs
            inputs = inputParser;
            inputs.addRequired('tel',@(x) isa(x,'telescope') );
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann'));
            inputs.addRequired('gs',@(x) isa(x,'source'));
            inputs.addParameter('noiseVar',0, @isnumeric );
            inputs.addParameter('misMat', [], @isnumeric );
            inputs.addParameter('track', [], @isnumeric );
            inputs.addParameter('Projector', [], @isnumeric );
            inputs.addParameter('ProjectorUnits', 1, @isnumeric );
            inputs.addParameter('flagTT', false, @islogical );
            inputs.addParameter('flagFocus', false, @islogical );
            inputs.addParameter('isXY', false, @islogical );
            inputs.addParameter('isAuto', false, @islogical );
            inputs.addParameter('isMap', false, @islogical );
            inputs.addParameter('isFFT', false, @islogical );
            inputs.addParameter('rmDiagonal', false, @islogical );
            inputs.addParameter('twoSteps', false, @islogical );
            inputs.addParameter('idxFitAtm', zeros(1,3*atm.nLayer), @isnumeric );
            inputs.addParameter('idxFitSys', zeros(1,7*length(gs)+3), @isnumeric );
            inputs.parse(tel,atm,wfs,gs,varargin{:});
            %Instantiation
            obj.tel      = tel;
            obj.atm      = atm;
            obj.wfs      = wfs;
            obj.gs       = gs;
            %Parameters
            obj.ns        = floor(2*sum(wfs.validLenslet(:)));
            obj.nwfs      = floor(length(gs));
            obj.nmeas     = floor(obj.ns*obj.nwfs);
            obj.pitch     = obj.tel.D/obj.wfs.lenslets.nLenslet;
            obj.flagTT    = inputs.Results.flagTT;
            obj.flagFocus = inputs.Results.flagFocus;
            obj.isXY      = inputs.Results.isXY;
            obj.isAuto    = inputs.Results.isAuto;
            obj.isMap     = inputs.Results.isMap;
            obj.isFFT     = inputs.Results.isFFT;
            obj.rmDiagonal= inputs.Results.rmDiagonal;
            obj.twoSteps  = inputs.Results.twoSteps;
            obj.idxFitAtm = inputs.Results.idxFitAtm;
            obj.idxFitSys = inputs.Results.idxFitSys;
            obj.noiseVar  = inputs.Results.noiseVar;
            %Mis-registration
            misMat = inputs.Results.misMat;
            if isempty(misMat) % set default values
                misMat            = zeros(5,obj.nwfs);
                misMat(1,:)       = 1;
                misMat(5,:)       = 1;
            end
            obj.G  = misMat(1,:);
            obj.th = misMat(2,:).*pi/180;
            obj.dX = misMat(3,:);
            obj.dY = misMat(4,:);
            obj.s  = misMat(5,:);
            %Additonal isoplanatic tip/tilt from telescope tracking error
            obj.track = inputs.Results.track;
            if isempty(obj.track) % set default values
                obj.track      = [0,0,0]; %xx - yy - xy common tip-tilt
            end
            %Sources location
            locSrc     = [obj.gs.directionVector];
            obj.xSrc   = locSrc(1,:);
            obj.ySrc   = locSrc(2,:);
            %Atmosphere
            obj.fr0    = [obj.atm.layer.fractionnalR0].*(atm.r0)^(-5/3.);
            obj.h      = [obj.atm.layer.altitude];
            obj.L0     = [obj.atm.layer.layeredL0];
            obj        = obj.atmToInit();
            %Projector
            P          = inputs.Results.Projector;
            if ~isempty(P)
                obj.Projector      = P;
                obj.ProjectorUnits = inputs.Results.ProjectorUnits;
            end
            
        end
        
        %% Get/Set noiseVar
        function out = get.noiseVar(obj)
            out = obj.p_noiseVar;
        end
        function set.noiseVar(obj,val)
            if isempty(obj.p_noiseVar)
                obj.p_noiseVar = zeros(2*obj.nwfs,1);
            end
            obj.p_noiseVar(:) = val;
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FEATURES FOR THE FITTING PROCEDURE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function covMat = covarianceModel(obj,x)
            if any(x)
                %Grabbing indexFit
                obj      = obj.indexFitToObj(x);
            end
            covMat = obj.getCovarianceMatrix();
        end
        function obj = atmToInit(obj)
            obj.initAtm = [obj.atm.r0^(-5/3.).*[obj.atm.layer.fractionnalR0], ...
                [obj.atm.layer.altitude],[obj.atm.layer.layeredL0]];
            obj.initSys = [obj.xSrc,obj.ySrc,obj.G,obj.th,obj.dX,obj.dY,obj.s,obj.track];
            obj.parameters = [obj.initAtm,obj.initSys];
        end
        function xdata = objToData(obj)
            xdata(1).value  = obj.fr0;
            xdata(2).value  = obj.h;
            xdata(3).value  = obj.L0;
            xdata(4).value  = obj.xSrc;
            xdata(5).value  = obj.ySrc;
            xdata(6).value  = obj.G;
            xdata(7).value  = obj.th;
            xdata(8).value  = obj.dX;
            xdata(9).value  = obj.dY;
            xdata(10).value = obj.s;
            xdata(11).value = obj.track;
            xdata(12).value = obj.flagTT;
        end
        function obj = dataToObj(obj,x)
            nL        = obj.atm.nLayer;
            nw        = obj.nwfs;
            obj.fr0   = x(1:nL);
            obj.h     = x(1+nL:2*nL);
            obj.L0    = x(1+2*nL:3*nL);
            obj.xSrc  = x(1+3*nL:3*nL+nw);
            obj.ySrc  = x(1+3*nL+nw:3*nL+2*nw);
            obj.G     = x(1+3*nL+2*nw:3*nL+3*nw);
            obj.th    = x(1+3*nL+3*nw:3*nL+4*nw);
            obj.dX    = x(1+3*nL+4*nw:3*nL+5*nw);
            obj.dY    = x(1+3*nL+5*nw:3*nL+6*nw);
            obj.s     = x(1+3*nL+6*nw:3*nL+7*nw);
            obj.track = x(1+3*nL+7*nw:3*nL+7*nw+3);
        end
        function obj = indexFitToObj(obj,x)
            %Common parameters
            nL       = obj.atm.nLayer;
            nw       = obj.nwfs;
            indexFit = [obj.idxFitAtm,obj.idxFitSys];
            %Split the indexFit array into sub-index vectors
            idxr0    = indexFit(1:nL);
            idxh     = indexFit(1+nL:2*nL);
            idxL0    = indexFit(1+2*nL:3*nL);
            idxxSrc  = indexFit(1+3*nL:3*nL+nw);
            idxySrc  = indexFit(1+3*nL+nw:3*nL+2*nw);
            idxG     = indexFit(1+3*nL+2*nw:3*nL+3*nw);
            idxth    = indexFit(1+3*nL+3*nw:3*nL+4*nw);
            idxdX    = indexFit(1+3*nL+4*nw:3*nL+5*nw);
            idxdY    = indexFit(1+3*nL+5*nw:3*nL+6*nw);
            idxs     = indexFit(1+3*nL+6*nw:3*nL+7*nw);
            idxtrack = indexFit(1+3*nL+7*nw:3*nL+7*nw+3);
            %Updating the obj class
            nr0=0;nh=0;nL0=0;nx=0;ny=0;nG=0;nth=0;ndX=0;ndY=0;ncg=0;nt=0;
            
            if any(idxr0)
                nr0 = sum(idxr0);
                tmp = x(1:nr0);
                obj.fr0(idxr0==1) = tmp;
            end
            if any(idxh)
                nh  = sum(idxh==1)+nr0;
                tmp = x(1+nr0:nh);
                obj.h(idxh==1) = tmp;
            end
            if any(idxL0)
                nL0 = sum(idxL0==1)+nh;
                tmp = x(1+nh:nL0);
                obj.L0(idxL0==1) = tmp;
            end
            if any(idxxSrc)
                nx  = sum(idxxSrc==1)+nL0;
                tmp = x(1+nL0:nx);
                obj.xSrc(idxxSrc==1) = tmp;
            end
            if any(idxySrc)
                ny  = sum(idxySrc==1)+nx;
                tmp = x(1+nx:ny);
                obj.ySrc(idxySrc==1) = tmp;
            end
            if any(idxG)
                nG  = sum(idxG==1)+ny;
                tmp = x(1+ny:nG);
                obj.G(idxG==1) = tmp;
            end
            if any(idxth)
                nth = sum(idxth==1)+nG;
                tmp = x(1+nG:nth);
                obj.th(idxth==1) = tmp;
            end
            if any(idxdX)
                ndX = sum(idxdX==1)+nth;
                tmp = x(1+nG:ndX);
                obj.dX(idxdX==1) = tmp;
            end
            if any(idxdY)
                ndY = sum(idxdY==1)+ndX;
                tmp = x(1+ndX:ndY);
                obj.dY(idxdY==1) = tmp;
            end
            if any(idxs)
                ncg  = sum(idxs==1) + ndY;
                tmp = x(1+ndY:ncg);
                obj.s(idxs==1) = tmp;
            end
            if any(idxtrack)
                nt  = sum(idxtrack==1)+ncg;
                tmp = x(1+ncg:nt);
                obj.track(idxtrack==1) = tmp;
            end
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MODEL DESCRIPTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function covMat = getCovarianceMatrix(obj)
            
            %Instantation
            covMat = zeros(obj.nmeas);
            nLayer = length(obj.h);
            
            % ---------------- Loop on atmospheric layer -----------------%
            for kLayer = 1:nLayer
                %Grabs a single altitude layer
                if ~obj.isFFT  % Hudgin Model
                    covMat   = covMat + obj.getCovarianceMatrixSingleLayer(kLayer);
                else % FFT Model
                    covMat   = covMat + obj.getCovarianceMatrixSingleLayerFFT(kLayer);
                end
            end
            
            % --------------------- Noise options ------------------------%
            if ~obj.rmDiagonal
                if any(obj.noiseVar)
                    covNoise = cell(2*obj.nwfs);
                    for k = 1:2*obj.nwfs
                        covNoise{k} = obj.noiseVar(k).*eye(obj.ns/2);
                    end
                    obj.Cnn  = blkdiag(covNoise{:});
                    covMat   = covMat + obj.Cnn;
                end
            end
            
            % ---------------- Tracking telescope error ------------------%
            if ~obj.flagTT
                covMat = addTrackingError(obj,covMat);
            end
            
            % --------------- Relative WFS sensitivities -----------------%
            covMat = addCentroidGainError(obj,covMat);
            
            % ------------------ Tip-tilt filtering ----------------------%
            if obj.flagTT
                covMat = obj.manageTipTilt(covMat);
            end
            
            % ------------------- 2-steps management ---------------------%
            if obj.twoSteps
                P      = obj.groundLayerFilter;
                covMat = P*covMat*P';
            end
            
            
            if obj.rmDiagonal
                covMat(1:obj.nmeas+1:end) = 0;
            end
            obj.covMat = covMat;
        end
        function covMat = addTrackingError(obj,covMat)
            nx    = obj.ns/2;
            
            for iGs = 1:obj.nwfs
                % Defining indexes
                xp = 1 + (iGs-1)*obj.ns :(iGs-1)*obj.ns + nx;
                yp = 1 + (iGs-1)*obj.ns + nx : iGs*obj.ns;
                for jGs = 1:iGs
                    xq = 1 + (jGs-1)*obj.ns :(jGs-1)*obj.ns + nx;
                    yq = 1 + (jGs-1)*obj.ns + nx : jGs*obj.ns;
                    %Adding tracking error
                    covMat(xp,xq) = covMat(xp,xq) + obj.track(1);
                    covMat(yp,yq) = covMat(yp,yq) + obj.track(2);
                    covMat(xp,yq) = covMat(xp,yq) + obj.track(3);
                    covMat(yp,xq) = covMat(yp,xq) + obj.track(3);
                    
                    if(jGs~=iGs)
                        covMat(xq,xp) = covMat(xq,xp) + obj.track(1);
                        covMat(yq,yp) = covMat(yq,yp) + obj.track(2);
                        covMat(xq,yp) = covMat(xq,yp) + obj.track(3);
                        covMat(yq,xp) = covMat(yq,xp) + obj.track(3);
                    end
                end
            end
        end
        function covMat = addCentroidGainError(obj,covMat)
            for iGs = 1:obj.nwfs
                % Defining indexes
                rp = 1 + (iGs-1)*obj.ns :iGs*obj.ns;
                for jGs = 1:iGs
                    rq = 1 + (jGs-1)*obj.ns :jGs*obj.ns;
                    %Multiplying by relative sensitivity
                    covMat(rp,rq) = covMat(rp,rq)*obj.s(iGs)*obj.s(jGs);
                    if jGs~=iGs
                        covMat(rq,rp) = covMat(rq,rp).*obj.s(iGs)^2;
                    end
                end
            end
        end
        
        function covMap = getCovarianceMapSinglePair(obj,iGs,jGs)
            %% getCovarianceMapSinglePair
            %
            % updated by Y.O. 15/1/2017
            nSub = obj.wfs(1).lenslets.nLenslet;
            n=nSub*2-1;
            z = meshgrid((1-nSub:nSub-1));
            z = z+z'*1i;
            p = ifft2(fft2(obj.wfs(1).validLenslet,n,n).*fft2(obj.wfs(1).validLenslet,n,n));
            p = p>.5;
            covxx = zeros(n);
            covyy = zeros(n);
            if obj.isXY
                covxy = zeros(n);
            end
            for kLayer = 1:length(obj.h)
                if obj.fr0(kLayer)<0.05
                    continue;
                end
                g1 = 1 - obj.h(kLayer)/obj.gs(iGs).height;
                %g2 = 1 - obj.h(kLayer)/obj.gs(jGs).height;
                d1 = obj.pitch*g1;
                dd1  = obj.h(kLayer)*tan([obj.xSrc(jGs) obj.ySrc(jGs)]);
                dd2  = obj.h(kLayer)*tan([obj.xSrc(iGs) obj.ySrc(iGs)]);
                dd = dd2 - dd1;
                zh = z*d1 - dd(1) - dd(2)*1i;
                Dphi = @(z_) obj.dphi(abs(z_),obj.L0(kLayer));
                A = Dphi(zh);
                covxx = covxx + obj.fr0(kLayer)*(-A+Dphi(zh+d1)+Dphi(zh-d1)-A);
                covyy = covyy + obj.fr0(kLayer)*(-A+Dphi(zh+d1*1i)+Dphi(zh-d1*1i)-A);
                if obj.isXY
                    covxy = covxy + obj.fr0(kLayer)*(Dphi(zh+d1/2*1i+d1/2)-Dphi(zh-d1/2*1i+d1/2)-...
                        Dphi(zh+d1/2*1i-d1/2)+Dphi(zh-d1/2*1i-d1/2));
                end
            end
            covxx = covxx.*p;
            covyy = covyy.*p;
            if obj.isXY
                covxy = covxy.*p;
            end
            if obj.flagTT
                mode = [2 3];
                if obj.flagFocus
                    mode = [2 3 4];
                end
                mask = obj.wfs.validLenslet;
                %covxx = tools.tipTiltFilter4covMap(covxx,nSub,mask);
                covxx = tools.slopeModeFilter4covMap(covxx,nSub,mask,mask,'x','x',mode);
                %covyy = tools.tipTiltFilter4covMap(covyy,nSub,mask);
                covyy = tools.slopeModeFilter4covMap(covyy,nSub,mask,mask,'y','y',mode);
                if obj.isXY
                    %covxy = tools.tipTiltFilter4covMap(covxy,nSub,mask);
                    covxy = tools.slopeModeFilter4covMap(covxy,nSub,mask,mask,'x','y',mode);
                end
            end
            covMap = [covxx; covyy];
            if obj.isXY
                covMap = [covMap; covxy];
            end
            lambda = obj.atm.wavelength;
            covMap = 0.5*covMap.*(constants.radian2arcsec*lambda/2/pi/obj.pitch)^2;
        end
        
        function covMat = getCovarianceMapSingleLayer(obj,kLayer)
            %% getCovarianceMapSingleLayer
            %
            % updated by Y.O. 15/1/2017
            nSub = obj.wfs(1).lenslets.nLenslet;
            n=nSub*2-1;
            z = meshgrid((1-nSub:nSub-1));
            z = z+z'*1i;
            covMat = zeros(2*n*obj.nwfs);
            p = ifft2(fft2(obj.wfs(1).validLenslet,n,n).*fft2(obj.wfs(1).validLenslet,n,n));
            p = p>.5;
            for iGs = 1:obj.nwfs
                g1 = 1 - obj.h(kLayer)/obj.gs(iGs).height;
                for jGs = 1:iGs
                    %g2 = 1 - obj.h(kLayer)/obj.gs(jGs).height;
                    d1 = obj.pitch*g1;
                    dd1  = obj.h(kLayer)*tan([obj.xSrc(jGs) obj.ySrc(jGs)]);
                    dd2  = obj.h(kLayer)*tan([obj.xSrc(iGs) obj.ySrc(iGs)]);
                    dd = dd2 - dd1;
                    zh = z*d1 - dd(1) - dd(2)*1i;
                    Dphi = @(z_) obj.dphi(abs(z_),obj.L0(kLayer));
                    A = Dphi(zh);
                    covxx = -A+Dphi(zh+d1)+Dphi(zh-d1)-A;
                    covyy = -A+Dphi(zh+d1*1i)+Dphi(zh-d1*1i)-A;
                    covxy = Dphi(zh+d1/2*1i+d1/2)-Dphi(zh-d1/2*1i+d1/2)-...
                        Dphi(zh+d1/2*1i-d1/2)+Dphi(zh-d1/2*1i-d1/2);
                    covxx = covxx.*p;
                    covyy = covyy.*p;
                    covxy = covxy.*p;
                    ki = (1:2*n)+(2*n)*(iGs-1);
                    kj = (1:2*n)+(2*n)*(jGs-1);
                    covMat(ki,kj) = [covxx covxy; covxy covyy];
                    if iGs~=jGs
                        covMat(kj,ki) = [rot90(covxx,2) rot90(covxy,2); rot90(covxy,2) rot90(covyy,2)];
                    end
                end
            end
            lambda = obj.atm.wavelength;
            ps     = (constants.radian2arcsec*lambda/2/pi/obj.pitch)^2;
            covMat = obj.fr0(kLayer)*0.5*covMat.*ps; %Covariance in arcsec^2
        end
        
        function covMat = getCovarianceMatrixSingleLayer(obj,kLayer)
            %Instantation
            covMat = zeros(obj.nmeas);
            
            %Loop on WFS
            for iGs = 1:obj.nwfs
                [xi,yi] = getSubapertureCoordinates(obj,kLayer,iGs);
                [Xi,Yi] = meshgrid(xi,yi);
                for jGs = 1:iGs
                    [xj,yj] = getSubapertureCoordinates(obj,kLayer,jGs);
                    [Xj,Yj] = meshgrid(xj,yj);
                    u       = Xj - Xi';
                    v       = Yj' - Yi;
                    %Covariance matrix computation
                    tmp = separationToCovarianceMatrix(obj,u,v,obj.fr0(kLayer),obj.h(kLayer),obj.L0(kLayer),iGs,jGs);
                    %Concatenating values
                    covMat(obj.slrange(iGs),obj.slrange(jGs)) = tmp;
                    if jGs ~= iGs
                        covMat(obj.slrange(jGs),obj.slrange(iGs)) = tmp';
                    end
                end
            end
            %Covariance in arcsec^2
            phaseToAngle = obj.atm.wavelength/2/pi/obj.pitch;
            covMat       = 0.5*covMat.*(constants.radian2arcsec*phaseToAngle)^2;
            
        end
        
        function covMap = getCovarianceMapSinglePairFFT(obj,iGs,jGs)
            %% getCovarianceMapSinglePairFFT
            %
            % updated by Y.O. 15/1/2017
            nSub = obj.wfs(1).lenslets.nLenslet;
            NF = obj.NF; % Fourier spectrum resolution
            sf = obj.sf;%obj.sf; % Fourier spectrum oversampling factor
            lambda = obj.atm.wavelength;
            [fx0,fy0] = freqspace(NF,'meshgrid'); % 0:2
            n=nSub*2-1;
            b0 = NF/2+1;
            b  = (1-nSub:nSub-1)*sf + b0;
            p = ifft2(fft2(obj.wfs(1).validLenslet,n,n).*fft2(obj.wfs(1).validLenslet,n,n));
            p = p>.5;
            covxx = zeros(NF);
            covyy = zeros(NF);
            if obj.isXY
                covxy = zeros(NF);
            end
            
            for kLayer = 1:length(obj.h)
                if obj.fr0(kLayer)<0.05
                    continue;
                end
                g1 = 1 - obj.h(kLayer)/obj.gs(iGs).height;
                %g2 = 1 - obj.h(kLayer)/obj.gs(jGs).height;
                dd1  = obj.h(kLayer)*tan([obj.xSrc(jGs) obj.ySrc(jGs)]);
                dd2  = obj.h(kLayer)*tan([obj.xSrc(iGs) obj.ySrc(iGs)]);
                fx = fx0*sf/2/(obj.pitch*g1);
                fy = fy0*sf/2/(obj.pitch*g1);
                f = hypot(fx,fy);
                dd = dd2 - dd1;
                delta = lambda*lambda/NF/NF*sf*sf/(obj.pitch*g1)/(obj.pitch*g1);
                phasor = exp( 2*1i*pi.*( dd(1)*fx + dd(2)*fy ) ); % shift in the Fourier space
                spectrum = @(u,v) delta*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
                    obj.spectrum(f,obj.L0(kLayer)).*...
                    (tools.sinc(obj.pitch*g1*fx).*tools.sinc(obj.pitch*g1*fy)).^2.*phasor;
                specxx = spectrum([1,0],[1,0]);
                specyy = spectrum([0,1],[0,1]);
                covxx = covxx + obj.fr0(kLayer)*real( fftshift( fft2( fftshift( specxx ) ) ) );
                covyy = covyy + obj.fr0(kLayer)*real( fftshift( fft2( fftshift( specyy ) ) ) );
                if obj.isXY
                    specxy = spectrum([1,0],[0,1]);
                    covxy = covxy + obj.fr0(kLayer)*real( fftshift( fft2( fftshift( specxy ) ) ) );
                end
            end
            covxx = covxx(b,b).*p;
            covyy = covyy(b,b).*p;
            if obj.isXY
                covxy = covxy(b,b).*p;
            end
            if obj.flagTT
                mode = [2 3];
                if obj.flagFocus
                    mode = [2 3 4];
                end
                mask = obj.wfs.validLenslet;
                %covxx = tools.tipTiltFilter4covMap(covxx,nSub,mask);
                covxx = tools.slopeModeFilter4covMap(covxx,nSub,mask,mask,'x','x',mode);
                %covyy = tools.tipTiltFilter4covMap(covyy,nSub,mask);
                covyy = tools.slopeModeFilter4covMap(covyy,nSub,mask,mask,'y','y',mode);
                if obj.isXY
                    %covxy = tools.tipTiltFilter4covMap(covxy,nSub,mask);
                    covxy = tools.slopeModeFilter4covMap(covxy,nSub,mask,mask,'x','y',mode);
                end
            end
            covMap = [covxx; covyy];
            if obj.isXY
                covMap = [covMap; covxy];
            end
            covMap = (constants.radian2arcsec)^2*covMap;
        end
        
        function covMat = getCovarianceMapSingleLayerFFT(obj,kLayer)
            %% getCovarianceMapSingleLayerFFT
            %
            % updated by Y.O. 15/1/2017
            nSub = obj.wfs(1).lenslets.nLenslet;
            NF = obj.NF; % Fourier spectrum resolution
            sf = obj.sf; % Fourier spectrum oversampling factor
            lambda = obj.atm.wavelength;
            [fx0,fy0] = freqspace(NF,'meshgrid'); % 0:2
            n=nSub*2-1;
            b0 = NF/2+1;
            b  = (1-nSub:nSub-1)*sf + b0;
            p = ifft2(fft2(obj.wfs(1).validLenslet,n,n).*fft2(obj.wfs(1).validLenslet,n,n));
            p = p>.5;
            covMat = zeros(2*n*obj.nwfs);
            for iGs = 1:obj.nwfs
                g1 = 1 - obj.h(kLayer)/obj.gs(iGs).height;
                for jGs = 1:iGs
                    %g2 = 1 - obj.h(kLayer)/obj.gs(jGs).height;
                    dd1  = obj.h(kLayer)*tan([obj.xSrc(jGs) obj.ySrc(jGs)]);
                    dd2  = obj.h(kLayer)*tan([obj.xSrc(iGs) obj.ySrc(iGs)]);
                    fx = fx0*sf/2/(obj.pitch*g1);
                    fy = fy0*sf/2/(obj.pitch*g1);
                    f = hypot(fx,fy);
                    dd = dd2 - dd1;
                    delta = lambda*lambda/NF/NF*sf*sf/(obj.pitch*g1)/(obj.pitch*g1);
                    phasor = exp( 2*1i*pi.*( dd(1)*fx + dd(2)*fy ) ); % shift in the Fourier space
                    spectrum = @(u,v) delta*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
                        obj.spectrum(f,obj.L0(kLayer)).*...
                        (tools.sinc(obj.pitch*g1*fx).*tools.sinc(obj.pitch*g1*fy)).^2.*phasor;
                    specxx = spectrum([1,0],[1,0]);
                    specyy = spectrum([0,1],[0,1]);
                    specxy = spectrum([1,0],[0,1]);
                    covxx = real( fftshift( fft2( fftshift( specxx ) ) ) );
                    covyy = real( fftshift( fft2( fftshift( specyy ) ) ) );
                    covxy = real( fftshift( fft2( fftshift( specxy ) ) ) );
                    ki = (1:2*n)+(2*n)*(iGs-1);
                    kj = (1:2*n)+(2*n)*(jGs-1);
                    cxx = covxx(b,b).*p;
                    cyy = covyy(b,b).*p;
                    cxy = covxy(b,b).*p;
                    covMat(ki,kj) = [cxx cxy; cxy cyy];
                    if iGs~=jGs
                        covMat(kj,ki) = [rot90(cxx,2) rot90(cxy,2); rot90(cxy,2) rot90(cyy,2)];
                    end
                end
            end
            covMat = obj.fr0(kLayer)*(constants.radian2arcsec)^2*covMat;
        end
        
        function covMat = getCovarianceMatrixSingleLayerFFT(obj,kLayer)
            %% getCovarianceMatrixSingleLayerFFT
            %
            % updated by Y.O. 15/1/2017
            covMap = getCovarianceMapSingleLayerFFT(obj,kLayer);
            nSub = obj.wfs(1).lenslets.nLenslet;
            nm=size(covMap,1)/(2*obj.nwfs);
            covMap=mat2cell(covMap,ones(2*obj.nwfs,1)*nm,ones(1,2*obj.nwfs)*nm);
            covMat = cellfun(@(x) tools.covMap2Matrix(x,nSub,nSub), covMap, 'UniformOutput', false);
            covMat = cell2mat(covMat);
            covMat(~repmat(obj.wfs(1).validLenslet(:),2*obj.nwfs,1),:)=[];
            covMat(:,~repmat(obj.wfs(1).validLenslet(:),2*obj.nwfs,1))=[];
            
        end
        
        function covMat = separationToCovarianceMatrix(obj,u,v,fr0l,hl,L0l,igs,jgs)
            %Edge size
            g1     = 1 - hl/obj.gs(igs).height;
            g2     = 1 - hl/obj.gs(jgs).height;
            size1  = g1*obj.pitch;
            size2  = g2*obj.pitch;
            
            %Mid-points distance
            d1     = size1/2;
            d2     = size2/2;
            ac     = d1 - d2;
            ad     = d1 + d2;
            bc     = -ad;
            bd     = -ac;
            
            %Covariance matrix computation
            dphi2  = @(x,y) obj.dphi(hypot(x,y),L0l);
            covxx =-dphi2(u+ac,v) + dphi2(u+ad,v) + dphi2(u+bc,v) - dphi2(u+bd,v);
            
            covyy =-dphi2(u,v+ac) + dphi2(u,v+ad) + dphi2(u,v+bc) - dphi2(u,v+bd);
            
            covxy =-dphi2(u+d1,v-d2) + dphi2(u+d1,v+d2) + dphi2(u-d1,v-d2)...
                - dphi2(u-d1,v+d2);
            %Concatenation
            covMat                                = zeros(obj.ns);
            covMat(1:obj.ns/2,1:obj.ns/2)         = covxx;
            covMat(obj.ns/2+1:end,1+obj.ns/2:end) = covyy;
            covMat(1:obj.ns/2,1+obj.ns/2:end)     = covxy;
            covMat(1+obj.ns/2:end,1:obj.ns/2)     = covxy;
            %Multiplying by r0^-5/3
            covMat = covMat.*fr0l;
        end
        
        function [X,Y] = getSubapertureCoordinates(obj,kLayer,iGs)
            %Defining the grid
            D     = obj.tel.D;
            nL    = obj.wfs.lenslets.nLenslet;
            d     = D/nL;
            x     = linspace(-D/2.+d/2.,D/2.-d/2., nL);
            [X,Y] = meshgrid(x,x);
            %Pupil magnification
            X     = obj.G(iGs).*X;
            Y     = obj.G(iGs).*Y;
            %Pupil rotation
            X     = X.*cos(obj.th(iGs)) + Y.*sin(obj.th(iGs));
            Y     = Y.*cos(obj.th(iGs)) - X.*sin(obj.th(iGs));
            %Pupils shifts
            X     = X + obj.dX(iGs);
            Y     = Y + obj.dY(iGs);
            %Keep valid lenslets
            idx   = obj.wfs.validLenslet;
            X     = X(idx);
            Y     = Y(idx);
            %Magnification factor from the cone effect
            g     = 1 - obj.h(kLayer)/obj.gs(1).height;
            %Separation
            X     = g.*X + obj.h(kLayer)*obj.xSrc(iGs);
            Y     = g.*Y + obj.h(kLayer)*obj.ySrc(iGs);
        end
        function [dH,Hmax]  = getTomographicResolution(obj)
            xS = obj.xSrc;
            yS = obj.ySrc;
            dH   = obj.tel.D/obj.wfs.lenslets.nLenslet/max(hypot(xS,yS));
            Hmax = obj.tel.D/min(hypot(xS,yS));
        end
        
        %% get covariance map
        function cov = getCovarianceMap(obj)
            
            cov = cell(1,nchoosek(obj.nwfs,2)+obj.isAuto);
            
            % auto-correlation
            if obj.isAuto
                if obj.isFFT
                    cov{1} = obj.getCovarianceMapSinglePairFFT(1,1);
                else
                    cov{1} = obj.getCovarianceMapSinglePair(1,1);
                end
            end
            
            % add noise
            %nSub = obj.wfs(1).lenslets.nLenslet;
            %if ~obj.rmDiagonal
            %    cov{1}(nSub,nSub) = cov{t}(nSub,nSub) + obj.noiseVar(2*iGs-1);
            %    cov{1}(3*nSub-1,nSub) = cov{t}(3*nSub-1,nSub) + obj.noiseVar(2*iGs);
            %end
            
            % cross-correlation
            t = obj.isAuto;
            for iGs = 1:obj.nwfs
                for jGs = iGs+1:obj.nwfs
                    t = t+1;
                    if obj.isFFT
                        cov{t} = obj.getCovarianceMapSinglePairFFT(iGs,jGs);
                    else
                        cov{t} = obj.getCovarianceMapSinglePair(iGs,jGs);
                    end
                end
            end
            cov = cell2mat(cov);
            
            % ground layer subtraction 
            if obj.twoSteps
                cov = obj.groundLayerFilter4Map(cov);
            end
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILTERING MATRICES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function covMat = manageTipTilt(obj,covMat)
            
            F  = obj.tiptiltFilter();
            n = obj.ns;
            for iGs = 1:obj.nwfs
                ri = 1 + (iGs-1)*n:iGs*n;
                for jGs = 1:iGs
                    rj  = 1 + (jGs-1)*n:jGs*n;
                    tmp = F*covMat(ri,rj)*F';
                    covMat(ri,rj) = tmp;
                    if jGs~=iGs
                        covMat(rj,ri) = tmp';
                    end
                end
            end
        end
        function out = tiptiltFilter(obj)
            out                            = -(2/obj.ns).*ones(obj.ns);
            out(1:obj.ns/2,1+obj.ns/2:end) = 0;
            out(1+obj.ns/2:end,1:obj.ns/2) = 0;
            out                            = eye(obj.ns) + out;
        end
        function out = avgSlopesMatrix(obj)
            out                            = (2/obj.ns).*ones(obj.ns);
            out(1:obj.ns/2,1+obj.ns/2:end) = 0;
            out(1+obj.ns/2:end,1:obj.ns/2) = 0;
        end
        function out = commonMatrix(obj)
            out = zeros(obj.nmeas);
            n  = obj.ns;
            
            for iGs = 1:obj.nwfs
                ri = 1 + (iGs-1)*n:iGs*n;
                for jGs = 1:iGs
                    rj  = 1 + (jGs-1)*n:jGs*n;
                    out(ri,rj) = eye(n);
                    if jGs~=iGs
                        out(rj,ri) = out(ri,rj);
                    end
                end
            end
            out = out./obj.nwfs;
        end
        function out = groundLayerFilter(obj)
            %Getting tip-tilt filter
            TTfilter = zeros(obj.nmeas);
            for iGs = 1:obj.nwfs
                ri = 1 + (iGs-1)*obj.ns:iGs*obj.ns;
                TTfilter(ri,ri) = obj.tiptiltFilter();
            end
            %%Getting average slopes
            GLkeeper = obj.commonMatrix();
            %%Filtering the HO ground layer
            out = eye(obj.nmeas) - GLkeeper*TTfilter;
        end
        function cov = groundLayerFilter4Map(obj,cov)
            % parameter
            nsub = obj.wfs(1).lenslets.nLenslet;
            n = 2*nsub-1;
            mask = obj.wfs(1).validLenslet;
            % number of map
            n1 = 2;
            if obj.isXY
                n1 = 3;
            end
            n2 = nchoosek(obj.nwfs,2)+1;
            % initialize ground component
            gl1 = cell(obj.nwfs,n1);
            gl2 = cell(obj.nwfs,n1);
            gl3 = cell(1,n1);
            for t=1:n1
                for iGs=1:obj.nwfs
                    gl1{iGs,t} = zeros(n);
                    gl2{iGs,t} = zeros(n);
                end
                gl3{1,t} = zeros(n);
            end
            % compute ground component
            cov = mat2cell(cov,ones(n1,1)*size(cov,1)/n1,ones(1,n2)*size(cov,2)/n2);
            k=1;
            for iGs = 1:obj.nwfs
                for jGs = iGs:obj.nwfs
                    if iGs == jGs %auto-correlation part
                        tmp = tools.tipTiltFilter4covMap(cov{t,1},nsub,mask);
                    else % cross-correlation part
                        k=k+1;
                        tmp = tools.tipTiltFilter4covMap(cov{t,k},nsub,mask);
                    end
                    for t=1:n1
                        gl1{iGs,t} = gl1{iGs,t} + tmp;
                        gl2{jGs,t} = gl2{jGs,t} + tmp;
                        gl3{1,t} = gl3{1,t} + tmp;
                        
                        if iGs~=jGs
                            tmp = rot90(tmp,2);
                            gl1{jGs,t} = gl1{jGs,t} + tmp;
                            gl2{iGs,t} = gl2{iGs,t} + tmp;
                            gl3{1,t} = gl3{1,t} + tmp;
                        end
                    end
                end
            end
            for t=1:n1
                for iGs=1:obj.nwfs
                    gl1{iGs,t} = gl1{iGs,t}/3;
                    gl2{iGs,t} = gl2{iGs,t}/3;
                end
                gl3{1,t} = gl3{1,t}/9;
            end
            
            % subtraction
            for t=1:n1
                cov{t,1} = cov{t,1} - gl1{1,t} - gl2{1,t} + gl3{1,t};
            end
            k=1;
            for iGs = 1:obj.nwfs
                for jGs = iGs+1:obj.nwfs
                    k=k+1;
                    for t=1:n1
                        cov{t,k} = cov{t,k} - gl1{iGs,t} - gl2{jGs,t} + gl3{1,t};
                    end
                end
            end
            
            cov = cell2mat(cov);
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TOMOGRAPHIC ERROR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function idx = getScienceIndex(obj)
            tmp = obj.xSrc == 0 & obj.ySrc == 0 & [obj.gs.height] ==Inf;
            isaSci = sum(tmp);
            
            if isaSci
                idx = find(tmp);
            else
                idx = 0;
            end
        end
        function out = slindex(obj,idx)
            if idx<0
                out = obj.slindex(abs(idx) + 1 ) - 1;
            else
                out = 1 + obj.ns*(idx-1);
            end
        end
        function out = slrange(obj,idx)
            out = obj.slindex(idx):obj.slindex(-idx);
        end
        function out = norange(obj,idx)
            all = linspace(1,obj.nmeas,obj.nmeas);
            
            %if isscalar(idx)
            %    idx = [idx];
            %end
            
            for i = 1:length(idx)
                k = idx(i);
                all(obj.slrange(k)) = 0;
            end
            
            nn = find(all);
            if isvector(nn)
                if sum(diff(nn)) == length(diff(nn))
                    out = nn(1):nn(end);
                else
                    out = nn;
                end
            else
                out = [];
            end
        end
        function R = getMMMSEreconstructor(obj,covMat)
            
            %Splitting
            idxSci = obj.getScienceIndex();
            if idxSci~=0
                tmpCsg  = covMat(obj.slrange(idxSci),obj.norange(idxSci));
                tmpCgg  = covMat(obj.norange(idxSci),obj.norange(idxSci));
                %Inverting the cross-covariance matrix
                cond = 30;
                mCgg = pinv(tmpCgg,obj.getInversionTolerance(tmpCgg,cond));
                %Computing the MMSE reconstructor
                R    = tmpCsg*mCgg';
            else
                R = 0;
            end
            obj.R = R;
        end
        function Cee = getErrorCovarianceMatrix(obj,covMat,R,idxSci)
            %Splitting
            if nargin <4
                idxSci = obj.getScienceIndex();
            end
            if idxSci~=0
                obj.Csg  = covMat(obj.slrange(idxSci),obj.norange(idxSci));
                obj.Cgs  = covMat(obj.norange(idxSci),obj.slrange(idxSci));
                obj.Cgg  = covMat(obj.norange(idxSci),obj.norange(idxSci));
                obj.Css  = covMat(obj.slrange(idxSci),obj.slrange(idxSci));
                %Computing the covariance of the tomographic error
                Cee   = obj.Css - obj.Csg*R' - R*obj.Cgs + R*obj.Cgg*R';
            else
                Cee = 0;
            end
            obj.Cee = Cee;
        end
        function wfeTomo = getWaveFrontError(obj,Cee)
            if isempty(obj.Projector)
                wfeTomo = 0;
            else
                Cpp         = obj.Projector*Cee*obj.Projector';
                wfeTomo     = sqrt(trace(Cpp)).*obj.ProjectorUnits;
            end
            obj.wfeTomo = wfeTomo;
        end
        function psf = getPSF(obj,R,modes,npsf,overSampling,method,idxSci)
            if nargin < 7
                idxSci = obj.getScienceIndex();
            end
            %Get the covariance of open-loop slopes
            openCov = obj.getCovarianceMatrix();
            % Grabbed the residual covariance           
            resCov = obj.getErrorCovarianceMatrix(openCov,R,idxSci);
            % Propagation through the projector
            Mc  = obj.Projector*obj.ProjectorUnits;
            Cvv = Mc*resCov*Mc';
            % Get the OTF
            otf = psfrTools.modes2Otf(Cvv,modes,obj.tel.pupil,npsf,overSampling,method);
            % Fourier transform the OTF to get the psf
            psf = fourierTools.otf2psf(otf);
            psf = psf/sum(psf(:));
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISPLAY TOOLS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function map = getDeconvolvedMap(obj,Mat,iGs,jGs)
            %Grabbing maps
            autoi = obj.displaySingleMap(Mat(obj.slrange(iGs),obj.slrange(iGs)));
            autoj = obj.displaySingleMap(Mat(obj.slrange(jGs),obj.slrange(jGs)));
            corij = obj.displaySingleMap(Mat(obj.slrange(iGs),obj.slrange(jGs)));
            %Taking the x-part only
            n     = 2*obj.wfs.lenslets.nLenslet - 1;
            idx   = 1:n;
            axx   = 0.5*autoi(idx,idx) + 0.5*autoj(idx,idx);
            cxx   = corij(idx,idx);
            %Deconvolving
            map   = fftshift(ifft2(fft2(cxx)/(fft2(axx))));
        end
        function SNR = getSNRProfile(obj,Mat)
            %Grabbing all maps
            map  = obj.displayAllMap(Mat);
            %Applying the transformation
            F    = obj.groundLayerFilter();
            Mat2 =  F*Mat*F';
            map2 = obj.displayAllMap(Mat2);
            
            %Stacking cross-covariance map
            %nx    = obj.ns/2;
            nL    = obj.wfs.lenslets.nLenslet;
            n     = 2*nL - 1;
            nw    = obj.nwfs;
            dp    = obj.tel.D/nL;
            % Getting GS directions
            d     = [obj.gs.directionVector];
            d     = d(1:2,:);
            %Coordinates
            x      = -nL+1:nL-1;
            [X,Y]  = meshgrid(x);
            Y      = -Y;
            %Allocation in memory
            nLayer = 2;
            z      = [0.,10].*1e3;%linspace(0,20000,nLayer);
            SNR    = zeros(nw,nLayer,2);
            k      = 0;
            for i = 1:nw-1
                for j = i+1:nw
                    %Indexes
                    idxi = 1 + (i-1)*2*n:(2*i-1)*n;
                    idxj = 1 + (j-1)*2*n:(2*j-1)*n;
                    tmp  = map(idxj,idxi);
                    tmp2 = map2(idxj,idxi);
                    sep  = hypot(d(1,i) - d(1,j),d(2,i) - d(2,j));
                    m    = d(:,i) - d(:,j);
                    theta = atan(m(2)/m(1));
                    k = k+1;
                    
                    for l = 1:nLayer
                        %Grabbing the polar position of the peak
                        r       = z(l)*sep/dp;
                        %Converting into Cartesian positions
                        x0      = r*sin(theta)*sign(theta);
                        y0      = r*cos(theta)*sign(theta);
                        dx      = abs(X-ceil(x0));
                        dy      = abs(Y-ceil(y0));
                        %Define more precisely the peak position
                        box     = tmp.*(dx<=3 & dy<=3);
                        [ap,bp] = find(box==max(box(box~=0)));
                        xp      = X(ap,bp);
                        yp      = Y(ap,bp);
                        %Defining signal and noise
                        dx      = abs(X-xp);
                        dy      = abs(Y-yp);
                        mskSig  = dx <=2 & dy<=1;
                        sig     = tmp.*mskSig;
                        sig2    = tmp2.*mskSig;
                        %Defining noise
                        mskN    = ~mskSig;%(tmp.*(abs(tmp)<=max(tmp(:))/4));
                        nn      = tmp.*mskN;
                        nn2     = tmp2.*mskN;
                        %COmputing the SNR
                        tmpSNR  = max(sig(:))^2/rms(nn(:))^2;
                        tmpSNR2 = max(sig2(:))^2/rms(nn2(:))^2;
                        SNR(k,l,:)  = [tmpSNR,tmpSNR2];
                    end
                end
            end
        end
        function map = displayAllMap(obj,Mat)
            
            % Modified by Y.O.
            % using tools.covMatrix2Map to convert matrix to map.
            nL    = obj.wfs.lenslets.nLenslet;
            mask = obj.wfs.validLenslet;
            nm=size(Mat,1)/(2*obj.nwfs);
            Mat=mat2cell(Mat,ones(2*obj.nwfs,1)*nm,ones(1,2*obj.nwfs)*nm);
            map=cellfun(@(x) tools.covMatrix2Map(x,nL,nL,'mask1',mask,'mask2',mask), Mat, 'UniformOutput', false);
            map=cell2mat(map);
            
            %{
            nx = obj.ns/2;
            nL = obj.wfs.lenslets.nLenslet;
            n  = 2*nL - 1;
            nV = obj.wfs.nValidLenslet;
            tx = size(Mat,1)/nV/2;
            ty = size(Mat,2)/nV/2;
            
            %Computing the covariance map
            map = zeros(2*n*tx,2*n*ty);
            for i = 1:tx
                for j = 1:ty
                    indi1 = (i-1)*n*2 + 1:i*n*2;
                    indj1 = (j-1)*n*2 + 1:j*n*2;
                    indi2 = (i-1)*nx*2 + 1:i*nx*2;
                    indj2 = (j-1)*nx*2 + 1:j*nx*2;
                    Mat2  = Mat(indi2,indj2);
                    map(indi1,indj1) = obj.displaySingleMap(Mat2);
                end
            end
            %}
        end
        function map = displaySingleMap(obj,Mat)
            
            % Modified by Y.O.
            % using tools.covMatrix2Map to convert matrix to map.
            nL    = obj.wfs.lenslets.nLenslet;
            mask = obj.wfs.validLenslet;
            nm=size(Mat,1)/2;
            Mat=mat2cell(Mat,[nm;nm],[nm nm]);
            map=cellfun(@(x) tools.covMatrix2Map(x,nL,nL,'mask1',mask,'mask2',mask), Mat, 'UniformOutput', false);
            map=cell2mat(map);
            
            %{
            nx    = obj.ns/2;
            nL    = obj.wfs.lenslets.nLenslet;
            n     = 2*nL - 1;
            %Splitting the covariance matrix into four parts
            Matxx = Mat(1:nx,1:nx);
            Matxy = Mat(1:nx,1+nx:end);
            Matyy = Mat(1+nx:end,1+nx:end);
            Matyx = Mat(1+nx:end,1:nx);
            
            %Computing the covariance map
            map                  = zeros(2*n);
            map(1:n,1:n)         = getCovarianceMap(obj,Matxx);
            map(1:n,1+n:end)     = getCovarianceMap(obj,Matxy);
            map(1+n:end,1+n:end) = getCovarianceMap(obj,Matyy);
            map(1+n:end,1:n)     = getCovarianceMap(obj,Matyx);
            %}
        end
        % Please use tools.covMatrix2Map
        %{
function map = getCovarianceMap(obj,Mat)
           
            %Defining the valid sub-aperture
            nL    = obj.wfs.lenslets.nLenslet;
            idx   = obj.wfs.validLenslet;
            %Shifts arrays
            xx      = linspace(1,nL,nL);
            [XX,YY] = meshgrid(xx);
            xx      = XX(idx);
            yy      = YY(idx);
            [XX,YY] = meshgrid(xx,yy);
            dx      = XX - XX' + nL; %from shift to index
            dy      = YY - YY' + nL;
            %Map
            map     = zeros(2*nL-1);
            div     = map;
            ind     = dx(:) + (2*nL-1).*(dy(:) - 1);
            linMat  = Mat(:);
            
            for i=1:length(ind)
               map(ind(i)) = map(ind(i)) + linMat(i);
               div(ind(i)) = div(ind(i)) + 1;
            end
            
            div(div == 0) = 1;
            map           = map./div;obj.sys.nGs+ 1

        end
        %}
    end
    
    methods (Static)
        function out = dphi(rho,L0,b)
            if nargin <4
                b = 11/3;
            end
            
            %Computing the variance
            if 1
                L0ratio    = (L0).^(5./3);
                v          = (24.*gamma(6./5)./5).^(5./6).* ...
                    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0ratio;
                %Computing the phase covariance
                cov        = ones(size(rho));
                index      = rho~=0;
                u          = 2.*pi.*rho(index)./L0;
                cov(index) = 2^(1./6)/gamma(5./6)*u.^(5./6).*besselk(5./6,u);
                %Computing the phase structure function
                out        = 2.*v.*(1 - cov);
            else
                gb  = 2^(b-1)*gamma(b/2+1)^2*gamma(b/2+2)/gamma(b/2)/gamma(b+1);
                out = gb*rho.^(b-2);
            end
        end
        function out = spectrum(f,L0)
            %% SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)));
            out = out.*(f.^2 + 1./L0.^2).^(-11./6);
        end
        function tol = getInversionTolerance(Cgg,cond)
            l    = svd(Cgg);
            tol  = max(l(:))/cond;
        end
    end
end
