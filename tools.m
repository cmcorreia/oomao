classdef tools < utilities & lamTools & psfrTools
    
    methods (Static)
        %% distances
        function rho = dist(l,u, nPts)
            [x,y] = meshgrid(linspace(l,u,nPts));
            z = complex(x,y);
            rho  = abs(bsxfun(@minus,z(:),z(:).'));
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- CN2H COMPRESSING ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate compressed profile using mean-weighted compression
        
        function [Cn2eq, altEq] = eqLayers(Cn2, altitudes, nEqLayers, power,method)
            %{
            Cn2         ::  The input Cn2 profile (vector)
            altitudes   ::  The input altitudes (vector)
            nEqLayers   ::  The number of output equivalent layers (scalar)
            power       ::  the exponent of the turbulence (default 5/3)
            
            See: Saxenhuber17: Comparison of methods for the reduction of
            reconstructed layers in atmospheric tomography, App Op, Vol. 56, No. 10 / April 1 2017
            %}
            
            nCn2 = numel(Cn2);
            nAltitudes = numel(altitudes);
            
            if nCn2 == nEqLayers
                Cn2eq = Cn2;
                altEq   = altitudes;
            else
                if nargin < 4
                    power = 5/3;
                end
                if nargin < 5
                    method = 'mean-weight';
                end
                
                nSlab = floor(round(nCn2)/fix(nEqLayers));
                
                if strcmp(method, 'mean-weight')
                    ppp = 1;
                    posSlab =  round((linspace(0, nEqLayers-1, nEqLayers))*nSlab)+1;
                    for iii = 1:nEqLayers-1
                        if posSlab(iii) >= posSlab(iii+1)
                            posSlab(iii+1) = posSlab(iii)+1;
                        end
                    end
                    posSlab = [posSlab, nAltitudes+1];
                    Cn2eq = zeros(1,nEqLayers);
                    altEq = zeros(1,nEqLayers);
                    for ii = 1:nEqLayers
                        Cn2eq(ii) = sum(Cn2(posSlab(ii):posSlab(ii+1)-1));
                        altEq(ii) =  (sum(altitudes(posSlab(ii):posSlab(ii+1)-1) .^ (power) .* Cn2(posSlab(ii):posSlab(ii+1)-1))/Cn2eq(ii)) .^ (1./power);
                    end
                end
            end
        end
        
        %%
        function [Cn2eq] = eqProfile(atmSys,atmModel,Cn2eq0)
            %% COMPUTE THE EQUIVALENT PROFILE BY STIPULATING THE ALTITUDES
            %{
            EXAMPLE
            r0 = 0.16;
            L0 = 50;
            atm         = atmosphere(photometry.V0,r0,L0,...
                'fractionnalR0',[0.5170    0.1190    0.0630    0.0610    0.1050    0.0810    0.0540],...
                'altitude',[0         500        1000        2000        4000        8000       16000]);
            fprintf('Theta_0 of the original atmosphere: %2.3f arcsec\n', atm.theta0InArcsec)
            
            
            layersAltNew = [0, 8, 16]*1e3;
            nNewLayers = length(layersAltNew);
            for i = 1:nNewLayers
                atmModel.layer(i).altitude = layersAltNew(i);
            end
            
            [Cn2eq, altEq] = tools.eqLayers([atm.layer.fractionnalR0], [atm.layer.altitude], nNewLayers);
            atmModel         = atmosphere(photometry.V0,r0,L0,...
                'fractionnalR0',Cn2eq,'altitude',altEq);
            tools.eqProfile(atm, atmModel, Cn2eq)
            fprintf('Theta_0 of the new atmosphere: %2.3f arcsec\n', atmModel.theta0InArcsec)
            
            fprintf('Fractional strength of new atmosphere: %2.3f\n', [atmModel.layer.fractionnalR0])
            %}
            theta0InArcsec = atmSys.theta0InArcsec;
            fun = @(Cn2Model) tools.anisoplanaticAngleMinimiser(theta0InArcsec, atmModel, Cn2Model);
            Cn2Eq = fminsearch(fun,Cn2eq0);
            
        end
        
        %%
        function out = anisoplanaticAngleMinimiser(theta0InArcsec, atmModel, Cn2Model)
            
            nLayers = length(Cn2Model);
            Cn2Model(Cn2Model < 0) = 0;
            Cn2Model = Cn2Model/sum(Cn2Model);
            for iLayer = 1:nLayers
                atmModel.layer(iLayer).fractionnalR0 = Cn2Model(iLayer);
            end
            out = abs(theta0InArcsec - atmModel.theta0InArcsec)^2;
        end
        
        %%
        function imagerot = rotateIm(image,degree)
            switch mod(degree, 360)
                % Special cases
                case 0
                    imagerot = image;
                case 90
                    imagerot = rot90(image);
                case 180
                    imagerot = image(end:-1:1, end:-1:1);
                case 270
                    imagerot = rot90(image(end:-1:1, end:-1:1));
                    
                    % General rotations
                otherwise
                    % Convert to radians and create transformation matrix
                    a = degree*pi/180;
                    R = [+cos(a) +sin(a); -sin(a) +cos(a)];
                    
                    % Figure out the size of the transformed image
                    [m,n,p] = size(image);
                    dest = round( [1 1; 1 n; m 1; m n]*R );
                    dest = bsxfun(@minus, dest, min(dest)) + 1;
                    imagerot = zeros([max(dest) p],class(image));
                    
                    % Map all pixels of the transformed image to the original image
                    for ii = 1:size(imagerot,1)
                        for jj = 1:size(imagerot,2)
                            source = ([ii jj]-dest(1,:))*R.';
                            if all(source >= 1) && all(source <= [m n])
                                
                                % Get all 4 surrounding pixels
                                C = ceil(source);
                                F = floor(source);
                                
                                % Compute the relative areas
                                A = [...
                                    ((C(2)-source(2))*(C(1)-source(1))),...
                                    ((source(2)-F(2))*(source(1)-F(1)));
                                    ((C(2)-source(2))*(source(1)-F(1))),...
                                    ((source(2)-F(2))*(C(1)-source(1)))];
                                
                                % Extract colors and re-scale them relative to area
                                cols = bsxfun(@times, A, double(image(F(1):C(1),F(2):C(2),:)));
                                
                                % Assign
                                imagerot(ii,jj,:) = sum(sum(cols),2);
                            end
                        end
                    end
                    imagerot = tools.crop(imagerot,size(image));
            end
        end
        
        %% Invert a sparse measurement noise covariance matrix to use in a tomographic reconstructor
        function iCn = invertMeasurementNoiseCovarianceMatrix(input)
            if iscell(input)
                nGs = size(input,1);
                out = zeros([size(input{1}), nGs]);
                for iGs = 1:nGs
                    out(:,:,iGs) = input{iGs};
                end
            else
                nGs = size(input,3);
                out = input;
            end
            iCn = cell(nGs,1);
            nSlopes = size(out,1);
            for iGs = 1:nGs
                Cn = out(:,:,iGs);
                
                % Extract diagonal and cross-terms -> create sparse matrix
                noiseCovarDiag = diag(Cn);
                noiseCovarDiagP1 = diag(Cn,size(Cn,1)/2);
                B = zeros(size(Cn,1),3);
                B(:,1) = noiseCovarDiag;
                B(1:1:end/2,2) = noiseCovarDiagP1;
                B(end/2+1:1:end,3) = noiseCovarDiagP1;
                CnE = spdiags(B,[0,-nSlopes/2,nSlopes/2],nSlopes,nSlopes);
                iCn{iGs} = pinv(full(CnE));
                %iCn{iGs}(abs(iCn{iGs})< 1e-10) = 0;
                
                % extract only meaningful values on the main and two
                % cross-covar diagonals
                B(:,1) = diag(iCn{iGs});
                B(1:1:end/2,2) = diag(iCn{iGs},size(Cn,1)/2);
                B(end/2+1:1:end,3) = diag(iCn{iGs},size(Cn,1)/2);
                iCn{iGs} = spdiags(B,[0,-nSlopes/2,nSlopes/2],nSlopes,nSlopes);
            end
            iCn = blkdiag(iCn{:});
        end
                %% type I OMGI optimiser
        function g = omgi(nu, PSD, samplingFreq, pureDelay, noiseVar,verbose)
            if ~exist('verbose','var')
                verbose = 0;
            end
            T = 1/samplingFreq;
            s = @(x) 2*1i*pi*x;
            
            if size(PSD,1) > size(PSD,2)
                PSD = PSD';
            end
            
            % TFs
            hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T);
            hDac = hWfs; % Different meaning, same value
            hLag = @(x) exp(-pureDelay*s(x));
            hInt = @(x) 1./(1-exp(-s(x)*T));
            hDm  = @(x) 1;
            G = @(x,g) hWfs(x) .* hDac(x) .* hLag(x) .* g .* hInt(x) .* hDm(x);
            Gn = @(x,g) hDac(x) .* hLag(x) .* g .* hInt(x) .* hDm(x);
            % rejection transfer function
            E = @(x,g) abs(1./(1+G(x,g)));
            %Noise transfer
            hN = @(x,g) Gn(x,g) .* E(x,g);%./hWfs(x);
            %myfun = @(g) sum(hN(nu)*noiseVar + RTF*PSD);
            g = fminsearch(@(g)lamTools.integrate(g,nu,E,hN,PSD,noiseVar,T), 0.1);
            
            if verbose
                gains = linspace(0.001,0.5,100);
                for kGain = 1:length(gains)
                    outN(kGain) = trapz(nu, abs(hN(nu,gains(kGain))).^2*noiseVar*2*T);
                    outS(kGain) = trapz(nu, abs(E(nu,gains(kGain))).^2.*PSD);
                end
                figure
                semilogy(gains, outS,'r')
                hold on
                semilogy(gains, outN,'b')
                semilogy(gains, outS + outN,'k')
                xlabel('gain')
                ylabel('residual')
                title('single integrator optimal modal gian')
                legend('signal residual','noise residual','total residual')
                outNopt = trapz(nu, abs(hN(nu,g)).^2*noiseVar*2*T);
                outSopt = trapz(nu, abs(E(nu,g)).^2.*PSD);
                plot(g,outNopt + outSopt,'ro')
            end
        end
        
        function out = integrate(g,nu,E,hN, PSD, noiseVar,T)
            outN = trapz(nu, abs(hN(nu,g)).^2*noiseVar*2*T);
            outS = trapz(nu, abs(E(nu,g)).^2.*PSD);
            out = outS + outN;
        end
        
        %%
        function out = doubleIntParamOptim(nu, PSD, samplingFreq, pureDelay, varNoise)
            %% DOUBLE INTEGRATOR PARAMETER OPTIMISATION
            
            % out = doubleIntParamOptim(nu, PSD, samplingFreq, pureDelay,
            % varNoise) computes the double integrator gain and lead-filter
            % parameters for a 45 phase margin stability.
            % nu            :: the temporal frequency vector
            % PSD           :: the temporal PSD vector in units^2
            % samplingFreq  :: the loop sampling frequency in Hz
            % pureDelay     :: the loop pure delay
            % varNoise      :: noise variance in units^2
            % Created       :: C. Correia, Dec'15
            % Comment: based on do_optim_typeII from HIA and TMT,
            % developed originally by JPVeran
            
            T = 1/samplingFreq;
            s = @(x) 2*1i*pi*x;
            
            g0 = 1e-3;
            varTotalOld = inf;
            
            % OPEN-LOOP TRANSFER FUNCTION
            
            %             G = @(x) ((1-exp(-s(x)*T))./(s(x)*T)).^2.*...   % hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T); hDac = hWfs; % Different meaning, same value
            %                 exp(-tau*s(x)).*...                         % hLag = @(x) exp(-tau*s(x));
            %                 (1./(1-exp(-s(x)*T)).^2).*...               % hInt = @(x) 1./(1-exp(-s(x)*T)); squared for the double integrator
            %                 1;                                          % hDm  = @(x) 1;
            
            hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T);
            hDac = hWfs; % Different meaning, same value
            hLag = @(x) exp(-pureDelay*s(x));
            hInt = @(x) 1./(1-exp(-s(x)*T));
            hDm  = @(x) 1;
            G = @(x) hWfs(x) .* hDac(x) .* hLag(x) .* hInt(x) .* hInt(x) .* hDm(x);
            gcc = [];
            
            while 1
                % Open loop gain without phase lead
                %hOl = @(x) g0*G(x);
                
                [phaseMargin, fcross] = calcPhiMargin(g0*G(nu), nu);
                
                % Phase lead needed
                additionalPhaseNeeded=pi/4-phaseMargin;
                
                % calculating the lead filter parameters to achieve the required phase lead
                a=(1-sin(additionalPhaseNeeded))/(1+sin(additionalPhaseNeeded));
                f0=fcross*sqrt(a); %fprintf('fs=%g, fcross=%f\n',fs,fcross);
                Tlead=1/(2*pi*f0);
                
                % gain 1/g created by phase lead filter
                g=sqrt(a);
                
                % complete Hol. g is here to adjust g0 to remove the scaling caused by lead filter.
                %hOl1 = @(x) g * hOl(x) .* (1+Tlead*s(x))./(1+a*Tlead*s(x));
                hOl1 = @(x) g * g0*G(x) .* (1+Tlead*s(x))./(1+a*Tlead*s(x));
                % Rejection transfer function
                E = @(x) abs(1./(1+hOl1(x)));
                
                % Closed-loop transfer function
                %Ecl = @(x) hOl1(x) .* E(x);
                %Noise transfer
                hN = @(x) hOl1(x) .* E(x);%./hWfs(x);
                
                % RESIDUAL SIGNAL VARIANCE
                %rms=trapz(nus, PSD.*abs(Hrej.*Hrej));
                %varSignal = 2*quadgk( @(nu) PSD.*abs(E(nu)).^2 , 0 , Inf);
                varSignal = trapz(nu, PSD.*abs(E(nu)).^2);
                % NOISE GAIN
                %noiseGain=(trapz(nus, abs(Hn.*Hn))./trapz(nus, ones(length(nus),1)));
                %noiseGain = 2*quadgk( @(nu) PSD.*abs(hN(nu)).^2 , 0 , Inf)
                noiseGain = trapz(nu, abs(hN(nu)).^2)/(1/2/T);
                %Tot Error.
                varPropNoise=noiseGain*varNoise;
                varTotal=varSignal+varPropNoise;
                
                if varTotal < varTotalOld
                    %increase the gain.
                    g0=g0*1.01;
                    gcc = [gcc g0];
                    varTotalOld = varTotal;
                    out{1} = g0*g;
                    out{2} = a;
                    out{3} = Tlead;
                    out{4} = varSignal;
                    out{5} = varPropNoise;
                    out{6} = noiseGain;
                    out{7} = G;
                    out{8} = hN;
                    
                else
                    break
                end
            end
            function [margin, fcross]=calcPhiMargin(Hol, nus)
                %Computes Phase margin of the open loop transfer function.
                %Defined as the phase of Hol when abs(Hol) crosses 1.
                ind=abs(Hol)>1;
                indl=find(ind);
                indl=indl(end);
                if indl==length(Hol)
                    ph2=angle(Hol(end));
                    fcross=nus(end);
                else
                    abs1=log10(abs(Hol(indl:indl+1)));%interpolate the logrithm which is close to linear.
                    frac1=abs1(1)/(abs1(1)-abs1(2));%linear interpolation.
                    ph1=(angle(Hol(indl:indl+1)));
                    %Reduce distance
                    diff=(ph1(2)-ph1(1))/(2*pi);
                    diff=diff-fix(diff);
                    ph1(2)=ph1(1)+2*pi*diff;
                    
                    ph2=ph1(1)-(ph1(1)-ph1(2))*frac1;
                    
                    nu1=nus(indl:indl+1);
                    fcross=nu1(1)-(nu1(1)-nu1(2))*frac1;
                end
                %bring [-pi;pi] to [-2pi;0];
                normph=ph2/(2*pi);
                normph=normph-floor(normph)-1;
                ph2=normph*2*pi;
                if ph2>0
                    ph2=ph2-2*pi;
                end
                margin=pi+ph2;
            end
            
        end

                %% create a DM with pre-defined actuator layout
        function [px,py] = actuatorCoordinates(tel,option,pitch,inPupil)
            %% actuatorCoordinates for different DM models
            % [px py] = actuatorCoordinates(tel,option,pitch,inPupil)
            % Computes the actuator coordinates for different DM models for a
            % given telescope as
            %    LBT
            %    square or Fried
            %    random
            %    hexagon
            %    triangle (scale downed M4)
            % a pitch in meters and either across the whole squared
            % computational domain set by tel.D or within the pupil
            % tel.pupil
            
            
            D = tel.D + 2*pitch;
            obscurationRatio = tel.obstructionRatio;
            if ~exist('inPupil','var') || isempty(inPupil)
                inPupil = 1;
            end
            switch option
                % LBT
                case 'LBT'
                    N = 12;     % Initial number of actuators in first ring
                    n = 3;      % Additional actuators per rings
                    nRings = 9; % Number of rings
                    xcentre = 0;
                    ycentre = 0;
                    px = 0;
                    py = 0;
                    
                    for ind=1:nRings
                        theta = 0:pi/(N/2 + ((ind-1)*n)):2*pi;
                        px = [px, pitch*ind*cos(theta(1:end-1)) + xcentre];
                        py = [py, pitch*ind*sin(theta(1:end-1)) + ycentre];
                    end
                    
                    if inPupil==1
                        % Select in Pupil
                        dist = sqrt(px.^2 + py.^2);
                        valid = dist <= D/2+pitch/2 & dist >= obscurationRatio * D/2  -pitch/2;
                        px = px(valid);
                        py = py(valid);
                    else
                        vx = find(abs(px) <= D/2);
                        vy = find(abs(py(vx)) <= D/2);
                        px = px(vx(vy));
                        py = py(vx(vy));
                    end
                    % Square
                case {'square', 'fried'}
                    [px, py] = meshgrid(-D/2:pitch:D/2);
                    px = reshape(px, 1, numel(px));
                    py= reshape(py, 1, numel(py));
                    
                    if inPupil==1
                        % Select in Pupil
                        dist = sqrt(px.^2 + py.^2);
                        valid = dist <= D/2+pitch/2 & dist >= obscurationRatio * D/2  -pitch/2;
                        px = px(valid);
                        py = py(valid);
                        %scatter(px, py)
                    else
                        vx = find(abs(px) <= D/2);
                        vy = find(abs(py(vx)) <= D/2);
                        px = px(vx(vy));
                        py = py(vx(vy));
                    end
                    
                    % Random
                case 'random'
                    oc = tel.obstructionRatio * tel.D/2;
                    pr = oc + (tel.D/2-oc).*rand(250,1);
                    ptheta = 2*pi*rand(250,1);
                    px = pr.*cos(ptheta);
                    py = pr.*sin(ptheta);
                    
                    % % Hexagon
                case {'hexagon', 'triangle'}
                    % Hexagon
                    if strcmp(option, 'hexagon')
                        x_hexagon=[-1 -0.5 0.5 1];
                        y_hexagon=[0 -sqrt(3)/2 -sqrt(3)/2 0];
                    end
                    
                    % Hexagon centre
                    if strcmp(option, 'triangle')
                        x_hexagon=[0, 1.5, -1, -0.5,      0.5,        1];
                        %y_hexagon=[0, -0.5,   0, -0.5, -0.5, 0];
                        y_hexagon=[0, -.9,   0, -.9, -.9, 0]; % for M4
                    end
                    
                    N = tel.D/pitch/2+2;
                    M = tel.D/(pitch*0.9)/2+2;
                    px = [];
                    py = [];
                    
                    for nn=0:N-1
                        for mm=0:M-1
                            px = [px, x_hexagon+1.5+3*nn;];
                            py = [py, y_hexagon+.9+.9*2*mm];
                        end
                    end
                    
                    % If number of actuators fitting into square is odd make sure y
                    % coordinate has centre at 0.  I.e. if pitch is exactly and EVEN
                    % division of D.
                    px = px - max(px(:))/2;
                    if (tel.D/pitch)/2==round((tel.D/pitch)/2)
                        py = py - max(py(:))/2 - 1/2;
                    else
                        py = py - max(py(:))/2;
                    end
                    px = px * pitch;
                    py = py * pitch;
                    
                    if inPupil==1
                        % Select in Pupil
                        dist = sqrt(px.^2 + py.^2);
                        valid = dist <= D/2+pitch/2 & dist >= obscurationRatio * D/2  -pitch/2;
                        px = px(valid);
                        py = py(valid);
                    else
                        vx = find(abs(px) <= D/2);
                        vy = find(abs(py(vx)) <= D/2);
                        px = px(vx(vy));
                        py = py(vx(vy));
                    end
            end
            % scatter plot
            %scatter(px,py)
            %box on
            %title('DM actuator locations')
            %xlabel('meters')
            %ylabel('meters')
            %axis equal tight;
            
        end
        
    end
end
