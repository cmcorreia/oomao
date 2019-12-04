classdef fittingMatrix < handle
    % fittingMatrix Create the fitting matrix for MCAO system from a
    % reconstructor MMSE and a dmStack (multiDM);
    properties
        %the reconstructor
        mmse;
        %the dms Stack
        multiDm;
        %dms modes low resolution
        modesLowRes;
        %dms resolution
        resolution;
        %Optimization stars for the MCAO projector
        optStars;
        %Optimization weight for the MCAO projector
        optWeight;
        %Propagators from science to pupil through atmosphere
        Hss;
        %Propagators from science to pupil through dm Stack
        Hdm;
        %Optimal projector for MCAO at DM altitudes (fitting step)
        %Hopt;
        %Influence function of the dms;
        dmInfFuncMatrix;
        %Inverse Influence function of the dms;
        iFCell;
        iDmInfFuncMatrix;
        %Command Matrix;
        commandMatrix;
        %Number of pixel for each Dm Layer
        layersNPixel;
        %Dm diameter
        D_m
    end
    
    methods
        %% Constructor
        function obj = fittingMatrix(multiDm,mmse,optStars,modesLowRes,varargin)
            p = inputParser;
            addRequired(p,'multiDm', @(x) isa(x,'multiDeformableMirror'));
            addRequired(p,'mmse', @(x) isa(x,'sparseLMMSE'));
            addRequired(p,'optStars', @(x) isa(x,'source'));
            addOptional(p,'modesLowRes', [], @(x) isnumeric(x) || ...
                (isa(x,'influenceFunction') || isa(x,'gaussianInfluenceFunction') ...
                || isa(x,'splineInfluenceFunction') || isa(x,'zernike') ...
                || isa(x,'hexagonalPistonTipTilt')) || isa(x,'differentialGaussianInfluenceFunction'));
            addParameter(p,'optWeight',[], @isnumeric);
            addParameter(p,'resolution',[], @isnumeric);
            parse(p,multiDm,mmse,optStars,modesLowRes,varargin{:});
            obj.mmse        = p.Results.mmse;
            obj.multiDm     = p.Results.multiDm;
            obj.modesLowRes = p.Results.modesLowRes;
            obj.optStars    = p.Results.optStars;
            nStar           = numel(obj.optStars);
            obj.optWeight   = p.Results.optWeight;
            if isempty(obj.optWeight)
                obj.optWeight = 1/nStar*ones(nStar,1);
            end
            obj.resolution = p.Results.resolution;
            if isempty(obj.resolution)
                obj.resolution = (obj.mmse.nSub+1)*ones(obj.multiDm.nDm,1);
            end
            
            nAtmLayer           = obj.mmse.atmModel.nLayer;
            nDmLayer            = obj.multiDm.nDm;
            Hss                 = cell(nStar,nAtmLayer);
            Hdm                 = cell(nStar,nDmLayer);
            obj.dmInfFuncMatrix = [];
            obj.iDmInfFuncMatrix= [];
            obj.iFCell          = cell(nDmLayer,1);
            obj.layersNPixel    = zeros(nDmLayer,1);
            obj.D_m             = zeros(nDmLayer,1);
            nValidActuatorTotal = 0;
            
            for kDmLayer = 1:nDmLayer
                dm              = obj.multiDm.dms{kDmLayer};
                actCoord        = dm.modes.actuatorCoord;
                obj.D_m(kDmLayer) = max(real(actCoord(:)))-min(real(actCoord(:)));
                do              = obj.mmse.dSub;
                obj.layersNPixel(kDmLayer) = round(obj.D_m(kDmLayer)./do)+1;
                if isempty(obj.modesLowRes(kDmLayer).modes) 
                    dmLowRes        = deformableMirror(dm.nActuator,'modes',obj.modesLowRes(kDmLayer),...
                        'resolution',obj.layersNPixel(kDmLayer),...
                        'validActuator',dm.validActuator);
                    F = 2*dmLowRes.modes.modes;
                else
                    F = obj.modesLowRes(kDmLayer).modes;
                end
                
                obj.dmInfFuncMatrix = blkdiag(obj.dmInfFuncMatrix,F);
                iF = pinv(full(F));
                obj.iFCell{kDmLayer,1} = iF;
                obj.iDmInfFuncMatrix = blkdiag(obj.iDmInfFuncMatrix,iF);
                nValidActuatorTotal = nValidActuatorTotal + dm.nValidActuator;
            end
            
            resTotal   = sum(obj.layersNPixel.^2);
            HdmMean    = sparse(nValidActuatorTotal,nValidActuatorTotal);
            HdmAtmMean = sparse(nValidActuatorTotal,length(mmse.Hss(1,:)));
            %PdmMean    = sparse(resTotal,resTotal);
            %PdmAtmMean = sparse(resTotal,length(mmse.Hss(1,:)));
            
            fprintf('___ BI-LINEAR INTERPOLATION OPERATOR ___\n')
            for kGs=1:nStar
                %Hss propagator
                for kAtmLayer = 1:nAtmLayer
                    pitchAtmLayer = obj.mmse.dSub/obj.mmse.overSampling(kAtmLayer);
                    height        = obj.mmse.atmModel.layer(kAtmLayer).altitude;
                    % pupil center in layer
                    beta  = obj.optStars(kGs).directionVector*height;
                    scale = 1-height/obj.optStars(kGs).height;
                    Hss{kGs,kAtmLayer} = p_bilinearSplineInterp(...
                        obj.mmse.atmGrid{kAtmLayer,1}(:),...
                        obj.mmse.atmGrid{kAtmLayer,2}(:),...
                        pitchAtmLayer,...
                        real(obj.mmse.outputPhaseGrid)*scale+beta(1),...
                        imag(obj.mmse.outputPhaseGrid)*scale+beta(2));
                end
                %Hdm propagator
                for kDmLayer = 1:nDmLayer
                    pitchDmLayer = obj.mmse.dSub;
                    height       = obj.multiDm.zLocations(kDmLayer);
                    % pupil center in layer
                    beta  = obj.optStars(kGs).directionVector*height;
                    scale = 1-height/obj.optStars(kGs).height;
                    actCoord = obj.multiDm.dms{kDmLayer}.modes.actuatorCoord;
                    dmin = min(real(actCoord(:)));
                    dmax = max(real(actCoord(:)));
                    Dx = (obj.layersNPixel(kDmLayer)-1)*pitchDmLayer;
                    sx = dmin-(Dx-(dmax-dmin))/2;
                    [x,y] = meshgrid(linspace(0,1,obj.layersNPixel(kDmLayer))*Dx+sx);
                    Hdm{kGs,kDmLayer} = bilinearSplineInterp(...
                        x,...
                        y,...
                        pitchDmLayer,...
                        real(obj.mmse.outputPhaseGrid)*scale+beta(1),...
                        imag(obj.mmse.outputPhaseGrid)*scale+beta(2));
                end
                intDM = [Hdm{kGs,:}];
                intL  = [Hss{kGs,:}];
                HdmMean    = HdmMean    + (intDM*obj.dmInfFuncMatrix)'*(intDM*obj.dmInfFuncMatrix);
                HdmAtmMean = HdmAtmMean +  (intDM*obj.dmInfFuncMatrix)'*intL;
                %                 HdmMean    = HdmMean    + obj.optWeight(kGs,1)*(intDM*obj.dmInfFuncMatrix)'*(intDM*obj.dmInfFuncMatrix);
                %                 HdmAtmMean = HdmAtmMean +  obj.optWeight(kGs,1)*(intDM*obj.dmInfFuncMatrix)'*intL;
                %PdmMean    = PdmMean    +  intDM'*intDM;
                %PdmAtmMean = PdmAtmMean +  intDM'*intL;
            end
            HdmMean    = HdmMean/nStar;
            HdmAtmMean = HdmAtmMean/nStar;
            %PdmMean    = PdmMean/nStar;
            %PdmAtmMean = PdmAtmMean/nStar;
            obj.Hss  = Hss;
            obj.Hdm  = Hdm;
            %obj.Hopt = pinv(full(PdmMean)+1e-3*eye(resTotal))*PdmAtmMean;
            obj.commandMatrix = pinv(full(HdmMean),1)*HdmAtmMean;
            
            % local function
            function P = p_bilinearSplineInterp(xo,yo,do,xi,yi)
                
                ni = length(xi);
                
                nxo = length(xo);
                nyo = length(yo);
                no = nxo*nyo;
                
                % remove the interaporated points out of the original grid
                mask = xi>=xo(1) & yi>=yo(1) & xi<=xo(end) & yi<=yo(end);
                
                % index for the inteporated grid
                index = (1:ni)';
                index = index(mask);
                
                % x & y index for the original grid
                ox = floor((xi-xo(1))/do)+1;
                ox = ox(mask);
                oy = floor((yi-yo(1))/do)+1;
                oy = oy(mask);
                
                % bi-linear inteporation value
                fxo = abs(xi(mask)-(xo(1)+do*(ox-1)))/do;
                fyo = abs(yi(mask)-(yo(1)+do*(oy-1)))/do;
                s1 = (1-fxo).*(1-fyo);
                s2 = fxo.*(1-fyo);
                s3 = (1-fxo).*fyo;
                s4 = fxo.*fyo;
                
                % vectoraized index for the original grid
                o1 = oy+nyo*(ox-1);
                o2 = oy+nyo*ox;
                o3 = oy+1+nyo*(ox-1);
                o4 = oy+1+nyo*ox;
                
                % masking
                o1(s1==0)=[];
                i1 = index(s1~=0);
                s1(s1==0)=[];
                o2(s2==0)=[];
                i2 = index(s2~=0);
                s2(s2==0)=[];
                o3(s3==0)=[];
                i3 = index(s3~=0);
                s3(s3==0)=[];
                o4(s4==0)=[];
                i4 = index(s4~=0);
                s4(s4==0)=[];
                
                % intepolation matrix
                P1 = sparse(i1,o1,s1,ni,no);
                P2 = sparse(i2,o2,s2,ni,no);
                P3 = sparse(i3,o3,s3,ni,no);
                P4 = sparse(i4,o4,s4,ni,no);
                
                P = P1+P2+P3+P4;
            end
        end
    end
end