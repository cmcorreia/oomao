%% Pupil Class
% Create a telescope pupil giving a segment, their coordinates,
% phasing errors, reflexivity, etc etc...

classdef pupilGenerator < matlab.mixin.Copyable %% handle class + deep copy available via b=copy(a)

    properties
        % segment
        segRef;  % segment class
        segList;
        nSegments; % number of segment
        
        % pupil
        segmentCoord; % Coordinates of centers of each segments
        nPixels; % matrix size in px
        definition; % desired definition in Pixels
        matrix; % double matrix containing shape of pupil
        obstructionRatio; % given by d/D , must be 0 <...< 1
        coeffPhaseModes; % in m
        coeffReflexion; % between 0 & 1
        coeffDistance;% To see separations between segments and/or to fix round(coord) errors
        pixelRatio;
        SMC; % Area taken by all the centers in m
        centerX; % center of matrix in px
        centerY; % center of matrix in px
        spiders; % structure
        D;
        Dpx;
        
        % other
        unit;
        %wavelength; %% WARNING : should be fixed at 1 so that phase (=OPD in this case) remains in m
        phaseModesGlobal;
        nPhaseModes;
        phaseModesType;
        flagNoGap; % 1st method with a "donut" mask
        verbose;
    end
    
    methods
        %% Constructor
        function obj = pupilGenerator(varargin)
            % CHECKING INPUTS
            inputs = inputParser;
            addParameter(inputs,'segRef',segment(inf,1,100),@(x) isa(x,'segment'));
            addParameter(inputs,'segCoord',[0 0],@isnumeric);
            %addParameter(inputs,'wavelength',1,@isnumeric); %% WARNING : should be fixed at 1 so that phase (=OPD in this case)remains in nm
            addParameter(inputs,'definition',[],@isnumeric);
            addParameter(inputs,'D',[],@isnumeric);
            addParameter(inputs,'coeffPhaseModes',[],@isnumeric);
            addParameter(inputs,'coeffReflexion',[],@isnumeric);
            addParameter(inputs,'coeffDistance',1,@isnumeric);
            addParameter(inputs,'nPhaseModes',25,@isnumeric);
            addParameter(inputs,'phaseModesType','zernike',@isstring);
            addParameter(inputs,'obstructionRatio',0,@isnumeric);
            addParameter(inputs,'unit','m',@ischar);
            addParameter(inputs,'flagNoGap',0,@isnumeric);
            addParameter(inputs,'verbose',false,@islogical);
            addParameter(inputs,'spiders',struct('n',0,'angle',[],'width',0,'usrDefined',0,'usrMask',[]),@isstruct);
            inputs.parse(varargin{:});
            
            % INSTANTIATION
            obj.segRef               = inputs.Results.segRef;
            obj.segmentCoord         = inputs.Results.segCoord;
            %obj.wavelength           = inputs.Results.wavelength;
            obj.D                    = inputs.Results.D;
            obj.definition           = inputs.Results.definition;
            obj.coeffDistance        = inputs.Results.coeffDistance;
            obj.unit                 = inputs.Results.unit;
            obj.nPhaseModes          = inputs.Results.nPhaseModes;
            obj.coeffPhaseModes      = inputs.Results.coeffPhaseModes;
            obj.phaseModesType       = inputs.Results.phaseModesType;
            obj.coeffReflexion       = inputs.Results.coeffReflexion;
            obj.obstructionRatio   = inputs.Results.obstructionRatio;
            obj.spiders              = inputs.Results.spiders;
            obj.flagNoGap            = inputs.Results.flagNoGap;
            obj.verbose              = inputs.Results.verbose;
            
            % CREATE PHYSICAL DATA
            
            obj.nSegments=size(obj.segmentCoord,1);
            obj.pixelRatio=obj.coeffDistance*obj.segRef.pixelRatio;
            
            if strcmp(obj.unit,'m')
                obj.SMC=max(max(obj.segmentCoord))-min(min(obj.segmentCoord)); % size of area taken by the center
                maxOffset=0;
                for k=1:length(obj.segList)
                    maxOffset=max(maxOffset,max(obj.segList(k).posXError,obj.segList(k).posYError));
                end
                obj.SMC=obj.SMC+2*maxOffset;
                obj.nPixels=ceil(obj.SMC*(obj.pixelRatio))+obj.segRef.nPxInterp;
                %                 if isempty(obj.D)
                %                     obj.D=(max(obj.segmentCoord(:,2))-min(obj.segmentCoord(:,2))+2*obj.segRef(1).radius);
                %                 end
            elseif strcmp(obj.unit,'px')
                obj.nPixels=max(obj.segmentCoord(:))+round(obj.segRef.nPxInterp/2);
                obj.SMC=(obj.nPixels-obj.segRef.nPxInterp)/obj.pixelRatio;
                %                 if isempty(obj.D)
                %                     obj.D=((max(obj.segmentCoord(:,2))-min(obj.segmentCoord(:,2))+2*obj.segRef(1).radius))/obj.pixelRatio;
                %                 end
            end
            
            if mod(obj.nPixels,2) %odd
                obj.centerX=(obj.nPixels-1)/2+1;
            else % even
                obj.centerX=obj.nPixels/2;
            end
            obj.centerY=obj.centerX;
            
            % CALL METHODS TO BUILD PUPIL
            mkSegList(obj);
            mkPupil(obj);
            
            if obj.flagNoGap
                obj.fillGap(obj.findBestSegmentforFillGap);
            end
            [indX,indY]=find(obj.matrix);
            obj.Dpx=max(max(indX)-min(indX)+1,max(indY)-min(indY)+1);
            if isempty(obj.D)
                obj.D=obj.Dpx/obj.pixelRatio;
            end
            applyReflexivity(obj,1:obj.nSegments,obj.coeffReflexion);
            applyPhaseError(obj,1:obj.nSegments,obj.coeffPhaseModes);
            mkCentObs(obj);
            mkSpiders(obj,obj.spiders);
        end
        
        %% Creating Pupil
        function mkSegList(obj)
            SLT=tic;
            obj.segList=segment.empty(0,obj.nSegments);
            for k=1:size(obj.segmentCoord,1)
                obj.segList(k)=obj.segRef;
                obj.segList(k).posX=obj.segmentCoord(k,1);
                obj.segList(k).posY=obj.segmentCoord(k,2);
                obj.segList(k).matrixBase=[]; % free space if keepHR on segRef was 1
            end
            if obj.verbose
                disp(['List of segments created in : ' num2str(toc(SLT))]);
            end
        end
        
        function mkPupil(obj)
            tic
            obj.matrix = zeros(obj.nPixels);
            % Filling pupil
            for k=1:obj.nSegments
                [posx,posy]           = obj.getPosSeg(k);
                obj.matrix(posx,posy) = obj.matrix(posx,posy) + obj.segList(k).matrix;
            end
            
            if obj.verbose
                disp(['Pupil done in ',num2str(toc),'s']);
            end
        end
        
        function gapMatrix=fillGap(obj,refSegNumber)
            %%
            % hypothesis :
            % all segments are regulary placed
            % non overlaying pixels.
            % all segments have the same geometry.
            [sx,~]      = find(obj.segRef.matrix);
            [posx,posy] = obj.getPosSeg(refSegNumber);
            tmp         = obj.matrix(posx,posy);
            gapMatrix   =~tmp;
            
            
            % Filling oblique sides (up left, up right)
            infX = min(sx);
            supX = round(size(obj.segRef.matrix,1)/2);
            infY = 1;
            supY = size(obj.segRef.matrix,2);
            
            obj.segRef.matrix(infX:supX, infY:supY)=...
                obj.segRef.matrix(infX:supX, infY:supY)+...
                gapMatrix(infX:supX, infY:supY);
            
            % Filling up side
            [~,Y]=find(obj.segRef.matrix(infX,:));
            obj.segRef.matrix(1:infX-1, min(Y):max(Y))=...
                obj.segRef.matrix(1:infX-1, min(Y):max(Y))+...
                gapMatrix(1:infX-1, min(Y):max(Y));
            
            obj.mkSegList;
            obj.mkPupil;
            
            if obj.verbose
                figure('name','Gap Matrix');
                imagesc(gapMatrix);
                axis equal
                
                figure('name','Upside Gap Matrix');
                imagesc(gapMatrix(infX:supX, infY:supY));
                axis equal
                
                figure('name','New Reference Segment')
                imagesc(obj.segRef.matrix);
                axis equal
            end
        end
        
        function segNumber=findBestSegmentforFillGap(obj)
            segNumber=1;
            k=1;
            flagFind=0;
            while ~flagFind && k<=obj.nSegments
                [posx,posy]=obj.getPosSeg(k);
                tmp=obj.matrix(posx,posy);
                if obj.verbose
                    figure(777);
                    imagesc(tmp);
                    pause(1/100);
                end
                
                if (tmp(1,1))&&(tmp(1,end))&&(tmp(1,round(size(tmp,2)/2)))
                    segNumber=k;
                    flagFind=1;
                end
                
                k=k+1;
            end
            
            if obj.verbose
                disp(['Segment n??' num2str(segNumber) 'selected to find gap']);
            end
            
        end
        
        %% Change segments
        
        function rmSegment(obj,indSeg)
            tic
            for k=1:length(indSeg)
                [posx,posy]=obj.getPosSeg(indSeg(k));
                obj.matrix(posx,posy)=obj.matrix(posx,posy).*(~abs(obj.segList(indSeg(k)).matrix));
                %obj.matrix(posx,posy)=obj.matrix(posx,posy)-obj.segList(indSeg(k)).matrix;
            end
            if obj.verbose
                disp(['Segments removed in ',num2str(toc),' s']);
            end
        end
        
        function applyPhaseError(obj,indSeg,coeffPhase)
            tic
            if ~isempty(coeffPhase)
                % i : number of segment to modify
                % coeffZern : i x Nmodes matrix
                
                [Ncoeffs,Nmodes]=size(coeffPhase);
                if Nmodes>obj.nPhaseModes               
                    obj.nPhaseModes=Nmodes;
                end
                
                if isempty(obj.phaseModesGlobal) || obj.nPhaseModes<Nmodes
                    obj.phaseModesGlobal=obj.computeModes(obj.phaseModesType,obj.nPhaseModes);
                end
                
                obj.rmSegment(indSeg);
                if Ncoeffs==length(indSeg)
                    for k=1:Ncoeffs
                        obj.segList(indSeg(k)).coeffPhaseModes=coeffPhase(k,:);
                        [obj.segList(indSeg(k)).matrix, obj.segList(indSeg(k)).phaseError]=obj.segList(indSeg(k)).applyPhase(obj.phaseModesGlobal);
                        [posx,posy]=obj.getPosSeg(indSeg(k));
                        obj.matrix(posx,posy)=obj.matrix(posx,posy)+obj.segList(indSeg(k)).matrix;
                    end
                else
                    error('Number of Phase coefficient doesn''t match number of segments to modify');
                end
                
                mkCentObs(obj);
                mkSpiders(obj,obj.spiders);
                
                t=toc;
                if obj.verbose
                    switch obj.phaseModesType
                        case 'zernike'
                            disp([num2str(Nmodes) ' zernike modes applied in ',num2str(t),'s']);
                        otherwise
                            disp([num2str(Nmodes) ' modes (phase) applied in ',num2str(t),'s']);
                    end
                end
            end
            
        end
        
        function applyReflexivity(obj,indSeg,coeffReflexion)
            if ~isempty(coeffReflexion)
                
                obj.rmSegment(indSeg);
                
                if length(indSeg)==length(coeffReflexion) || (length(indSeg)==size(coeffReflexion,3)) % "cube" of reflexion is possible
                    for k=1:length(indSeg)
                        [obj.segList(indSeg(k)).matrix, obj.segList(indSeg(k)).reflexivity]= obj.segList(indSeg(k)).applyReflexivity(coeffReflexion(k));
                        [posx,posy]=obj.getPosSeg(indSeg(k));
                        obj.matrix(posx,posy)=obj.matrix(posx,posy)+obj.segList(indSeg(k)).matrix;
                    end                    
                    
                    mkCentObs(obj);
                    mkSpiders(obj,obj.spiders);
                    
                    if obj.verbose
                        disp(['Reflexion applied in ',num2str(toc),'s']);
                    end
                else
                    error('Number of Reflexion coefficient doesn''t match number of segments');
                end
            end
        end
        
        function shiftSegment(obj,indSeg,x,y)
            % Move a segment giving its number, X offset and Y offset
            tic
            obj.rmSegment(indSeg);
            for k=1:length(indSeg)
                obj.segList(indSeg(k)).posXError=obj.segList(indSeg(k)).posXError+x(k);
                obj.segList(indSeg(k)).posYError=obj.segList(indSeg(k)).posYError+y(k);
                [posx,posy]=obj.getPosSeg(indSeg(k));
                obj.matrix(posx,posy)=obj.matrix(posx,posy)+obj.segList(indSeg(k)).matrix;
            end
            
            mkCentObs(obj);
            mkSpiders(obj,obj.spiders);
            if obj.verbose
                disp(['Segments moved in ' num2str(toc) ' s']);
            end
        end

        function rotateSegment(obj,indSeg,rValue)
            obj.rmSegment(indSeg);
            for k=1:length(indSeg)
                obj.segList(indSeg(k)).angleError=rValue(k);
                obj.segList(indSeg(k)).matrix=obj.segList(indSeg(k)).mkSegment;
                % old reflexivity ans phase applied while re-making matrix
                [posx,posy]=obj.getPosSeg(indSeg(k));
                obj.matrix(posx,posy)=obj.matrix(posx,posy)+obj.segList(indSeg(k)).matrix;
            end
            
            mkCentObs(obj);
            mkSpiders(obj,obj.spiders);
        end
        
        function shrinkSegment(obj,indSeg,rdValue)
            % reduce size of a segment giving a value between 0 & 1
            obj.rmSegment(indSeg);
            for k=1:length(indSeg)
                obj.segList(indSeg(k)).sizeError=rdValue(k);
                obj.segList(indSeg(k)).matrix=obj.segList(indSeg(k)).mkSegment;
                % old reflexivity ans phase applied while re-making matrix
                [posx,posy]=obj.getPosSeg(indSeg(k));
                obj.matrix(posx,posy)=obj.matrix(posx,posy)+obj.segList(indSeg(k)).matrix;
            end
            
            mkCentObs(obj);
            mkSpiders(obj,obj.spiders);
        end
        
        function [posx,posy]=getPosSeg(obj,indSeg)
            switch obj.unit
                case 'm'
                    [sx,sy]=size(obj.segList(indSeg).matrix);
                    if mod(sx,2) %impair
                        Ax=obj.centerX+(-(sx-1)/2:(sx-1)/2);
                    else %pair
                        Ax=obj.centerX+(-sx/2+1:sx/2);
                    end
                    if mod(sy,2) %impair
                        Ay=obj.centerY+(-(sy-1)/2:(sy-1)/2);
                    else %pair
                        Ay=obj.centerY+(-sy/2+1:sy/2);
                    end
                    posy=Ax+round(obj.pixelRatio*(obj.segList(indSeg).posX+obj.segList(indSeg).posXError));
                    %%
                    % Same convention as ESO doc,
                    % but when y decreases, the segment goes up (a bit weird)
                    posx=Ay+round(obj.pixelRatio*(obj.segList(indSeg).posY+obj.segList(indSeg).posYError));
                    
                    %%
                    % For a logical visualisation : if y decrease, the segment goes down
                    % posx=Ay-round(obj.pixelRatio*(obj.segList(i).posY+obj.segList(i).posYError));
                case 'px'
                    [sx,sy]=size(obj.segList(indSeg).matrix);
                    if mod(sx,2) %impair
                        Ax=(-(sx-1)/2:(sx-1)/2);
                    else %pair
                        Ax=(-sx/2+1:sx/2);
                    end
                    if mod(sy,2) %impair
                        Ay=(-(sy-1)/2:(sy-1)/2);
                    else %pair
                        Ay=(-sy/2+1:sy/2);
                    end
                    posy=obj.segList(indSeg).posX+obj.segList(indSeg).posXError+Ax;
                    posx=obj.segList(indSeg).posY+obj.segList(indSeg).posYError+Ay;
                otherwise
                    posx=1:size(obj.segList(indSeg).matrix,1);
                    posy=1:size(obj.segList(indSeg).matrix,2);
            end
        end
        
        
        %% Change Pupil
        
        function mkCentObs(obj)
            tic
            if obj.obstructionRatio~=0
                [rr ,cc] = meshgrid(linspace(-1,1,obj.nPixels));
                C=(sqrt((rr).^2+(cc).^2)<=obj.obstructionRatio);
                obj.matrix(C)=0;
                if obj.verbose
                    disp(['Central Obstruction done in ',num2str(toc),'s']);
                end
            end
            
        end
        
        function mkSpiders(obj,spiderStruct)
            
            nSpiders         = spiderStruct.n;
            spidersAngle     = spiderStruct.angle;
            spidersWidth     = spiderStruct.width;
            spiderUsrDefined = spiderStruct.usrDefined;
            s                = obj.nPixels;
            tic
            if (nSpiders~=0) && ~spiderUsrDefined
                              
                phi     = atan(spidersWidth/(2*obj.D));
                theta   = [-phi,phi,-phi+pi,+phi+pi];
                [X,Y]   = meshgrid(linspace(-1,1,ceil(s/2)));
                spiders = zeros(ceil(s/2));
                
                for k = 1:nSpiders
                    xv      = 2*cos(theta+spidersAngle(k));
                    yv      = 2*sin(theta+spidersAngle(k));
                    [in,on] = inpolygon(X,Y,xv,yv);
                    spiders = spiders|in|on;
                end
                spidersMatrix = logical(tools.interpolate(double(spiders),s));
                spidersMatrix=~spidersMatrix;
                     
%                     % If IPT is available
%                     midy=round(s/2);
%                     
%                     if strcmp(obj.unit,'m')
%                         halfWidthPx=(spidersWidth/2)*obj.pixelRatio;
%                     elseif strcmp(obj.unit,'px')
%                         halfWidthPx=round(spidersWidth/2);
%                     else
%                         halfWidthPx=0;
%                     end
%                     X               = [s s 1 1];
%                     Y               = [midy-halfWidthPx midy+halfWidthPx midy+halfWidthPx midy-halfWidthPx];
%                     horSpider       = poly2mask(X,Y,s,s);
%                     spidersMatrix   = zeros(s);
%                     
%                     for k=1:nSpiders
%                         spidersMatrix = spidersMatrix|imrotate(horSpider,spidersAngle(k)*180/pi,'bilinear','crop');
%                     end
%                     spidersMatrix=~spidersMatrix;
                                
            elseif spiderUsrDefined
                spidersMatrix=spiderStruct.usrMask;
            else
                spidersMatrix=ones(size(obj.matrix));
            end
            obj.matrix=obj.matrix.*spidersMatrix;
            t=toc;
            if obj.verbose
                disp(['Spiders done in ',num2str(t),'s']);
            end            
        end
        
        function out=zeroPad(obj,ratio,flagReplace)
            if nargin<3
                flagReplace=1;
            end
            tic
            out=tools.enlargePupil(obj.matrix,ratio);
            if flagReplace
                obj.matrix=out;
                obj.nPixels=size(obj.matrix,1);
                oldCenterX=obj.centerX;
                oldCenterY=obj.centerY;
                if mod(obj.nPixels,2) %odd
                    obj.centerX=(obj.nPixels-1)/2+1;
                else % even
                    obj.centerX=obj.nPixels/2;
                end
                obj.centerY=obj.centerX;
                
                if strcmp('px',obj.unit)
                    for k=1:obj.nSegments
                        obj.segList(k).posX=obj.segList(k).posX-(oldCenterX-obj.centerX);
                        obj.segList(k).posY=obj.segList(k).posY-(oldCenterY-obj.centerY);
                    end
                end
            end
            if obj.verbose
                disp(['Pupil zero padded in ',num2str(toc),' s']);
            end
        end
        
        function out=rmZeroBorder(obj,flagReplace)
            disp('removing border full of zeros');
            if nargin<2
                flagReplace=1;
            end
            [i,j]=find(obj.matrix);
            Z1=min(min(i),min(j));
            Z2=max(max(i),max(j));
            out=obj.matrix(Z1:Z2,Z1:Z2);
            
            if flagReplace==1
                obj.matrix=out;
                obj.nPixels=size(obj.matrix,1);
                
                oldCenterX=obj.centerX;
                oldCenterY=obj.centerY;
                if mod(obj.nPixels,2) %odd
                    obj.centerX=(obj.nPixels-1)/2+1;
                else % even
                    obj.centerX=obj.nPixels/2;
                end
                obj.centerY=obj.centerX;
                
                if strcmp('px',obj.unit)
                    for k=1:obj.nSegments
                        obj.segList(k).posX=obj.segList(k).posX-(oldCenterX-obj.centerX);
                        obj.segList(k).posY=obj.segList(k).posY-(oldCenterY-obj.centerY);
                    end
                end                
            end
        end
        
        function out=resize(obj,newNPixels,flagReplace,method)
            if nargin<3
                flagReplace=0;
            end
            if nargin<4
                method='nearest';
            end
            
            xi       = linspace(-1,1,length(obj.matrix));
            [Xi,Yi]  = meshgrid(xi);
            xo       = linspace(-1,1,newNPixels);
            [Xo,Yo]  = meshgrid(xo);
            out      = interp2(Xi,Yi,obj.matrix,Xo,Yo,method);
            
            if flagReplace
                obj.matrix=out;
                %obj.segRef.matrix=tools.interpolate(obj.segRef.matrix,round(obj.segRef.nPixels*newNPixels/obj.nPixels));
                %obj.segRef=segment(obj.segRef.nSides,obj.segRef.radius,round(obj.segRef.nPixels*newNPixels/obj.nPixels));
                %obj.segRef.nPixels=length(obj.segRef.matrix);
                obj.pixelRatio=obj.pixelRatio*newNPixels/obj.nPixels;
                obj.nPixels=size(obj.matrix,1);
                if mod(obj.nPixels,2) %odd
                    obj.centerX=(obj.nPixels-1)/2+1;
                else % even
                    obj.centerX=obj.nPixels/2;
                end
                obj.centerY=obj.centerX;
                
            end
        end
        
        function out=computeModes(obj,type,nModes)
            [sx,sy]=size(obj.segRef.matrix);
            switch type
                case 'zernike'
                    zern = zernike(1:nModes,'resolution',max(sx,sy),'logging',false);
                    out=reshape(zern.modes,sx,sy,nModes);
                otherwise
                    out=zeros(sx,sy);
            end
        end
        
        %% Display and Disp overload
        
%         function disp(obj,what)
%             if nargin<2
%                 what='both';
%             end
%             objName=inputname(1);
%             
%             figWidth = 1250; % pixels
%             figHeight = 1025;
%             rect = [700 0 figWidth figHeight];            
%             
%             switch what
%                 case 'r'
%                     figure('name',[objName '  reflexiviy'],'NumberTitle','off','OuterPosition',rect);
%                     colormap(parula(1024));
%                     imagesc(abs(obj.matrix));
%                     title('Reflexivity');
%                     axis equal
%                     colorbar
%                     
%                 case 'p'
%                     figure('name',[objName '  phase map'],'NumberTitle','off','OuterPosition',rect);
%                     colormap(jet(1024));
%                     tmp=obj.matrix;
%                     tmp((tmp==0))=NaN;  % better contrast
%                     tmp=angle(tmp);
%                     imagesc(tmp);
%                     title('Phase');
%                     axis equal
%                     colorbar
%                 case 'both'
%                     figure('name',[objName ' : reflexiviy'],'NumberTitle','off','OuterPosition',rect);
%                     colormap(parula(1024));
%                     imagesc(abs(obj.matrix));
%                     title('Reflexivity');
%                     axis equal
%                     colorbar
%                     figure('name',[objName ' : phase map'],'NumberTitle','off','OuterPosition',rect);
%                     colormap(jet(1024));
%                     tmp=obj.matrix;
%                     tmp((tmp==0))=NaN; % better contrast
%                     tmp=angle(tmp);
%                     imagesc(tmp);
%                     title('Phase');
%                     axis equal
%                     colorbar;
%                 otherwise
%                     disp('not valid argument');
%             end
%         end
        
        %         function display(obj)
        %             obj.segRef
        %             disp('___PUPIL___');
        %             disp(['Made of ', num2str(obj.nSegments),' segments']);
        %             disp(['Size of the matrix : ', num2str(obj.nPixels(1)),'x',num2str(obj.nPixels(2))]);
        %             disp(['One meter is : ',num2str(obj.pixelRatio),' pixels']);
        %             disp(['Outer Radius : ', num2str(obj.radius), ' m']);
        %             disp(['Inner Radius : ', num2str(obj.inRadius), ' m']);
        %             disp(['All Glass Radius : ', num2str(obj.allGlassRadius), ' m']);
        %         end
        
        
        
    end
    
    methods (Static)
        
        %% Demo
        function demoKeckMeters
            %close all;
            %Parameters
            %lambda=500e-9; %m
            segmentResolution=200; %px
            segmentRadius=0.9;  % m
            segmentNSides=6;
            
            %segment creation
            S=segment(segmentNSides,segmentRadius,segmentResolution);
            
            % definition of the segments' centers coordinates, in meters
            SV(:,1)=[1.35000000000000;0;-1.35000000000000;-1.35000000000000;0;1.35000000000000;2.70000000000000;1.35000000000000;0;-1.35000000000000;-2.70000000000000;-2.70000000000000;-2.70000000000000;-1.35000000000000;0;1.35000000000000;2.70000000000000;2.70000000000000;4.05000000000000;2.70000000000000;1.35000000000000;0;-1.35000000000000;-2.70000000000000;-4.05000000000000;-4.05000000000000;-4.05000000000000;-4.05000000000000;-2.70000000000000;-1.35000000000000;0;1.35000000000000;2.70000000000000;4.05000000000000;4.05000000000000;4.05000000000000];
            SV(:,2)=[-0.783000000000000;-1.56600000000000;-0.783000000000000;0.783000000000000;1.56600000000000;0.783000000000000;-1.56600000000000;-2.34900000000000;-3.13200000000000;-2.34900000000000;-1.56600000000000;0;1.56600000000000;2.34900000000000;3.13200000000000;2.34900000000000;1.56600000000000;0;-2.34900000000000;-3.13200000000000;-3.91500000000000;-4.69800000000000;-3.91500000000000;-3.13200000000000;-2.34900000000000;-0.783000000000000;0.783000000000000;2.34900000000000;3.13200000000000;3.91500000000000;4.69800000000000;3.91500000000000;3.13200000000000;2.34900000000000;0.783000000000000;-0.783000000000000];
            
            % Random errors
            R=0.9+rand(1,length(SV))*(1-0.9); % random reflexivity btw .95 and 1
            piston=rand(1,length(SV))*10e-9; % random piston btw 0 and 10 nm
            tip=rand(1,length(SV))*5e-9; % random tip btw 0 and 5 nm
            tilt=rand(1,length(SV))*5e-9; % random tilt btw 0 and 5 nm
            Z=[piston' tip' tilt'];
            
            % Spider Structure
            SPIDER.n=3;
            SPIDER.angle=[0 pi/3 2*pi/3];
            SPIDER.width=0.0254; % in m
            SPIDER.usrDefined=0;
            SPIDER.usrMask=[];
            
            % pupil creation
            % if coordinates are good enough (see documentation for more details)
            % it is possible to fill the gaps between the segments
            % with the flagNoGap option
            P=pupilGenerator('segRef',S,'segCoord',SV,'flagNoGap',1,...
                'coeffDistance',1,'coeffReflexion',R,'coeffPhaseModes',Z);
            
            % P.disp % if no argument : display reflexivity AND phase
            %P.disp('r'); % display reflexivity only
            %P.disp('p'); % display phase only
            
            % change pupil after creation : remove, rotate, move, change
            % reflexivity, change phase
            P.rmSegment([7 17 21]); % remove segment n??7, n??17 and n??21
            P.applyPhaseError(18,[0 0 0 1e-10 0 0 0 0 10e-9]); % add coma and trefoil
            P.shiftSegment(30,-3,0.55); % move segment n??36 of -3m on x and 55cm on y
            P.rotateSegment(30,pi/6);
            P.applyReflexivity([1 30],[ 0.5 0.75]);
            P.shrinkSegment(2,0.5);
            P.disp;
            
            % Making spiders and central obstruction is possible
            P.mkSpiders(SPIDER);
            P.obstructionRatio=0.2375;
            P.mkCentObs;
            %P.disp('r');
            % also available at creation via :
            % P=pupilGenerator(S,SV,'spiders',SPIDER,'obstructionRatio',0.2375);
            
            % once all modif done, one can resize the pupil giving new size in px
            % at 04/05/2017 : modif methods above don't work  well (or at all...) if used after
            % the resize...
            P.resize(300,1); % 300 pixels. the "1" means that the matrix will be replaced in the pupil
            P.zeroPad(2,1);
            %P.disp('p');
        end
        
        function demoKeckPixels
            %close all;
            %Parameters
            %lambda=500e-9; %m
            segmentResolution=200; %px
            segmentRadius=0.9;  % m
            segmentNSides=6;
            
            %segment creation
            S=segment(segmentNSides,segmentRadius,segmentResolution);
            
            % definition of the segments' centers coordinates, in pixels
            SV(:,1)=[772;622;472;472;622;772;922;772;622;472;322;322;322;472;622;772;922;922;1072;922;772;622;472;322;172;172;172;172;322;472;622;772;922;1072;1072;1072];
            SV(:,2)=[535;448;535;709;796;709;448;361;274;361;448;622;796;883;970;883;796;622;361;274;187;100;187;274;361;535;709;883;970;1057;1144;1057;970;883;709;535];
            
            % Random errors
            R=0.9+rand(1,length(SV))*(1-0.9); % random reflexivity btw .95 and 1
            piston=rand(1,length(SV))*10e-9; % random piston btw 0 and 10 nm
            tip=rand(1,length(SV))*5e-9; % random tip btw 0 and 5 nm
            tilt=rand(1,length(SV))*5e-9; % random tilt btw 0 and 5 nm
            Z=[piston' tip' tilt'];
            
            % Spider Structure
            SPIDER.n=3;
            SPIDER.angle=[0 pi/3 2*pi/3];
            SPIDER.width=3; % in px
            SPIDER.usrDefined=0;
            SPIDER.usrMask=[];
            % pupil creation
            % if coordinates are good enough (see documentation for more details)
            % it is possible to fill the gaps between the segments
            % with the flagNoGap option
            P=pupilGenerator('segRef',S,'segCoord',SV,'flagNoGap',1,...
                'unit','px','coeffReflexion',R,'coeffPhaseModes',Z);
            
            % P.disp % if no argument : display reflexivity AND phase
            P.disp('r'); % display reflexivity only
            P.disp('p'); % display phase only
            
            % change pupil after creation : remove, rotate, move, change
            % reflexivity, change phase
            P.rmSegment([7 17 21]); % remove segment n??7, n??17 and n??21
            P.applyPhaseError(18,[0 0 0 1e-10 0 0 0 0 10e-9]); % add coma and trefoil
            P.shiftSegment(30,-315,25); % move segment n??36 of -315 px on x and -25px on y
            P.rotateSegment(30,pi/6);
            P.applyReflexivity([1 30],[0.5 0.75]);
            P.shrinkSegment(2,0.5);
            P.disp;
            
            % Making spiders and central obstruction is possible
            P.mkSpiders(SPIDER);
            P.obstructionRatio=0.2375;
            P.mkCentObs;
            P.disp('r');
            % also available at creation via :
            % P=pupilGenerator(S,SV,'unit','px','spiders',SPIDER,'obstructionRatio',0.2375);
            
            % once all modif done, one can resize the pupil giving new size in px
            % at 04/05/2017 : modif methods above don't work  well (or at all...) if used after
            % the resize...
            P.resize(300,1); % 300 pixels. the "1" means that the matrix will be replaced in the pupil
            P.zeroPad(2,1);
            P.disp('p');
        end
        
        function out = demoELTMeters
            %close all;
            %Parameters
            %lambda=500e-9; %m
            segmentResolution=50; %px
            segmentRadius=1.3/2;  % m
            segmentNSides=6;
            
            %segment creation
            S=segment(segmentNSides,segmentRadius,segmentResolution);
            
            % definition of the segments' centers coordinates, in meters
            vertices = load('/data/HARMONI/SCAO/ESO_ELT_R1.3/TEL-INS-IF/Segmentation/SegmentVertices.txt');
            SV(:,1) = vertices(:,1);
            SV(:,2) = vertices(:,2);
            
            % Spider Structure
            SPIDER.n=6;
            SPIDER.angle=[0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3];
            SPIDER.width=0.5; % in m
            SPIDER.usrDefined=0;
            SPIDER.usrMask=[];
            
            % pupil creation
            % if coordinates are good enough (see documentation for more details)
            % it is possible to fill the gaps between the segments
            % with the flagNoGap option
            P=pupilGenerator('segRef',S,'segCoord',SV,'flagNoGap',0,...
                'coeffDistance',1);
            
            % P.disp % if no argument : display reflexivity AND phase
            P.disp('r'); % display reflexivity only
            
            % change pupil after creation : remove, rotate, move, change
            % reflexivity, change phase
            
            % Making spiders and central obstruction is possible
%             P.mkSpiders(SPIDER);
%             P.obstructionRatio=0.3;
%             P.mkCentObs;
%             P.disp('r');
            % also available at creation via :
%             P=pupilGenerator(S,SV,'spiders',SPIDER,'obstructionRatio',0.2375);
%             P.disp('r');
            out = P;
        end
        
    end
end