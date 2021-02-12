classdef geomGrad < handle
    properties (SetAccess = public)
        validLenslet;
        resolution;
        pupil;
        edgePix;
        nSlope;
        lenslets; % compatibility with the diffractive shackHartmann
        % wavefront sensor tag
        % sets the display ON or OFF (useful for cluster execution)
        graphicalDisplay = 0;
        
        tag = 'G-SHACK-HARTMANN';
    end
    properties (SetObservable=true)
        % measurements
        slopes=0; %units of Shannon-sampled pixels, i.e. \lambda/2/d (\lambda - WFS wavelength, d-size of the sub-aperture in m)
        % gradient matrix
        gradOp;
    end
    properties (Dependent)
        % offsets in fraction of the telescope diameter
        offset;
        % rotation in radians
        rotation
    end
    
    properties (Access=private)
        % vector of offsets
        p_offset = [0;0;0];
    end
    
    methods
        %% Constructor
        function obj = geomGrad(validLenslet, pupil, varargin)
            %error(nargchk(1,3,nargin));
            p = inputParser;
            p.addRequired('validLenslet', @(x) isnumeric(x) || islogical(x));
            p.addRequired('pupil', @(x) isnumeric(x) || islogical(x));
            p.addOptional('edgePix', 0, @isnumeric);
            p.parse(validLenslet, pupil,varargin{:});
            
            obj.validLenslet    = validLenslet;
            obj.resolution      = size(pupil,1);
            obj.pupil           = pupil;
            obj.edgePix         = p.Results.edgePix;
            obj.nSlope          = sum(validLenslet(:))*2;
            obj.lenslets.nLenslet= sum(validLenslet(:));
            computeGradOp(obj);
        end
        
        %% Get and Set rotation
        function set.rotation(obj,val)
            if isempty(val)
                obj.p_offset(3,:) = 0;
            else
                obj.p_offset(3,length(val)) = 0;
                obj.p_offset(3,:) = val;
            end
            %if length(obj.rotation) == 1 && length(obj.rotation) < obj.nSrc
            %    obj.rotation = repmat(obj.rotation,1,obj.nSrc);
            %end
        end
        
        function out = get.rotation(obj)
            out = obj.p_offset(3,:);
        end
        %% Get and Set Offsets
        function set.offset(obj,val)
            % check offset's magnitude
            %             if any(abs(val)) > 1
            %                 warning('Offsets larger than a sub-aperture not supported. Truncating the offsets to 1')
            %                 val(val > 1) = 1;
            %                 val(val < -1) = -1;
            %             end
            % check offset's size consistency
            if any(size(val)) == 1
                if size(val,1) >  size(val,2)
                    val = val';
                end
            end
            if isscalar(val) % create a 2x1 vector
                offset_ = [1;1]*val;
            elseif size(val,2) >  size(val,1)
                if size(val,2) == 2 % if a 1x2 vector, transpose
                    offset_ =val';
                elseif size(val,1) == 1 % if a 1 x N vector, replicate along the other direction
                    offset_ = repmat(val, 2,1);
                else
                    offset_ = val;
                end
            end
            obj.p_offset = offset_;
            obj.p_offset(3,size(obj.p_offset,2)) = 0;
            %if size(obj.offset,2) == 1 && size(obj.offset,2) < obj.nSrc
            %    obj.offset = repmat(obj.offset,1,obj.nSrc);
            %end
        end
        function out = get.offset(obj)
            out = obj.p_offset(1:2,:);
        end
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(obj, src)
            nSrc = length([src.nSrc]);
            obj.slopes = [];
            ndimsPhase = ndims(src(1).phase);
            %nRes = size(src(1).phase,1);
            [nx, ny, nz] = size(src(1).phase);
            %             if ndimsPhase == 3
            %                 nRes  = size(src(1).phase,1);
            %                 nAct = size(src(1).phase,3);
            %             end
            
            for iSrc = 1:nSrc
                % apply any offsets
                if (~isempty(obj.offset) && any(obj.offset(:) ~= 0) || ~isempty(obj.rotation) && any(obj.rotation(:) ~= 0) ) && iSrc <= size(obj.offset,2)
                    R = 0.5;
                    [x0,y0] = meshgrid(linspace(-R,R,nx), linspace(-R,R,ny));
                    phasei = zeros(nx, nx, nz);
                    xi = x0 - obj.offset(1,iSrc);
                    yi = y0 - obj.offset(2,iSrc);
                    ri = obj.rotation(iSrc);
                    for kWave = 1:nz
                        if ri ~= 0
                            phasei(:,:,kWave) = rotate(src(iSrc).phaseUnmasked(:,:,kWave),ri*180/pi);
                            %phasei(:,:,kWave) = tools.rotateIm(src(iSrc).phaseUnmasked(:,:,kWave),ri*180/pi);
                            phasei(:,:,kWave) = interp2(x0,y0,phasei(:,:,kWave),xi,yi,'cubic',0);
                        else
                            phasei(:,:,kWave) = interp2(x0,y0,src(iSrc).phaseUnmasked(:,:,kWave),xi,yi,'cubic',0);
                        end
                    end
                    
                    ph = reshape(phasei,nx*nx,nz);
                    obj.slopes(:,:,iSrc) = obj.gradOp*ph(obj.pupil(:),:);
                    
                else
                    if ndimsPhase == 3
                        ph = reshape(src(iSrc).phase,nx*nx,nz);
                        obj.slopes(:,:,iSrc) = obj.gradOp*ph(obj.pupil(:),:);
                    else
                        obj.slopes(:,iSrc) = obj.gradOp*src(iSrc).phase(obj.pupil);
                    end
                end
            end
        end
        %% Compute gradients
        function obj = computeGradOp(obj)
            
            nSlopes   = 2*nnz(obj.validLenslet);
            G         = zeros(nSlopes,nnz(obj.pupil)); iL = 1;
            nLenslet  = size(obj.validLenslet,1);
            nPxSubap = obj.resolution/nLenslet;
            nPxSubap2= nPxSubap^2;
            for j = 1:nLenslet
                for i = 1:nLenslet
                    % for a valid lenslet
                    if obj.validLenslet(i,j) == 1
                        
                        % valid sub-aperture pixels
                        subapMask = zeros(obj.resolution);
                        subapMask((i-1)*nPxSubap+1:i*nPxSubap,(j-1)*nPxSubap+1:j*nPxSubap) = 1;
                        subapMask = subapMask.*obj.pupil;
                        
                        % X gradients
                        
                        maskLet = zeros(obj.resolution); nPixUsed = 0;
                        if nnz(subapMask) == nPxSubap2 % fully illuminated sub-aperture
                            maskLet((i-1)*nPxSubap+1:i*nPxSubap,(j-1)*nPxSubap+1) = -1;
                            maskLet((i-1)*nPxSubap+1:i*nPxSubap,j*nPxSubap) = 1;
                            if obj.edgePix
                                nPixUsed = nPxSubap2;
                            else
                                nPixUsed = nPxSubap2-nPxSubap;
                            end
                        else % partially illuminated sub-aperture
                            for ik = (i-1)*nPxSubap+1:i*nPxSubap
                                lig = subapMask(ik,:);
                                if nnz(lig) >= 2
                                    jl = find(lig,1,'first');
                                    jr = find(lig,1,'last');
                                    nPixUsed = nPixUsed + (jr-jl+obj.edgePix);
                                    maskLet(ik,jl) = -1;% left pixel
                                    maskLet(ik,jr) = 1;% right pixel
                                end
                            end
                            
                        end
                        maskLet = maskLet*nPxSubap/nPixUsed/pi;
                        % rasterize
                        M = maskLet(:);
                        % keep only valid pixels
                        Mpup = M(obj.pupil,1);
                        % populate gradient matrix
                        if obj.edgePix
                            G(iL,:) = Mpup'*(nPxSubap/(nPxSubap-1));
                        else
                            G(iL,:) = Mpup';
                        end
                        
                        % Y gradients
                        maskLet = zeros(obj.resolution); nPixUsed = 0;
                        if nnz(subapMask) == nPxSubap2 % fully illuminated sub-aperture
                            maskLet((i-1)*nPxSubap+1,(j-1)*nPxSubap+1:j*nPxSubap) = -1;
                            maskLet(i*nPxSubap, (j-1)*nPxSubap+1:j*nPxSubap) = 1;
                            if obj.edgePix
                                nPixUsed = nPxSubap2;
                            else
                                nPixUsed = nPxSubap2-nPxSubap;
                            end
                        else
                            for jk = (j-1)*nPxSubap+1:j*nPxSubap
                                col = subapMask(:,jk);
                                if nnz(col) >= 2
                                    ih = find(col,1,'first');
                                    il = find(col,1,'last');
                                    nPixUsed = nPixUsed + (il-ih+obj.edgePix);
                                    maskLet(ih,jk) = -1;% pixel du HAUT
                                    maskLet(il,jk) = 1;% pixel du BAS
                                end
                            end
                        end
                        maskLet = maskLet*nPxSubap/nPixUsed/pi;
                        % rasterize
                        M = maskLet(:);
                        % keep only valid pixels
                        Mpup = M(obj.pupil,1);
                        % populate gradient matrix
                        if obj.edgePix
                            G(iL+nSlopes/2,:) = Mpup'*(nPxSubap/(nPxSubap-1));
                        else
                            G(iL+nSlopes/2,:) = Mpup';
                        end
                        iL = iL+1;
                        
                    end
                end
            end
            [i,j,s] = find(G);
            obj.gradOp = sparse(i,j,s,nSlopes,nnz(obj.pupil));
        end
    end
    
end
