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
        slopes=0;
        % gradient matrix
        gradOp;
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
        
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(obj, src)
            nSrc = length([src.nSrc]);
            obj.slopes = [];
            nRes = size(src(1).phase,1);
            if length(size(src(1).phase)) == 3
                nRes  = size(src(1).phase,1);
                nActu = size(src(1).phase,3);
            end
            for iSrc = 1:nSrc
                if length(size(src(1).phase)) == 3
                    ph = reshape(src(iSrc).phase,nRes*nRes,nActu);                    
                    obj.slopes(:,:,iSrc) = obj.gradOp*ph(obj.pupil(:),:);
                else
                    obj.slopes(:,iSrc) = obj.gradOp*src(iSrc).phase(obj.pupil);
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
