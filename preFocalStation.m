classdef preFocalStation < handle
    
    % preFocalStation
    % Creates an object pfs that offsets and rotates the elctric field from
    % the source
    %
    % pfs = preFocalStation(src,offset); creates a prefocalStation
    %
    
    
    properties
        % offsets in meters
        offset;
        % rotation in radians
        rotation
        % telescope diameter for normalisation
        D
        % number of sources
        nSrc
        % display info
        tag = 'PREFOCAL STATION';
    end
    
    
    properties (Dependent)%, SetAccess=private)

    end
    
    properties (Access=private)
    end
    
    methods
        
        %% Constructor
        function obj = preFocalStation(tel,src,offset,varargin)
            p = inputParser;
            p.addRequired('tel', @(x) isa(x,'telescopeAbstract') );
            p.addRequired('src', @(x) isa(x,'source') );
            p.addRequired('offset', @isnumeric);
            p.addOptional('rotation', [], @isnumeric);
            p.parse(tel,src,offset, varargin{:})
            
            obj.D   = p.Results.tel.D;
            obj.nSrc = src(1).nSrc;
            obj.offset = p.Results.offset;
            obj.rotation = p.Results.rotation;
            
        end
        
        %% Destructor
        function delete(obj)
        end
        
        %% get and set d
        function set.D(obj,val)
            obj.D = val;
        end
        function out = get.D(obj)
            out = obj.D;
        end
        
        %% Get and Set rotation
        function set.rotation(obj,val)
            if isempty(val)
                obj.rotation = 0;
            else
                obj.rotation = val;
            end
            if length(obj.rotation) == 1 && length(obj.rotation) < obj.nSrc
                obj.rotation = repmat(obj.rotation,1,obj.nSrc);
            end
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
                obj.offset = [1;1]*val;
            elseif size(val,2) >  size(val,1)
                if size(val,2) == 2 % if a 1x2 vector, transpose
                    obj.offset =val';
                elseif size(val,1) == 1 % if a 1 x N vector, replicate along the other direction
                    obj.offset = repmat(val, 2,1);
                else
                    obj.offset = val;
                end
            end
            if size(obj.offset,2) == 1 && size(obj.offset,2) < obj.nSrc
                obj.offset = repmat(obj.offset,1,obj.nSrc);
            end
        end
        
        %%
        function display(obj)
            %% DISPLAY Display text object information
            
            fprintf('___ %s ___\n',obj.tag)
            fprintf('Number of sources: %d  \n', obj.nSrc)
            fprintf('----------------------------------------------------\n')
            
        end
        
       
        
        
        %% Relay preFocalStation to source relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(obj, src)% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            
            % apply any offsets
            if ~isempty(obj.offset) && any(obj.offset(:) ~= 0)
                [nx, ny, nz] = size(src(1).phaseUnmasked);
                R = obj.D/2;
                [x0,y0] = meshgrid(linspace(-R,R,nx), linspace(-R,R,ny));
                nSrc_ = numel(src);
                phasei = zeros(nx, ny, nz);
                for iSrc = 1:nSrc_
                    xi = x0 - obj.offset(1,iSrc);
                    yi = y0 - obj.offset(2,iSrc);
                    ri = obj.rotation(iSrc);
                    for kWave = 1:nz
                        if ri ~= 0
                            phasei(:,:,kWave) = rotate(src(iSrc).phaseUnmasked(:,:,kWave),ri);
                            phasei(:,:,kWave) = interp2(x0,y0,phasei(:,:,kWave),xi,yi,'cubic',0);
                        else
                            phasei(:,:,kWave) = interp2(x0,y0,src(iSrc).phaseUnmasked(:,:,kWave),xi,yi,'cubic',0);
                        end
                    end
                    src(iSrc).resetPhase(phasei);
                end
            end
        end
  
    end
    
    % ---------------------------------------------------------------------
    methods (Static)
        
        
        
        
        
    end
    
    % ---------------------------------------------------------------------
    methods (Access=private)
        %% CALCSPATIALFILTERMASK computes the mask transmission array
        function [spatialFilterMask] = calcSpatialFilterMask(obj)
            
            % Create empty mask
            spatialFilterMask = zeros(obj.resolution, obj.resolution);
            % Fill mask with ones where it is transmissive
            ech = 2;
            spatialFilterMask(obj.resolution/2 - obj.N * ech/2 +1 : obj.resolution/2 + obj.N * ech/2, obj.resolution/2 - obj.N * ech/2 +1: obj.resolution/2 + obj.N * ech/2) = 1.;
            
        end
        
        %% APPLYSPATIALFILTER propagates the electric field through the mask
        function applySpatialFilter(obj, src)
            
            phase0 = src.phase;
            
            % Computes original electric field in pupil plane
            srcPupilPlane  = padarray(src.amplitude .* exp(1i*mod(src.phase, 2*pi)), [obj.resolution/2,obj.resolution/2], 'both');
            %             srcPupilPlane  = padarray(src.wave, [obj.resolution/2,obj.resolution/2], 'both');
            
            % Computes electric field in focal plane
            srcFocalPlane  = fftshift(fft2(fftshift(srcPupilPlane)));
            
            % Computes  electric field after filtering
            srcFocalPlaneF = srcFocalPlane .* padarray(obj.spatialFilterMask, [obj.resolution/2,obj.resolution/2],0,'both');
            
            % Computes  electric field in pupil plane after filtering
            pupil = logical(abs(src.wave));
            srcPupilPlaneF = fftshift(ifft2(fftshift(srcFocalPlaneF)));
            srcPupilPlaneF = pupil .* srcPupilPlaneF(obj.resolution/2+1:3*obj.resolution/2, obj.resolution/2+1:3*obj.resolution/2);
            
            % Computes final electric field in pupil plane
            src.amplitude  = abs(srcPupilPlaneF);
            %src.phase      = -1*src.phase;
            %src.phase      =  my_phase_unwrap_2D(angle(srcPupilPlaneF), src.amplitude, 400,400);
            src.phase      =  angle(srcPupilPlaneF) - phase0; % Due to " src.phase = " operation, that adds the new phase to previous one
        end
        
        
    end
    
end
