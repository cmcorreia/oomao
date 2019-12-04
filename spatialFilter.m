classdef spatialFilter < handle

    % spatialFilter
    % Creates an object spatialFilter, this object applies to 
    % the electric field (hence to the source object) a filtering 
    % in the focal plane, whose shape is square, of clear opening size N
    %
    % sf = spatialFilter(5,nPx); creates a spatial filter with 5 times the
    % full telescope diffraction on nPx support
    %
    
    
    properties
        % size of the spatial filter opening, in unit of full telescope pupil diffraction 
        N;
        % size of the pupil plane field support, in pixels
        resolution;
        % mask of the spatial filter
        spatialFilterMask;
        % display info
        tag = 'SPATIAL FILTER';
    end
    
    
    properties (Dependent)%, SetAccess=private)
    end
    
    properties (Access=private)
    end
    
    methods
        
        %% Constructor
        function obj = spatialFilter(N,resolution)
            p = inputParser;
            p.addRequired('N', @isnumeric);
            p.addRequired('resolution', @isnumeric);
            p.parse(N,resolution)
            
            obj.N = p.Results.N;
            obj.resolution = p.Results.resolution;

	        [spatialFilterMask] = calcSpatialFilterMask(obj);
            obj.spatialFilterMask = spatialFilterMask;
        end
        
        %% Destructor
        function delete(obj)
        end

        function display(obj)
        %% DISPLAY Display text object information

            fprintf('___ %s ___\n',obj.tag)
            fprintf('Size of spatial filter : %d Diffraction size \n',obj.N)
            fprintf('Size of field support  : %d x %d pixels \n',obj.resolution,obj.resolution)
            fprintf('----------------------------------------------------\n')

        end
        
        
        %% IMAGESC displays the mask function 
        function imagesc(obj)
            imagesc(obj.spatialFilterMask);
            title(obj.tag);
            axis image;
            colorbar;
        end
        
        
        %% Relay spatialFilter to source relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(obj, src)% wfs=pyramid(nLenslet,nPx,'modulation',modul);
            applySpatialFilter(obj,src);
        
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
