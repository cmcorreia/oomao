classdef fourierModes < handle
    %% fourierModes
    % Creates object containg fourierModes
    
    properties
        % maximum resolution of mode (i.e. related to the maximum spatial
        % frequency)
        N;
        resolution;
        % # of modes
        nMode;
        %  tag
        % Modes generated
        modes;
        % x spatial frequency
        l;
        % y spatial frequency
        k;
        % cos or sin mode
        CS;
        tag = 'FOURIER MODES';
        
    end
    
    
    properties (Dependent)%, SetAccess=private)
        
        
        
    end
    
    properties (Access=private)
        %p_p;
    end
    
    methods
        
        %% Constructor
        function obj = fourierModes(N,resolution)
            p = inputParser;
            p.addRequired('N', @isnumeric);
            p.addRequired('resolution', @isnumeric);
            p.parse(N,resolution)
            
            obj.N = p.Results.N;
            obj.resolution = p.Results.resolution;
            obj.nMode = obj.N^2;
            
            [modes,l,k,CS] = calcModes(obj);
            obj.modes = reshape(modes,obj.resolution^2,obj.nMode);
            obj.l = l;
            obj.k = k;
            obj.CS = CS;
            
            
        end
        
        %% Destructor
        function delete(obj)
            
        end
        
        
        %% Get modes
        %         function modes = get.modes(obj)
        %             [modes,l,k,CS] = calcModes(obj);
        %             obj.modes = reshape(modes,obj.resolution^2,obj.nMode);
        %             obj.l = l;
        %             obj.k = k;
        %             obj.CS = CS;
        %         end
        %         function set.modes(obj,val)
        %             obj.modes = val;
        %         end
        function out = mtimes(obj,c)
            %% MTIMES Zernike multiplication
            %
            % v = obj*c multiply the Zernike modes in obj by the
            % vector c
            
            out = obj.modes*c;
        end
        
        
        %         function display(obj)
        %             %% DISPLAY Display object information
        %             %
        %             % display(obj) prints information about the atmosphere+telescope object
        %
        %             fprintf('___ %s ___\n',obj.tag)
        %             if obj.nMode>10
        %                 fprintf(' . %d modes: [ %d , %d ]\n',obj.nMode,obj.n(1),obj.n(end))
        %             elseif obj.nMode>1
        %                 fprintf(' . %d modes: ',obj.nMode)
        %                 fprintf('%d,',obj.n)
        %                 fprintf('\b\n')
        %             else
        %                 fprintf(' . mode: %d \n',obj.n(1))
        %             end
        %             fprintf('----------------------------------------------------\n')
        %
        %         end
        
        
        
        function imagesc(obj,N)
            imagesc(reshape(obj.modes(:,N),sqrt(size(obj.modes,1)), sqrt(size(obj.modes,1))));
        end
        
        
    end
    
    % ---------------------------------------------------------------------
    methods (Static)
        
        
        
        
        
    end
    
    % ---------------------------------------------------------------------
    methods (Access=private)
        
        function [modes,l,k,CS] = calcModes(obj)
            
            nPx = obj.resolution;
            N = obj.N;
            
            % Strore modes
            modes = zeros(N,N,N^2);
            
            % x/y values
            [n,m] = meshgrid(0:nPx-1,0:nPx-1);
            modes = zeros(nPx,nPx,N*N);
            
            % Spatial frequencies
            ll = -floor(obj.N/2):floor((obj.N-1)/2);
            kk = -floor(obj.N/2):floor((obj.N-1)/2);
            
            n0 = ceil((N+1)/2);
            modeN = 0;
            l = zeros(1,N^2);
            k = zeros(1,N^2);
            CS = zeros(1,N^2);
            
            for i=1:N
                for j=1:N
                    
                    modeN = modeN+1;
                    if ceil(N/2)==N/2 && i==1 && j>n0
                        fmode = sin(2*pi/nPx*(-ll(i)*n-kk(j)*m));
                        CS(modeN) = -1;
                    elseif i<n0 || (i==n0 && j<=n0)
                        fmode = cos(2*pi/nPx*(-ll(i)*n-kk(j)*m));
                        CS(modeN) = 1;
                    else
                        fmode = sin(2*pi/nPx*(-ll(i)*n-kk(j)*m));
                        CS(modeN) = -1;
                    end
                    
                    modes(:,:,modeN) = fmode;
                    l(modeN) = ll(i);
                    k(modeN) = kk(j);
                end
            end
            
        end
        
    end
    
end