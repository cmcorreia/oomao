classdef multiCalibrationVault < handle
    %% MULTICALIBRATIONVAULT Create a multiCalibrationVault object
    
    properties
        % Cell with calibrationVault objects
        calibVault
        % number of guide-stars
        nGs
        % number of DMs
        nDm
    end
    
    methods
        
        %% Constructor
        function obj = multiCalibrationVault(calibVault,varargin)
            obj.calibVault = calibVault;
            obj.nGs = size(calibVault,2);
            obj.nDm = size(calibVault,1);
        end
        
        function slopesPsolMat = mtimes(obj,dmStack)
            slopesPsolCell = cell(obj.nGs, obj.nDm);
            for kStar = 1: obj.nGs
                for kDm = 1:obj.nDm
                    if kDm == 1
                        slopesPsolCell{kDm,kStar} = obj.calibVault{kDm,1}.D*dmStack.dms{kDm}.coefs;
                    else
                        slopesPsolCell{kDm,kStar} = obj.calibVault{kDm,kStar}.D*dmStack.dms{kDm}.coefs;
                    end
                end
                slopesPsolMat(:,kStar) = sum([slopesPsolCell{:,kStar}],2);
            end
        end
    end
end
