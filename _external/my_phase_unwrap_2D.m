        function [phase_est_unwrap] = my_phase_unwrap_2D(phase_est, mask_phase, pupwidthX, imwidthX)
        
        % -------------------------------------
        % Phase unwrapping
        % -------------------------------------
        if nargin == 1
            mask_phase = ones(size(phase_est));
        end
        Npsim = length(mask_phase);
       
        
        % -------------------------------------
        % Phase unwrapping
        % -------------------------------------
        [MGx,MGy]       = build_HudginDiscGradient(Npsim, 0, 1, mask_phase);
        
        init = ((imwidthX/pupwidthX)/2-1+1/2)*pupwidthX + 1;
        phase_est = phase_est(init:init+Npsim, init:init+Npsim);
        grad_x_phase_est = MGx*phase_est(:);
        grad_y_phase_est = MGy*phase_est(:);
        
        
        grad_x_phase_est = rem(grad_x_phase_est+pi, 2*pi)-pi;
        grad_x_phase_est = rem(grad_x_phase_est-pi, 2*pi)+pi;
        grad_y_phase_est = rem(grad_y_phase_est+pi, 2*pi)-pi;
        grad_y_phase_est = rem(grad_y_phase_est-pi, 2*pi)+pi;
        
        phase_est_unwrap = (MGx'*MGx + MGy'*MGy)\([MGx' MGy']*[grad_x_phase_est;grad_y_phase_est]);
        phase_est_unwrap = reshape(phase_est_unwrap, Npsim+1,Npsim+1);
        phase_est_unwrap = phase_est_unwrap(1:Npsim, 1:Npsim);
        id_pup = find(mask_phase);
        phase_est_unwrap = phase_est_unwrap - mean(phase_est_unwrap(id_pup));
   