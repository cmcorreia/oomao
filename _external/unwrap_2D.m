function im_unwrapped = unwrap_2D(im, im_mask)
% Objective:
%   Compute the PSF of an AO image
% INPUT:
%   phase       ::  A matrix with the phase [pupwidth x pupwidth]
%
% OUTPUT:
%  phase_out    ::  2-pi period unwrapped phase in 2 D
%Created by     ::  Carlos Correia 
%Creation date  ::  22 march 12
%Change Record: ::  
% ------------HEADER END-----------------

%% REPLACE WITH YOUR imAGES
if nargin < 2
    im_mask=ones(size(im));                     %Mask (if applicable)
end

%%
im_mag=abs(im);                             %Magnitude image
im_phase=angle(im);                         %Phase image
im_unwrapped=zeros(size(im));               %Zero starting matrix for unwrapped phase
adjoin=zeros(size(im));                     %Zero starting matrix for adjoin matrix
unwrapped_binary=zeros(size(im));           %Binary image to mark unwrapped pixels

%% Calculate phase quality map
im_phase_quality=PhaseDerivativeVariance(im_phase);   

%% Identify starting seed point on a phase quality map
minp=im_phase_quality(2:end-1, 2:end-1); minp=min(minp(:));
maxp=im_phase_quality(2:end-1, 2:end-1); maxp=max(maxp(:));
%figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off; title('Phase quality map'); 
%uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
%[xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
[xpoint,ypoint] = find(im_phase_quality == minp);

%% Unwrap
colref=round(xpoint); rowref=round(ypoint);
im_unwrapped(rowref,colref)=im_phase(rowref,colref);                        %Save the unwrapped values
unwrapped_binary(rowref,colref,1)=1;
if im_mask(rowref-1, colref, 1)==1 adjoin(rowref-1, colref, 1)=1; end       %Mark the pixels adjoining the selected point
if im_mask(rowref+1, colref, 1)==1 adjoin(rowref+1, colref, 1)=1; end
if im_mask(rowref, colref-1, 1)==1 adjoin(rowref, colref-1, 1)=1; end
if im_mask(rowref, colref+1, 1)==1 adjoin(rowref, colref+1, 1)=1; end
im_unwrapped=GuidedFloodFill(im_phase, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

% figure; imagesc(im_mag), colormap(gray), axis square, axis off; title('Magnitude image'); 
% figure; imagesc(im_phase), colormap(gray), axis square, axis off; title('Wrapped phase'); 
% figure; imagesc(im_unwrapped), colormap(gray), axis square, axis off; title('Unwrapped phase'); 


% if nargin < 2
%     im_mask=ones(size(im));                     %Mask (if applicable)
% end%%
% 
% im_mag=abs(im);                             %Magnitude image
% im_phase=angle(im);                         %Phase image
% 
% 
% %%  Set parameters
% max_box_radius=4;                           %Maximum search box radius (pixels)
% threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image
% 
% %% Unwrap
% residue_charge=PhaseResidues(im_phase, im_mask);                            %Calculate phase residues
% branch_cuts=BranchCuts(residue_charge, max_box_radius, im_mask);            %Place branch cuts
% [im_unwrapped, rowref, colref]=FloodFill(im_phase, branch_cuts, im_mask);   %Flood fill phase unwrapping
% 
% %% Display results
% figure; imagesc(residue_charge), colormap(gray), axis square, axis off, title('Phase residues (charged)');
% figure; imagesc(branch_cuts), colormap(gray), axis square, axis off, title('Branch cuts');
% figure; imagesc(immultiply(im_phase,im_mask)), colormap(gray), axis square, axis off, title('Wrapped phase');
% tempmin=min(min(im_unwrapped));          %This bit is just done to create a pleasing display when a mask is used
% temp=(im_unwrapped==0);
% temp_im=im_unwrapped;
% temp_im(temp)=tempmin;
% figure; imagesc(temp_im), colormap(gray), axis square, axis off, title('Unwrapped phase');
% 
