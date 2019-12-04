function [MGx,MGy,Gxx,Gyy] = build_HudginDiscGradient(Nsubap, sub_samp, dx, mask)
% ------------HEADER-----------------
% Objective         ::  Build discrete 2-D gradient operator Gamma based on
%                       Hudgin geometry WFS. Based on function
%                       build_GammaSub.m from Qiang and Vogel under ../fdpcg/
%  input parameters ::
%                    1.Nsubap:      number of sub-apertures across the pupil
%                    2.sub_sampl:  subsampling of sensor grid. 0 if no
%                                  subsampling, 1 for 2-point subsampling,
%                                  2 for 3-point subsampling and so forth.
%                                  See Change record, 02/10/2008
%                    3. dx:        grid spacing, in meter
%                    4. mask:      sensor mask
%
%  output parameters ::
%                    1. Gamma_x:   2-D gradient operator with sensor masked
%                    2. Gamma_y:   2-D gradient operator with sensor masked
%                    3. Gxx:       periodic 1-D gradient operator - not any
%                                  longer
%                    4. Gyy:       periodic 1-D gradient operator
%                    5. srIdx:     sensor position where signals are measured - no longer available to
%                                  comply with sub-sampling procedure
%                    6. srMsk:     sensor mask - no longer available to
%                                  comply with sub-sampling procedure
%
%Created by         ::  Carlos Correia
%Creation date      ::  16/02/2007
%Change Record:     ::  
%                       ccorreia 02/10/2008
%                        - Update of the sub-sampled case. Now if there is
%                       sub-sampling, the 'sub_samp' argin means how many
%                       points there are between the corners of the
%                       sub-aperture. For example: If sub_samp == 1 there
%                       is a 1-point subsampling, with the sub-aperture
%                       making 3x3 pixels size. If sub_samp == 0 then each
%                       sub-aperture is a 2x2 pixels size.
%
%                       ccorreia 22/08/2008
%                       - Now I take care of the fact there are Nx(N-1) slopes in X and (N-1)xN in Y. This
%                         same result is in the work of Ren and Dekany 2004.
% ------------HEADER END----------------

if sub_samp > 0
    % add a zeros band rightewards and downwards
    mask = [mask zeros(size(mask,1),1)];
    mask = [mask;zeros(1,size(mask,2))];
    mask = subsample_mask(mask,sub_samp);
    % remove the added band
    mask = mask(1:end-1,1:end-1);
end

twoN = (Nsubap+1)+ (Nsubap)*sub_samp;
twoN_sq = twoN^2;
nskip = 2^sub_samp;
n_act = Nsubap+1;

% build 1-D periodic gradient operators
onevec = ones(twoN_sq,1);
Gyy = spdiags([onevec -onevec], [0 sub_samp+1], twoN,twoN);
Gyy(twoN,1) = -1;

   


Gyy1d = Gyy;
%cc, 10/10/06, I added the "floor" in the following cycles
% for i = 1 : twoN
%     if mod(floor(i-1-nskip/2), nskip) ~= 0
%         Gyy1d(i,:) = Gyy1d(i,:) * 0;
%     end
% end


% kron the 1-D gradient operators to retrieve 2-D gradient operators
%ccorreia: although it becomes a 2-D gradient operator, do not confound
%with the Laplacian (which is a curvature operator (2nd-order
%derivative)
% see Fast wave-front reconstruction by solving the Sylvester equation with the alternating direction implicit method
% Hongwu Ren and Richard Dekany, 2004. Explanations therein.
% cc, 21/08/2008 these next two lines have been updated. Now I take
% care of the fact there are Nx(N-1) slopes in X and (N-1)xN in Y. This
% same result is in the work of Ren and Dekany
Gy = kron(Gyy1d(1:length(Gyy1d),:), eye(length(Gyy1d)-1,length(Gyy1d)));
Gx = kron(eye(size(Gyy1d)),Gyy1d(1:length(Gyy1d)-1,:));

%Gy = kron(Gyy1d, eye(size(Gyy1d)));
%Gx = kron(eye(size(Gyy1d)),Gyy1d);


%cc : These "for" cycles could be replaced by a dyaddown and dyadup procedure:
%test=dyaddown(dyaddown(mask_est)');
%test1=dyadup(dyadup(test)');
%..has many times as the subsampling requires (once for no subsampling,
%twice for 2-element subsampling and so forth
%cc, 10/10/06, I added the floor (2x) in the following cycles
%clear srIdx;
% idxCnt = 0;
% for i = 1 : n_act-1
%     m = floor((i-1)*nskip + nskip/2 + 1);
%     for j = 1 : n_act-1
%         n = floor((j-1)*nskip + nskip/2 + 1);
%         if mask_extended(m, n) > 0
%             idxCnt = idxCnt + 1;
%             srIdx(idxCnt) = (m-1)*twoN + n;
%         end
%     end
% end

if nargin > 3
    srIdx = find(mask);
else
    srIdx = 1:twoN_sq-twoN;
end


%srIdx = 1:nnz(mask);
srMask = sparse(twoN,twoN);
%Positions pointed by srIdx are assigned to one, The remaining are
%zero-valued
srMask(srIdx) = ones(size(srIdx));
M_sr = spdiags(srMask(:), 0, twoN_sq-twoN, twoN_sq-twoN);
%M_sr = spdiags(srMask(:), 0, twoN_sq,twoN_sq);

MGx = M_sr * Gx;
MGy = M_sr * Gy;
MGx = MGx(srIdx,:);
MGy = MGy(srIdx,:);

