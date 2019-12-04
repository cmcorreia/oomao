function [ControlMatrix] = tsvd(G, no_svd_to_truncate, U, S, V)
% ------------HEADER-----------------
% Objective         ::  Compute generalised inverse using the truncated-SVD
%                       method
%  input paramsters:
%                    1.G:          Matrix to invert  
%                    2.no_svd_to_truncate:  number of svd values to truncate 
%                    3. U, S, V:   the [U S V] = svd(G) are the svd
%                                  decomposition matrices
%  output parameters:
%                    1. ControlMatrix: the generalised inverse matrix  
%
%Created by         ::  Carlos Correia 
%Creation date      ::  16/02/2007
%Change Record:     ::  
% ------------HEADER END----------------

if nargin < 3
    if issparse(G)
        [U, S, V] = svds(G);
    else
        [U, S, V] = svd(G);
    end
end
if nargin == 1
    no_svd_to_truncate = 0;
end
keep = size(diag(S),1) - no_svd_to_truncate;
Ured = U(:,1:keep);
Vred = V(:,1:keep);
S = S(1:keep,1:keep);
invS =diag(1./diag(S));
invS = invS(1:keep,1:keep);
ControlMatrix = Vred * invS * Ured';


