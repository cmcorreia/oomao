function newArray=padArray2(array,n,m)

% function newArray=padArray2(array,n,m)
%
% This function takes a 2D matrix and extends its size n and m times by
% copying the rows and columns like so:
% [1 2 3] becomes [1 1 1 2 2 2 3 3 3]
% [4 5 6] [1 1 1 2 2 2 3 3 3]
% [4 4 4 5 5 5 6 6 6]
% [4 4 4 5 5 5 6 6 6]
% for an n = 2, m = 3.

[r,c]=size(array);
rowinds=reshape(repmat([1:r],n,1),r*n,1);
colinds=reshape(repmat([1:c],m,1),1,c*m);
indies=sub2ind([r,c],repmat(rowinds,1,c*m),repmat(colinds,r*n,1));
newArray=array(indies);
end

