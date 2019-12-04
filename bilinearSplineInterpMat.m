function varargout = bilinearSplineInterpMat(gs,atm,tel,groundGrid,tau,overSamp)
%% BILINEARSPLINEINTERP Bilinear interpolation matrix
%
% [H,mask] = bilinearSplineInterpMat(gs,atm,tel,groundGrid) computes the
% bilinear interpolation matrices for a given system made of source,
% atmosphere and telescope objects and from the pupil sampling grid at the
% ground

if nargin > 4 
    frameTime   = tau;
else
    frameTime = 0;
end

if nargin < 6
    overSamp   = ones(1,atm.nLayer);
end

nLayer      = atm.nLayer;
nGs         = length(gs);
nPxGround   = length(groundGrid);

% Ground sampling
% [xGround,yGround] ...
%            = utilities.cartAndPol(nPxGround,tel.R);
[xGround,yGround] ...
           = meshgrid(linspace(-1,1,nPxGround)*tel.R);
xGround    = xGround(groundGrid);
yGround    = yGround(groundGrid);
groundGrid = groundGrid(:)';

pitchGround = tel.D/(nPxGround-1);

% cell of bi-linear interpolation operators
H    = cell(nGs , nLayer);
% cell of phase layer mask where the bilinear spline are non-zeros
mask = cell(1   , nLayer);

fprintf('___ BI-LINEAR INTERPOLATION OPERATOR ___\n')


for kLayer = 1:nLayer
    
    pitchLayer = pitchGround/overSamp(kLayer);
    
    % Layer sampling
    D = atm.layer(kLayer).D;
%     nPxLayer        = ceil(D/pitchGround) + 1;
%    nPxLayer        = floor(D/pitchGround) + 1;
    nPxLayer        = floor(D/pitchLayer) + 1;
%    newD = (nPxLayer-1)*pitchGround;
    newD = (nPxLayer-1)*pitchLayer;
    while newD<D
        nPxLayer        = nPxLayer + 2;
%        newD = (nPxLayer-1)*pitchGround;
        newD = (nPxLayer-1)*pitchLayer;
    end       
    D = newD;
    ws = atm.layer(kLayer).windSpeed;
    theta = atm.layer(kLayer).windDirection;
    deltax = ws*cos(theta)*frameTime;
    deltay = ws*sin(theta)*frameTime;
%    [xLayer,yLayer] = utilities.cartAndPol(nPxLayer,D/2);
     [xLayer,yLayer] = meshgrid(linspace(-1,1,nPxLayer)*D/2);

     xLayer = xLayer - deltax;
     yLayer = yLayer - deltay;

    mask{kLayer} = false;
    
    for kGs=1:nGs
        
        fprintf(' [%d,%d]',kGs,kLayer)
        titleString = sprintf('guide star %d,layer %d',kGs,kLayer);
        height = atm.layer(kLayer).altitude;
        % pupil center in layer
        beta   = gs(kGs).directionVector*height;
        scale = 1-height/gs(kGs).height;
%        if height==0
%            mask{kLayer} = mask{kLayer} | groundGrid;
%            nH = numel(xGround);
%            mH = nPxGround^2;
%            iH = 1:nH;
%            jH = find(groundGrid);
%            sH = ones(1,length(jH));
%             figure
%             title(titleString)
%             hold
%             H{kGs,kLayer} = sparse(iH,jH,sH,nH,mH);
%             plot(xLayer(:),yLayer(:),'.')
%                   plot(xGround*scale+beta(1),yGround*scale+beta(2),'r.')

%         elseif kGs == nGs
%            % figure
%            % title(titleString)
%            % hold
%             H{kGs,kLayer} = bilinearSplineInterp(xLayer,yLayer,pitchGround,xGround*scale + beta(1)-deltax,yGround*scale + beta(2)-deltay);
%             mask{kLayer} = mask{kLayer} | ( ~all(H{kGs,kLayer}==0) );
%        else
           % figure
           % title(titleString)
           % hold
           %% CC added :: sign + before beta instead of - :: makes the calculation HCphiH' equal to Cxx the spatio-angular covariance matrix
            %H{kGs,kLayer} = bilinearSplineInterp(xLayer,yLayer,pitchGround,xGround*scale + beta(1),yGround*scale + beta(2));
            H{kGs,kLayer} = bilinearSplineInterp(xLayer,yLayer,pitchLayer,xGround*scale + beta(1),yGround*scale + beta(2));
            mask{kLayer} = mask{kLayer} | ( ~all(H{kGs,kLayer}==0) );

                %  plot(xLayer(:),yLayer(:),'.')
                %  plot(xGround*scale+beta(1),yGround*scale+beta(2),'r.')

%             
%        end
        
    end
    
    fprintf('\n')
    
end

fprintf('----------------------------\n')

varargout{1} = H;
varargout{2} = mask;

end

function H = bilinearSplineInterp(xo,yo,do,xi,yi)
%% BILINEARSPLINEINTERP Bilinear interpolation
%
% H = bilinearSplineInterp(xo,yo,do,xi,yi,d); computes the sparse matrix H
% to perform the bilinear interpolation zi = H*zo, where zo and zi are
% defined on the meshes [xo;yo] and [xi;yi], respectively. do is the mesh
% step size of [xo;yo]

xo = xo(:)';
yo = yo(:)';
xi = xi(:);
yi = yi(:);

ni = length(xi);
if ni~=length(yi)
    error('xi and yi must have the same length.')
end
no = length(xo);
if no~=length(yo)
    error('xo and yo must have the same length.')
end

u = bsxfun(@minus,xi,xo)/do;
v = bsxfun(@minus,yi,yo)/do;

H = linearSpline(u).*linearSpline(v);

% [i,j,s] = find(H);
% [m,n]   = size(H);
% as      = abs(s);
% index   = as > eps(max(as)); % round-off error filter
%
% H = sparse(i(index),j(index),s(index),m,n);


    function y = linearSpline(x)
        %% LINEARSPLINE Linear spline function
        %
        % y = linearSpline(x) computes the function y = 1 - |x| for |x|<1 and y = 0
        % elsewhere
        
        [m,n] = size(x);
        x     = abs(x);
        index = x < 1;
        [i,j] = find(index);
        s     = 1 - x(index);
        y = sparse(i,j,s,m,n);
        
    end

end
