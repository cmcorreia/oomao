%{
------------HEADER-----------------
Objective         ::  Compute analytically the modal projection of Zernike poluynomials using Noll's 
                      indexing onto a smaller, x and y displaced and
                      rotated pupil. Based on the function TransformC by
                      Lundstrum, JOSA A 2008
                      
INPUT_VARS
zernWFS           ::  a zernike object
tel               ::  telescope
gs                ::  the guide-star asterism
atm               ::  atmosphere

OUTPUT VARS
Proj              ::  Projection matrix with [nGs x Nz] by [Nz x Nlayers]
                      dimensions
Projcell          ::  Cell version of Proj

Created by        ::  Carlos Correia
Creation date     ::  June 12
Change Record:    ::  
------------HEADER END----------------

%}

function [Proj Projcell] = analyticalSmallFootprintExpansion(zernWfs,tel,gs,val,flagPistonRemoval)

if ~exist('flagPistonRemoval','var')
    flagPistonRemoval = 1;
end

nMode = zernWfs.nMode;
A2N   =  ANSI2Noll(nMode);
nMode = length(A2N)-1;
if isa(val,'atmosphere')
    altitudes = [val.layer.altitude];
    nL = val.nLayer;
else
    if isa(val,'double')
        altitudes = val;
        nL = length(val);
    end
end

nGs = length(gs);


Proj = zeros((nMode+1), (nMode+1),  nGs, nL);
ProjPistonRM = zeros((nMode), (nMode),  nGs, nL);
for kLayer = 1:nL
    fprintf('@(Projection)> ')
    for kPd = 1:nGs
        fprintf('gs#%d/Layer#%d - \n',kPd,kLayer)
        src = gs(kPd);
        D0      = tel.diameterAt(altitudes(kLayer))*1000;
        Dh = tel.D*1000;
        delta = altitudes(kLayer).*tan(src.zenith).*...
            [cos(src.azimuth),sin(src.azimuth)];
        tx = delta(1)*1000;
        ty = delta(2)*1000;
        for kMode = 2:nMode+1
            vec0 = zeros(1,nMode+1);
            vec0(kMode) = 1;
%             C2 = TransformC([D0 (vec0')'], Dh, tx, ty,0);
%             tmp = A2N*C2(2:nMode+2);
            
            C21 = TransformC([D0 (vec0(1:kMode)')'], Dh, tx, ty,0);
            C2x = zeros(nMode+1,1);
            C2x(1:kMode) = C21(2:kMode+1);
            tmp = A2N*C2x;
            
            if kLayer > 0
                ww(kMode) = find(tmp, 1, 'last')-1;
                %Proj(:,kMode,kPd, kLayer) = A2N*C2(2:end);
            end
            Proj(:,kMode,kPd, kLayer) = tmp;
            
        end
    end
end


[~,I]    = sort(ww+1);
Projcell = cell(nGs, nL);
for kLayer = 1:nL
    %fprintf('@(Projection)> ')
    for kPd = 1:nGs
        %fprintf('pd#%d/dm#%d - ',kPd,kLayer)
        Proj(:,:,kPd, kLayer) = Proj(:,I,kPd, kLayer);
        Projcell{kPd, kLayer} = Proj(:,:,kPd, kLayer);
        if flagPistonRemoval == 1
            ProjPistonRM(:,:,kPd, kLayer) = Proj(2:end,2:end,kPd, kLayer);
            Projcell{kPd, kLayer} = ProjPistonRM(:,:,kPd, kLayer);
        end
        %Proj{kPd, kLayer} = Proj{kPd, kLayer}(:, I);
    end
end


Proj = cell2mat(Projcell);
%Proj = reshape(Proj, nL*(nMode+1), nL*(nMode+1));
%idx = [2:45 47:90];
%Proj1 = Proj(:,idx);


   
