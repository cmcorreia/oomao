classdef psfrTools
    % Gathering useful functions for PSF-Reconstruction/QSO/CG classes
    
    methods (Static)
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- GENERAL RECEPIES ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [m,n] = minmax(x)
            m = min(x(:));
            n = max(x(:));
        end
        
        function out = sr2wfe(SR,lambda)
            out = 1e9*sqrt(-log(SR))*lambda/2/pi;
        end
        
        function out = wfe2sr(wfe,lambda)
            out = 100.*exp(-(2*pi*wfe*1e-9/lambda)^2);
        end
        
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- IMAGES PROCESSING ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
        function [img,flux] = processImage(frame,bg,flat,badPixMap,overSampling,varargin)
            %{
            im              :: raw image
            N               :: nPixOut number of pixels on the cropped, output PSF
            badPiMap        :: Bad pixels map calibration
            rebin           :: interpolation factor used when recentring the PSF
            nFit            :: number of sample points on the OTF used in a
                               polynomial fit to estimate the noise level at the OTF(0)
                               location -- ONLY RELEVANT WITH SPATIALLY
                               UNCORRELATED NOISE
            thresholding    :: threshold image to user-defined value
                                pixels
            flagCCD         :: flag to deconvolve from detector pixel response
            peakLocation    :: default 'middlePixel' meaning the the max of
                               the PSF is centred at one single pixel. Otherwise it's at the
                               intersection of 4 pixels
            %}
            
            inputs = inputParser;
            inputs.addParameter('fovInPixel',0,@isnumeric);
            inputs.addParameter('rebin',0,@isnumeric);
            inputs.addParameter('masking',false,@islogical);
            inputs.addParameter('tfccd',false,@islogical);
            inputs.addParameter('nfit',0,@isnumeric);
            inputs.addParameter('thresholding',-Inf,@isnumeric);            
            inputs.addParameter('peakLoc','middlePixel',@ischar);
            inputs.parse(varargin{:});
            
            fovInPixel   = inputs.Results.fovInPixel;
            rebin        = inputs.Results.rebin;
            masking      = inputs.Results.masking;
            tfccd        = inputs.Results.tfccd;
            nfit         = inputs.Results.nfit;
            thresholding = inputs.Results.thresholding;
            peakLoc      = inputs.Results.peakLoc;
            
            %01. ------ Unbiasing from background and flat                 
                        
            % Background map            
            if any(bg(:))
                rawImg = frame - median(bg,3);      
            else
                rawImg = frame;
            end
            
            % Bad pixel map
            rawImg = ~badPixMap.*rawImg;
                                    
            % Flat field map
            rawImg = rawImg./flat;
            rawImg(rawImg~=rawImg) = 0;     
            
            if median(frame(:)) ~=0
                rawImg = rawImg - median(rawImg(:));
            end
            
            %02. ------- Take care of bad pixels by interpolating from the avalaible neighbors
            if any(badPixMap(:))
                supBad = tools.createDeadPixFrame(badPixMap);
                img    = tools.corrDeadPixFrame(supBad, rawImg );
            else
                img   = rawImg;
            end                      
                                                
            %03. ------- Find the PSF peak
            % Perform a median filtering to get rid of hot pixels
            n    = 1;
            border = 10;
            img2 = img(border+1:end-border,border+1:end-border);     
            for i=1:n
               img2(:) = medfilt1(img2(:));
            end
            
            [x0,y0] = find(img2 == max(img2(:)));
            x0 = x0 + border;
            y0 = y0 + border;
            
            % Explore  non-filtered value to really pinpoint the max value
            nLoc = length(x0);
            pVal = zeros(1,nLoc);
            for i=1:nLoc
                pVal(i) = img(x0(i),y0(i));
            end            
            x0 = x0(find(pVal == max(pVal)));
            y0 = y0(find(pVal == max(pVal)));
            
            %04. ------- Crop PSF to the wished field of view
            if fovInPixel
                nImg    = size(img,2);
                if fovInPixel < 2*min([nImg-x0 nImg-y0])
                    idx = round(x0 - fovInPixel/2 + 1):round(x0+fovInPixel/2);
                    idy = round(y0 - fovInPixel/2 + 1):round(y0+fovInPixel/2);
                    img = img(idx,idy);
                end
            end           
            mx = max(img(:));
                  
            
            %05. -------- Finely recenter the image on the midlle pixel
            if rebin
                [~,otf] = tools.recenterPsf(img,rebin);
            else
                otf = tools.psf2otf(img);
            end
                      
            %06. ------ Remove numerical noise using frequencies values higher than D/lambda
            if masking
                notf        = size(otf,2);
                u1D         = (-1+1/2/notf:2/notf:1-1/2/notf)*overSampling;
                [u2Dx,u2Dy] = meshgrid(u1D);
                msk         = hypot(u2Dx,u2Dy)<=1;
                otf         = otf.*msk;
            end
            
            %07. ------ Compensate for the CCD pixel transfer function
            if tfccd
                if ~exist('u2Dx','var')
                    notf        = size(otf,2);
                    u1D         = (-1+1/2/notf:2/notf:1-1/2/notf)*overSampling;
                    [u2Dx,u2Dy] = meshgrid(u1D);
                end
                ftccd    = abs(sinc(u2Dx/pi/2)).*abs(sinc(u2Dy/pi/2));
                ftccd    = ftccd/max(ftccd(:));
                msk      = ftccd > 1e-3;
                otf(msk) = otf(msk)./ftccd(msk);
            end
            
            %08 ..... Remove the background by interpoling the OTF peak
            if nfit
                otfRad = radial(otf);
                uRad   = linspace(0,1,nfit+2);
                fit    = polyfit(uRad(2:end),otfRad(2:nfit+2),nfit);
                yfit   = 0*uRad;
                for k =1:nfit+1
                    yfit   = yfit + fit(k)*uRad.^(nfit-k+1);
                end
                otf(end/2+1,end/2+1) = yfit(1);
            end
            
            %09. ......  Back to the angular separation domain
            if rebin | nfit | masking | tfccd                
                img = tools.otf2psf(otf,peakLoc);
                img = img/max(img(:))*mx;
            end
                         
            
            %10. ...... Thresholding            
            img(img<thresholding) = 0;
            
            %11. ......  Estimate the Aperture flux
            flux  = tools.getFlux(img);  
            
        end
        
        
        function [imCor,otf_lr] = recenterPsf(psf,overSampling)
            flux          = sum(psf(:));
            [npsfx,npsfy] = size(psf);
            npsfx2        = npsfx*overSampling;
            npsfy2        = npsfy*overSampling;
            
            % Get the high-resolution PSF
            if overSampling > 1
                psf_hr = tools.interpolateOtf(psf,npsfx2);
            else
                psf_hr = psf;                
            end
            % Get the max value
            mx        = max(psf_hr(:));
            [idx,idy] = find(psf_hr == mx);
            dx        = floor(npsfx2/2-idx)+1;
            dy        = floor(npsfy2/2-idy)+1;
            if (dx~=0) | (dy~=0)
                % Get the OTF
                otf_hr = tools.psf2otf(psf_hr);
                % Apply the Phasor
                [u,v]     = freqspace(length(otf_hr),'meshgrid');
                fftPhasor = exp(-1i.*pi.*(u*dy+v*dx));
                otf_hr    = otf_hr.*fftPhasor;
                % Get the PSF low-resolution
                imCor  = tools.otf2psf(otf_hr);
                imCor  = tools.interpolateOtf(imCor,npsfx);
                imCor  = flux*imCor/sum(imCor(:));
                otf_lr = tools.psf2otf(imCor);
            else
                imCor = psf;
                otf_lr = tools.psf2otf(imCor);
            end
        end
        
        function frame = createDeadPixFrame(badPixelMap)
            %dpframe = createDeadPixFrame(badPixelMap)
            %badPixelMap is the map of dead pixels
            %frame  is the image to be corrected
            
            %The dead pixel is replaced by a weighted average of the neighbours,
            %1 2 1
            %2 X 2
            %1 2 1
            %when they are "available". "Available" means that the sum of the
            %weights of neighbouring pixels must exceeds 4.
            
            %If no neighbouring pixel is available, the dead pixel is not
            %corrected, but a new "dead pixel map" is created, and the function is
            %called once again (recursive calls).
            
            % Get the number of dead pixels
            [sx,sy]      = size(badPixelMap);
            npixnoncorr  = 0;
            [nnx,nny]    = find(badPixelMap);
            nn1D         = find(badPixelMap(:));
            nDeadPix     = length(nn1D);
            %Instantiation
            tmp          = badPixelMap*0;
            frame        = zeros(nDeadPix,10,2); %2nd row: #pixel (one pixel + 8 neighbors)
            frame(:,:,1) = 1;                    %3rd row: adresses
            
            %loop on Pixel
            for i=1:nDeadPix
                nb = 2;
                frame(i,1,1) = nn1D(i);  % 1st row = bad pixel
                frame(i,2,1) = 0;        % number of used neighbour pixel for correction
                x            = nnx(i);
                y            = nny(i);
                wcum         = 0;
                
                % Edges neighbours
                if x>0 && x<=sx && y+1>0 && y+1<=sy
                    if ~badPixelMap(x,y+1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i) + sx;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                
                if x>0 && x<=sx && y-1>0 && y-1<=sy
                    if ~badPixelMap(x,y-1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-sx;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                if x+1>0 && x+1<=sx && y>0 && y<=sy
                    if ~badPixelMap(x+1,y)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)+1;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                if x-1>0 && x-1<=sx && y>0 && y<=sy
                    if ~badPixelMap(x-1,y)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-1;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                
                % Diagonal neighbours
                if x+1>0 && x+1<=sx && y+1>0 && y+1<=sy
                    if ~badPixelMap(x+1,y+1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)+1+sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                if x-1>0 && x-1<=sx && y+1>0 && y+1<=sy
                    if ~badPixelMap(x-1,y+1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-1+sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                if x+1>0 && x+1<=sx && y-1>0 && y-1<=sy
                    if~badPixelMap(x+1,y-1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)+1-sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                if x-1>0 && x-1<=sx && y-1>0 && y-1<=sy
                    if~badPixelMap(x-1,y-1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-1-sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                
                % Take decision regarding the number of avalaible neighbours
                if wcum<4   %not enough neigbours
                    npixnoncorr      = npixnoncorr + 1;
                    tmp(x,y)         = tmp(x,y) + 1;
                    frame(i,3:end,1) = 1;    % pixel adresses set to 1
                    frame(i,:,2)     = 0;    % weights set to 0
                else
                    frame(i,2,1)     = nb;    %number of correcting pixels
                end
            end
            
            if npixnoncorr
                frame_suppl = tools.createDeadPixFrame(tmp);
                nSup        = size(frame_suppl,1);
                %Frame concatenation
                new_frame                     = zeros(nDeadPix+nSup,10,2);
                new_frame(1:nDeadPix,:,:)     = frame;
                new_frame(1+nDeadPix:end,:,:) = frame_suppl;
                frame                         = new_frame;
            end
            
        end
        
        function imCor = corrDeadPixFrame( frame, im )
            % imCor = corrDeadPixFrame(frame,im)
            %Correcting the bad pixel on a image from the bad pixel map
            
            npixd = size(frame,1);
            imCor = im;
            for i=1:npixd
                w =  sum(frame(i,3:end,2));
                if w~=0
                    imCor(frame(i,1,1)) = sum( im(frame(i,3:end,1)) .* frame(i,3:end,2)) / double(w);
                end
            end
        end
        
        function imagerot = rotateIm(image,degree)                        
            switch mod(degree, 360)
                % Special cases
                case 0
                    imagerot = image;
                case 90
                    imagerot = rot90(image);
                case 180
                    imagerot = image(end:-1:1, end:-1:1);
                case 270
                    imagerot = rot90(image(end:-1:1, end:-1:1));

                    % General rotations
                otherwise                    
                    % Convert to radians and create transformation matrix
                    a = degree*pi/180;
                    R = [+cos(a) +sin(a); -sin(a) +cos(a)];
                    
                    % Figure out the size of the transformed image
                    [m,n,p] = size(image);
                    dest = round( [1 1; 1 n; m 1; m n]*R );
                    dest = bsxfun(@minus, dest, min(dest)) + 1;
                    imagerot = zeros([max(dest) p],class(image));
                    
                    % Map all pixels of the transformed image to the original image
                    for ii = 1:size(imagerot,1)
                        for jj = 1:size(imagerot,2)
                            source = ([ii jj]-dest(1,:))*R.';
                            if all(source >= 1) && all(source <= [m n])
                                
                                % Get all 4 surrounding pixels
                                C = ceil(source);
                                F = floor(source);
                                
                                % Compute the relative areas
                                A = [...
                                    ((C(2)-source(2))*(C(1)-source(1))),...
                                    ((source(2)-F(2))*(source(1)-F(1)));
                                    ((C(2)-source(2))*(source(1)-F(1))),...
                                    ((source(2)-F(2))*(C(1)-source(1)))];
                                
                                % Extract colors and re-scale them relative to area
                                cols = bsxfun(@times, A, double(image(F(1):C(1),F(2):C(2),:)));
                                
                                % Assign
                                imagerot(ii,jj,:) = sum(sum(cols),2);                                
                            end
                        end
                    end
                    imagerot = tools.crop(imagerot,size(image));
            end
        end                               
  
                        
        function out = enlargePupil(P,n)
            % TO BE MERGED WITH enlargeOtf
            [nx,ny,nf] = size(P);
            %if n*nx ~= round(n*nx)
            %    fprintf('Warning: non integer support, ratio adjusted to %f\n',round(n*nx)/nx);
            %end
            xi      = floor((nx*n-nx)/2+1);
            xf      = floor((nx*n+nx)/2);
            yi      = floor((ny*n-ny)/2+1);
            yf      = floor((ny*n+ny)/2);
            if length(size(P)) == 2
                out  = zeros(round(n*nx),round(ny*n));
                out(xi:xf,yi:yf) = P;
            elseif length(size(P)) == 3
                out = zeros(round(n*nx),round(ny*n),nf);
                out(xi:xf,yi:yf,:) = P;
            end
        end
        
        function out = enlargeOtf(otf,n)
            if true
                out = padarray(otf,floor((n-1)*size(otf)/2),'both');
            else
                % Otf sizes
                nx  = size(otf,1);
                nx2 = round(n*nx);
                out = zeros(nx2);
                
                % Zero-padding
                if ~iseven(nx2)
                    idx = floor(0.5*(nx2-nx) + 1): floor(0.5*(nx2 + nx));
                else
                    idx = floor(0.5*(nx2-nx)+1: 0.5*(nx2+nx));
                end
                idx = floor(nx2/2 + 1 - nx/2): floor(nx2/2 + nx/2);
                out(idx,idx) = otf;
            end
        end
        
        function out = interpolateOtf(otf,nRes,method)
            if nargin < 3
                method = 'spline';
            end
            % Define angular frequencies vectors
            notf = size(otf,1);
            
            if iseven(notf)
                u1D = (-notf/2:1:notf/2-1)*2/notf;
            else
                u1D = (-floor(notf/2):1:floor(notf/2))*2/notf;
            end               
            if iseven(nRes)
                u1D2 = (-nRes/2:1:nRes/2-1)*2/nRes;
            else
                u1D2 = (-floor(nRes/2):1:floor(nRes/2))*2/nRes;
            end          
           % Create meshgrids
           [Uxi,Uyi]= meshgrid(u1D);
           [Uxo,Uyo]= meshgrid(u1D2);
            
            % Interpolation
            out = interp2(Uxi,Uyi,otf,Uxo,Uyo,method);
        end
        
        function out = interpolate(P,nRes,method)
            if nargin<3
                method='linear';
            end
            
            nP       = length(P);
            xi       = linspace(-1,1,nP);
            xo       = linspace(-1,1,nRes);
            [Xi,Yi]  = meshgrid(xi);
            [Xo,Yo]  = meshgrid(xo);
            out      = interp2(Xi,Yi,P,Xo,Yo,method);
        end
               
        
        function im_ = translateImage(image,dX,val)
            if nargin < 3
                val = 0;
            end
            
            % Zero-pad the suppot to avoid fft buffer circulation effect
            [nx,ny] = size(image);
            % nx -> rows
            % ny -> columns
            im_ = padarray(image,floor([nx,ny]/2),'both');            
            % Derives the fourier phasor
            dx = dX(1);
            dy = dX(2);
            [u,v] = freqspace(size(im_),'meshgrid');
            phasor = exp(-1i*pi*(u*dy+v*dx));
            % Translates the image
            otf = tools.psf2otf(im_);
            otf = otf/max(otf(:));
            im_ = tools.otf2psf(otf.*phasor);
            im_ = im_/sum(im_(:))*sum(image(:));
            % Empty zone filling     
            if any(size(im_)>size(image))
                ixi = round(1:nx/2 + dx);
                ixf = round(3*nx/2+dx+1):size(im_,1);
                iyi = round(1:ny/2 + dy);
                iyf = round(3*ny/2+dy+1):size(im_,2);
                im_(ixi,:) = val;
                im_(ixf,:) = val;
                im_(:,iyi) = val;
                im_(:,iyf) = val;
                % Cropping
                im_ = tools.crop(im_,[nx,ny]); 
            end
            
        end
                                                               
         function [imCor,otf_lr] = recenterPSF(psf,overSampling)
            flux          = sum(psf(:));
            [npsfx,npsfy] = size(psf);
            npsfx2        = npsfx*overSampling;
            npsfy2        = npsfy*overSampling;
            
            % Get the high-resolution PSF
            if overSampling > 1
                psf_hr = tools.interpolateOtf(psf,npsfx2);
            else
                psf_hr = psf;                
            end
            % Get the max value
            mx        = max(psf_hr(:));
            [idx,idy] = find(psf_hr == mx);
            idx = idx(1);
            idy = idy(1);
            dx        = floor(npsfx2/2-idx)+1;
            dy        = floor(npsfy2/2-idy)+1;
            if (dx~=0) | (dy~=0)
                % Get the OTF
                otf_hr = tools.psf2otf(psf_hr);
                % Apply the Phasor
                [u,v]     = freqspace(length(otf_hr),'meshgrid');
                fftPhasor = exp(-1i.*pi.*(u*dy+v*dx));
                otf_hr    = otf_hr.*fftPhasor;
                % Get the PSF low-resolution
                imCor  = tools.otf2psf(otf_hr);
                imCor  = tools.interpolateOtf(imCor,npsfx);
                imCor  = flux*imCor/sum(imCor(:));
                otf_lr = tools.psf2otf(imCor);
            else
                imCor = psf;
                otf_lr = tools.psf2otf(imCor);
            end
        end
           
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- FOURIER RECIPIES ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
       
        function out = fftCorrel(x,y)
            nPts = length(x(:));
            out  = ifft2(fft2(x).*conj(fft2(y)))./nPts;
        end
        
        function out = rolledFFT(x)
            nPts = length(x);
            if ~iseven(nPts)
                out  = fftshift(fft2(fftshift(x)))./nPts^2;
            else
                out  = ifftshift(fft2(ifftshift(x)))./nPts^2;
            end
        end
        
       function out = convolve(object,PSF)
           out = fftshift(ifft2(fft2(object).*fft2(PSF)));
           out = sum(object(:)).*out./sum(out(:));
       end             
                
       function out = telescopeOtf(pupil,overSampling)
           extendedPup  = tools.enlargeOtf(pupil,overSampling);
           out = fftshift(tools.fftCorrel(extendedPup,extendedPup));
       end
       
       function out = telescopePsf(pupil,overSampling,peakLocation)
            if nargin <3
               peakLocation = 'middlePixel';
            end
            nSize = size(pupil,1);
            
            if overSampling >1
                otf = tools.telescopeOtf(pupil,overSampling);
                otf = tools.interpolateOtf(otf,nSize);
                out = tools.otf2psf(otf,peakLocation);
            else
                
                otf = tools.telescopeOtf(pupil,2);
                otf = tools.interpolateOtf(otf,nSize/overSampling);
                out = tools.otf2psf(otf);
                out = tools.interpolateOtf(out,nSize);
           end
           
       end
       
       function out = pupil2psf(pupil,phase,overSampling)
           if nargin <2
               phase = 0.*pupil;
           end
           if nargin < 3
               overSampling = 1;
           end
           otf = tools.pupil2otf(pupil,phase,overSampling);
           out = tools.otf2psf(otf);
       end
       
       function out = pupil2otf(pupil,phase,overSampling)
           if nargin <2
               phase = 0.*pupil;
           end
           if nargin < 3
               overSampling = 1;
           end
           P   = tools.enlargeOtf(pupil,2*overSampling); 
           phi = tools.enlargeOtf(phase,2*overSampling);
           E   =  P.*exp(1i.*phi);
           out = fftshift(tools.fftCorrel(E,E));
       end
       
       function out = psd2cov(psd,pixelScale)
           % CAUTION !! THE FFT2 ALGORTIHM IS COMPUTING SPATIAL FREQUENCIES
           % ON AN EVEN SIZE GRID. THE PSD MUST BE PROVIDED ON SUCH A SAME
           % WAY TO CONSERVE THE SYMMETRY
           out = fft2(fftshift(psd));
           out = out.*pixelScale^2;
       end
       
       function out  = cov2sf(cov)
           out = 2*max(cov(:)) - cov - conj(cov);
       end
       
       function out = sf2otf(sf)
           out  = exp(-0.5 * sf);
       end
       
       function out = otf2psf(otf,peakLocation)
           if nargin <2
               peakLocation = 'middlePixel';
           end
           
           nSize     = size(otf,1);
           [u,v]     = freqspace(nSize,'meshgrid');
           
           if strcmp(peakLocation,'middlePixel')
               % ensures to put the psf max on the middle pixel
               fftPhasor = 1;
           else
               fftPhasor = exp(1i.*pi.*(u+v).*0.5);
               %Max spead out on the four middle pixels of the image
           end
           if  iseven(nSize)
               out       = real(fftshift(ifft2(fftshift(otf.*fftPhasor))));
           else
               out       = real(fftshift(ifft2(ifftshift(otf.*fftPhasor))));
           end
           out       = out./sum(out(:));
       end
       
       function out = otfShannon2psf(otfShannon,Samp,fov)
           
           if Samp > 1
               %     % Zero-pad the OTF to get the good PSF pixel scale
               otf    = padarray(otfShannon,round((Samp-1)*size(otfShannon)/2),'both');
               % Interpolate the OTF to crop the PSF to the desired FOV
               otf    = tools.interpolateOtf(otf,fov);
               otf    = otf/max(otf(:));
               out    = tools.otf2psf(otf);
           elseif Samp ==1
               otf    = tools.interpolateOtf(otf,fov);
               otf    = otf/max(otf(:));
               out    = tools.otf2psf(otf);
           else
               % OTF are derived for a Nyquist-sampled PSF
               otf        = tools.interpolateOtf(otfShannon,fov/Samp);
               otf(otf<0) = 0;
               otf        = otf/max(otf(:));
               out        = tools.otf2psf(otf);
               out        = tools.interpolateOtf(out,obj.fov);
           end           
       end
       
       function out = psf2otf(psf)
           out = (tools.rolledFFT(psf))/sum(psf(:));
       end
       
       function out = psd2otf(psd,pixelScale)
           cov = tools.psd2cov(psd,pixelScale);
           sf  = tools.cov2sf(cov);
           out = tools.sf2otf(sf);
       end
       
       function out = psd2psf(psd,pixelScale)
           otf = fftshift(tools.psd2otf(psd,pixelScale));
           out = tools.otf2psf(otf);
       end
       
       
      
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%---------------- TELEMETRY PROCESSING  ------------------------
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
          function out = getNoiseVariance(s,varargin)
            
            inputs = inputParser;
            inputs.addRequired('s',@isnumeric);
            inputs.addParameter('nfit',0,@isnumeric);
            inputs.addParameter('nshift',0,@isnumeric);
            inputs.addParameter('psd',false,@islogical);
            inputs.addParameter('RTF',[],@isnumeric);
            inputs.parse(s,varargin{:});
            
            nfit   = inputs.Results.nfit;
            nshift = inputs.Results.nshift;
            psd    = inputs.Results.psd;
            RTF    = inputs.Results.RTF;
            
            % Mean removal
            s       = squeeze(s);
            s       = bsxfun(@minus,s,mean(s,2));
            
            [nS,nF] = size(s);
            out     = zeros(nS);
            
            if nfit
                % Polynomial fitting procedure
                delay   = linspace(0,1,nfit+2);
                for i=1:nS
                    g      = (ifft(fft(s(i,:)).*conj(fft(s(i,:))))/nF);
                    mx     = max(g(:));
                    fit    = polyfit(delay(2:end),g(2:nfit+2),nfit);
                    yfit   = 0*delay;
                    for k =1:nfit+1
                        yfit   = yfit + fit(k)*delay.^(nfit-k+1);
                    end
                    out(i,i) = mx - yfit(1);
                end
                out = diag(out);
            elseif nshift
                ds  = s - circshift(s,[0,-nshift]);
                out = s*ds'/nF;
                %equal to the difference between the auto-corr of s
                %and the nshift frame-delay cross-cor of s taken at zero
            elseif psd
                if isempty(RTF)
                    RTF = 1;
                end
                
                fftS = fft(s,[],2)/nF;
                fftS = fftS(:,1:floor(end/2))./RTF;
                cfftS= conj(fftS);
                nn = round(nF/10);
                
                for i=1:nS
                    % Get the cross PSD
                    crossPSD  = abs(bsxfun(@times,fftS(i,:),cfftS(i:end,:)));
                    % Estimate the noise plateau
                    out(i,i:end)  = mean(crossPSD(:,end-nn:end),2);
                end
                out = transpose(out) + out - diag(diag(out));
            end
        end
        
        function out = getOLslopes(s,u,MI,dt)
            s = squeeze(s);
            u = squeeze(u);
            out = s + MI*(dt.*circshift(u,-2,2) + (1-dt).*circshift(u,-1,2));
        end
        
        
       function varargout = getr0L0fromDMcommands(dmModes,waveFront,pupil,D,varargin)
           inputs = inputParser;
           inputs.addRequired('dmModes',@isnumeric);       
           inputs.addRequired('waveFront',@isnumeric );                        
           inputs.addRequired('pupil',@(x) isnumeric(x) || islogical(x) ); 
           inputs.addRequired('D',@isnumeric); 
           inputs.addParameter('nMin',4,@islogical);               
           inputs.addParameter('nMax',100,@islogical);  
           inputs.addParameter('mskPhase',true(1,size(dmModes,2)),@islogical);  
           inputs.addParameter('fitL0',false,@islogical);  
           inputs.addParameter('best',false,@islogical);  
           inputs.addParameter('covNoise',[],@isnumeric);  
           inputs.parse(dmModes,waveFront,pupil,D,varargin{:});
           
           nMin     = inputs.Results.nMin;
           nMax     = inputs.Results.nMax;
           mskPhase = inputs.Results.mskPhase;
           fitL0    = inputs.Results.fitL0;
           best     = inputs.Results.best;
           covNoise = inputs.Results.covNoise;
           
           % Derive the zernike reconstructor
           nRes   = sqrt(size(dmModes,1));%tel.resolution;
           zern   = zernike(nMin:nMax,'resolution',nRes);
           zernP  = zernike(nMin:nMax,'resolution',nRes,'pupil',pupil);
           % Project the modes onto the telescope pupil
           zModes = zern.modes;
           zpModes= zernP.modes;
           proj   = (zpModes' * zpModes)\zpModes'*zModes;
           % DM commands to zernike
           u2z    = proj*pinv(full(zModes))*dmModes(:,mskPhase);
           z      = u2z*waveFront(mskPhase,:);
           % Noise estimation
           if isempty(covNoise)
               covNoise = tools.getNoiseVariance(waveFront(mskPhase,:),'nshift',1);
           end
           %Zernike variance distribution
           c0   = std(z,[],2).^2 - diag(u2z*covNoise*u2z');
           
           
           if ~best
               [r0,L0,dr0,dL0] = tools.fitZernikeVariance(D,nMin,nMax,c0,fitL0);
               iBest           = 1;
               varargout{5}    = nMin;
           else
               nM   = nMin:nMax/4;
               iN   = length(nM);
               r0   = zeros(1,iN);
               L0   = zeros(1,iN);
               dr0  = zeros(1,iN);
               dL0  = zeros(1,iN);
               for i=1:iN
                   [r0(i),L0(i),dr0(i),dL0(i)] = tools.fitZernikeVariance(D,...
                       nM(i),nMax,c0(i:end),fitL0);
               end
               eps   = hypot(dr0./r0,dL0./L0);
               iBest = find(eps == min(eps));
               r0    = r0(iBest);
               L0    = L0(iBest);
               dr0   = dr0(iBest);
               dL0   = dL0(iBest);
               varargout{5} = nM(iBest);
           end
           
           varargout{1} = r0;
           varargout{2} = L0;
           varargout{3} = dr0;
           varargout{4} = dL0;
           
           if nargout > 5
               varargout{6} = covNoise;
           end
       end
       function varargout = fitZernikeVariance(D,nMin,nMax,c0,fitL0)
           %Fitting procedure
           FUN = @(x,xdata) (tools.zernikeVarianceModel(x,xdata));
           X0  = D./[0.2,50];
           ub  = D./[0.01,D];
           lb  = D./[1,100];
           
           if ~fitL0
               X0 = X0(1);
               ub = ub(1);
               lb = lb(1);
           end
                                 
           opt = optimoptions(@lsqcurvefit,'MaxIter',1e2,'TolFun',1e-10,...
               'TolX',1e-10,'MaxFunEvals',3e2);
           
           [beta,~,R,~,~,~,J] = lsqcurvefit(FUN,X0,nMin:nMax,c0,lb,ub,opt);
           
           %[beta,R,J] = nlinfit(nMin:nMax,c0,FUN,X0);
           varargout{1}  = abs(D/beta(1));
           if ~fitL0
               varargout{2}  = Inf;
           else
               varargout{2}  = abs(D/beta(2));
           end
           
           if isreal(beta)
               tmp        = diff(nlparci(beta,R,'jacobian',J),1,2);
               % Outputs
               varargout{3} = D*tmp(1)/beta(1)^2;
               if ~fitL0
                   varargout{4} = 0;
               else
                   varargout{4} = D*tmp(2)/beta(2)^2;
               end
           else
               varargout{3} = Inf;
               varargout{4} = Inf;
           end
       end
       function out = zernikeVarianceModel(x,xdata)
           %% ZERNIKEVARIANCE Zernike coefficients variance
           %
           % out = variance(modes,atmosphere) computes the
           % variance of Zernike coefficients from the modes and the
           % atmosphere object
           %
           % out = variance(zernike,atmosphere) computes the
           % variance of Zernike coefficients from the Zernike polynomials
           % object and the atmosphere object
           %
           % See also zernike, atmosphere
           
           % SAME FUNCITON AS PHASESTATS.ZERNIKEVARIANCE -> TO BE MERGED
           
           dr0     = abs(x(1));
           if length(x) == 1
               dL0  = 0;
           else
               dL0 = abs(x(2));
           end
           jv      = xdata;
           [nv,mv] = tools.nmOrder(jv);
           nv0     = nv;
           index   = diff(nv)~=0;
           jv      = [jv(index) jv(end)];
           mv      = [mv(index) mv(end)];
           nv      = [nv(index) nv(end)];
           nf      = length(nv);
           out     = zeros(length(jv),1);
           
           for cpt = 1:nf
               
               j = jv(cpt);
               n = nv(cpt);
               m = mv(cpt);
               
               out(nv0==n,1) = zernCovCoef(dr0,dL0,j,j,n,m,n,m);
               
           end
           function out = zernCovCoef(dr0,dL0,i,j,ni,mi,nj,mj)
               if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                   if dL0==0
                       if i==1 && j==1
                           out = Inf;
                       else
                           out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                               (dr0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                               newGamma(-5./6+(ni+nj)./2,...
                               [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                       end
                   else
                       out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                           (dr0./dL0).^(5./3)./dL0.^2.*...
                           sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                           UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*dL0);
                   end
               else
                   out = 0;
               end
               function out = newGamma(a,b)
                   % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                   % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                   % out = newGamma(a,b)
                   
                   out = prod(gamma(a))./prod(gamma(b));
               end
               function out = UnParamEx4q2(mu,alpha,beta,p,a)
                   % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                   % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                   % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                   % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                   % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                   
                   a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                   b1 = [1+alpha+beta 1+alpha 1+beta];
                   a2 = [(1-mu)./2+p 1+p p];
                   b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                   
                   out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                       newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                       pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                       newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                       pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                   function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                       % POCHAMMERSERIES Computes power series in Pochammer notation
                       % pochammerSeries(p,q,a,b,z)
                       % pochammerSeries(p,q,a,b,z,tol)
                       % pochammerSeries(p,q,a,b,z,[],nmax)
                       % pochammerSeries(p,q,a,b,z,tol,nmax)
                       
                       if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                           
                           if p==length(a) && q==length(b)
                               
                               switch nargin
                                   case 6
                                       nmax = 1e3;
                                   case 7
                                       if isempty(tol)
                                           tol = 1e-6;
                                       end
                                   otherwise
                                       tol = 1e-6;
                                       nmax = 1e3;
                               end
                               
                               out = zeros(size(z));
                               
                               indz = find(z==0);
                               if ~isempty(indz)
                                   out(indz) = 1;
                               end
                               
                               indnz = find(z~=0);
                               if ~isempty(indnz)
                                   z = z(indnz);
                                   ck = 1;
                                   step = Inf;
                                   k = 0;
                                   som = ck;
                                   while (k<=nmax) && (step>tol)
                                       ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                       step = abs(abs(ck)-abs(ckp1));
                                       som = som + ckp1;
                                       k = k+1;
                                       ck = ckp1;
                                   end
                                   if step>tol
                                       warning('pochammerSeries','Maximum iteration reached before convergence')
                                   end
                                   out(indnz) = som;
                               end
                               
                           else
                               error('p and q must be the same length than vectors a and b, respectively')
                               
                           end
                           
                       else
                           error('This generalized hypergeometric function doesn''t converge')
                       end
                   end
               end
               
           end
       end       
       function [n,m] = nmOrder(i)
           % [n,m] = nmOrder(i,)
           %Give the radial and azimutal order of the ith zernike
           n = floor((-1.+sqrt(8*(i-1)+1))/2);
           p = (i-(n.*(n+1))/2);
           k = mod(n,2);
           m = floor((p+k)/2)*2 - k;
       end
                                                                                     
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- IMAGES STATISTICS ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = getAreaFromPupil(pupil,D)
            ps  = D./size(pupil);
            out = sum(pupil(:))*ps(1)*ps(2);
        end
        
        function X = getFocalGrid(resolution,pixelScale)
            xx     = (-resolution/2:1:resolution/2-1)*pixelScale;
            [y,x]  = meshgrid(xx);
            X      = [x,y];
        end
        
        function [S,Srad] = getSurfaceBrightness(psf,zeroPoint,ps)
            %Get magnitude
            [Flux,bg] = tools.getFlux(psf);
            psf       = psf - bg;
            m         = -2.5*log10(Flux/zeroPoint);
            %Get surface brightness
            X         = tools.getFocalGrid(size(psf,2),ps);
            r         = hypot(X(:,1:end/2),X(:,end/2+1:end));
            A         = pi*r.^2;
            S         = -2.5*log10(psf/zeroPoint);
            Srad      = radial(S);
        end
        
        function [Flux,bg,ron,msk] = getFlux(psf)
            %Define the inner circle
            npx      = length(psf);
            x        = linspace(-1,1,npx);
            [X,Y]    = meshgrid(x);
            r        = hypot(X,Y);
            msk      = r>1;
            % Computing the residual background
            psfNoise = psf .*msk;
            bg       = median(psfNoise(:));
            % Computing the read-out noise
            ron      = std(psfNoise(:));
            %Computing the normalized flux
            psfFlux  = psf.*(r<=1);
            Flux     = sum(psfFlux(:) -bg);
        end
        
        function [FWHM,dFWHM,aRatio,theta,beta] = getFWHM(psf,pixelScale,rebin,method)
            
            % Gaussian and Moffat fitting are not really efficient on
            % anisoplanatic PSF. Prefer the coutour function in such a
            % case. The cutting method is not compliant to PSF not oriented
            % along x or y-axis.
            
            if nargin < 3
                rebin = 4;
            end            
            if nargin < 4
                method = 'contour';
            end
            %Interpolation   100*tools.getStrehl(camImg,pupMat,overSampling)         
            im2     = tools.interpolateOtf(psf,rebin*size(psf,1));
            if strcmp(method,'cutting')
                % Brutal approach when the PSF is centered and aligned
                % x-axis FWHM
                imx     = im2(:,floor(end/2+1));
                idx     = imx >= (max(imx(:))/2.);
                w       = find(idx==1);
                FWHMx   = (max(w) - min(w))/rebin*pixelScale;%sqrt(4.*sum(idx)/pi)*pixelScale/rebin;
                % y-axis FWHM
                imy     = im2(floor(end/2+1),:);
                idx     = imy >= (max(imy(:))/2.);
                w       = find(idx==1);
                FWHMy   = (max(w) - min(w))/rebin*pixelScale;%sqrt(4.*sum(idx)/pi)*pixelScale/rebin;                
                theta   = 0;               
            elseif strcmp(method,'contour')
                % Contour approach~: something wrong about the ellipse
                % orientation
                C       = contourc(im2,max(im2(:))*[0.5 0.5]);
                if ~isempty(C)
                    % centering the ellispe
                    mx      = [max(C(1,2:end)),max(C(2,2:end))];
                    mn      = [min(C(1,2:end)),min(C(2,2:end))];
                    cent    = (mx+mn)/2;
                    wx      = C(1,2:end) - cent(1);
                    wy      = C(2,2:end) - cent(2);
                    % Get the module
                    r       = hypot(wx,wy)/rebin*pixelScale;
                    % Getting the FWHM
                    FWHMx   = 2*max(r);
                    FWHMy   = 2*min(r);
                    % Getting the ellipse orientation
                    xm      = wx(r == max(r));
                    ym      = wy(r == max(r));
                    theta   = mean(abs(cart2pol(xm,ym)*180/pi));%mean(180*atan(ym./xm)/pi);
                    % Angle are counted positively in the reverse clockwise
                else
                    FWHMx = 0;
                    FWHMy = 0;
                end
                % direction.                
            elseif strcmp(method,'Gaussian')
                xdata   = tools.getFocalGrid(size(psf),pixelScale);
                x0      = [max(psf(:)),20,20,1,0,0,0];
                f       = @(x,xdata) tools.gaussianModel(x,xdata);
                beta    = lsqcurvefit(f,x0,xdata,psf);
                FWHMx   = beta(2)*2*sqrt(2*log(2));
                FWHMy   = beta(3)*2*sqrt(2*log(2));
                theta   = beta(4);
            elseif strcmp(method,'Moffat')
                xdata   = tools.getFocalGrid(size(psf),pixelScale);
                x0      = [max(psf(:)),20,20,1,0,0,0];
                f       = @(x,xdata) tools.moffatModel(x,xdata);
                beta    = lsqcurvefit(f,x0,xdata,psf);
                FWHMx   = 2*beta(2)*sqrt(2^(1./beta(4))-1);
                FWHMy   = 2*beta(3)*sqrt(2^(1./beta(4))-1);
                theta   = beta(5);
            end
            % Get Ellipticity
            FWHM   = max(FWHMx,FWHMy);
            dFWHM  = pixelScale/rebin/2;
            aRatio = max(FWHMx/FWHMy,FWHMy/FWHMx);
        end
        
        function [SR,dSR] = getStrehl(psf0,pupil,nyquistSampling)
            
            [~,~,ron] = tools.getFlux(psf0);
            psf     = tools.recenterPSF(psf0,4);            
            % Get the OTF
            otf     = tools.psf2otf(psf);
            otf     = otf/max(otf(:));
            notf    = size(otf,2);
            % Get the Diffraction-limit OTF
            otfDL   = tools.telescopeOtf(pupil,2*nyquistSampling);
            otfDL   = tools.interpolate(otfDL,notf,'spline');
            otfDL   = otfDL/max(otfDL(:));
            % Get the Strehl
            u       = -1+1/2/notf:2/notf:1-1/2/notf;
            SR      = abs((trapz(u,trapz(u,otf))/trapz(u,trapz(u,otfDL))));
            % Get the uncertainty
            dSR     = ron/sum(otfDL(:));
        end
        
        function ee = getEnsquaredEnergy(eeWidthInDiff,psf,D,nyquistSampling)
            
            % Get the otf and the strehl
            otf    = real(tools.psf2otf(psf));
            if nyquistSampling>=1
                otf    = tools.crop(otf,round(size(otf,1)/nyquistSampling));
            end
            otf    = otf/max(otf(:));
            
            % entrapped energy
            a        = eeWidthInDiff/D;%(wvl/D*constants.radian2arcsec)/D;
            nOtf     = length(otf);
            u        = linspace(-1,1,nOtf).*D;
            [x,y]    = meshgrid(u);
            eeFilter = a^2*sinc(x.*a).*sinc(y.*a);
            ee       = trapz(u,trapz(u,otf.*eeFilter));
            
        end
        
        function out = getFVU(xtrue,xest,nbox)
            if nargin > 2
                n   = length(xtrue);
                idx = floor(n/2+1-nbox/2):floor(n/2+nbox/2);
                xest = xest(idx,idx);
                xtrue= xtrue(idx,idx);
            end
            MSE = sum(sum((xest-xtrue).^2));
            VarX= sum(sum((xtrue - mean(xtrue(:))).^2));
            out = MSE/VarX;            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- MODEL FUNCTIONS ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
        function [out,imFit,beta] = findStellarParameters(im,psfModel,initParam)
            
            
            %1. Model
            nIm          = size(im);
            stellarModel = @(x,xdata) tools.stellarFieldModel(x(1:end-1),xdata,nIm) + x(end);
            normFactor   = sum(im(:));
            
            %2. Fitting options and initial guess
            nS    = (length(initParam)-1)/3;
            lX    = initParam(1:2*nS) - 2;
            uX    = initParam(1:2*nS) + 2;
            lF    = zeros(1,nS);
            uF    = 10*ones(1,nS);
            lb    = [lX,lF,-5*rms(im(:))];
            ub    = [uX,uF,5*rms(im(:))];
            
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',3e2,...
                'InitDamping',1,'Display','iter');
            
            %3. Fitting procedure
            [beta,~,R,~,~,~,J] = lsqcurvefit(stellarModel,initParam,psfModel,...
                im/normFactor,lb,ub,opt);
            
            %4. Unpacking results
            xS = beta(1:nS);
            yS = beta(nS+1:2*nS);
            fS = beta(2*nS+1:3*nS)*normFactor;
            
            %5. Measurements uncertainties
            dbeta = diff(nlparci(beta,R,'jacobian',J),1,2);           
            dX    = dbeta(1:nS);
            dY    = dbeta(nS+1:2*nS);
            dF    = dbeta(2*nS+1:3*nS)*normFactor;
            
            % Concatenating results            
            
            out = zeros(6,nS);
            for iS=1:nS
                out(1,iS) = xS(iS);
                out(2,iS) = yS(iS);
                out(3,iS) = fS(iS);
                out(4,iS) = dX(iS);
                out(5,iS) = dY(iS);
                out(6,iS) = dF(iS);
            end            
            
            % Fitted model
            imFit = stellarModel([xS,yS,fS,beta(end)],psfModel);
        end
        
        
        function out = stellarFieldModel(pStars,psfModel,nIm)
            
            nS = length(pStars)/3;
            xS = pStars(1:nS);
            yS = pStars(1+nS:2*nS);
            fS = pStars(1+2*nS:3*nS);
                
            out = 0*psfModel;
             for iS = 1:nS
                % Translate the PSF
                psf_i = tools.translateImage(psfModel,[xS(iS),yS(iS)]);
                % Update the image and Flux scaling:
                out = out + psf_i*fS(iS);                
             end
             out = tools.crop(out,nIm);
        end
        
                                        
      
        function out = multipleImage(x,im)
            %Grab locations
            n       = length(x)/3;
            [nx,ny] = size(im);
            xloc    = x(1:n);
            yloc    = x(n+1:2*n);
            flux    = x(2*n+1:end);
            
            if any(abs(xloc) > nx) || any(abs(yloc) > ny)
                id = abs(xloc) > nx | abs(yloc) > ny;
                nn = sum(id);
                if nn == 1
                    fprintf('Warning: 1 PSF are out the field\n');
                else
                    fprintf('Warning: %d PSFs are out the field\n',nn);
                end
                flux(id) = 0;
            end
            
            % For xloc=yloc=0 and flux=1, the procedure does modify the PSF
            % while it mustn't. To be check.
            %nx      = 2*nx;
            %ny      = 2*ny;
            out     = zeros(nx,ny);
            otf     = tools.psf2otf(im);
            %otf     = tools.interpolateOtf(otf,nx);
            otf     = otf/max(otf(:));
            [u,v]   = freqspace([nx ny],'meshgrid');
            for i=1:n
                yi    = yloc(i);
                xi    = xloc(i);
                fftPhasor = exp(-1i*pi*(u*xi + v*yi));
                map_i = tools.otf2psf(otf.*fftPhasor);
                % Adding images
                out = out + flux(i)*map_i/sum(map_i(:));
            end
        end
        
        
        function out = gaussianModel(x,xdata)          
            % ------- Grabbing parameters ---------%
            I0 = x(1);          %Amplitude
            ax = x(2);          %x spreading
            ay = x(3);          %y-spreading
            th = x(4)*pi/180.;  %rotation
            x0 = x(5);          %x-shift
            y0 = x(6);          %y-shift
            
            % ------- Including shifts ---------
            n     = size(xdata,2);
            X     = xdata(:,1:n/2);
            Y     = xdata(:,n/2+1:end);
            %Shifts
            X     = X - x0;
            Y     = Y - y0;
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Gaussian expression
            out = I0.*exp(-0.5*((Xr./ax).^2 + (Yr./ay).^2) );           
        end
        
        function out = multiGaussian(x,xdata)
            
            nG  = length(x)/6;
            out = zeros(size(xdata,1));
            for i=1:nG
                out = out + tools.gaussianModel(x(1+(i-1)*6:i*6),xdata);
            end
        end
        
        function out = moffatModel(x,xdata)
            
            % ------- Grabbing parameters ---------%
            I0 = x(1);          %Amplitude
            ax = x(2);          %x spreading
            ay = x(3);          %y-spreading
            be = x(4);          %center slope
            th = x(5);          %rotation
            x0 = x(6);          %x-shift
            y0 = x(7);          %y-shift
            
            % ------- Including shifts ---------
            n     = size(xdata,2);
            X     = xdata(:,1:n/2);
            Y     = xdata(:,n/2+1:end);
            %Shifts
            X     = X - x0;
            Y     = Y - y0;
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Moffat expression
            out = I0.*(1. + (Xr./ax).^2 + (Yr./ay).^2).^(-be);
            %Normalization
            out = out*(be-1)/pi/ax/ay;
        end
       
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- ANISOPLANATISM------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [out,ATF] = anisoPsfToPsfFilter(im0,tel,atm,gs,src,nqSampling,npt,flux,varargin)
            inputs = inputParser;
            inputs.addRequired('im0',@isnumeric);
            inputs.addRequired('tel',@(x) isa(x,'telescope'));
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('gs',@(x) isa(x,'source'));
            inputs.addRequired('src',@(x) isa(x,'source'));
            inputs.addRequired('nqSampling',@isnumeric);
            inputs.addRequired('npt',@isnumeric);
            inputs.addRequired('flux',@isnumeric);
            inputs.addParameter('filter','none',@ischar);
            inputs.addParameter('otfFlag',false,@islogical);
            inputs.addParameter('flagPupilAverage',false,@islogical);
            inputs.parse(im0,tel,atm,gs,src,nqSampling,npt,flux,varargin{:});
            
            filter  = inputs.Results.filter;
            otfFlag = inputs.Results.otfFlag;
            flagP   = inputs.Results.flagPupilAverage;
            if otfFlag
                otf0 = im0;
                npsf = size(im0,2);
                if nqSampling > 1
                    nc = 2*round(npsf/nqSampling/2);
                elseif nqSampling ==1
                    nc = npsf;
                else
                    nc = round(npsf/nqSampling); %TBD
                end
                    
            else
                if ~any(im0(:))
                    % Define the diffraction-limited PSF;
                    P    = tools.enlargePupil(tel.pupil,2*nqSampling);
                    im0  = tools.pupil2psf(P);
                    im0  = tools.crop(im0,size(im0,1)/2/nqSampling);
                end
                npsf  = size(im0,2);
                F0  = sum(im0(:));
                im0  = im0/F0;
                %On-axis OTF
                otf0  = tools.psf2otf(im0);
                if nqSampling > 1
                    nc = 2*round(npsf/nqSampling/2);
                    otf0  = tools.crop(otf0,nc);
                    otf0  = otf0/max(otf0(:));
                elseif nqSampling ==1
                    nc = npsf;
                    otf0  = otf0/max(otf0(:));
                else
                    nc = round(npsf/nqSampling); %TBD
                end
                
            end
            
            % ---- Computing ATFs
            
            if numel(gs) > 1
                ATFho = tools.anisoplatismFilter(tel,atm,gs(1),src,npt,nc,'filter','ttout','flagPupilAverage',flagP);
                ATFtt = tools.anisoplatismFilter(tel,atm,gs(2),src,npt,nc,'filter','ttonly','flagPupilAverage',flagP);
                ATF   = ATFtt.*ATFho;
            else
                ATF = tools.anisoplatismFilter(tel,atm,gs,src,npt,nc,'filter',filter,'flagPupilAverage',flagP);
            end
            
            %Applying ATF
            nSrc = numel(src);
            out  = zeros(npsf,npsf*nSrc);
            for iSrc=1:nSrc
                idx   = 1 + (iSrc-1)*npsf:iSrc*npsf;
                if nqSampling > 1
                    tmp   = padarray(otf0.*ATF(:,:,iSrc),[(npsf-nc)/2 (npsf-nc)/2],0,'both');
                elseif nqSampling ==1
                    tmp   = otf0.*ATF;
                else
                    tmp  = otf0.*ATF; %TBD
                end
                out(:,idx) = tools.otf2psf(tmp).*flux(iSrc);
            end
        end
        
        function [ATF,Cani] = anisoplatismFilter(tel,atm,gs,src,npt,npsf,varargin)
            
            inputs = inputParser;
            inputs.addRequired('tel',@(x) isa(x,'telescope'));
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addRequired('gs',@(x) isa(x,'source'));
            inputs.addRequired('src',@(x) isa(x,'source'));
            inputs.addRequired('npt',@isnumeric);
            inputs.addRequired('npsf',@isnumeric);
            inputs.addParameter('filter','none',@ischar);
            inputs.addParameter('iHdm',1,@isnumeric);
            inputs.addParameter('flagPupilAverage',false,@islogical);
            inputs.parse(tel,atm,gs,src,npt,npsf,varargin{:});
            
            filter           = inputs.Results.filter;
            iHdm             = inputs.Results.iHdm;
            flagPupilAverage = inputs.Results.flagPupilAverage;
            
            
            ATF  = ones(npsf,npsf);
            Cani = zeros(npt^2,npt^2);
            go   = sum(sum(abs(bsxfun(@minus,[src.directionVector],gs.directionVector)))) +  any(bsxfun(@minus,[src.height],gs.height));
            
            if go
                % Defining a low resolution telescope for computation
                tel            = tools.duplicateTelescope(tel);
                tel.resolution = npt;
                tel.pupil      = ones(npt);
                % Instantiating things
                npt2 = sum(tel.pupil(:));
                d    = tel.D/(npt-1);
                msk  = logical(tel.pupil);
                nSrc = numel(src);
                ATF  = ones(npsf,npsf,nSrc);
                Cani = zeros(npt2,npt2,nSrc);
                % modified the atmosphere with respect to the src
                % wavelength and telescope elevation
                wvl0               = atm.wavelength;
                atm.wavelength     = src(1).wavelength;
                airmass            = 1/cos(tel.elevation*pi/180);
                atm.r0             = atm.r0/airmass^(3/5);
                for i=1:atm.nLayer
                    atm.layer(i).altitude = atm.layer(i).altitude*airmass;
                end
                
                gs.height  = gs.height*airmass;
                for iSrc = 1:numel(src)
                    src(iSrc).height = src(iSrc).height*airmass;
                end
                % Diffraction-limited OTF
                otfDL = tools.zonalCovarianceToOtf(zeros(npt2),npsf,tel.D,d,msk);
                mskOtf= otfDL > 1e-3;
                
                % Tip-tilt filter if any
                F = 1;
                if strcmp(filter,'ttonly') || strcmp(filter,'ttout')
                    
                    %[X,Y]  = meshgrid((1:npt),(1:npt));
                    %TT     = [X(:),Y(:)];
                    %F      = TT*pinv(TT);
                    
                    x      = (-1+1/npt:2/npt:1-1/npt);
                    [X,Y]  = meshgrid(x,x);
                    TT     = [X(:),Y(:)];
                    F      = TT*pinv(TT)/sqrt(sqrt(2));
                    if strcmp(filter,'ttout')
                        zern   = zernike(2:3,'resolution',tel.resolution, 'pupil',tel.pupil);
                        TT     = zern.modes(tel.pupilLogical,:);
                        F  = eye(npt2) -  TT*pinv(TT);
                    end
                end
                % Grabbing covariance matrices
                Cgg = phaseStats.spatioAngularCovarianceMatrix(npt,tel.D,atm,gs,'mask',msk);
                
                for iSrc=1:nSrc
                    if any(src(iSrc).directionVector ~= gs.directionVector) || gs.height ~= src(iSrc).height
                        
                        Css = phaseStats.spatioAngularCovarianceMatrix(npt,tel.D,atm,src(iSrc),'mask',msk);
                        Csg = phaseStats.spatioAngularCovarianceMatrix(npt,tel.D,atm,src(iSrc),'srcCC',gs,'mask',msk);
                        
                        % Anisoplanatic matrix
                        Cani(:,:,iSrc)= F*(Css + Cgg - Csg{1} - Csg{1}')*F';
                        
                        % DM Filtering
                        Cani(:,:,iSrc)= iHdm*Cani(:,:,iSrc)*iHdm';
                        
                        %Get OTF
                        if flagPupilAverage
                            covMap        = utilities.covMatrix2Map(Cani(:,:,iSrc),npt,npt);
                            tmp          = exp(covMap - covMap(floor(end/2+1),floor(end/2+1)));
                            ATF(:,:,iSrc)= tools.interpolateOtf(tmp,npsf);
                        else
                            tmp           = tools.zonalCovarianceToOtf(Cani(:,:,iSrc),npsf,tel.D,d,msk);
                            tmp(mskOtf)   = tmp(mskOtf)./otfDL(mskOtf);
                            ATF(:,:,iSrc) = tmp/max(tmp(:));
                        end
                    end
                end
                
                
                % Back to the initial atmosphere
                atm.r0             = atm.r0*airmass^(3/5);
                for i=1:atm.nLayer
                    atm.layer(i).altitude = atm.layer(i).altitude/airmass;
                end
                atm.wavelength     = wvl0;                
            end
        end
        
         
        function out = mcDonald(x)
            out = x.^(5/6.).*besselk(5./6,x)./(2^(5/6.)*gamma(11/6.)) ;
            out(find(x == 0)) = 3/5.;
        end
              
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- PSF-R ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function phiPara = dmFilter(phiAtmo,dmSq,iF)
            %Computes coefs
            dmSq.coefs = iF*phiAtmo(:);
            phiPara  = dmSq.surface;
        end
               
         function [otf,dphi] = modes2Otf(Cmm,modes,pupil,npsf,overSampling,varargin)
            inputs = inputParser;
            inputs.addRequired('Cmm',@isnumeric);
            inputs.addRequired('modes',@isnumeric);
            inputs.addRequired('pupil',@isnumeric);
            inputs.addRequired('npsf',@isnumeric);
            inputs.addRequired('overSampling',@isnumeric);
            inputs.addParameter('method','Vii',@ischar);
            inputs.parse(Cmm,modes,pupil,npsf,overSampling,varargin{:});
           
            method       = inputs.Results.method;
            
            % Autocorrelation of the pupil expressed in pupil
            nPx        = sqrt(size(modes,1));
            pupExtended= tools.enlargePupil(double(pupil),2*overSampling);
            fftPup     = fft2(pupExtended);
            conjPupFft = conj(fftPup);
            G          = fftshift(real(fft2(fftPup.*conjPupFft)));
            % Defining the inverse
            den        = zeros(size(G));
            msk        = G./max(G(:)) > 1e-7;
            den(msk)   = 1./G(msk);
            
            
            if any(Cmm(:)) && strcmp(method,'Vii')
                % Diagonalizing the Cvv matrix
                [U,S]   = svd(Cmm);
                s       = diag(S);
                nModes  = length(s);
                M       = modes * U;
                %loop on actuators
                
                buf = zeros(size(pupExtended));
                
                for k=1:nModes
                    Mk   = reshape(M(:,k),nPx,nPx);
                    Mk   = tools.enlargePupil(Mk,2*overSampling);
                    % Vii computation
                    Vk   = real(fft2(Mk.^2.*pupExtended).*conjPupFft) - abs(fft2(Mk.*pupExtended)).^2;
                    % Summing modes into dm basis
                    buf  = buf + s(k) .* Vk;
                end
                
                dphi     = den.*fftshift(real(fft2(2.*buf)));
                otf      = G.*exp(-0.5*dphi);
                
                
            elseif any(Cmm(:)) && strcmp(method,'Uij')
                nm   = size(modes,2);
                dphi = 0*pupExtended;
                
                %Double loops on modes
                for i=1:nm
                    Mi = reshape(modes(:,i),nPx,nPx);
                    Mi = tools.enlargePupil(Mi,2*overSampling);
                    for j=1:i
                        %Getting modes + interpolation to twice resolution
                        Mj    = reshape(modes(:,j),nPx,nPx);
                        Mj   = tools.enlargePupil(Mj,2*overSampling);
                        term1 = real(fft2(Mi.*Mj.*pupExtended).*conjPupFft);
                        term2 = real(fft2(Mi.*pupExtended).*conj(fft2(Mj.*pupExtended)));
                        % Uij computation
                        Uij   = real(ifft2(term1-term2));
                        %Summing terms
                        fact = double(i~=j) + 1;
                        dphi = dphi + fact*Cmm(i,j).*Uij;
                    end
                end
                dphi = fftshift(2*dphi).*den.*msk;
                otf  = G.*exp(-0.5*dphi);
            elseif any(Cmm(:)) && strcmp(method,'singleIF')
                % NOT TESTED !!!
                %Get the covariance map
                N      = size(pupil,1);
                Cvvmap = tools.covMatrix2Map(Cmm,sqrt(N),sqrt(N));
                ncmap  = N/2+1;                
                %defining the new map in the real domain
                start = int(ncmap - D/N);
                stop  = int(ncmap + D/N);
                step  = d;
                rr    = start:step:stop;                
                cv        = zeros(N);
                cv(rr,rr) = Cvvmap;                
                %TF of the influence function autocorrelation                
                A   = abs(fft2(modes)).^2;                
                %computing the correlation function
                C  = real(fft2( A * fft2(fftshift(cv)) )) / length(fi);
                %getting the phase structure function
                dphi  = 2*(max(C) - C);                               
                %computing the OTF
                otf  = G.*exp(-0.5*dphi);               
            else
                % Diffraction-limit case
                G    = G./max(G(:));
                otf  = G;
                dphi = 0*G;
            end
            
            % Interpolation of the OTF => determination of the PSF fov
            otf = otf.*(G>1e-5);
            otf = tools.interpolateOtf(otf,npsf);
            otf = otf./max(otf(:));
            dphi= tools.interpolateOtf(dphi,npsf);
        end
        
        function otf = zonalCovarianceToOtf(Cphi,npsf,D,dp,idxValid)
            
            % Grabbing the valid actuators positions in meters
            loc  = tools.pointWiseLocation(D,dp,idxValid);
            %OTF sampling
            nPh  = round(D/dp+1);
            %nU1d = 2*nPh;
            u1D  = (-nPh:1:nPh-1)*dp;
            nU1d = length(u1D);
            % Determining couples of point with the same separation
            [shiftX,shiftY] = tools.mkotf_indpts(nU1d,nPh,u1D,loc,dp);
            % WF Amplitude
            amp0 = ones(nnz(idxValid),1);
            
            %% Long-exposure OTF
            otf = tools.mkotf(shiftX,shiftY,nU1d,amp0,dp,-0.5*Cphi);
            %Interpolation
            otf      = tools.interpolateOtf(otf,npsf);
            otf      = otf./max(otf(:));
        end
        
        function out = pointWiseLocation(D,dp,idxValid)
            % Defining the point-wise locations
            xloc                 = -D/2:dp:D/2;
            [actuLocX, actuLocY] = meshgrid(xloc);
            actuLocX             = actuLocX(idxValid);
            actuLocY             = actuLocY(idxValid);
            out                  = [actuLocX(:), actuLocY(:)];
        end
        
        function otf = mkotf(indptsc,indptsc2,nU1d,ampl,dp,C_phi)
            
            %Instantiation
            otf         = zeros(nU1d);
            C_phip_diag = exp(diag(C_phi));
            C_phipi     = exp(-2*C_phi);
            C_phi_size2 = size(C_phi,2);
            
            for iu=1:nU1d
                for ju=1:nU1d
                    indpts  = indptsc{iu, ju};
                    indpts2 = indptsc2{iu, ju};
                    
                    if isempty(indpts)
                        otf(iu,ju) = 0;
                    else
                        myarg      = C_phip_diag(indpts2).*C_phip_diag(indpts)...
                            .*C_phipi(C_phi_size2*(indpts-1) + indpts2);
                        kernel     = (conj(ampl(indpts2)) .* ampl(indpts))' * myarg;
                        otf(iu,ju) = kernel*dp^2;
                    end % if isempty(mypts)
                end %ju
            end %iu
            
            dc  = sum(abs(ampl).^2) * dp^2;
            otf = otf/dc;
        end
        
        function [indptsc,indptsc2] = mkotf_indpts(nU1d,nPh,u1D,loc,dp)
            
            % index pts in a 3x bigger array
            locInPitch = loc/dp;
            minLoc     = min(locInPitch);
            ninloc     = size(loc,1);
            nc         = 3*nPh;
            n          = nc-nPh;
            n1         = (n-mod(n,2))/2 + 1 + mod(nPh,2);
            minLoc2    = minLoc -(n1-1);
            loc2       = round(locInPitch - ones(ninloc,1)*minLoc2+1);
            %embed the loc2 inside ncxnc array.
            indx_emb       = loc2(:,1)+(loc2(:,2)-1)*nc;
            mask           = zeros(nc);
            mask(indx_emb) = 1;
            indptsc        = cell(nU1d);
            indptsc2       = cell(nU1d);
            
            for iu=1:nU1d
                for ju=1:nU1d
                    u2D        = [u1D(iu)  u1D(ju)];
                    uInPitch   = u2D/dp;
                    iniloc_sh2 = loc2 + ones(ninloc,1)*round(uInPitch);
                    indxsh     = iniloc_sh2(:,1) + (iniloc_sh2(:,2)-1)*nc;
                    
                    %index of points in iniloc_sh2 that are intersect with iniloc_sh
                    %indpts is the overlapping points in the shifted array
                    indptsc{ju,iu}  = find(mask(indxsh));
                    mask2           = zeros(nc);
                    mask2(indxsh)   = 1;
                    indptsc2{ju,iu} = find(mask2(indx_emb));
                end
            end
        end
        
       
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- DISPLAY TOOLS ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function displayContours(psfClass,nRes)
            
            % Get meta-image and contours
            [im,x2]      = tools.psfGridAssembly(psfClass,nRes);
            res          = 10;
            [SRim,Xo,Yo] = tools.getContours(psfClass,'Strehl',res);
            FWHMim       = tools.getContours(psfClass,'FWHM',res);
            eim          = tools.getContours(psfClass,'Ellipticity',res);
            
            % ---------------- PSF grid ---------------------
            figure
            imagesc(x2,x2,im);
            xlabel('PSF position in x (arcsec)','interpreter','latex');
            ylabel('PSF position in y (arcsec)','interpreter','latex');
            set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            
            % ---------------- PSF properties contours ---------------------
            nC = 8;
            figure
            contour(Xo,Yo,SRim',nC,'ShowText','on','fill','on');
            xlabel('PSF position in x (arcsec)','interpreter','latex');
            ylabel('PSF position in y (arcsec)','interpreter','latex');
            title('Strehl-ratio contour','interpreter','latex');
            set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            
            figure
            contour(Xo,Yo,FWHMim',nC,'ShowText','on','fill','on');
            xlabel('PSF position in x (arcsec)','interpreter','latex');
            ylabel('PSF position in y (arcsec)','interpreter','latex');
            title('FWHM contour','interpreter','latex');
            set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            
            if any(eim-1)
                figure
                contour(Xo,Yo,eim,nC,'ShowText','on','fill','on');
                xlabel('PSF position in x (arcsec)','interpreter','latex');
                ylabel('PSF position in y (arcsec)','interpreter','latex');
                title('Ellipticity contour','interpreter','latex');
                set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            end
            % ---------------- PSF properties plots ---------------------
            r        = hypot([psfClass.xSrc],[psfClass.ySrc]);
            [rr,idx] = sort(r);
            SR       = [psfClass.Strehl];
            FWHM     = [psfClass.FWHM];
            E        = [psfClass.Ellipticity];
            
            figure
            plot(rr,SR(idx),'bo','MarkerSize',5,'MarkerFaceColor','b');
            xlabel('Distance to on-axis (arcsec)','interpreter','latex');
            ylabel('Strehl ratio (\%)','interpreter','latex');
            set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            
            figure
            plot(rr,FWHM(idx),'bo','MarkerSize',5,'MarkerFaceColor','b');
            ylabel('FWHM (mas)','interpreter','latex');
            xlabel('Distance to on-axis (arcsec)','interpreter','latex');
            set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            
            if any(eim-1)
                figure
                plot(rr,E(idx),'bo','MarkerSize',5,'MarkerFaceColor','b');
                ylabel('Aspect ratio','interpreter','latex');
                xlabel('Distance to on-axis (arcsec)','interpreter','latex');
                set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
            end
        end
        
        function [out,Xo,Yo] = getContours(psf,properties,res)
            
            if strcmp(properties,'Strehl')
                val = [psf.Strehl];
            elseif strcmp(properties,'FWHM')
                val = [psf.FWHM];
            elseif strcmp(properties,'Ellipticity')
                val = [psf.Ellipticity];
            elseif strcmp(properties,'FluxRatio')
                val = [psf.FluxInFWHM];
            end
            
            %% Geometry
            x    = [psf.xSrc];
            y    = [psf.ySrc];
            nimx = size(find(floor(diff(sort(x)))~=0),2) + 1;
            nimy = size(find(floor(diff(sort(y)))~=0),2) + 1;
            Nmxx = max(x(:));
            Nmnx = min(x(:));
            Nmxy = max(y(:));
            Nmny = min(y(:));
            
            %% INTERPOLATION
            xi       = linspace(Nmnx,Nmxx,nimx);
            yi       = linspace(Nmny,Nmxy,nimy);
            [Xi,Yi]  = meshgrid(xi, yi);
            %Interpolation for x-axis
            xo = [];
            if nimx>1
                for i=1:nimx-1
                    tmp = linspace(xi(i),xi(i+1),res);
                    if i==1
                        xo = tmp;
                    else
                        xo = [xo tmp(2:end)];
                    end
                end
            else
                xo = xi(1).*ones(res,1)';
            end
            %Interpolation for y-axis
            yo = [];
            if nimy>1
                for i=1:nimy-1
                    tmp = linspace(yi(i),yi(i+1),res);
                    if i==1
                        yo = tmp;
                    else
                        yo = [yo tmp(2:end)];
                    end
                end
            else
                yo = yi(1).*ones(res,1)';
            end
            
            [Xo,Yo]  = meshgrid(xo, yo);
            
            %% Reshaping
            valim     = reshape(val,nimx,nimy)';
            if nimx>1 && nimy>1
                out    = interp2(Xi,Yi,valim,Xo,Yo);
            elseif nimx>1 && nimy==1
                out    = interp1(xi,valim,xo);
                out    = reshape(repelem(out,res),res,res);
            elseif nimy>1 && nimx==1
                out    = interp1(yi,valim,yo);s
                out    = reshape(repelem(out',res),res,res);
            elseif nimx==1 && nimy == 1
                out = valim;
            end% Instantiate model for non-windowed case

        end
        
        function [m_psf,x2] = psfGridAssembly(psf,nRes,varargin)
            
            
            inputs = inputParser;
            inputs.addParameter('x',0,@isnumeric);
            inputs.addParameter('y',0,@isnumeric);
            inputs.addParameter('psInMas',10,@isnumeric);
            inputs.parse(varargin{:});
            x       = inputs.Results.x;
            y       = inputs.Results.y;
            ps      = inputs.Results.psInMas/1e3;
            
            if isa(psf,'psfStats')
                x   = [psf.xSrc]; %in arcsec
                y   = [psf.ySrc];
                ps  = psf(1).pixelScaleInMas/1e3;
                psf = [psf.image];
            end
            
            
            %1. Reshaping
            nP0   = size(psf,1);
            if length(size(psf)) ==3
                nIm   = size(psf,3);
            else
                nIm   = size(psf,2)/nP0;
                psf   = reshape(psf,nP0,nP0,nIm);
            end
            %2. Crop if required
            if nRes < nP0
                xi  = round((nP0 - nRes)/2 + 1);
                xf  = round((nP0 + nRes)/2);
                psf = psf(xi:xf,xi:xf,:);
                nP0 = nRes;
            end
            %3. Define the meta image
            if length(x) >1
                x = x(:)';
            end
            if length(y) > 1
                y = y(:)';
            end
            nimx  = size(find(floor(diff(sort(x)))~=0),2)+1;
            nimy  = size(find(floor(diff(sort(y)))~=0),2)+1;
            fov   = 2*max([max(abs(x(:))) max(abs(y(:)))]);
            %nPall = ceil(fov/ps) + nP0;%
            nPall = max(nimx,nimy)*nP0;
            ps    = fov/(nPall - nP0);
            m_psf = zeros(nPall);
            % on-axis reference
            nRefX = nPall/2+1;
            nRefY = nPall/2+1;
            
            %4. Fill up the meta-image
            for k=1:nIm
                % Grab indexes
                x0  = x(k)/ps + nRefX;
                y0  = y(k)/ps + nRefY;
                idx = round(x0-nP0/2+1:x0+nP0/2);
                idy = round(y0-nP0/2+1:y0+nP0/2);
                % Fill the meta-image
                m_psf(idy,idx) = psf(:,:,k);
            end
            
            %5. Compute the new geometric grid in physical units
            x2 = linspace(-1,1,nPall)*ps*nPall/2;
        end
        
        function makeAxisSquared(n)
            for i=1:numel(n)
                axesHandles = findobj(get(figure(n(i)),'Children'), 'flat','Type','axes');
                % Set the axis property to square
                axis(axesHandles,'square');
            end
        end
                 
      
    end
end


