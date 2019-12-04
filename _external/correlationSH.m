function [dX, dY] = correlationSH(image, f_ref, algo)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    algo = 1;
end

dxdy = size(image);
if sum(image(:)) == 0 % unilluminated subap (before SH initialization)
    dX = 0;
    dY = 0;
else
    f_image = fftshift(fft2(ifftshift(image)));
    corrFunc = abs2(fftshift(ifft2(ifftshift(f_image.*conj(f_ref)))));% centered around [N/2, N/2]
    if algo == 2
        % Poyneer interpolation
        [yPix, xPix]= ind2sub(size(image), find(corrFunc == max(corrFunc(:))));
        if (xPix<2) || (yPix<2) || (xPix>dxdy(1)-1) || (yPix>dxdy(1)-1)
            dX=0;
            dY=0;
        else
            dX = 0.5*(corrFunc(yPix, xPix-1) - corrFunc(yPix, xPix+1)) / (corrFunc(yPix, xPix+1) + corrFunc(yPix, xPix-1) -2* corrFunc(yPix, xPix)) + xPix - (dxdy(1)/2+1);
            dY = 0.5*(corrFunc(yPix-1, xPix) - corrFunc(yPix+1, xPix)) / (corrFunc(yPix+1, xPix) + corrFunc(yPix-1, xPix) - 2 * corrFunc(yPix, xPix)) + yPix - (dxdy(1)/2+1);
        end
    else
        % corrFunc = corrFunc - 0.0;
        corrFunc(corrFunc<0) = 0;
        [X Y] = cog(corrFunc);
        dX = X - (dxdy(1)/2+1);
        dY = Y - (dxdy(1)/2+1);
    end
end

end