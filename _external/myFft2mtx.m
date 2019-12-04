function out = myFft2mtx(in,outPx, sampling,dir)
if ~exist('dir','var')
    dir = 1;
end
if ~exist('c','var')
    c = 0.5;
end
cst = 2*1i*pi;

%sampling = outPx/sampling;
inPx = size(in,1);
x = ((0:inPx-1)+c - inPx/2)/inPx;
xo = ((0:outPx-1)+c - outPx/2);
if numel(sampling) > 1
    if isa(in,'gpuArray')
        e1 = gpuArray(zeros(outPx, inPx, numel(sampling)));
        e2 = gpuArray(zeros(inPx, outPx, numel(sampling)));
    else
        e1 = zeros(outPx, inPx, numel(sampling));
        e2 = zeros(inPx, outPx, numel(sampling));
    end
    for k = 1:numel(sampling)
        e1(:,:,k) = exp(dir*cst*x'*xo/sampling(k));
        e2(:,:,k) = exp(dir*cst*xo'/sampling(k)*x);
    end
else
    if isa(in,'gpuArray')
        e1 = gpuArray(exp(dir*cst*x'*xo/sampling));
        e2 = gpuArray(exp(dir*cst*xo'/sampling*x));
    else
        e1 = exp(dir*cst*x'*xo/sampling);
        e2 = exp(dir*cst*xo'/sampling*x);
    end
end

%if inPx == outPx
%    fun = @(A,B) A*B;
%    out = bsxfun(fun, bsxfun(fun, e2, in), e1);
%else
if numel(sampling) > 1
    if isa(in,'gpuArray')
        out = gpuArray(zeros(outPx,outPx,numel(sampling)));
    end
    for k = 1:numel(sampling)
        out(:,:,k) = fftshift(e2(:,:,k)*in(:,:,k)*e1(:,:,k));
    end
else
    out = fftshift(e2*in*e1);
end
%end

% ; Grilles de coordonnï¿½es dans l'espace rï¿½el
% X = transpose(k - Na/2)/Na
% Y = X
% 
% ; Grilles de coordonnï¿½es dans l'espace de Fourier
% U = transpose(l - Nb/2)*param/Nb
% V = U
% 
% IF keyword_set(idl_sign) THEN BEGIN
%    IF keyword_set(inverse) THEN dir = 1 ELSE dir = -1
% ENDIF ELSE BEGIN
%    IF keyword_set(inverse) THEN dir = -1 ELSE dir = 1
% ENDELSE
% 
% 
% 
% ; Calcul de la transformï¿½e de Fourier. La normalisation est telle que total(|TF|^2) = 1
% TF = (param/(Na*Nb))*exp(dir*2.*ii*pi*U##transpose(X))##$
%       image##exp(dir*2.*ii*pi*Y##transpose(V))
   