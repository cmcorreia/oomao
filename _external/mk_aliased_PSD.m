function [W,Sx, Sy] = mk_aliased_PSD(f, m, d, r0, L0, pistRem, T, td, vlx, vly)
%Scope: Create aliased PSD

%   INPUT PARAMETERS
% f                          ::    a vector with central frequencies
% m                          ::    factor of max freq to be aliased
% d                          ::    sampling period, max(f) = 1/2/d
% r0, L0                     ::    Fried parameter and outer scale of atmosphere
% pistRem                   ::    flag to remove or not piston; if
%                                   pistRem > 0 then telescopeDiameter =
%                                   pistRem
% vlx, vly                   ::    The wind velocity in x and y
% T                          ::    The WFS integration time
% td                         ::    pure delay if any

%
%   OUTPUT PARAMETERS
% W                          ::    the phase aliased PSD
% WSx, WSy                   ::    SH gradients aliased PSD
%
%   CREATION                Carlos Correia, March 2014
%
%   RECORD CHANGE
%
if ~exist('vlx','var')
    vlx = 0;
end
if ~exist('vly','var')
    vly = 0;
end
if ~exist('T','var')
    T = 0;
end
if ~exist('td','var')
    td = 0;
end
if ~exist('pistRem','var')
    flagRemovePiston = 1;
else
    flagRemovePiston = pistRem;
end
D = pistRem;
if length(r0) > 1
    fractionnalR0 = r0(2:end);
    r0 = r0(1);
else
    fractionnalR0 = 1;
end
nLayer = length(fractionnalR0);
if isvector(f)
    [fx, fy] = meshgrid(f,f);
else
    fx = f;
    fy = f';
    f = fx(1,:);
end
f0 = 1/L0;
W = 0;
Sx = 0;
Sy = 0;


for mi = -m:m
    for ni = -m:m
        if 1 %mi~=0 || ni~=0
            fm = fx - mi/d;
            fn = fy - ni/d;
            Int = 0;
            for kLayer = 1 : nLayer
                Int             = Int + fractionnalR0(kLayer)*exp(1i*pi*(fm*vlx(kLayer)*td+fn*vly(kLayer)*td)).*...
                    sinc(vlx(kLayer)*T*fm).*sinc(vly(kLayer)*T*fn);
            end
            %Int = 1;
            Av = sinc(d*fm).*sinc(d*fn).*exp(1i*pi*d*(fm+fn)).*Int;
            Avsq = abs(Av).^2;
            
            if flagRemovePiston
                % --- piston removal ---
                %D              = 8;
                besselargx     = pi*D/d*f * d; % here I should not divide by dx because the frequency vector is already in propor units
                [BX, BY] = meshgrid(besselargx,besselargx);
                fltr   = 2*besselj(1,sqrt(BX.^2+BY.^2))./(sqrt(BX.^2+BY.^2));
                fltr(isnan(fltr)) = 1;
                PR = 1-abs(fltr).^2;
                
                W_mn = (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6).*PR;
            else
                W_mn = (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6);
            end
            W = W + W_mn;

            if nargout > 1
                Sx = Sx + W_mn.*abs(1i*2*pi*d*fm).^2.*Avsq;
                
                Sy = Sy + W_mn.*abs(1i*2*pi*d*fn).^2.*Avsq;
            end
            
        end
    end
end

W = W*0.0229*r0^(-5/3);
if nargout > 1
    Sx = Sx*0.0229*r0^(-5/3);
    Sy = Sy*0.0229*r0^(-5/3);
end
