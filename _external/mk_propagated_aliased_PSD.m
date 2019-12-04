function [PSD_ALIASING, PSD_SERVOLAG] = mk_propagated_aliased_PSD(f, m, d, r0, L0, Rx, Ry, vlx, vly, T, td, pistRem, g, filter, f_filter)
%Scope: Create aliased PSD

%   INPUT PARAMETERS
% f                          ::    a matrix with central frequencies
% m                          ::    factor of max freq to be aliased
% d                          ::    sampling period, max(f) = 1/2/d
% r0, L0                     ::    Fried parameter and outer scale of atmosphere
% Rx, Ry                     ::    The reconstructor in the spatial frequencyy domain, x and y
% vlx, vly                   ::    The wind velocity in x and y
% T                          ::    The WFS integration time
% td                         ::    pure delay if any
% pistRem                    ::    flag to remove or not piston, if
%                                   pistRem > 0 then telescopeDiameter =
%                                   pistRem
% g, leak                    ::    loop gain and leak integrator factor
%
%   OUTPUT PARAMETERS
% W                          ::    the phase aliased PSD
% WSx, WSy                   ::    SH gradients aliased PSD
%
%   CREATION                Carlos Correia, March 2014
%
%   RECORD CHANGE
%
if ~exist('pistRem','var') || isempty(pistRem)
    flag_remove_piston = 0;
else
    flag_remove_piston = 1;
end
D = pistRem;

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
if ~exist('g','var')
    g = [];
end
if isempty(g) || (isscalar(g) && g == 0)
    CL = 0;
else
    CL = 1;
end
if ~exist('leak','var')
    leak = 1;
end

if ~exist('filter','var')
    filter = 'SH';
end

if ~exist('f_filter','var')
    f_filter = [];
end

if isvector(f)
    [fx, fy] = meshgrid(f,f);
else
    fx = f;
    fy = f';
end
if length(r0) > 1
    fractionnalR0 = r0(2:end);
    r0 = r0(1);
else
    fractionnalR0 = 1;
end
nLayer = length(fractionnalR0);
f0 = 1/L0;
Phi_alias = 0;

w = 2*1i*pi*d;
for mi = -m:m
    for ni = -m:m
        if mi~=0 || ni~=0
            fm = fx - mi/d;
            fn = fy - ni/d;
            if flag_remove_piston
                % --- piston removal ---
                %D              = 8;
                PR = pistonRemoval(D,d,f,mi,ni);
            end
            %--------------------------------------------------------------
            % CLOSED-LOOP CASE
            %--------------------------------------------------------------
            if CL
                
                if flag_remove_piston
                    W_mn = (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6).*PR;
                else
                    W_mn = (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6);
                end
                
                if ~isempty(f_filter)
                    W_mn = W_mn.*(abs(fm) < f_filter).* (abs(fn) < f_filter);
                end
                
                %q = (Rx.*w.*fm + Ry.*w.*fn) .* (sinc(d*fm).*sinc(d*fn)) .* ...
                %    (sinc(fm*vlx*T).* sinc(fn*vly*T).*exp(1i*2*pi*fm*vlx*td).*exp(1i*2*pi*fn*vly*td).*g./(1-2*(leak-g)*cos(blx).*cos(bly) + (leak-g)^2));
                if strcmp(filter,'pyr') || strcmp(filter,'pyramid') || strcmp(filter,'Conan') || strcmp(filter,'Fauvarque') 
                    Q = (Rx.*sign(fm) + Ry.*sign(fn)) .* (sinc(d*fm).*sinc(d*fn));
                else
                    Q = (Rx.*w.*fm + Ry.*w.*fn) .* (sinc(d*fm).*sinc(d*fn));
                end
                avr = 0;
                for kLayer = 1:nLayer
                    avr = avr + fractionnalR0(kLayer)*...
                        (sinc(fm*vlx(kLayer)*T).* sinc(fn*vly(kLayer)*T).*...
                        exp(1i*2*pi*fm*vlx(kLayer)*td).*exp(1i*2*pi*fn*vly(kLayer)*td).*g);
                end
                q =  Q.* avr;
                Phi_alias = Phi_alias + W_mn .* (q.*conj(q));
                
                % CODE USED FOR THE JOSA-A ANTI-ALIASING WIENER RECONSTRUCTION
                %                 Phi_alias = Phi_alias + W_mn ...
                %                     .*((Rx.*fm) + (Ry.*fn)).^2 ...
                %                     .* (sinc(d*fm).*sinc(d*fn)).^2 ...
                %                     .* (sinc(fm*vlx*T).* sinc(fn*vly*T).*exp(1i*2*pi*fm*vlx*td).*exp(1i*2*pi*fn*vly*td).*g./(1-2*(leak-g)*cos(blx).*cos(bly) + (leak-g)^2)).^2;
                % EVEN OLDER CODE
                %                 Phi_alias = Phi_alias + piston_rem .* (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6).*(Rx.*fm).^2.*sinc(dx*fm).^2 .* sinc(fm*vlx*T).^2 ...
                %                     .*(Ry.*fn).^2.*sinc(dx*fn).^2 .* sinc(fm*vly*T).^2 ...
                %                     *g^2 .* 1./(1-2*(1-g)*cos(blx).*cos(bly) + (1-g)^2);
            else
                %--------------------------------------------------------------
                % OPEN-LOOP CASE
                %--------------------------------------------------------------
                if flag_remove_piston
                    W_mn = (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6).*PR;
                else
                    W_mn = (abs((fm.^2 + fn.^2)) + f0^2).^(-11/6);
                end
                
                if ~isempty(f_filter)
                    W_mn = W_mn.*(abs(fm) < f_filter).* (abs(fn) < f_filter);
                end
                
                if strcmp(filter,'pyr') || strcmp(filter,'pyramid')
                    Q = (Rx.*sign(fm) + Ry.*sign(fn)) .* (sinc(d*fm).*sinc(d*fn));
                else
                    Q = (Rx.*w.*fm + Ry.*w.*fn) .* (sinc(d*fm).*sinc(d*fn));
                end
                %Q = (Rx.*fm + Ry.*fn) .* (sinc(d*fm).*sinc(d*fn));
                %Q = Rx*sign(mi) + Ry*sign(ni);
                avr = 0;
                for kLayer = 1:nLayer
                    avr = avr + fractionnalR0(kLayer)*...
                        (sinc(fm*vlx(kLayer)*T).* sinc(fn*vly(kLayer)*T).*...
                        exp(1i*2*pi*fm*vlx(kLayer)*td).*exp(1i*2*pi*fn*vly(kLayer)*td));
                end
                q =  Q.* avr;
                Phi_alias = Phi_alias + W_mn .* (q.*conj(q));
                
                % CODE USED FOR THE JOSA-A ANTI-ALIASING WIENER RECONSTRUCTION
                %                 Phi_alias = Phi_alias + W_mn...
                %                     .*(abs(Rx.*w.*fm + Ry.*w.*fn).^2) ...
                %                     .* (sinc(d*fm).*sinc(d*fn)).^2 ...
                %                     .* (sinc(fm*vlx*T).* sinc(fn*vly*T).*exp(1i*2*pi*fm*vlx*td).*exp(1i*2*pi*fn*vly*td)).^2;
                %                 imagesc(abs(Phi_alias)), %pause
            end
            
        end
    end
end

Phi_alias(isnan(Phi_alias)) = 0;
if flag_remove_piston
PSD_ALIASING = r0^(-5/3)*0.0229.*Phi_alias.*pistonRemoval(D,d,f,0,0);
else
PSD_ALIASING = r0^(-5/3)*0.0229.*Phi_alias;
end    
PSD_ALIASING = real(PSD_ALIASING);
if nargout > 1
    % ---------------------------------------------------
    % SERVO-LAG
    % ---------------------------------------------------
    if isempty(g)
        g = 0;
    end
    
    if CL
        blx = 2*pi*fx*vlx*T;
        bly = 2*pi*fy*vly*T;
        
        % Eq. (34) Flicker's "Analytical evaluations of closed-loop adaptive optics spatial power spectral densities"
        %         PSD_SERVOLAG = r0^(-5/3)*0.0229.*(abs((fx.^2 + fy.^2)) + f0^2).^(-11/6).*...
        %             abs((Rx.*fy + Ry.*fy).^2 .* ...
        %             (sinc(d*fx).*sinc(d*fy)).^2 .*...
        %             (1-sinc(fx*vlx*T).*sinc(fy*vly*T).*exp(1i*2*pi*fx*vlx*td).*exp(1i*2*pi*fy*vly*td).*g.*(exp(1i*(blx+bly))-(leak-g))./(1-2*(leak-g)*cos(blx).*cos(bly) + (leak-g)^2))).^2;
        
        
        % using arbitrary transfer functions
        ft = blx + bly;
        Z = exp(-1i*ft);
        
        % 1/ INTEGRATOR
        Hint = 1./(1-Z);
        RTF = 1./(1+g.*Z.^(-2).*Hint);
        %NTF = Z.^(-1).*g.*Hint.*RTF;
        
        q = (Rx.*fx + Ry.*fy) .*  sinc(d*fx).*sinc(d*fy) .*  (RTF);
        PSD_SERVOLAG = r0^(-5/3)*0.0229.*(abs((fx.^2 + fy.^2)) + f0^2).^(-11/6).*(q.*conj(q));
        
        %     PSD_SERVOLAG = r0^(-5/3)*0.0229.*(abs((fx.^2 + fy.^2)) + f0^2).^(-11/6).*...
        %         abs((Rx.*fy + Ry.*fy).^2 .* ...
        %         (sinc(d*fx).*sinc(d*fy)).^2 .*...
        %         (RTF)).^2;
        %     %a =  abs((Rx.*fy + Ry.*fy) .* (sinc(d*fx).*sinc(d*fy)) .* RTF).^2;
        % 2/ LQG
        
    else
        q = (Rx.*fx + Ry.*fy) .* sinc(d*fx).*sinc(d*fy) .* (1-sinc(fx*vlx*T).*sinc(fy*vly*T).*exp(1i*2*pi*fx*vlx*td).*exp(1i*2*pi*fy*vly*td));
        PSD_SERVOLAG = r0^(-5/3)*0.0229.*(abs((fx.^2 + fy.^2)) + f0^2).^(-11/6).*(q.*conj(q));
        %     PSD_SERVOLAG = r0^(-5/3)*0.0229.*(abs((fx.^2 + fy.^2)) + f0^2).^(-11/6).*...
        %         abs((Rx.*fy + Ry.*fy).^2 .* ...
        %         (sinc(d*fx).*sinc(d*fy)).^2 .*...
        %         (1-sinc(fx*vlx*T).*sinc(fy*vly*T).*exp(1i*2*pi*fx*vlx*td).*exp(1i*2*pi*fy*vly*td))).^2;
    end
end

function PR = pistonRemoval(D,d,f,mi,ni)
if ~isvector(f)
    f = f(1,:);
end
besselargx     = pi*D/d*(f -mi/d) * d; % here I should not divide by dx because the frequency vector is already in proper units
besselargy     = pi*D/d*(f -ni/d) * d;
[BX, BY]       = meshgrid(besselargx,besselargy);
fltr           = 2*besselj(1,sqrt(BX.^2+BY.^2))./(sqrt(BX.^2+BY.^2));
fltr(isnan(fltr)) = 1;
PR = 1-abs(fltr).^2;

