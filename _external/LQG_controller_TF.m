%{
------------HEADER-----------------
Objective         ::  Compute the LQG transfer functions

Comments:

INPUT VARS
   A,B,C           :: State-space matrices A, B and C
   d               :: the delay (positive); if negative then matrix D is
                      assumed null D=0; 
                      if d is complex, then a two-sided TF is computed by mirroring the frequency variable f = [-f f];
   M               :: Kalman-filter estimator gain (L = AxM)
   K               :: LQG control gain
   T               :: Sampling-interval
   f               :: frequency-vector
   type            :: (string) estimator or predictor
   causal          :: if causal==0 case where u = -Kinf x_{k+1|k}, otherwise in the causal case u = -Kinf x_{k|k-1} 
OUTPUT VARS
   Gi              :: the complex magnitude
   f               :: The frequency vector
   RTF             :: Rejection TF
   NTF             :: Noise rejection TF
   NCPTF           :: Non-common path rejection TF
Created by         :: C. Correia
Creation date      :: 01/12/2010

Change Record:     :: Jun15: revised the predictor and estimator TFs
                             - added the general case when D \neq 0 and B
                             \neq 0 with arbitrary delays
------------HEADER END----------------
%}

% %%
% % EXAMPLE :: COMPARE Z-TRANSFORM MATLAB BUILT-IN TFs TO MY LQG TFs
% %% Kalman filter design - example with CL TFs
% 
% % 1 step delay 
% Ts = 1/100;
% z = tf('z',Ts);
% 
% A = [0.99 0;1 0];
% B = [0;0];
% C = [0 1];
% D = -1;
% 
% sys = ss(A, [B [1;0]], C, [D 0], Ts);
% 
% QN = 0.5;
% RN = 0.5;
% [KEST,Linf,P,Minf] = kalman(sys,QN,RN,0);
% 
% d = 1;
% 
% 
% W = (eye(2)-1/z*A*(eye(2)-Minf*C));
% 
% Kinf = -[1 0];
% 
% Hcl = -Kinf/W*1/z*A*Minf;
% 
% D = (1/z)^d;
% 
% 
% %% 2 step delay
% % Ts = 1/100;
% % z = tf('z',Ts);
% % 
% % A = [1.996004 -0.996008 0;1 0 0; 0 1 0];
% % B = [0;0;0];
% % C = [0 0 1];
% % D = -1;
% % 
% % sys = ss(A, [B [1;0;0]], C, [D 0], Ts);
% % 
% % QN = 0.5;
% % RN = 0.5;
% % [KEST,Linf,P,Minf] = kalman(sys,QN,RN,0);
% % 
% % d = 2;
% % 
% % 
% % W = (eye(size(A))-1/z*A*(eye(size(A))-Minf*C));
% % 
% % Kinf = -[1 0 0];
% % 
% % Hcl = -Kinf/W*1/z*A*Minf;
% % 
% % D = (1/z)^d;
% 
% %--------------------------------------------------------------------------
% % 1.1 PREDICTOR, CAUSAL CASE
% Hlqg = -1/(1+1/z*Kinf/W*A*Minf*D)*1/z*Kinf/W*A*Minf; % This is the correct predictor LQG TF for s_k = C x_k - D u_k :: % revised Jun15 with an extra delay
% 
% RTFz       = 1/(1+z^(-d)*Hlqg);
% [Gi,floc, RTF, NTF] = LQG_controller_TF(A,B,C,d,Minf,Kinf,Ts,[],'predictor',1);
% 
% MAG = bode(RTFz, 2*pi*floc);
% 
% loglog(floc, abs(RTF).^2,'m')
% hold on
% loglog(floc, squeeze(MAG).^2,'k--')
% 
% MAGn = bode(z*Hlqg*RTFz, 2*pi*floc);
% loglog(floc, squeeze(MAGn).^2,'k:')
% 
% %--------------------------------------------------------------------------
% % 1.2/ PREDICTOR CASE
% Hlqg = -1/(1+Kinf/W*(A*Minf*D+B))*Kinf/W*A*Minf; % This is the correct predictor LQG TF for s_k = C x_k - D u_k
% 
% 
% RTFz       = 1/(1+z^(-d)*Hlqg);
% [Gi,floc, RTF, NTF] = LQG_controller_TF(A,B,C,d,Minf,Kinf,Ts,[],'predictor',0);
% 
% MAG = bode(RTFz, 2*pi*floc);
% 
% loglog(floc, abs(RTF).^2,'k')
% hold on
% loglog(floc, squeeze(MAG).^2,'r--')
% 
% MAGn = bode(z*Hlqg*RTFz, 2*pi*floc);
% loglog(floc, squeeze(MAGn).^2,'r:')
% 
% %--------------------------------------------------------------------------
% % 3/ ESTIMATOR CASE
% W = (eye(size(A))-1/z*(eye(size(A))-Minf*C)*A);
% 
% Hlqg = -1/(1+Kinf/W*Minf*D)*Kinf/W*Minf; % This is the correct estimator LQG TF for s_k = C x_k - D u_k
% RTFz       = 1/(1+z^(-d)*Hlqg);
% [Gi,floc, RTF, NTF] = LQG_controller_TF(A,B,C,d,Minf,Kinf,Ts,[],'estimator');
% 
% MAG = bode(RTFz, 2*pi*floc);
% 
% 
% loglog(floc, abs(RTF).^2,'c')
% hold on
% loglog(floc, squeeze(MAG).^2,'g--')
% 
% MAGn = bode(z*Hlqg*RTFz, 2*pi*floc);
% loglog(floc, squeeze(MAGn).^2,'g:')
% 
% 
% title('Rejection transfer functions')
% xlabel('Frequency, [Hz]')
% ylabel('[dB]')
% axis tight
% legend('Predictor, causal, numerical',...
%         'Predictor, causal, Bode',...
%     'Predictor, non-causal, numerical',...
%         'Predictor, non-causal, Bode',...
% 'Estimator, causal, numerical',...
%         'Estimator, causal, Bode')
% 
% %% If matrix D is set to zero and the delay accounted for in the state transition matrix, then the following example provides the same RTF as previously
% 
% A = [0.99 0 0;1 0 0; 0 0 0];
% B = [0;0;1];
% C = [0 1 -1];
% D = 0;
% 
% sys = ss(A, [B [1;0;0]], C, [D 0], Ts);
% 
% QN = 0.5;
% RN = 0.5;
% [KEST,Linf,P,Minf] = kalman(sys,QN,RN,0);
% 
% %--------------------------------------------------------------------------
% % 1/ PREDICTOR CAUSAL CASE
% W = (eye(3)-1/z*A*(eye(3)-Minf*C));
% Kinf = -[1 0 0];
% 
% Hcl = -Kinf/W*1/z*A*Minf;
% 
% D = 0;
% Hlqg = -1/(1+1/z*Kinf/W*(A*Minf*D+B))*1/z*Kinf/W*A*Minf; % This is the correct predictor LQG TF for s_k = C x_k - D u_k :: % revised Jun15 with an extra delay
% 
% 
% RTFz       = 1/(1+z^(-d)*Hlqg);
% [Gi,floc, RTF, NTF] = LQG_controller_TF(A,B,C,-d,Minf,Kinf,Ts,[],'predictor',1);
% 
% MAG = bode(RTFz, 2*pi*floc);
% 
% loglog(floc, abs(RTF).^2,'k.')
% hold on
% loglog(floc, squeeze(MAG).^2,'r.')
% 
% %--------------------------------------------------------------------------
% % 3/ ESTIMATOR CASE
% I = eye(3);
% W = (eye(3)-1/z*(eye(3)-Minf*C)*A);
% D = 0;
% Hlqg = -(1+Kinf/W*(1/z*(I-Minf*C)*B + (1/z)^d*Minf*D))\Kinf/W*Minf; % This is the correct estimator LQG TF for s_k = C x_k - D u_k
% %Hlqg = -(1+Kinf/W*(1/z*(I+Minf*C)*B ))\(Kinf/W*Minf); % This is the correct estimator LQG TF for s_k = C x_k - D u_k
% 
% RTFz       = 1/(1+z^(-d)*Hlqg);
% [Gi,floc, RTF, NTF] = LQG_controller_TF(A,B,C,-d,Minf,Kinf,Ts,[],'estimator',1);
% 
% MAG = bode(RTFz, 2*pi*floc);
% 
% loglog(floc, abs(RTF).^2,'k.')
% hold on
% loglog(floc, squeeze(MAG).^2,'r.')
% 
% title('Rejection transfer functions')
% xlabel('Frequency, [Hz]')
% ylabel('[dB]')
% axis tight
% legend('Predictor, causal, numerical',...
%         'Predictor, causal, Bode',...
%         'Estimator, causal, numerical',...
%         'Estimator, causal, Bode')
% 
%     % CONTROLLER SYNTHESIS
% rlqg = lqgreg(KEST,Kinf)
% z = tf('z',Ts); 
% 
% % FEEDBACK CONTROLLER
% sysCL = feedback(1/z, rlqg)
% [MAG,  ~, W] = bode(sysCL,2*pi*floc);
% loglog(W/2/pi, squeeze(MAG).^2);



function [Gi,f, RTF, NTF, olRTF, NCPTF] = LQG_controller_TF(A,B,C,d,M,K,T,f,type,causal)

imath    = sqrt(-1);
I        = eye(size(A));
Np       = 500-1; %
if isempty(f)
    f        = linspace(1e-4,1/(2*T),Np); f = logspace(-3,log10(1/(2*T)),Np);
    if ~isreal(d)
        f = [-f f];
    end
end
if ~exist('d','var') || isempty(d)
    d = 2; % default command delay is 2
end
d = real(d);

if ~exist('causal','var')
    causal = 1; % case where u = -Kinf x_{k|k-1}, otherwise in the non-causal case u = -Kinf x_{k+1|k}
end
%--------------------------------------------
% --- numerical TF frequency-by-frequency ---
%--------------------------------------------
% --- >>> SEE NOTE ENTITLED "LQG transfer functions" filename LQG_TFs.pdf
% and Correia10, PhD thesis, pg 104

zi       = exp(-2*pi*imath*f*T);
QQi      = zeros(size(A,1), size(A,2),length(f));
Gi       = zeros(size(f));

L = A*M; % input M must be Minf and not Linf
%M = pinv(A)*L;
D = 1;
warning off all
if d < 0
    Dss = 0;
else
    Dss = 1;
end
Dss = repmat(Dss,size(C,1),1); % make sure it's the same left-hand size as C
% --- this method uses state-space A,B,C,(D neq 0): no command states in
% the main state
if strcmp(type,'predictor')
    if causal
        for fi = 1:length(f)
            invzi = zi(fi)^(-1);
            D = Dss*invzi^d; % assume measurement model has d step delay
            QQi(:,:,fi) = (I - invzi*(A - L*C));
            %Gi(:,fi) = -(1 + K/QQi(:,:,fi)*invzi*L*D)\K/QQi(:,:,fi)*L;%*invzi;
            Gi(:,fi) = -(1 + invzi*K/QQi(:,:,fi)*(L*D + B))\invzi*K/QQi(:,:,fi)*L;
        end
    else
        for fi = 1:length(f)
            invzi = zi(fi)^(-1);
            D = Dss*invzi^d; % assume measurement model has d step delay
            QQi(:,:,fi) = (I - invzi*(A - L*C));
            %Gi(:,fi) = -(1 + K/QQi(:,:,fi)*invzi*L*D)\K/QQi(:,:,fi)*L;%*invzi;
            Gi(:,fi) = -(1 + K/QQi(:,:,fi)*(L*D + B))\K/QQi(:,:,fi)*L;
            
        end
    end
else
    for fi = 1:length(f)
        invzi = zi(fi)^(-1);
        D = Dss*invzi^d;% assume measurement model has d step delay
        %QQi(:,:,fi) = (I - invzi*(I - M*C)*A);
        %Gi(:,fi) = -(1 + K/QQi(:,:,fi)*(M*D + invzi*(I - M*C)*B))\K/QQi(:,:,fi)*M;
        QQi = (I - invzi*(I - M*C)*A);
        Gi(:,fi) = -(1 + K/QQi*(M*D + invzi*(I - M*C)*B))\K/QQi*M;
        
        %D = invzi;
        %QQi(:,:,fi) = (I - invzi*A + invzi*M*C*A);
        %Gi(:,fi) = -(1 + K/QQi(:,:,fi)*M*D)\K/QQi(:,:,fi)*M;
    end
end


% --- Averaging effect of the WFS in Fourier space: sinc(0:1/2)
sinci         = sinc(f*T);
sinci =1;
% --- residual phase and noise tranfer functions ---
if d > 0
    RTF       = (1+sinci.*zi.^(-d).*(Gi)).^(-1);
    olRTF     = (Gi)./(1+sinci.*zi.^(-d).*Gi);
    NTF       = zi.^(-1).*(Gi)./(1+sinci.*zi.^(-d).*Gi);
    NCPTF     = zi.^(-2).*(Gi)./(1+sinci.*zi.^(-d).*Gi);
else
    RTF       = (1+sinci.*zi.^(d).*(Gi)).^(-1);
    olRTF     = (Gi)./(1+sinci.*zi.^(d).*Gi);
    NTF       = zi.^(-1).*(Gi)./(1+sinci.*zi.^(d).*Gi);
    NCPTF     = zi.^(-2).*(Gi)./(1+sinci.*zi.^(d).*Gi);
end




