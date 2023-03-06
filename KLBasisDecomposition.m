function [kl_basis, M2C, eig, coef] = KLBasisDecomposition(tel,atm,dm,oversampling,vectorised,dmLowRes)
%[kl_basis, Bp, eig, coef] = KLBasisDecomposition(influenceFunctions,D,r0,L0,oversampling)
% This function computes the Kharunen Loeve Basis in the DM space applying a doule digaonalization (geometrical and stastical).
% This code is adapted for OOMAO from an ESO code based on the work develloped by E. Gendron (1995).
%   _ tel is a telescope object
%   _ atm is an atmosphere object
%   _ dm is a deformable mirror object
%   _ oversampling is the oversampling factor used to compute the FFT of
%     the Influence Functions.
%   By default, the Piston and Tip Tilt are removed from the Inf. Functions
%   this can be changed easily.

if ~exist('vectorised','var') || isempty(vectorised)
    vectorised = 0;
end

useNewCode = 1; % use this version of the code with updates on how to generate a vector space orthogonal to the original space spanned by the IFs
% number of pixels in the pupil
Spup=tel.pixelArea;

% reshape influence functions
IF=reshape(full(dm.modes.modes),tel.resolution,tel.resolution,dm.nValidActuator);
nAct=dm.nValidActuator;
nRes = size(IF,1);

% truncate influence functions with Telescope Pupil.
% I was having edge effect without doing it.
IF = bsxfun(@times, IF, tel.pupil);
%% REMOVE MODES FROM INFLUENCE FUNCTIONS
% compute modes to remove
Z=zernike(tel,1:3);
% normalize with RMS
rm_modes=Z.modes;

rm_modes=rm_modes*diag(1./rms(rm_modes(tel.pupilLogical,:)));


rm_modes2D=reshape(rm_modes,nRes,nRes,Z.nMode);
% Keep pure Tip & Tilt
TT=reshape(rm_modes2D(:,:,2:3),nRes*nRes,2);
% Remove piston and TT from influence functions
number_of_modesToBeRemoved = size(rm_modes2D,3);
coef = zeros(nAct, number_of_modesToBeRemoved);

fprintf('Removing Piston and TT...')
if ~vectorised
    n_rm_modes=size(rm_modes2D,3);
    tStart = tic;
    for cpt=1:dm.nValidActuator
        cmpt=zeros(nRes,nRes);
        for kj=1:n_rm_modes
            mode=squeeze(rm_modes2D(:,:,kj));
            coef=sum(sum(IF(:,:,cpt).*mode))/sum(sum(mode.*mode));
            cmpt=cmpt+mode.*coef;
        end
        IF(:,:,cpt)=IF(:,:,cpt)-cmpt;
    end
    fprintf('Low-order removal took: %2.2f s\n',toc(tStart))
else
    if ~useNewCode
        % Equivalent code :: vectorised :: 10x faster
        tStart = tic;
        'Using vectorised option'
        normCoeff = diag(rm_modes'*rm_modes);
        proj = rm_modes'*reshape(IF, nRes^2,nAct);
        coef = bsxfun(@times, proj, 1./normCoeff);
        IF = reshape(IF, nRes^2,nAct) - rm_modes*coef;
        IF = reshape(IF, nRes, nRes, []);
        fprintf('Low-order removal took: %2.2f s\n',toc(tStart))
    else
        % NEW NEW Equivalent code :: vectorised :: 10x faster
        tStart = tic;
        'Using vectorised option'
        IF = reshape(IF, nRes^2,nAct);
        Delta = IF'*IF;
        tau = pinv(Delta)*IF'*rm_modes;
        %IF = IF*(eye(324) - tau*pinv(tau));
        IF = IF*(eye(nAct) - tau*pinv(tau'*Delta*tau)*tau'*Delta); % this
        % is the formulation in Gendron/Ferreira18 and Verinaud that seems
        % to lead to consistent results as previous line. However, the
        % generated KL modes are much more similar to Zernike's in the
        % low-order range than the previous. This implementation is
        % therefore preferred
        
        IF = reshape(IF, nRes, nRes, []);
        fprintf('Low-order removal took: %2.2f s\n',toc(tStart))
    end
end
disp('done!')
%%
% figure(108)
% subplot(2,2,1)
% plot(IF0'*IForig*tau)
% title('Original code')
% colorbar
% 
% subplot(2,2,2)
% plot(IF0'*IF1*tau)
% title('Vectorised Original code')
% colorbar
% 
% subplot(2,2,3)
% plot(IF0'*IF2*tau)
% title('Orthogonal in \tau code')
% colorbar
% 
% subplot(2,2,4)
% plot(IF0'*IF3*tau)
% title('Gendron/Ferreira/Verinaud code')
% colorbar

%% GEOMETRIC COVARAIANCE MATRIX
IFperp = reshape(IF,nRes^2,nAct);
IFCovariance = IFperp'*IFperp/Spup;
    
% Compute diagonalisation of matrix IFCovariance.
[U1,S1]=svd(IFCovariance);
S1=diag(real(S1));

% Filter out the modes that were removed from the influence functions
S1((end-number_of_modesToBeRemoved+1):end) = inf;

% store values in eig
eig{1}=S1;

U1=real(U1);
M = zeros(nAct);

% Normalization of M with squared root of the eigen values
if ~vectorised
    for x = 1:nAct
        M(:,x) = U1(:,x)/sqrt(S1(x));
    end
else
    M = U1*diag(1./sqrt(S1));
end

%% STATISTICAL COVARIANCE MATRIX (FOURIER SPACE)
% 1) Fourier transform of the influence functions:
if ~vectorised
    FT_IF = complex(zeros(oversampling*nRes, oversampling*nRes, nAct));
    for x = 1:nAct
        support = zeros(oversampling*nRes, oversampling*nRes);
        support(1:nRes, 1:nRes) = IF(:,:,x);
        FTsupport = fft2(support);
        FT_IF(:,:,x)=FTsupport;
    end
else
    FT_IF = complex(zeros(oversampling*nRes, oversampling*nRes, nAct));
    
    FT_IF= fft2(IF, oversampling*nRes,  oversampling*nRes);
end
clear IF % Last use of IF variable - JFS and MTa 2019-10-29

% 2) Generation of Phase Spectrum in Fourier Space
sp_freq        = generateDistanceGrid(oversampling*nRes)/(oversampling*tel.D);
phase_spectrum = phaseStats.spectrum(sp_freq, atm);

% normalization factor
norma=(Spup*Spup)*(oversampling*tel.D)*(oversampling*tel.D);

fprintf('Statistical Covariance Matrix Computation...')

IF0 = 0*FT_IF;
if ~vectorised
    for cpt = 1:nAct
        IF0(:,:,cpt) = FT_IF(:,:,cpt).*phase_spectrum;
    end
else
    IF0 = bsxfun(@times, FT_IF, phase_spectrum);
end
IF2 = reshape(IF0,(oversampling*nRes)^2,nAct)'; % add transpose, remove from line below
IF3 = conj(reshape(FT_IF,(oversampling*nRes)^2,nAct));

clear IF0     % Last use of IF0 variable - JFS and MTa 2019-10-29
clear FT_IF  % Last use of FT_IF variable - JFS and MTa 2019-10-29

IFFTCovariance = (real(IF2)*real(IF3) + imag(IF2)*imag(IF3))/norma; % remove transpose, add to line above
disp('done!')

% For some reason IF2'*IF3 was not giving the correct result, so I do the
% multiplication separately between real and image parts. I realized that
% the imaginary part of very small, so I do not compute it anymore. the
% full computation would be:
% IFFTCovariance = (real(IF2)*real(IF3) + imag(IF2)*imag(IF3))/norma + ...
%    1i*((real(IF2)*imag(IF3) + imag(IF2)*real(IF3))/norma);

if 0 %ccorreia 9/12/2022: working on a single EIGEN decomposition instead...
    if exist('dmLowRes','var') & ~isempty(dmLowRes)
    [x,y] = meshgrid(linspace(-tel.D/2, tel.D/2,sqrt(size(dmLowRes.modes.modes,1))));
    COVMAT = phaseStats.covarianceMatrix(x+1i*y,atm);
    iIFLowRes = pinv(full(dmLowRes.modes.modes));
    COVU = iIFLowRes*COVMAT*iIFLowRes';
    [U2,S2]=svd(COVU);
    else
    
    COVMAT = phaseStats.covarianceMatrix(dm.modes.actuatorCoord,atm);
    iIFCovariance = pinv(IFCovariance);
    iIFperp = iIFCovariance*IFperp';
    
    Hp=COVMAT;%IFFTCovariance;
    
    %[U2,S2]=svd((iIFperp*IFperp)'*Hp*(iIFperp*IFperp)); % = eig(Hp'*Hp); in theory. In practice I do not find that...
    
    [U2,S2]=svd(iIFCovariance*IFCovariance*Hp*(iIFCovariance*IFCovariance)');
    end
   
    
    S2=diag(real(S2));
    eig{2}=S2;
    
    U2=real(U2);
    
    % Modes To Commands matrix
    Bp0=U2;
    
    % KL Modes excluding TT
    kl_basis0SingleDiag = IFperp*Bp0;
    
    % KL Basis in the phase space, including TT (The 3 last modes are removed)
    kl_basisSingleDiag=[TT kl_basis0SingleDiag(:,1:end-3)];
    
    %PLOT
    % modes have a sound shape, 
    for i = 1:25, subplot(5,5,i), imagesc(reshape(kl_basis0SingleDiag(:,i),480,480)),colorbar, end 
    %... but are not orthogonal wrt phase (although close to the S1
    %decomposition...
    plot(diag(kl_basis0SingleDiag'*kl_basis0SingleDiag))
    hold on
    plot(S1)
    set(gca, 'YScale','log')
end

% double diagonalisationb
Hp=M'*IFFTCovariance*M; % should be in the general case Hp = inv(M'*IFCovariance*M) * IFFTCovariance * inv(M'*IFCovariance*M)'; Since (M'*IFCovariance*M) = I, the equation simplifies

[U2,S2]=svd(Hp);

S2=diag(real(S2));
eig{2}=S2;

U2=real(U2);

% Modes To Commands matrix
Bp0=M*U2;
% KL Modes excluding TT
kl_basis0 = IFperp*Bp0;

% KL Basis in the phase space, including TT (The 3 last modes are removed)
kl_basis=[TT kl_basis0(:,1:end-3)];

% Mode To Command matrix - Project on full span of the DM IF space [correia 27/02/2022]
M2C=pinv(full(dm.modes.modes(tel.pupilLogical,:)))*kl_basis(tel.pupilLogical,:)/2;

end



function R = generateDistanceGrid(N,M)
%
%	Returns an (N,M) floating array in which:
%
%	R(i,j) = SQRT(F(i)^2 + G(j)^2)   where:
%		 F(i) = i  IF 0 <= i <= n/2
%		      = n-i  IF i > n/2
%		 G(i) = i  IF 0 <= i <= m/2
%              = m-i  IF i > m/2

if (nargin<2)
    M=N;
end

R=zeros(N,M);
for x=1:N
    if (x-1<=N/2)
        f=(x-1)^2;
    else
        f=(N-x+1)^2;
    end
    for y=1:M
        if (y-1<=M/2)
            g=(y-1)^2;
        else
            g=(M-y+1)^2;
        end
        R(x,y)=sqrt(f+g);
    end
end
end