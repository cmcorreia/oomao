function [kl_basis, M2C, eig, coef] = KLBasisDecomposition(tel,atm,dm,oversampling)
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

% number of pixels in the pupil
Spup=tel.pixelArea;

% reshape influence functions
IF=reshape(full(dm.modes.modes),tel.resolution,tel.resolution,dm.nValidActuator);
nAct=dm.nValidActuator;
nRes = size(IF,1);

% truncate influence functions with Telescope Pupil.
% I was having edge effect without doing it.
IF=IF.*repmat(tel.pupil,1,1,nAct);
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
n_rm_modes=size(rm_modes2D,3);

for cpt=1:dm.nValidActuator
    cmpt=zeros(nRes,nRes);
    for kj=1:n_rm_modes
        mode=squeeze(rm_modes2D(:,:,kj));
        coef=sum(sum(IF(:,:,cpt).*mode))/sum(sum(mode.*mode));
        cmpt=cmpt+mode.*coef;
    end
    IF(:,:,cpt)=IF(:,:,cpt)-cmpt;
end

disp('done!')

%% GEOMETRIC COVARAIANCE MATRIX
IF1 = reshape(IF,nRes^2,nAct);
IFCovariance = IF1'*IF1/Spup;

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
for x = 1:nAct
    M(:,x) = U1(:,x)/sqrt(S1(x));
end

%% STATISTICAL COVARIANCE MATRIX (FOURIER SPACE)
% 1) Fourier transform of the influence functions:
FT_IF = complex(zeros(oversampling*nRes, oversampling*nRes, nAct));
for x = 1:nAct
    support = zeros(oversampling*nRes, oversampling*nRes);
    support(1:nRes, 1:nRes) = IF(:,:,x);
    FTsupport = fft2(support);
    FT_IF(:,:,x)=FTsupport;
end

clear IF % Last use of IF variable - JFS and MTa 2019-10-29

% 2) Generation of Phase Spectrum in Fourier Space 
sp_freq        = generateDistanceGrid(oversampling*nRes)/(oversampling*tel.D);
phase_spectrum = phaseStats.spectrum(sp_freq, atm);

% normalization factor
norma=(Spup*Spup)*(oversampling*tel.D)*(oversampling*tel.D);

fprintf('Statistical Covariance Matrix Computation...')

IF0 = 0*FT_IF;
for cpt = 1:nAct
    IF0(:,:,cpt) = FT_IF(:,:,cpt).*phase_spectrum;
end
IF2 = reshape(IF0,(oversampling*nRes)^2,nAct)'; % add transpose, remove from line 106
IF3 = conj(reshape(FT_IF,(oversampling*nRes)^2,nAct));

clear IF0     % Last use of IF0 variable - JFS and MTa 2019-10-29
clear FT_IF  % Last use of FT_IF variable - JFS and MTa 2019-10-29

IFFTCovariance = (real(IF2)*real(IF3) + imag(IF2)*imag(IF3))/norma; % remove transpose, add to line 101
disp('done!')

% For some reason IF2'*IF3 was not giving the correct result, so I do the
% multiplication separately between real and image parts. I realized that
% the imaginary part of very small, so I do not compute it anymore. the
% full computation would be:
% IFFTCovariance = (real(IF2')*real(IF3) + imag(IF2')*imag(IF3))/nrm + ...
%     1i*((real(IF2')*imag(IF3) + imag(IF2')*real(IF3))/nrm);

% double diagonalisation
Hp=M'*IFFTCovariance*M;

[U2,S2]=svd(Hp);

S2=diag(real(S2));
eig{2}=S2;

U2=real(U2);

% Modes To Commands matrix
Bp0=M*U2;
% KL Modes excluding TT
kl_basis0 = IF1*Bp0;

% KL Basis in the phase space, including TT (The 3 last modes are removed)
kl_basis=[TT kl_basis0(:,1:end-3)];

% Mode To Command matrix
M2C=pinv(IF1(tel.pupilLogical,:))*kl_basis(tel.pupilLogical,:);

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
