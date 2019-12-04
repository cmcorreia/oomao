%> @file TransformC.m
%> @brief TransformC returns transformed Zernike coefficient set, C2, from the original set, C1, 
%> both in standard ANSI order, with the pupil diameter in mm as the first term. Scaling and translation is performed first and then rotation.
%> @author Linda Lundstrom and Peter Unsbo
%> @date   March 2007 (from the original article)
%> 
%> @par TransformC - transforming Zernike coeffs for shifted and rotated pupils
%> The @b TransformC function was proposed in the article ``Transformation of 
%> Zernike coefficients: scaled, translated, and rotated wavefronts with circular 
%> and elliptical pupils'' by Linda Lundstrom and Peter Unsbo, published in 
%> Vol. 24, No. 3, March 2007,  J. Opt. Soc. Am. A (pages 569-577).
%> 
%> The algorithm presents the means to transform Zernike coefficients analytically 
%> with regard to concentric scaling, translation of pupil center, and rotation. 
%> The transformations are described both for circular and elliptical pupils. 
%> 
%> TransformC returns transformed Zernike coefficient set, C2, from the original 
%> set, C1, both in standard ANSI order, with the pupil diameter in mm as the first term. 
%> Scaling and translation is performed first and then rotation.
%======================================================================
%> @param C1     = the original set of Zernike coefficients [vector], with the pupil diameter in mm as the first term.
%> @param dia2	 = new (desired) pupil diameter [mm] for concentric scaling
%> @param tx	 = Cartesian translation coordinates [mm]
%> @param ty	 = Cartesian translation coordinates [mm]
%> @param thetaR = angle of rotation [degrees]
%> @retval C2	 = transformed Zernike coefficient set in standard ANSI order with the pupil diameter in mm as the first term.
% ======================================================================
function C2 = TransformC (C1,dia2,tx,ty,thetaR)
% ''TransformC'' returns transformed Zernike coefficient set, C2, from the original set, C1,
% both in standard ANSI order, with the pupil diameter in mm as the first term.
% dia2 - new pupil diameter [mm]
% tx, ty - Cartesian translation coordinates [mm]
% thetaR - angle of rotation [degrees]
% Scaling and translation is performed first and then rotation.

dia1=C1(1); % Original pupil diameter
C1=C1(2:end);
etaS=dia2/dia1; % Scaling factor
etaT=2*sqrt(tx^2+ty^2) / dia1; % Translation in Cartesian coordinates
thetaT=atan2(ty, tx);
thetaR=thetaR*pi/180; % Rotation in radians
jnm=length(C1)-1; 

nmax=ceil((-3+sqrt(9+8*jnm))/2); 
jmax=nmax*(nmax+3)/2;

S=zeros(jmax+1,1); S(1:length(C1))=C1; C1=S; clear S
P=zeros(jmax+1); % Matrix P transforms from standard to Campbell order
N=zeros(jmax+1); % Matrix N contains the normalization coefficients
R=zeros(jmax+1); % Matrix R is the coefficients of the radial polynomials

CC1=zeros(jmax+1,1); % CC1 is a complex representation of C1
counter=1;

for m=-nmax:nmax % Meridional indexes
    for n=abs(m):2:nmax % Radial indexes
        jnm=(m+n*(n+2))/2;
        P(counter,jnm+1)=1;
        N(counter,counter)=sqrt(n+1);
        for s=0:(n-abs(m))/2
            R(counter-s,counter)=(-1)^s*factorial(n-s)/(factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s));
        end
        if m<0, CC1(jnm+1)=(C1((-m+n*(n+2))/2+1)+1i*C1(jnm+1))/sqrt(2);
        elseif m==0, CC1(jnm+1)=C1(jnm+1);
        else
            CC1(jnm+1)=(C1(jnm+1)-1i*C1((-m+n*(n+2))/2+1))/sqrt(2);
        end
        counter=counter+1;
    end
end

ETA=[]; % Coordinate-transfer matrix

for m=-nmax:nmax
    for n=abs(m):2:nmax
        ETA=[ETA P*(transform(n,m,jmax,etaS,etaT,thetaT,thetaR))];
    end
end

C=inv(P)*inv(N)*inv(R)*ETA*R*N*P;

CC2=C*CC1;
C2=zeros(jmax+1,1); % C2 is formed from the complex Zernike coefficients, CC2

for m=-nmax:nmax
    for n=abs(m):2:nmax
        jnm=(m+n*(n+2))/2;
        if m<0, C2(jnm+1)=imag(CC2(jnm+1)-CC2((-m+n*(n+2))/2+1))/sqrt(2);
        elseif m==0, C2(jnm+1)=real(CC2(jnm+1));
        else
            C2(jnm+1)=real(CC2(jnm+1)+CC2((-m+n*(n+2))/2+1))/sqrt(2);
        end 
    end 
end
C2=[dia2;C2];


%======================================================================
%> @brief The sub-function @b transform returns 
%> @param n      = Meridional indexes of Zernike coefficients
%> @param m      = Radial indexes of Zernike coefficients
%> @param jmax	 = maximum Meridional index
%> @param etaS	 = Scaling factor, etaS=dia2/dia1; 
%> @param etaT	 = Translation in Cartesian coordinates
%> @param thetaT = angle of rotation [degrees]
%> @param thetaR = angle of rotation [degrees]
%> @retval Eta	 =  The matrix \f$\eta\f$ is formed successively, column by column, with the function @b transform, which
%> performs the calculations of Eqs. (20), (23), and (25) from the article ``Transformation of 
%> Zernike coefficients: scaled, translated, and rotated wavefronts with circular 
%> and elliptical pupils'' by Linda Lundstrom and Peter Unsbo, published in 
%> Vol. 24, No. 3, March 2007,  J. Opt. Soc. Am. A (pages 569-577). The
%> implemented code uses complex Zernike coefficients for
%> the calculations and includes Eqs. (7) and (8) to convert
%> between complex and real coefficients. The matrix  \f$\eta \f$ will be a diagonal matrix with each element equal
%> to \f$\eta_s^n \f$, where \f$n\f$ is the exponent of the corresponding \f$\rho\f$ term in \f$ \langle \rho M|\f$.
% ======================================================================
function Eta=transform(n,m,jmax,etaS,etaT,thetaT,thetaR)
% Returns coefficients for transforming a ro^n * exp(i*m*theta)-term (from Eq. 20) into '-terms
Eta=zeros(jmax+1,1);

for p=0:(n+m)/2
    for q=0:(n-m)/2
        nnew=n-p-q; 
        mnew=m-p+q;
        jnm=(mnew+nnew*(nnew+2))/2;
        Eta(floor(jnm+1))=Eta(floor(jnm+1))+nchoosek((n+m)/2,p)*nchoosek((n-m)/2,q)...
            * etaS^(n-p-q)*etaT^(p+q)*exp(1i*((p-q)*(thetaT-thetaR)+m*thetaR));
        %%% nchoosek returns the binomial coefficient, defined as n!/((n–k)! k!). 
    end 
end