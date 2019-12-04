%> @file ANSI2Noll.m
%> @brief Convert ANSI ordering to Noll ordering.
%> @author Carlos Correia
%> @date   18 June 2012
%> 
%> @par ANSI2Noll - Convert ANSI ordering to Noll ordering
%> There are many conventions on the Zernike expansions, which differ in the 
%> ordering of modes.
%> 
%> The Zernike ordering scheme used in the telescope optics domain (e.g. ZEMAX 
%> standard Zernike coefficients), based on @b Noll's @b concept, is shown below 
%> for the first 21 terms. Here, the ordering number j is determined by ordering 
%> polynomial with lower radial order first, and for given radial order with odd 
%> number for sine function and even number for the cosine.
%> 
%> @image html zernike_noll.png
%> 
%> @b Noll's @b scheme also has a vertical expansion symmetry similar to Wyant's, 
%> but since it directly uses radial order n, it unfolds somewhat differently. So 
%> for n=6, the term for m=0 is secondary spherical aberration (j=22), for m=2 it 
%> is tertiary astigmatism (j=23 and 24 for the sine and cosine function, 
%> respectively), for m=4 it is secondary quadrafoil (j=25 and 26) and for m=6 the 
%> hexafoil (j=27 and 28).
%> 
%> Zernike scheme commonly used in ophthalmology - the @b ANSI @b standard @b 
%> Zernike @b expansion - is different from those used in assessing telescope 
%> aberrations. The ANSI standard single-index scheme (Zernike pyramid) routinely 
%> used in ophthalmology for evaluating eye aberrations. It is graphically 
%> presented as a pyramid resulting from Zernike term expansion as a function of 
%> radial order n and angular frequency m, with the latter being numerically 
%> positive for cosine function of \f$\theta\f$, and negative for the sine 
%> function 
%> 
%> @image html zernike_pyramide.png
%> 
%> Zernike terms expansion pyramid is a function of term's radial degree (or order) 
%> n and azimuthal frequency m. It is the basis for classifying aberrations as 
%> lower \f$(n\leq 2)\f$ and higher-order \f$(n>2)\f$ in ophthalmology. Shown are 
%> the top 20 terms.
%======================================================================
%> @param nModes     = The number of modes
%> @retval ANSIt	 = The nMode x nModes matrix 
% ======================================================================
function ANSIt = ANSI2Noll(nModes)
[Nj, Mj] = make_zernike_index(nModes); %Output the radial order and azimuthal frequency indexes for up to a given number of Nmodes

ANSIt = sparse(numel(Nj), numel(Nj));

for n = 0:max(Nj)
    Nollj = find(Nj == n);
    nn = numel(Nollj);

    if iseven(nn)

        if iseven(Nollj(min(numel(Nollj),2)))
            coss = Nollj(2:2:end);
            sins = Nollj(3:2:end);
        else
            sins = Nollj(2:2:end);
            coss = Nollj(3:2:end);
        end

        ANSIj = [sins(end:-1:1) Nollj(1) coss];

    else

        if iseven(Nollj(min(numel(Nollj),2)))
            coss = Nollj(2:2:end);
            sins = Nollj(3:2:end);
        else
            sins = Nollj(2:2:end);
            coss = Nollj(3:2:end);
        end
        
        ANSIj = [sins(end:-1:1) Nollj(1) coss];

    end
    
    for rr = 1:numel(Nollj)
        ANSIt(ANSIj(rr), Nollj(rr)) = 1;
    end
    
end


% ------------HEADER-----------------
% Objective:
%   Convert ANSI ordering to Noll ordering
%INPUT:  
%   Nmodes       ::  The number of modes
% 
%OUTPUT:
%   MAT          ::  The nMode x nModes matrix 
%
%Created by      ::  Carlos Correia 
%Creation date   ::  18 June 12
%Change Record:  ::  
% ------------HEADER END-----------------