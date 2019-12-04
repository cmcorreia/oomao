function x = myLevenbergMarquardt(fun,init,obs,e,varargin)
%% MYLEVENBERGMARQUARDT developed by Y.O.
% This function solves non-linear inverse problem with the
% Levenberg-Marquardt method. This code may be faster than lsqnonlin in
% MATLAB but less stable than the lsqnonlin.
%
%
%


%Check inputs
inputs = inputParser;
inputs.addRequired( 'fun'   , @(x) isa(x,'function_handle') );
inputs.addRequired( 'init'  , @isnumeric );
inputs.addRequired( 'obs'   , @isnumeric );
inputs.addRequired( 'e'     , @isnumeric );
inputs.addParameter( 'nMAX'  , 100  , @isnumeric );
inputs.addParameter( 'l0'    , 100  , @isnumeric );
inputs.addParameter( 'lup'   , 5   , @isnumeric );
inputs.addParameter( 'ldown' , 5   , @isnumeric );
inputs.addParameter( 'delta' , 1e-3 , @isnumeric );

inputs.parse(fun,init,obs,e,varargin{:});
nMAX  = inputs.Results.nMAX;
l     = inputs.Results.l0;
lup   = inputs.Results.lup;
ldown = inputs.Results.ldown;
delta = inputs.Results.delta;

n = length(init);
m = length(obs(:));
f = m-n-1;

% Initialization
obs = obs(:);
x = log(init(:)); % nonnegative
y = fun(exp(x)); y = y(:);
dy = obs-y;
chi = sum(dy.^2); % kai2
fprintf('  Iter         chi2             gradient       parameter        lambda\n');
fprintf('  %03d    %10.3e\n',0,chi/f);
if chi/f<e
    x = exp(x);
    return;
end
J = numericalJacobianUpdate(fun,x,delta);
JJ = J'*J;
K = diag(diag(JJ));
L = JJ+l*K;
R = J'*dy;

y_old = y;
dy_old = dy;
chi_old = chi;

% Iteration
for k=1:nMAX
    
    % update paramter
    h = L\R;
    
    % Evaluation
    y = fun(exp(x+h)); y = y(:);
    dy = obs-y;
    chi = sum(dy.^2);
    
    fprintf('  %03d    %10.3e    %10.3e    %10.3e    %10.3e\n',k,chi/f,max(abs(R)),max(abs(h./x)),l);
    if chi/f<e
        break;
    end
    
    p = (chi_old-chi)/(h'*(l*K*h+R));
    if p>0
        
        x = x+h;
        
        if mod(k,2*n)==0 || chi>chi_old
            J = numericalJacobianUpdate(fun,x,delta);
        else
            J = rank1JacobianUpdate(J,y_old,y,h);
        end
        
        y_old = y;
        dy_old = dy;
        chi_old = chi;
        
        l = max(l/ldown,1e-5);
                
        JJ = J'*J;
        K = diag(diag(JJ));
        L = JJ+l*K;
        R = J'*dy;
        
    else
        
        l = min(l*lup,1e+5);
        
        J = numericalJacobianUpdate(fun,x,delta);
        
        JJ = J'*J;
        K = diag(diag(JJ));
        L = JJ+l*K;
        R = J'*dy_old;
    end
end

x = exp(x);

end

%% numerical update of Jacobian
function J = numericalJacobianUpdate(fun, x, delta)

x = exp(x);
nx = length(x);
dx = delta*(1+abs(x));

y1 = fun(x);

J=cell(1,nx);
for k=1:nx
    newx = x;
    newx(k) = newx(k) + dx(k);
    y2 = fun(newx);
    dfun = (y2(:)-y1(:));
    J{1,k} = dfun(:)./abs(dx(k));
end
J=cell2mat(J);
end

%% rank-1 update of Jacobian
function newJ = rank1JacobianUpdate(J,y1,y2,h)
newJ=J+((y2-y1-J*h)*h')./(h'*h);
end
