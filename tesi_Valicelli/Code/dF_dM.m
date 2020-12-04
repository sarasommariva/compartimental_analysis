% der_M = dF_dM(k1,M,Ci,t0,C0,t,alpha) computes the derivative 
% with respect to M of F(t) = alpha * C(t):
% - C(t) is the analitic solution of the compartmental ODE's system 
%   in terms of concentrations for the FDG model with nc compartments and
%   known input function Ci;
% - alpha vector of weights to sum the compartment concentrations.
%..........................................................................
% - ODE's system:
% dC/dt = M*C + k1*Ci*e1 C(0)=[0;0]
% - Analytic solution:
% C(t) = k1 * int_t0^t exp(M(t-u)) * Ci(u) * ei du + exp(M(t-t0)) * C0 
%..........................................................................
% F(t) = alpha * C(t) = alpha * [ k1 * int_t0^t exp(M(t-u)) * Ci(u) * ei du + exp(M(t-t0)) * C0 ]
% dF_dM(t) = alpha * [(dC/dM(t)).H] 
%          = alpha * [int_t0^t kron(C(u).',exp((t-u)M) du] * H(:)
%..........................................................................
% - INPUT:
% k1 is a scalar.
% M is a nc x nc matrix.
% Ci is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t0 is a scalar.
% C0 is a vector of length nc.
% t is the time vector.
% alpha is a vector of length nc.
% - OUTPUT:
% dC_dM is an array of size nc^2 x nc^2 x length(t).
% Each nc^2 x nc^2 matrix gives the values of the derivative at times t.
% der_M is a matrix of size nc^2 x length(t).

function der_M = dF_dM(k1,M,Ci,t0,C0,t,alpha)
  
    ngl = 8; [x,w] = gauss_legendre(ngl);

    nc = length(M);
    nt = length(t);
        
    dC_dM = zeros(nc,nc^2,nt);
        
    C = concentration(k1,M,Ci,t0,C0,t);

    Cu = @(u)( concentration(k1,M,Ci,t0,C0,u) );
    f = @(u)( kron((Cu(u)).',expm((t(1)-u)*M)) );
    
    dC_dM(:,:,1) = quadglv(f,t0,t(1),x,w);
    
    for n=2:nt
        
        Cu = @(u)( concentration(k1,M,Ci,t(n-1),C(:,n-1),u) );
        f = @(u)( kron((Cu(u)).',expm((t(n)-u)*M)) );
        
        dC_dM(:,:,n) = expm((t(n)-t(n-1))*M) * dC_dM(:,:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
   
    der_M = reshape(alpha*reshape(dC_dM,nc,nc^2*nt),nc^2,nt);
     
end    