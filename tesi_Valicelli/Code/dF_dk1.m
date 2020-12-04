% der_k1 = dF_dk1(M,Ci,t0,t,alpha) computes the derivative with respect
% to k1 of F(t) = alpha * C(t):
% - C(t) is the analitic solution of the compartmental ODE's system 
%   in terms of concentrations for the FDG model with nc compartments and
%   known input function Ci;
% - alpha vector of weights to sum the compartment concentrations.
%..........................................................................
% ODE's system:
% dC/dt = M*C + k1*Ci*e1 C(t0)=C0
% analytic solution:
% C(t) = k1 * int_t0^t exp(M(t-u)) * Ci(u) * e1 du + exp(M(t-t0)) * C0 
%..........................................................................
% F(t) = alpha * C(t) = alpha * [ k1 * int_t0^t exp(M(t-u)) * Ci(u) * e1 du + exp(M(t-t0)) * C0 ]
% dF_dk1(t) = alpha * dC_dk1(t) = alpha * [ int_t0^t exp(M(t-u)) * Ci(u) * e1 du ]
%..........................................................................
% - INPUT:
% M is a nc x nc matrix.
% Ci is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t0 is a scalar.
% t is the time vector.
% alpha is a vector of length nc.
% - OUTPUT:
% dC_dk1 is a matrix of sixe nc x length(t)
% der_k1 is a vector of size 1 x length(t)

function der_k1 = dF_dk1(M,Ci,t0,t,alpha) 

    ngl = 8; [x,w] = gauss_legendre(ngl);

    nc = length(M);
    nt = length(t);
    
    dC_dk1 = zeros(nc,nt);

    f = @(u)( Ci(u) * getcols(expm((t(1)-u)*M),1) );
    dC_dk1(:,1) = quadglv(f,t0,t(1),x,w);

    for n=2:nt 
        
        f = @(u)( Ci(u) * getcols(expm((t(n)-u)*M),1) );
        dC_dk1(:,n) = expm((t(n)-t(n-1))*M) * dC_dk1(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
    der_k1 = (alpha * dC_dk1);
    
end