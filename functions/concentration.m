% C = concentration(k1,M,Ci,t0,C0,t) computes the analytic solution 
% of the compartmental ODE's system in terms of concentrations 
% for the FDG model with nc compartments and known input function Ci. 
%..........................................................................
% - ODE's system:
% dC/dt = M*C + k1*Ci*e1 C(t0)=C0
% - Analytic solution:
% C(t) = k1 * int_t0^t exp(M(t-u)) * Ci(u) * ei du + exp(M(t-t0)) * C0 
%..........................................................................
% - INPUT:
% k1 is a scalar.
% M is a nc x nc matrix.
% Ci is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t0 is a scalar.
% C0 is a vector of length nc.
% t is the time vector.
% - OUTPUT:
% C is a matrix of size nc x length(t)

% NB: allow directly the computation of the concentration starting from t0 != 0

function C = concentration(k1,M,Ci,t0,C0,t)

    ngl = 8; [x,w] = gauss_legendre(ngl);

    nc = length(M);
    nt = length(t);
    
    C = zeros(nc,nt);

    f = @(u)( k1 * Ci(u) * getcols(expm((t(1)-u)*M),1) );
    C(:,1) = expm((t(1)-t0)*M) * C0(:) + quadglv(f,t0,t(1),x,w);

    for n=2:nt 
        
        f = @(u)( k1 * Ci(u) * getcols(expm((t(n)-u)*M),1) );
        C(:,n) = expm((t(n)-t(n-1))*M) * C(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
end