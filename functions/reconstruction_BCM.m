function [k1x,k2x,k3x,k5x,k6x,Cx,Cxdata,relerr,iter] = ...
    reconstruction_BCM(Cdata,Ca, Vb, Vi, v, t,t0,C0,k1,k2,k3,k5,k6) 

%% Input
% Cdata : n_time x 1 array of double
%   Total concentration
% Ca : function handle
%   Input function
% Vb : double
%   Blood volume function
% Vi : double
%   Interstizial volume fraction
% v : double
%   v = V_er/V_cyt ratio of the intracellular volumes of ER and cytosol
% t : 1xn_time array of double
%   Time instants
% t0 : double
%   Time for intial state of the Cauchy problem
% C0 : 1x3 array of double
%   Initial state of the Cauchy problem
% k1, k2, k3, k5, k6: doubles
%   Initial value of the rate constants for the iterative method.

narg_1 = 8;
narg_opt = 13;

Vr = v/(1+v);

% To sum the compartment concentrations
alpha = [Vi+(1-Vr)*(1-Vb-Vi),(1-Vr)*(1-Vb-Vi),Vr*(1-Vb-Vi)];

% Initialization of kinetic parameters
switch nargin
    case narg_1 % real data
        k1x = rand(1);
        k2x = rand(1);
        k3x = rand(1);
        k5x = 100*rand(1);
        k6x = 0;        
    case narg_opt % simulation
        Kx = num2cell([k1,k2,k3,k5,k6]+(rand(1,5)-4/5).*[k1,k2,k3,k5,k6]);
        [k1x,k2x,k3x,k5x,k6x] = deal(Kx{:});
end
x = [k1x;k2x;k3x;k5x;k6x]; 

% Solve DIRECT PROBLEM (with initial parameters)
Mx = [[-(k2x+k3x);k3x;0],[0;-k5x;k5x],[k6x;0;-k6x]];
Cx = concentration(k1x,Mx,Ca,t0,C0,t);
Cxdata = (alpha*Cx + Vb*Ca(t))';

% Relative error
relerr = norm(Cdata-Cxdata)/norm(Cdata);

%% GAUSS - NEWTON Regularized method 

% Iteration numbers
nit = 0;
iter = 0;
cont = 0;
nit_max = 30; 

% Initializations
vect_relerr = zeros(3*nit_max,1);
cell_K = cell(1,3*nit_max);

% Lower bound for the relative error
toll = sqrt(sum(Cdata))/norm(Cdata)*5e-1; 

while relerr>toll
    
    % Count iteration's number
    nit = nit+1; % nit = number of iteration for a single choise of the initial parameters 
    iter = iter+1; % iter = numer of iteration for the whole algorithm
    
    % Derivative with respect to the parameters 
    aux1 = dF_dk1(Mx,Ca,t0,t,alpha)'; 
    auxM = dF_dM(k1x,Mx,Ca,t0,C0,t,alpha)';
    D = [aux1 auxM];
    D = [D(:,1),-D(:,2),D(:,3)-D(:,2),D(:,7)-D(:,6),D(:,8)-D(:,10)];
    
    % Regularization parameter --------------------------------------------
    % values for lambda to compute V_lambda
    switch nargin
        case narg_1 % real data
            vect_lambda = 5e4:5e4:5e6;
        case narg_opt % simulation
            vect_lambda = 5e2:5e2:5e4;
    end
    vect_lambda = vect_lambda';    
    % Generalized Cross Validation
    [r,~] = GCV(D,Cdata-Cxdata,vect_lambda);
    %----------------------------------------------------------------------
    
    % '\' solve the linear system Dh=Z with Z=Cdata-Cxdata h=increment
    % ==> Tychonov regularization: (rI+D'*D)h=D'*Z
    h = (r*diag([1,1,1,1,1])+D.'*D)\(D.'*(Cdata-Cxdata));
    
    % Check for positive coefficients
    xph = x+h;
    if any(xph<=0)
        h(xph<=0)=0;
    end
    x = x+h; 
    
    % Refresh the parameters
    k1x = x(1); k2x = x(2); k3x = x(3); k5x = x(4); k6x = x(5);
    
    % Solve DIRECT PROBLEM  --> refresh data
    Mx = [[-(k2x+k3x);k3x;0],[0;-k5x;k5x],[k6x;0;-k6x]];
    Cx = concentration(k1x,Mx,Ca,t0,C0,t);
    Cxdata = (alpha*Cx + Vb*Ca(t))';
    
    % Relative error
    [relerrprec,relerr] = deal(relerr,norm(Cdata-Cxdata)/norm(Cdata));
    
    % Store the relative error 'relerr' and the solution 'x' for each iteration 'iter'
    vect_relerr(iter) = relerr;
    cell_K{iter} = x;
    
    %.....................................................................%
    if  ( nit>=15 && relerr>=0.3 ) || ( relerr>=0.3 && abs(relerr-relerrprec)<1e-3 ) || ( nit>=25 && abs(relerr-relerrprec)<1e-4 ) || (nit>=nit_max)        
        
        nit = 0;
        cont = cont+1;
        
        if cont == 3
            
            vect_relerr(vect_relerr==0) = [];
            x = cell_K{vect_relerr==min(vect_relerr)};
            k1x = x(1); k2x = x(2); k3x = x(3); k5x = x(4); k6x = x(5);
            
            % Solve DIRECT PROBLEM (with the parameters leading to the smallest relerr)
            Mx = [[-(k2x+k3x);k3x;0],[0;-k5x;k5x],[k6x;0;-k6x]];
            Cx = concentration(k1x,Mx,Ca,t0,C0,t);
            Cxdata = (alpha*Cx + Vb*Ca(t))';
            
            relerr = min(vect_relerr);
            
            break
        else
            
            % Initialization of kinetic parameters
            switch nargin
                case narg_1 % real data
                    k1x = rand(1);
                    k2x = rand(1);
                    k3x = rand(1);
                    k5x = 100*rand(1);
                    k6x = 0;
                case narg_opt % simulation
                    Kx = num2cell([k1,k2,k3,k5,k6]+(rand(1,5)-4/5).*[k1,k2,k3,k5,k6]);
                    [k1x,k2x,k3x,k5x,k6x] = deal(Kx{:});
            end
            x = [k1x;k2x;k3x;k5x;k6x];
            
            % Solve DIRECT PROBLEM (with initial parameters)
            Mx = [[-(k2x+k3x);k3x;0],[0;-k5x;k5x],[k6x;0;-k6x]];
            Cx = concentration(k1x,Mx,Ca,t0,C0,t);
            Cxdata = (alpha*Cx + Vb*Ca(t))';
            
            % Relative error
            relerr = norm(Cdata-Cxdata)/norm(Cdata);

            continue
            
        end
        
    end
    
%.........................................................................%
%     figure(1)
%     plot(t,Cdata,'b')
%     hold on
%     plot(t,Cxdata,'c--')
% 
%     pause
%.........................................................................%

end

end
