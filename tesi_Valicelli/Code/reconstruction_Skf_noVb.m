function [k1x,k2x,k3x,k4x,Cf0x,Cp0x,Cx,Cxdata,relerr,iter,cont,nit] = ...
    reconstruction_Skf_noVb(Cdata,Ca,Vb,Vi,t,t0,C0,k1,k2,k3,k4,Cf0,Cp0)

%% Input
% Cdata : n_time x 1 array of double
%   Total concentration
% Ca : function handle
%   Input function
% Vb : double
%   Blood volume function
% Vi : double
%   Interstizial volume fraction
% t : 1xn_time array of double
%   Time instants
% t0 : double
%   Time for intial state of the Cauchy problem
% C0 : 1x2 array of double
%   Initial state of the Cauchy problem
% k1, k2, k3, k4: doubles
%   Initial value of the rate constants for the iterative method.

%% Output
% k1x,k2x,k3x,k4x : double
%   Estimated value for the rate constants
% Cx : 2 x n_time array of double
%   Estimated concentration for each compartment
% Cxdata : n_time x 1
%   Estimated total concentration
% relerr : double
%   Relative error for the total concentration
% iter : int
%   Number of iterations

%%
narg_1 = 7;
narg_opt = 13;

% To sum the compartment concentrations
alpha = [1-Vb,1-Vb-Vi];

% Initialization of kinetic parameters
% switch nargin
%     case narg_1 % real data
%         k1x = rand(1);
%         k2x = rand(1);
%         k3x = rand(1);
%         k4x = 0;
%         Cf0x =0; %rand(100);
%         Cp0x = 0;%rand(100);
%     case narg_opt % simulation
%         Kx = num2cell([k1,k2,k3,k4,Cf0,Cp0]+(rand(1,6)-4/5).*[k1,k2,k3,k4,Cf0,Cp0]);%[(rand(1,4)-4/5).*[k1,k2,k3,k4],0.1*Cf0,0.1*Cp0]);
%         [k1x,k2x,k3x,k4x,Cf0x,Cp0x] = deal(Kx{:});
% end
            k1x = k1;
            k2x = k2;
            k3x = k3;
            k4x = k4;
            Cf0x =Cf0;
            Cp0x = Cp0;
         
x = [k1x;k2x;k3x;k4x;Cf0x;Cp0x];

% Solve DIRECT PROBLEM (with initial parameters)
Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
Cx = concentration(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t);
Cxdata = (alpha*Cx + Vb*Ca(t))';

% Relative error
relerr = norm(Cdata-Cxdata)/norm(Cdata);

%% GAUSS - NEWTON Regularized method

% Iteration numbers
nit = 0;
iter = 0;
cont = 0;
nit_max = 100;%100;%30;

% Initializations
vect_relerr = zeros(3*nit_max,1);
cell_K = cell(1,3*nit_max);

% Lower bound for the relative error
toll = 5e-3;%sqrt(sum(Cdata))/norm(Cdata)*5e-1;

while relerr>toll
    
    % Count iteration's number
    nit = nit+1; % nit = number of iteration for a single choise of the initial parameters
    iter = iter+1; % iter = numer of iteration for the whole algorithm
    
    % Derivative with respect to the parameters
    aux1 = dF_dk1(Mx,Ca,t0,t,alpha)';
    auxM = dF_dM(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t,alpha)';
    D = [aux1 auxM];
    D = [D(:,1),-D(:,2),D(:,3)-D(:,2),D(:,4)-D(:,5)]; %qui aggiungo altri parametri
    auxcf0 = dF_dCf(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t,alpha)';
    auxcp0 = dF_dCp(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t,alpha)';
%     D = [D, alpha*expm((t-t0)*Mx)*[1,0]']; %Cf0
%     D = [D, alpha*expm((t-t0)*Mx)*[0,1]']; %Cp0
    D = [D, auxcf0 ,auxcp0];

    
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
    h = (r*diag([1,1,1,1,1,1])+D.'*D)\(D.'*(Cdata-Cxdata)); %Da tenere controllata questa riga
    
    % Check for positive coefficients
    xph = x+h;
    if any(xph<=0)
        h(xph<=0)=0;
    end
    x = x+h;
    
    % Refresh the parameters
    k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4); Cf0x = x(5); Cp0x = x(6);
    
    % Solve DIRECT PROBLEM  --> refresh data
    Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
    Cx = concentration(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t);
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
            k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4); Cf0x = x(5); Cp0x = x(6);
            
            % Solve DIRECT PROBLEM (with the parameters leading to the smallest relerr)
            Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
            Cx = concentration(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t);
            Cxdata = (alpha*Cx + Vb*Ca(t))';
            
            relerr = min(vect_relerr);
            
            break
        else
            
            % Initialization of kinetic parameters
%             switch nargin
%                 case narg_1 % real data
%                     k1x = rand(1);
%                     k2x = rand(1);
%                     k3x = rand(1);
%                     k4x = 0;
%                     Cf0x = 0;
%                     Cp0x = 0;
%                    
%                 case narg_opt % simulation
%                     Kx = num2cell([k1,k2,k3,k4,Cf0,Cp0]+(rand(1,6)-4/5).*[k1,k2,k3,k4,Cf0,Cp0]);%[(rand(1,4)-4/5).*[k1,k2,k3,k4],0.1*Cf0,0.1*Cp0]);%num2cell([k1,k2,k3,k4,Cf0,Cp0,Vb]+(rand(1,4)-4/5).*[k1,k2,k3,k4,Cf0,Cp0,Vb]);
%                     [k1x,k2x,k3x,k4x,Cf0x,Cp0x] = deal(Kx{:});
%             end
            k1x = k1;
            k2x = k2;
            k3x = k3;
            k4x = k4;
            Cf0x =Cf0;
            Cp0x = Cp0;
         
            x = [k1x;k2x;k3x;k4x;Cf0x;Cp0x];
            
            % Solve DIRECT PROBLEM (with initial parameters)
            Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
            Cx = concentration(k1x,Mx,Ca,t0,[Cf0x,Cp0x],t);
            Cxdata = (alpha*Cx + Vb*Ca(t))';
            
            % Relative error
            relerr = norm(Cdata-Cxdata)/norm(Cdata);
            
            continue
            
        end
        
    end
  
% .........................................................................%
%     figure(1)
%     plot(t,Cdata,'b')
%     hold on
%     plot(t,Cxdata,'c--')
%     
%     pause
% .........................................................................%
%     
end

end
