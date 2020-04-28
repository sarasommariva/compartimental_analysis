%% Simulation to test the reliability of the BCM model
% with respect to the Skf model

% Scheme BCM --> Skf:
% 0) given time t, blood input function Ca (simulated data)
% 1) select ground truth values for the parameters (k1,k2,k3,k5,k6)
% 2) solve the BCM direct problem with
%  - blood volume fraction (tumor) Vb = 0.15
%  - interstizial volume fraction (tumor) Vi = 0.3
%  - ratio of the intracellular volumes of ER and cytosol v=V_er/V_cyt=0.17
%  - alpha = [Vi+(1-Vr)*(1-Vb-Vi),(1-Vr)*(1-Vb-Vi),Vr*(1-Vb-Vi)]
%    with Vr = V_er/(V_er+V_cyt) = v/(1+v)
%  and compute the Ct datum
% 3) for 50 run of the algorithm:
%    3a) add white gaussian noise on Ct
%    3b) reconstruct by the use of the Sokoloff model with
%  - blood volume fraction (tumor) Vb = 0.15
%  - interstizial volume fraction (tumor) Vi = 0.3
%  - alpha = [1-Vb,1-Vb-Vi]
%    and obtain the parameters (k1,k2,k3,k4)
%    3c) reconstruct by the use of the BCM model to test GN

%% 0) synthetic data

% time vector
time = [7.5 22.50262 37.50008 52.50036 67.50043 82.50053 97.50076 ...
    112.5008 127.501005 142.501 165.0011 195.0016 225.0018 255.002 ...
    285.0021 334.1693 394.1695 454.1707 514.171015 574.1711 685.1639 ...
    835.1641 1065.2923 1370.788 1675.5025 1980.81945 2285.738];
t = (time./60)';
ntime = length(t);

% simulated Input Function
Ca = @(tt)(6.6991e+06*(tt.^5.4135e+00).*exp(-tt/9.65e-02)+3.611738e+02*exp(-tt/2.85452e+01)).';
%--> Ca [arterial] approximated by sum/product of exponential functions

%% 1) ground truth values for the BCM parameters
k1_gt_BCM = 0.4;
k2_gt_BCM = 0.2;
k3_gt_BCM = 0.7;
k5_gt_BCM = 0.5;
k6_gt_BCM = 0.03;

K_gt_BCM = [k1_gt_BCM; k2_gt_BCM; k3_gt_BCM; k5_gt_BCM; k6_gt_BCM];

%% 2) BCM direct problem

% Volume fractions
Vb = 0.15; % blood volume fraction (tumor)
Vi = 0.3; % interstizial volume fraction (tumor)
v = 0.17; Vr = v/(1+v);

% To sum the compartment concentrations
alpha_BCM = [Vi+(1-Vr)*(1-Vb-Vi),(1-Vr)*(1-Vb-Vi),Vr*(1-Vb-Vi)];

% Solve DIRECT PROBLEM (with initial parameters)
M_gt_BCM = [[-(k2_gt_BCM+k3_gt_BCM);k3_gt_BCM;0],[0;-k5_gt_BCM;k5_gt_BCM],[k6_gt_BCM;0;-k6_gt_BCM]];
C_BCM = concentration(k1_gt_BCM,M_gt_BCM,Ca,0,[0;0;0],t);
Ct_BCM = (alpha_BCM*C_BCM + Vb*Ca(t))';

%% 3) inverse problem

% arrays reconstructed solutions
tent = 50;

% Sokoloff
k1_Skf_vec = zeros(tent,1); k2_Skf_vec = zeros(tent,1);
k3_Skf_vec = zeros(tent,1); k4_Skf_vec = zeros(tent,1);
relerr_Skf_vec = zeros(tent,1); iter_Skf_vec = zeros(tent,1);
Cx_Skf_vec = cell(tent,1);

% BCM
k1_BCM_vec = zeros(tent,1); k2_BCM_vec = zeros(tent,1); 
k3_BCM_vec = zeros(tent,1); k5_BCM_vec = zeros(tent,1); k6_BCM_vec = zeros(tent,1);
relerr_BCM_vec = zeros(tent,1); iter_BCM_vec = zeros(tent,1);
Cx_BCM_vec = cell(tent,1);

for n = 1:tent
    
    disp(' '); disp(['n = ',num2str(n)]); 
    
    % add white gaussian noise
    Ct_noisy = awgn(Ct_BCM,30,'measured');
    
    % reconstruction by the Sokoloff model
    disp('reconstruction Sokoloff...')
    [k1_Skf_vec(n),k2_Skf_vec(n),k3_Skf_vec(n),k4_Skf_vec(n),...
        Cx_Skf_vec{n},Cxdata_Skf,relerr_Skf_vec(n),iter_Skf_vec(n)] = ...
        reconstruction_Skf(Ct_noisy,Ca,t,0,[0;0],k1_gt_BCM,k2_gt_BCM,k3_gt_BCM*0.2,k6_gt_BCM);
    
    % reconstruction by the BCM model
    disp('reconstruction BCM...')
    [k1_BCM_vec(n),k2_BCM_vec(n),k3_BCM_vec(n),k5_BCM_vec(n),k6_BCM_vec(n),...
        Cx_BCM_vec{n},Cxdata_BCM,relerr_BCM_vec(n),iter_BCM_vec(n)] = ...
        reconstruction_BCM(Ct_noisy,Ca,t,0,[0;0;0],k1_gt_BCM,k2_gt_BCM,k3_gt_BCM,k5_gt_BCM,k6_gt_BCM);
    
end

% Sokoloff --
% mean values of the reconstructed parameters
Km_Skf = [mean(k1_Skf_vec); mean(k2_Skf_vec); mean(k3_Skf_vec); mean(k4_Skf_vec)];
% standard deviations of the reconstructed parameters 
Kstd_Skf = [std(k1_Skf_vec); std(k2_Skf_vec); std(k3_Skf_vec); std(k4_Skf_vec)];

% solve the direct problem to obtain the 'mean' compartment concentrations
M_Skf = [[-(Km_Skf(2)+Km_Skf(3));Km_Skf(3)],[Km_Skf(4);-Km_Skf(4)]];
Cxm_Skf = concentration(Km_Skf(1),M_Skf,Ca,0,[0;0],t);

% store the results
K_Skf.k1 = k1_Skf_vec; K_Skf.k2 = k2_Skf_vec; K_Skf.k3 = k3_Skf_vec; K_Skf.k4 = k4_Skf_vec;
K_Skf.relerr = relerr_Skf_vec; K_Skf.iter = iter_Skf_vec; K_Skf.comp = Cx_Skf_vec;
K_Skf.mean = Km_Skf; K_Skf.std = Kstd_Skf; K_Skf.comp_mean = Cxm_Skf;

% BCM --
% mean values of the reconstructed parameters
Km_BCM = [mean(k1_BCM_vec); mean(k2_BCM_vec); mean(k3_BCM_vec); mean(k5_BCM_vec); mean(k6_BCM_vec)];
% standard deviations of the reconstructed parameters
Kstd_BCM = [std(k1_BCM_vec); std(k2_BCM_vec); std(k3_BCM_vec); std(k5_BCM_vec); std(k6_BCM_vec)];

% solve the direct problem to obtain the 'mean' compartment concentrations
M_BCM = [[-(Km_BCM(2)+Km_BCM(3));Km_BCM(3);0],[0;-Km_BCM(4);Km_BCM(4)],[Km_BCM(5);0;-Km_BCM(5)]];
Cxm_BCM = concentration(Km_BCM(1),M_BCM,Ca,0,[0;0;0],t);

% store the results
K_BCM.k1 = k1_BCM_vec; K_BCM.k2 = k2_BCM_vec; K_BCM.k3 = k3_BCM_vec; K_BCM.k5 = k5_BCM_vec; K_BCM.k6 = k6_BCM_vec;
K_BCM.relerr = relerr_BCM_vec; K_BCM.iter = iter_BCM_vec; K_BCM.comp = Cx_BCM_vec;
K_BCM.mean = Km_BCM; K_BCM.std = Kstd_BCM; K_BCM.comp_mean = Cxm_BCM;

%% save 

save(strcat('simulation_BCM_Skf.mat'),'t','Ca','K_gt_BCM','C_BCM','Ct_BCM','K_Skf','K_BCM');


















