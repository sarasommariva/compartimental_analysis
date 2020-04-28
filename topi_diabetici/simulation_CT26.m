%% Analysis of CT26 tumor data by the use of
% - the Sokoloff model
% - the BCM model accounting for the ER
% for 50 run of the algorithm

%% EXPERIMENTS

path_data_CT26 = '/Users/mara.scussolini/Documents/Dottorato/COMPARTIMENTALE/TOPI cancer/CT26 (VARIE TAC - controlli)';

mice_CT26 = {'#0 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#1 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#2 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#4 CT26 NO STS 2012.11.14 PRIMA PET';...
    '#1 CT26 NO STS 2012.11.20 SECONDA PET';...
    '#3 CT26 NO STS 2012.11.21 SECONDA PET'};

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

%% for each mouse
for i=1:length(mice_CT26)

disp(' '); disp(['mouse = ',mice_CT26{i}]); 
mouse = mice_CT26{i};

% data import
if i <= 4
    folder = strcat('/PRIMA PET/',mouse,'/');
else
    folder = strcat('/SECONDA PET/',mouse,'/');
end  
[t,Ca,Ct] = import_voistat(path_data_CT26,folder);
    
% data correction
if i == 1 || i == 3
    Ca(21) = (Ca(20)+Ca(22))/2;
    Ct(21) = (Ct(20)+Ct(22))/2;
end
t = t'; Ca = @(tt)(interp1([0 t],[0 Ca'],tt,'linear',0));
    
%% Sokoloff model
disp(' '); disp('reconstruction Sokoloff...');
for n=1:tent
    
    disp(['n = ',num2str(n)]); 

    [k1_Skf_vec(n),k2_Skf_vec(n),k3_Skf_vec(n),k4_Skf_vec(n),...
        Cx_Skf_vec{n},Cxdata_Skf,relerr_Skf_vec(n),iter_Skf_vec(n)] = ...
        reconstruction_Skf(Ct,Ca,t,0,[0;0]);
    
end

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

%% BCM model
disp(' '); disp('reconstruction BCM...');
for n=1:tent
    
    disp(['n = ',num2str(n)]); 
    
    [k1_BCM_vec(n),k2_BCM_vec(n),k3_BCM_vec(n),k5_BCM_vec(n),k6_BCM_vec(n),...
        Cx_BCM_vec{n},Cxdata_BCM,relerr_BCM_vec(n),iter_BCM_vec(n)] = ...
        reconstruction_BCM(Ct,Ca,t,0,[0;0;0]);

end

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

save(strcat('K_',mouse,'.mat'),'mouse','t','Ca','Ct','K_Skf','K_BCM');

end












