
clc
clear 
close all
folder_results = fullfile('.', 'data');
load(fullfile(folder_results, 'data_mice.mat'))

% 3.1. Solve forward model with the selected parameters
time = mice.time;
c_t = mice.c_t;
c_b = mice.c_b;

% Original values
Vb_0 = mice.Vb;
Vi_0 = mice.Vi; 
k1_0 = mice.k1; 
k2_0 = mice.k2; 
k3_0 = mice.k3; 
k4_0 = mice.k4;

% M = [[-(k2_0+k3_0);k3_0],[k4_0;-k4_0]];
% alpha = [1-Vb_0,1-Vb_0-Vi_0];
T = 16;
t_0 = time(T);
c_b = @(tt)(interp1([0 time],[0 c_b'], tt,'linear',0));

% Compute original initial conditions C_0_skf
[ct_rec_0, c_comp_rec_0] = forward_Skf(c_b, Vb_0, Vi_0, time, 0, [0,0], k1_0, k2_0, k3_0, k4_0);
C_0_skf_0 = c_comp_rec_0(:,T);
% Compute solution with reference values
[ct_rec, c_comp_rec] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, C_0_skf_0, k1_0, k2_0, k3_0, k4_0);


% Modified values: random values in a range of [k-Pk,k+Pk]
P =1;
% Search random in an interval [k-P*k,k*Pk] 
% if a = k-P*k; b = k+P*k; then 
% K_mod = a+(b-a).*rand(1); with (b-a) = 2*P*k

Vb = 2*P*mice.Vb.*rand(1)+mice.Vb-P*mice.Vb;
Vi = 2*P*mice.Vi.*rand(1)+mice.Vi-P*mice.Vi; 
k1 = 2*P*mice.k1.*rand(1)+mice.k1-P*mice.k1; 
k2 = 2*P*mice.k2.*rand(1)+mice.k2-P*mice.k2; 
k3 = 2*P*mice.k3.*rand(1)+mice.k3-P*mice.k3; 
k4 = 2*P*mice.k4.*rand(1)+mice.k4-P*mice.k4;
% C_0_skf=2*P*C_0_skf_0.*rand(1)+C_0_skf_0-P*C_0_skf_0;
 Cp0= 2*P*C_0_skf_0(2).*rand(1)+C_0_skf_0(2)-P*C_0_skf_0(2);
 Cf0 = ct_rec(1)-Cp0;
 C_0_skf = [Cf0,Cp0];



% Modified  values
% Vb = mice.Vb+P*mice.Vb;
% Vi = mice.Vi+P*mice.Vi; 
% k1 = mice.k1+P*mice.k1;
% k2 = mice.k2+P*mice.k2;
% k3 = mice.k3+P*mice.k3;
% k4 = mice.k4+P*mice.k4;
% %C_0_skf =C_0_skf_0+P*C_0_skf_0; 
% Cp0 = C_0_skf_0(2)+P*C_0_skf_0(2);
% Cf0 = ct_rec(1)-Cp0;
% C_0_skf = [Cf0,Cp0];

% Specific values
%  k1 = 1.0865;%0.7195;%0.578;%0.5230;%0.7141;%0.6438;%
%  k2 = 0.5383;%1.1592;%1.4980;%1.2523;%0.8596;%0.7772;%
%  k3 = 0.0246;%0.0432;%0.0871;%0.0874;%0.0525;%0.0562;%
%  k4 = 0.0308;%0.0621;%0.0947;%0.0185;%0.0511;%0.0491;%
%  Cp0 =65.8560;% 70.6265;
%  Cf0 = ct_rec(1)-Cp0;
% C_0_skf =[Cf0,Cp0];%[195.3751,73.8798];%[ct_rec(1)-55.1799,55.1799];%[220.0824,58.0805];%[207.1280,54.6618];% 
%   Vb = 0.0661;%0.0532;%0.064;%0.0028;%0.0396;%0.0434;%
%  Vi = 0;




Ref = [k1_0,k2_0,k3_0,k4_0,C_0_skf_0(1),C_0_skf_0(2),Vb_0]';
Initialguess = [k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb]';



[k1x,k2x,k3x,k4x,Cx,Cxdata,relerr,iter] = reconstruction_Skf(ct_rec,c_b, Vb_0, Vi_0, time(T:end),t_0,C_0_skf_0,k1,k2,k3,k4);
Param = [k1x,k2x,k3x,k4x]';

[k1x_new,k2x_new,k3x_new,k4x_new,Cf0x_new,Cp0x_new,Vbx_new,Cx_new,Cxdata_new,relerr_new,iter_new,cont] =reconstruction_Skf_New(ct_rec,c_b, Vi, time(T:end),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb);
Param_new = [k1x_new,k2x_new,k3x_new,k4x_new,Cf0x_new,Cp0x_new,Vbx_new]';

[k1x_noVb,k2x_noVb,k3x_noVb,k4x_noVb,Cf0x_noVb,Cp0x_noVb,Cx_noVb,Cxdata_noVb,relerr_noVb,iter_noVb] = reconstruction_Skf_noVb(ct_rec,c_b, Vb_0, Vi_0, time(T:end),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2));
Param_noVb = [k1x_noVb,k2x_noVb,k3x_noVb,k4x_noVb,Cf0x_noVb,Cp0x_noVb]';   



% Table with all the parameters
AA = zeros(7,5);
AA(1:end,1) = Ref;
AA(1:end,2) = Initialguess;
AA(1:4,3) = Param;
AA(1:end,4) = Param_new;
AA(1:end-1,5) = Param_noVb;

% Compute C_t with the calculated parameters
[ct_reck, c_comp_reck] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, C_0_skf_0, k1x, k2x, k3x, k4x);
[ct_rec_new, c_comp_rec_new] = forward_Skf(c_b, Vbx_new, Vi, time(T:end), t_0, [Cf0x_new,Cp0x_new], k1x_new, k2x_new, k3x_new, k4x_new);
[ct_rec_noVb, c_comp_rec_noVb] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, [Cf0x_noVb,Cp0x_noVb], k1x_noVb, k2x_noVb, k3x_noVb, k4x_noVb);



figure(1)
plot(time(T:end),ct_rec,'Linewidth',1)
hold on
% plot(time(T:end),c_t(T:end),'--k','Linewidth',1)%(T:end)
% hold on
plot(time(T:end),ct_reck,'Linewidth',1)
hold on
plot(time(T:end),ct_rec_new,'Linewidth',1)
hold on
plot(time(T:end),ct_rec_noVb,'Linewidth',1)
legend('C_T simulated','C_T k_1-k_4 known','C_T all unknown','C_T V_b unknown','Location','southeast')
%legend('C_T simulated','C_T data','C_T k_1-k_4 known','C_T all unknown','C_T V_b unknown','Location','southeast')
axis([0,time(end),min([min(ct_rec),min(ct_reck),min(ct_rec_new),min(ct_rec_noVb)]),max([max(ct_rec),max(ct_reck),max(ct_rec_new),max(ct_rec_noVb)])])
xlabel('Time')
ylabel('C_T')
title('Comparison between simulated data')


% Errors
ERR = [norm(ct_rec-ct_reck,2)/norm(ct_rec,2),norm(ct_rec-ct_rec_new,2)/norm(ct_rec,2),norm(ct_rec-ct_rec_noVb,2)/norm(ct_rec,2)];
AA(8,3:end) = ERR;
AA = round(AA,4);

% Relative errors bewteen the parameters
Err2(1,:) = [abs(k1x-k1_0)/k1_0,abs(k1x_new-k1_0)/k1_0,abs(k1x_noVb-k1_0)/k1_0];
Err2(2,:) = [abs(k2x-k2_0)/k2_0,abs(k2x_new-k2_0)/k2_0,abs(k2x_noVb-k2_0)/k2_0];
Err2(3,:) = [abs(k3x-k3_0)/k3_0,abs(k3x_new-k3_0)/k3_0,abs(k3x_noVb-k3_0)/k3_0];
Err2(4,:) = [abs(k4x-k4_0)/k4_0,abs(k4x_new-k4_0)/k4_0,abs(k4x_noVb-k4_0)/k4_0];
Err2(5,:) = [0,abs(Cf0x_new-C_0_skf_0(1))/C_0_skf_0(1),abs(Cf0x_noVb-C_0_skf_0(1))/C_0_skf_0(1)];
Err2(6,:) = [0,abs(Cp0x_new-C_0_skf_0(2))/C_0_skf_0(2),abs(Cp0x_noVb-C_0_skf_0(2))/C_0_skf_0(2)];
Err2(7,:) = [0,abs(Vbx_new-Vb_0)/Vb_0,0];
Err2 = round(Err2,4); 

