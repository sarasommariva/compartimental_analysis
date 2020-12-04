clc
clear 
close all
folder_results = fullfile('.', 'data');
load(fullfile(folder_results, 'data_mice.mat'))
tic
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
AA = zeros(7,30);

% Modified values: random values in a range of [k-Pk,k+Pk]
P =1;
% Search random in an interval [k-P*k,k*Pk] 
% if a = k-P*k; b = k+P*k; then 
% K_mod = a+(b-a).*rand(1); with (b-a) = 2*P*k


% N = number of tests
N = 20;
par_loguni_rho = [-2,0];
 %k1 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                %+par_loguni_rho(1));
for i = 1:N
%Vb = 2*P*mice.Vb.*rand(1)+mice.Vb-P*mice.Vb;
%Vi = 2*P*mice.Vi.*rand(1)+mice.Vi-P*mice.Vi; 
k1 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                +par_loguni_rho(1));%2*P*mice.k1.*rand(1)+mice.k1-P*mice.k1; 
k2 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                +par_loguni_rho(1));%2*P*mice.k2.*rand(1)+mice.k2-P*mice.k2; 
k3 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                +par_loguni_rho(1));%2*P*mice.k3.*rand(1)+mice.k3-P*mice.k3; 

k4 = k4_0;
C_0_skf =C_0_skf_0;
Vb =Vb_0;
Vi = 0;


Ref = [k1_0,k2_0,k3_0,k4_0,C_0_skf_0(1),C_0_skf_0(2),Vb_0]';
Initialguess = [k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb]';
IG1(i,:) = Initialguess;

[k1x_new,k2x_new,k3x_new,k4x_new,Cf0x_new,Cp0x_new,Vbx_new,Cx_new,Cxdata_new,relerr_new,iter_new, cont_new,nit_new] =reconstruction_Skf_New(ct_rec,c_b, Vi, time(T:end),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb);
k1x_newvett(i) = k1x_new;
k2x_newvett(i) = k2x_new;
k3x_newvett(i)=k3x_new;
k4x_newvett(i)= k4x_new;
Cf0x_newvett(i) = Cf0x_new;
Cp0x_newvett(i)= Cp0x_new;
Vbx_newvett(i)= Vbx_new;
iternew1(i) = iter_new;
contnew1(i) = cont_new;
nitnew1(i) = nit_new;
Param_new1(i,:) = [k1x_newvett(i),k2x_newvett(i),k3x_newvett(i),k4x_newvett(i),Cf0x_newvett(i),Cp0x_newvett(i),Vbx_newvett(i)]';


% Table with all the parameters
 
 j =2*i; 
 AA(1:7,j-1) = Initialguess;
AA(1:7,j) = Param_new1(i,:);

BB(1:7,i) = Param_new1(i,:);

% Compute C_t with the calculated parameters
%[ct_reck, c_comp_reck] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, C_0_skf_0, k1x, k2x, k3x, k4x);
%[ct_rec_new, c_comp_rec_new] = forward_Skf(c_b, Vbx_new, Vi, time(T:end), t_0, [Cf0x_new,Cp0x_new], k1x_new, k2x_new, k3x_new, k4x_new);
%[ct_rec_noVb, c_comp_rec_noVb] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, [Cf0x_noVb,Cp0x_noVb], k1x_noVb, k2x_noVb, k3x_noVb, k4x_noVb);
% figure(i)
% plot(time(T:end),ct_rec,'Linewidth',1)
% hold on
% plot(time(T:end),ct_rec_new,'Linewidth',1)
% hold on
% legend('C_T simulated','C_T all unknown','Location','southeast')
% axis([0,time(end),min([min(ct_rec),min(ct_rec_new)]),max([max(ct_rec),max(ct_rec_new)])])
% xlabel('Time')
% ylabel('C_T')
% title('Comparison between simulated data')
AA = round(AA,4);

end

Mink1 = min(BB(1,:));
Mink2 = min(BB(2,:));
Mink3 = min(BB(3,:));

Maxk1 = max(BB(1,:));
Maxk2 = max(BB(2,:));
Maxk3 = max(BB(3,:));

%% Second Step k1 k2 k3 in a smaller range. Cf0 and Cp0 random. k4,Vb fixed
for i = 1:N
% Fixed in smaller intervals
k1 = Mink1+(Maxk1-Mink1).*rand(1); 
k2 = Mink2+(Maxk2-Mink2).*rand(1); 
k3 = Mink3+(Maxk3-Mink3).*rand(1); 

k4 = k4_0; 
Vb = Vb_0;

Cp0= 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                +par_loguni_rho(1));%2*P*C_0_skf_0(2).*rand(1)+C_0_skf_0(2)-P*C_0_skf_0(2);
Cf0 = (ct_rec(1)-Vb*c_b(t_0))/(1-Vb)-Cp0;
if Cf0<0
    break 
end 
C_0_skf = [Cf0,Cp0];

Initialguess3 = [k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb]';
IG2(i,:) = Initialguess3;
[k1x_new,k2x_new,k3x_new,k4x_new,Cf0x_new,Cp0x_new,Vbx_new,Cx_new,Cxdata_new,relerr_new,iter_new,contnew,nitnew] =reconstruction_Skf_New(ct_rec(1:7),c_b, Vi, time(T:22),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb);
k1x_newvett(i) = k1x_new;
k2x_newvett(i) = k2x_new;
k3x_newvett(i)= k3x_new;
k4x_newvett(i) = k4x_new;
Cf0x_newvett(i) = Cf0x_new;
Cp0x_newvett(i) = Cp0x_new;
Vbx_newvett(i) = Vbx_new;
iternew2(i) = iter_new;
contnew2(i) = cont_new;
nitnew2(i) = nitnew;
Param_new2(i,:) = [k1x_newvett(i),k2x_newvett(i),k3x_newvett(i),k4x_newvett(i),Cf0x_newvett(i),Cp0x_newvett(i),Vbx_newvett(i)]';


% Table with all the parameters
 
j =2*i; 
AA2(1:7,j-1) = Initialguess3;
AA2(1:7,j) = Param_new2(i,:);

BB2(1:7,i) = Param_new2(i,:);

% Compute C_t with the calculated parameters
% [ct_rec_new, c_comp_rec_new] = forward_Skf(c_b, Vbx_new, Vi, time(T:end), t_0, [Cf0x_new,Cp0x_new], k1x_new, k2x_new, k3x_new, k4x_new);
% figure(i+5)
% plot(time(T:end),ct_rec,'Linewidth',1)
% hold on
% plot(time(T:end),ct_rec_new,'Linewidth',1)
% hold on
% legend('C_T simulated','C_T all unknown','Location','southeast')
% axis([0,time(end),min([min(ct_rec),min(ct_rec_new)]),max([max(ct_rec),max(ct_rec_new)])])
% xlabel('Time')
% ylabel('C_T')
% title('Comparison between simulated data')

AA2 = round(AA2,4);


end


%MinCf = min(BB2(5,:));
MinCp = min(BB2(6,:));

%MaxCf = max(BB2(5,:));
MaxCp = max(BB2(6,:));


%% Third V_b and k4
for i = 1:N
Vb = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                +par_loguni_rho(1));%2*P*mice.Vb.*rand(1)+mice.Vb-P*mice.Vb;
%Vi = 2*P*mice.Vi.*rand(1)+mice.Vi-P*mice.Vi; 
k4 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
                +par_loguni_rho(1));%2*P*mice.k4.*rand(1)+mice.k4-P*mice.k4;

% Fixed in smaller intervals
k1 = Mink1+(Maxk1-Mink1).*rand(1); 
k2 = Mink2+(Maxk2-Mink2).*rand(1); 
k3 = Mink3+(Maxk3-Mink3).*rand(1); 
%Cf0 = MinCf+(MaxCf-MinCf).*rand(1);
Cp0 = MinCp+(MaxCp-MinCp).*rand(1);
%Cf0 = ct_rec(1)-Cp0;
Cf0 = (ct_rec(1)-Vb*c_b(t_0))/(1-Vb)-Cp0;
C_0_skf =[Cf0,Cp0];%C_0_skf_0;

Initialguess2 = [k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb]';
IG3(i,:) = Initialguess2;

[k1x_new,k2x_new,k3x_new,k4x_new,Cf0x_new,Cp0x_new,Vbx_new,Cx_new,Cxdata_new,relerr_new,iter_new, cont_new,nitnew] =reconstruction_Skf_New(ct_rec,c_b, Vi, time(T:end),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb);
k1x_newvett(i) = k1x_new;
k2x_newvett(i) = k2x_new;
k3x_newvett(i)=k3x_new;
k4x_newvett(i)= k4x_new;
Cf0x_newvett(i) = Cf0x_new;
Cp0x_newvett(i)= Cp0x_new;
Vbx_newvett(i)= Vbx_new;
iternew3(i) = iter_new;
contnew3(i) = cont_new;
nitnew3(i) = nitnew;
Param_new3(i,:) = [k1x_newvett(i),k2x_newvett(i),k3x_newvett(i),k4x_newvett(i),Cf0x_newvett(i),Cp0x_newvett(i),Vbx_newvett(i)]';



% Table with all the parameters
 
 j =2*i; 
AA3(1:7,j-1) = Initialguess2;
AA3(1:7,j) = Param_new3(i,:);
BB3(1:7,i) = Param_new3(i,:);

% Compute C_t with the calculated parameters
%[ct_reck, c_comp_reck] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, C_0_skf_0, k1x, k2x, k3x, k4x);
[ct_rec_new, c_comp_rec_new] = forward_Skf(c_b, Vbx_new, Vi, time(T:end), t_0, [Cf0x_new,Cp0x_new], k1x_new, k2x_new, k3x_new, k4x_new);
%[ct_rec_noVb, c_comp_rec_noVb] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, [Cf0x_noVb,Cp0x_noVb], k1x_noVb, k2x_noVb, k3x_noVb, k4x_noVb);

AA3 = round(AA3,4);

figure(i+10)
plot(time(T:end),ct_rec,'Linewidth',1)
hold on
plot(time(T:end),ct_rec_new,'Linewidth',1)
hold on
legend('C_T simulated','C_T all unknown','Location','southeast')
axis([0,time(end),min([min(ct_rec),min(ct_rec_new)]),max([max(ct_rec),max(ct_rec_new)])])
xlabel('Time')
ylabel('C_T')
title('Comparison between simulated data')

% Errors
ERR1 = norm(ct_rec-ct_rec_new,2)/norm(ct_rec,2);
ERR(i,:) = ERR1;
AA3(8,j) = ERR(i,:);
AA3 = round(AA3,4);

% Relative errors bewteen the parameters
Err2(i,1) = abs(k1x_new-k1_0)/k1_0;
Err2(i,2) = abs(k2x_new-k2_0)/k2_0;
Err2(i,3) = abs(k3x_new-k3_0)/k3_0;
Err2(i,4) = abs(k4x_new-k4_0)/k4_0;
Err2(i,5) = abs(Cf0x_new-C_0_skf_0(1))/C_0_skf_0(1);
Err2(i,6) = abs(Cp0x_new-C_0_skf_0(2))/C_0_skf_0(2);
Err2(i,7) = abs(Vbx_new-Vb_0)/Vb_0;
Err2 = round(Err2,4); 


Errinit(i,1) = abs(k1-k1_0)/k1_0;
Errinit(i,2) = abs(k2-k2_0)/k2_0;
Errinit(i,3) = abs(k3-k3_0)/k3_0;
Errinit(i,4) = abs(k4-k4_0)/k4_0;
Errinit(i,5) = abs(Cf0-C_0_skf_0(1))/C_0_skf_0(1);
Errinit(i,6) = abs(Cp0-C_0_skf_0(2))/C_0_skf_0(2);
Errinit(i,7) = abs(Vb-Vb_0)/Vb_0;
Errinit= round(Errinit,4); 

end

internit = [nitnew1;iternew1;contnew1]';
internit2 = [nitnew2;iternew2;contnew2]';
internit3 = [nitnew3;iternew3;contnew3]';
save('New_20_1','AA','AA2','AA3','internit','internit2','internit3','Err2','Param_new1','Param_new2','Param_new3','Errinit','IG1','IG2','IG3')

%% Boxplot
figure(100)
boxplot(Err2, {'k1','k2','k3','k4','Cf0','Cp0','Vb'})
%ylim([0,5])

figure(101)
boxplot(Errinit, {'k1','k2','k3','k4','Cf0','Cp0','Vb'})


Prova = Err2;
Prova(find(contnew3 == 3),:)=[];

figure(200)
boxplot(Prova, {'k1','k2','k3','k4','Cf0','Cp0','Vb'})

toc