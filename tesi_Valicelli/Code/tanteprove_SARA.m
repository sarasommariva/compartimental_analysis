clc
clear 
close all

%% Step 1 Load data
% folder_results = fullfile('.', 'data');
folder_results = './';
load(fullfile(folder_results, 'data_mice.mat'))

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


%%
AA = zeros(7,30);

% Modified values: random values in a range of [k-Pk,k+Pk]
P =1;
% Search random in an interval [k-P*k,k*Pk] 
% if a = k-P*k; b = k+P*k; then 
% K_mod = a+(b-a).*rand(1); with (b-a) = 2*P*k

N = 2; % Tests number
for i = 1:N
Vb = 2*P*mice.Vb.*rand(1)+mice.Vb-P*mice.Vb;
Vi = 2*P*mice.Vi.*rand(1)+mice.Vi-P*mice.Vi; 
k1 = 2*P*mice.k1.*rand(1)+mice.k1-P*mice.k1; 
k2 = 2*P*mice.k2.*rand(1)+mice.k2-P*mice.k2; 
k3 = 2*P*mice.k3.*rand(1)+mice.k3-P*mice.k3; 
k4 = 2*P*mice.k4.*rand(1)+mice.k4-P*mice.k4;


% Random tra 0 e 1
% Vi = 0;
% k1 = rand(1);
% k2 = rand(1); 
% k3 =rand(1);
% k4 = rand(1);


% log uniforme 
% par_loguni_rho = [-2,0];
%  k1 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
%                 +par_loguni_rho(1));
% Vi = 0; 
% k2 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
%                 +par_loguni_rho(1));
% k3 =  10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
%                 +par_loguni_rho(1));
% k4 = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(1) ...
%                 +par_loguni_rho(1));

% % C_0_skf=2*P*C_0_skf_0.*rand(1)+C_0_skf_0-P*C_0_skf_0;
 %Cp0 = C_0_skf_0(2);
 Cp0 =C_0_skf_0(2);% 2*P*C_0_skf_0(2).*rand(1)+C_0_skf_0(2)-P*C_0_skf_0(2);
 Cf0 = (ct_rec(1)-Vb*c_b(t_0))/(1-Vb)-Cp0;
 if Cf0<0
     break
 end
 Cf0vett(i) = Cf0;
 Cp0vett(i) = Cp0;
 C_0_skf = [Cf0,Cp0];


Ref = [k1_0,k2_0,k3_0,k4_0,C_0_skf_0(1),C_0_skf_0(2),Vb_0]';
Initialguess = [k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb]';

IG(i,:) = Initialguess;

[k1x,k2x,k3x,k4x,Cx,Cxdata,relerr,iter] = reconstruction_Skf(ct_rec,c_b, Vb_0, Vi_0, time(T:end),t_0,C_0_skf_0,k1,k2,k3,k4);
Param = [k1x,k2x,k3x,k4x]';

[k1x_new,k2x_new,k3x_new,k4x_new,Cf0x_new,Cp0x_new,Vbx_new,Cx_new,Cxdata_new,relerr_new,iter_new,cont_new,nit_new] =reconstruction_Skf_New(ct_rec,c_b, Vi, time(T:end),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2),Vb);
k1x_newvett(i) = k1x_new;
disp(k1x_new)
k2x_newvett(i) = k2x_new;
k3x_newvett(i)=k3x_new;
k4x_newvett(i)= k4x_new;
Cf0x_newvett(i) = Cf0x_new;
Cp0x_newvett(i)= Cp0x_new;
Vbx_newvett(i)= Vbx_new;
iternew(i) = iter_new;
contnew(i) = cont_new;
nitnew(i) = nit_new;

Param_new(i,:) = [k1x_newvett(i),k2x_newvett(i),k3x_newvett(i),k4x_newvett(i),Cf0x_newvett(i),Cp0x_newvett(i),Vbx_newvett(i)]';

[k1x_noVb,k2x_noVb,k3x_noVb,k4x_noVb,Cf0x_noVb,Cp0x_noVb,Cx_noVb,Cxdata_noVb,relerr_noVb,iter_noVb,cont_noVb,nit_noVb] = reconstruction_Skf_noVb(ct_rec,c_b, Vb_0, Vi_0, time(T:end),t_0,C_0_skf,k1,k2,k3,k4,C_0_skf(1),C_0_skf(2));
k1x_noVbvett(i)=k1x_noVb;
k2x_noVbvett(i)=k2x_noVb;
k3x_noVbvett(i)=k3x_noVb;
k4x_noVbvett(i)=k4x_noVb;
Cf0x_noVbvett(i)=Cf0x_noVb;
Cp0x_noVbvett(i)=Cp0x_noVb;
iternoVb(i) = iter_noVb;
contnoVb(i) = cont_noVb;
nitnoVb(i) = nit_noVb;

Param_noVb(i,:) = [k1x_noVbvett(i),k2x_noVbvett(i),k3x_noVbvett(i),k4x_noVbvett(i),Cf0x_noVbvett(i),Cp0x_noVbvett(i)]';   



% Table with all the parameters
 
 j =3*i; 
% AA(1:end,1) = Ref;
 AA(1:7,j-2) = Initialguess;
%AA(1:4,i) = Param;
AA(1:7,j-1) = Param_new(i,:);
AA(1:6,j) = Param_noVb(i,:);

% Compute C_t with the calculated parameters
%[ct_reck, c_comp_reck] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, C_0_skf_0, k1x, k2x, k3x, k4x);
[ct_rec_new, c_comp_rec_new] = forward_Skf(c_b, Vbx_new, Vi, time(T:end), t_0, [Cf0x_new,Cp0x_new], k1x_new, k2x_new, k3x_new, k4x_new);
[ct_rec_noVb, c_comp_rec_noVb] = forward_Skf(c_b, Vb_0, Vi_0, time(T:end), t_0, [Cf0x_noVb,Cp0x_noVb], k1x_noVb, k2x_noVb, k3x_noVb, k4x_noVb);



figure(i)
plot(time(T:end),ct_rec,'Linewidth',1)
hold on
% plot(time(T:end),c_t(T:end),'--k','Linewidth',1)%(T:end)
% hold on
% plot(time(T:end),ct_reck,'Linewidth',1)
% hold on
plot(time(T:end),ct_rec_new,'Linewidth',1)
hold on
plot(time(T:end),ct_rec_noVb,'Linewidth',1)
legend('C_T simulated','C_T all unknown','C_T V_b unknown','Location','southeast')
%legend('C_T simulated','C_T data','C_T k_1-k_4 known','C_T all unknown','C_T V_b unknown','Location','southeast')
axis([0,time(end),min([min(ct_rec),min(ct_rec_new),min(ct_rec_noVb)]),max([max(ct_rec),max(ct_rec_new),max(ct_rec_noVb)])])
xlabel('Time')
ylabel('C_T')
title('Comparison between simulated data')
saveas(gcf,'prove.m')
% Errors
ERR1 = [norm(ct_rec-ct_rec_new,2)/norm(ct_rec,2),norm(ct_rec-ct_rec_noVb,2)/norm(ct_rec,2)];
ERR(i,:) = ERR1;
AA(8,j-1:j) = ERR(i,:);
AA = round(AA,4);


% Relative errors bewteen the parameters
Err2(i,1:2) = [abs(k1x_new-k1_0)/k1_0,abs(k1x_noVb-k1_0)/k1_0];
Err2(i,3:4) = [abs(k2x_new-k2_0)/k2_0,abs(k2x_noVb-k2_0)/k2_0];
Err2(i,5:6) = [abs(k3x_new-k3_0)/k3_0,abs(k3x_noVb-k3_0)/k3_0];
Err2(i,7:8) = [abs(k4x_new-k4_0)/k4_0,abs(k4x_noVb-k4_0)/k4_0];
Err2(i,9:10) = [abs(Cf0x_new-C_0_skf_0(1))/C_0_skf_0(1),abs(Cf0x_noVb-C_0_skf_0(1))/C_0_skf_0(1)];
Err2(i,11:12) = [abs(Cp0x_new-C_0_skf_0(2))/C_0_skf_0(2),abs(Cp0x_noVb-C_0_skf_0(2))/C_0_skf_0(2)];
Err2(i,13:14) = [abs(Vbx_new-Vb_0)/Vb_0,0];
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

iternew = iternew';
iternoVb = iternoVb';
contnew = contnew';
contnoVb = contnoVb';
itercont = [iternew,contnew,iternoVb,contnoVb];

nitnew = nitnew';
nitnoVb = nitnoVb';
internit = [nitnew,iternew,contnew,nitnoVb,iternoVb,contnoVb];
%save('20prove_3','AA','contnew','iternew','contnoVb','iternoVb','Err2','Param_new','Param_noVb','IG','Errinit','itercont','internit')

%% Boxplot errori
figure(100)
PROVA = [Err2(:,1),Err2(:,3),Err2(:,5),Err2(:,7),Err2(:,9),Err2(:,11),Err2(:,13)];
boxplot(PROVA, {'k1','k2','k3','k4','Cf0','Cp0','Vb'})
%saveas(gcf,'boxplot.fig')

figure(101)
PROVA2 = [Err2(:,2),Err2(:,4),Err2(:,6),Err2(:,8),Err2(:,10),Err2(:,12),Err2(:,14)];
boxplot(PROVA2, {'k1','k2','k3','k4','Cf0','Cp0','Vb'})
%saveas(gcf,'boxplot_noVb.fig')

% Errori initial guess
figure(102)
boxplot(Errinit, {'k1','k2','k3','k4','Cf0','Cp0','Vb'})
%saveas(gcf,'boxplot_IG.fig')

toc

mnew = mean(Param_new,1);
mnew = round(mnew,4);
dnew = std(Param_new,1);
dnew = round(dnew,4);

mnoVb = mean(Param_noVb,1);
mnoVb = round(mnoVb,4);
dnoVb = std(Param_noVb,1);
dnoVb = round(dnoVb,4);

%% Boxplot senza cont = 3
Prova = PROVA;
Prova(find(contnew == 3),:)=[];
figure(200)
boxplot(Prova)
figure(201) 
Prova2 = PROVA2;
Prova2(find(contnoVb == 3),:)=[];
boxplot(Prova2)
