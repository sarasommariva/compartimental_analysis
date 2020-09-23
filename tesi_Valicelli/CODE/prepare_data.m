clc
clear
close all

%% Step 1. Define folders and files
im = 16; % We'll select the im-th control mouse.

folder_data_original = fullfile('..', '..', 'topi_diabetici', 'results');
folder_results = fullfile('.', 'data');

file_data_original = fullfile(folder_data_original, 'diab_vs_control_preproc.mat');
file_result = fullfile(folder_data_original, sprintf('result_mouse_%d_temp.mat', im));

%% Step 2. Read and save variables
% 2.a. Load files
load(file_data_original)
load(file_result)

% 2.b. Read data
mice_or = mice_group_control(im);
time_or = time'; 

mice.name = 'toy_data';
mice.c_b = mice_or.if;
mice.c_t = mice_or.post;
mice.time = time_or(1:numel(mice.c_b));

% 2.c. Read solution inverse modeling
[~, ik] = min(ris.post.relerr_skf);

mice.k1 = ris.post.k1x_skf(ik);
mice.k2 = ris.post.k2x_skf(ik);
mice.k3 = ris.post.k3x_skf(ik);
mice.k4 = ris.post.k4x_skf(ik);
mice.Vb = 0.04;
mice.Vi = 0;

save(fullfile(folder_results, 'data_mice.mat'), 'mice');

%% Step 3. Some checks
clc
clear 
close all

folder_results = fullfile('.', 'data');
load(fullfile(folder_results, 'data_mice.mat'))

% 3.1. Solve forward model with the selected parameters
time = mice.time;
c_b = mice.c_b;
Vb = mice.Vb;
Vi = mice.Vi; 
k1 = mice.k1; k2 = mice.k2; k3 = mice.k3; k4 = mice.k4;
t_0 = 0;
C_0_skf = [0, 0];

c_b = @(tt)(interp1([0 time],[0 c_b'], tt,'linear',0));

ct_rec = forward_Skf(c_b, Vb, Vi, time, t_0, C_0_skf, k1, k2, k3, k4);

figure
subplot(2, 1, 1)
plot(mice.time, mice.c_b, 'linewidth', 2)
subplot(2, 1, 2)
plot(mice.time, mice.c_t, 'linewidth', 2)
hold on
plot(mice.time, ct_rec, 'r', 'linewidth', 2)







